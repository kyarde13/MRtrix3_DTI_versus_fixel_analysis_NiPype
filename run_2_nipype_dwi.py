#!/usr/bin/env python3

import os.path
import os
import glob
import nipype.interfaces.mrtrix3 as mrt

from nipype.interfaces.utility import IdentityInterface, Function
from nipype.interfaces.io import SelectFiles, DataSink
from nipype.pipeline.engine import Workflow, Node, MapNode

import argparse

import mrcat_function as mrcatfunc
import preproc_function as preprocfunc
import preprocess as preprocess
import fod2fixel_function as fod2fixelfunc
import fixel2peaks_function as fixel2peaksfunc
import mrcalc_function as mrcalcfunc
import utils as utils
import tensor2metric_function as tensor2metricfunc

#import fixel2peaks_func as fixel2peaksfunc
#import mrcalc_func as mrcalcfunc
#import tensor2metric_func as tensor2metricfunc

def create_DWI_workflow(
    subject_list,
    bids_dir,
    work_dir,
    out_dir,
    bids_templates,
):

    # create initial workflow
    wf = Workflow(name='DWI', base_dir=work_dir)

    # use infosource to iterate workflow across subject list
    n_infosource = Node(
        interface=IdentityInterface(
            fields=['subject_id']
        ),
        name="subject_source"
        # input: 'subject_id'
        # output: 'subject_id'
    )
    # runs the node with subject_id = each element in subject_list
    n_infosource.iterables = ('subject_id', subject_list)

    # select matching files from bids_dir
    n_selectfiles = Node(
        interface=SelectFiles(
            templates=bids_templates,
            base_directory=bids_dir
        ),
        name='get_subject_data'
    )
    wf.connect([
        (n_infosource, n_selectfiles, [('subject_id', 'subject_id_p')])
    ])

  
# DWIDenoise
# https://nipype.readthedocs.io/en/latest/api/generated/nipype.interfaces.mrtrix3.preprocess.html
    n_denoise = Node(
        interface=mrt.DWIDenoise(),
        name='n_denoise'
    )
    wf.connect([
        (n_selectfiles, n_denoise, [('DWI_all', 'in_file')])
    ])

    # datasink
    n_datasink = Node(
        interface=DataSink(base_directory=out_dir),
        name='datasink'
    )

    wf.connect([
        (n_selectfiles, n_datasink, [('all_b0_PA', 'all_b0_PA_unchanged')]),
        (n_denoise, n_datasink, [('out_file', 'all_b0_PA_denoised')])
    ])

########## I'VE ADDED IN ##########################################################################
    # MRDeGibbs
    # https://nipype.readthedocs.io/en/latest/api/generated/nipype.interfaces.mrtrix3.preprocess.html
    n_degibbs = Node(
        interface=mrt.MRDeGibbs(),
        name='n_degibbs'
    )
    wf.connect([
        (n_denoise, n_degibbs, [('out_file', 'in_file')])
    ])

    wf.connect([
        (n_degibbs, n_datasink, [('out_file', 'all_b0_PA_degibbs')])
    ])

   # DWI Extract
    n_dwiextract = Node(
        interface=mrt.DWIExtract(
            bzero=True,
            out_file='b0vols.mif'
        ),
        name='n_dwiextract'
    )

    wf.connect([
        (n_degibbs, n_dwiextract, [('out_file', 'in_file')])
    ])

    wf.connect([
        (n_dwiextract, n_datasink, [('out_file', 'noddi_b0_degibbs')])
    ])

    # MRcat
    n_mrcat = Node(
        interface=mrcatfunc.MRCat(
            #axis=3,
            out_file = 'b0s.mif'
        ),
        name='n_mrcat'
    )

    # Connect DTI_B0_PA to mrcat node
    wf.connect([
        (n_selectfiles, n_mrcat, [('DTI_B0_PA', 'in_file1')])
    ])

    wf.connect([
        (n_dwiextract, n_mrcat, [('out_file', 'in_file2')])
    ])

    # Output the mrcat file into file 'noddi_and_PA_b0s.mif'
    wf.connect([
        (n_mrcat, n_datasink, [('out_file', 'noddi_and_PA_b0s.mif')])
    ])

    # DWIfslpreproc
    n_dwifslpreproc = Node(
        interface=preprocfunc.DWIFslPreProc(
            out_file = 'preprocessedDWIs.mif',
            use_header = True
        ),
        name='n_dwifslpreproc'
    )

    # Connect output of degibbs to dwifslpreproc node
    wf.connect([
        (n_degibbs, n_dwifslpreproc, [('out_file', 'in_file')])
    ])
    # Connect output of mrcat to se_epi input
    wf.connect([
        (n_mrcat, n_dwifslpreproc, [('out_file', 'se_epi_file')])
    ])
    # Put output of dwifslpreproc into 'preprocessedDWIs.mif'
    wf.connect([
        (n_dwifslpreproc, n_datasink, [('out_file', 'preprocessedDWIs.mif')])
    ])

    # DWI bias correct
    n_dwibiascorrect = Node(
        interface = preprocess.DWIBiasCorrect(
            #out_file = 'ANTSpreprocessedDWIs.mif',
            use_ants = True
        ),
        name = 'n_dwibiascorrect',
        #use_ants = mrt.DWIBiasCorrect().inputs.use_ants,
        #use_fsl = mrt.DWIBiasCorrect().inputs.use_fsl,
    )

    wf.connect([
        (n_dwifslpreproc, n_dwibiascorrect, [('out_file', 'in_file')])
    ])
    wf.connect([
        (n_dwibiascorrect, n_datasink, [('out_file', 'ANTSpreprocessedDWIs.mif')])
    ]) 

    #DWI2mask
    n_dwi2mask = Node(
        interface=mrt.BrainMask(
            out_file = 'mask.mif'
        ),
        name='n_dwi2mask'
    )
    wf.connect([
        (n_dwibiascorrect, n_dwi2mask, [('out_file', 'in_file')])
    ])
    wf.connect([
        (n_dwi2mask, n_datasink, [('out_file', 'mask.mif')])
    ]) 

    ## A) Fixel-based analysis
    #DWI2response
    n_dwi2response = Node(
        interface=mrt.ResponseSD(
            algorithm = 'dhollander',
            wm_file = 'wm_res.txt',
            gm_file = 'gm_res.txt',
            csf_file = 'csf_res.txt'
        ),
        name='n_dwi2response'
    )

    wf.connect([
        (n_dwibiascorrect, n_dwi2response, [('out_file', 'in_file')])
    ])
    wf.connect([
        (n_dwi2response, n_datasink, [('wm_file', 'wm_res.txt')])
    ]) 
    wf.connect([
        (n_dwi2response, n_datasink, [('gm_file', 'gm_res.txt')])
    ]) 
    wf.connect([
        (n_dwi2response, n_datasink, [('csf_file', 'csf_res.txt')])
    ]) 

    #DWI2fod
    n_dwi2fod = Node(
        interface=mrt.ConstrainedSphericalDeconvolution(
            algorithm = 'msmt_csd',
            wm_odf = 'wmfod.mif',
            gm_odf = 'gmfod.mif',
            csf_odf = 'csffod.mif'
        ),
        name='n_dwi2fod'
    )
    # connect outputs of dwi2fod into dwi2response
    wf.connect([
        (n_dwibiascorrect, n_dwi2fod, [('out_file', 'in_file')])
    ])
    wf.connect([
        (n_dwi2response, n_dwi2fod, [('wm_file', 'wm_txt')])
    ])   
    wf.connect([
        (n_dwi2response, n_dwi2fod, [('gm_file', 'gm_txt')])
    ])  
    wf.connect([
        (n_dwi2response, n_dwi2fod, [('csf_file', 'csf_txt')])
    ])  
    # output wmfod file from dwi2fod
    wf.connect([
        (n_dwi2fod, n_datasink, [('wm_odf', 'wmfod.mif')])
    ])
    wf.connect([
        (n_dwi2fod, n_datasink, [('gm_odf', 'gmfod.mif')])
    ])
    wf.connect([
        (n_dwi2fod, n_datasink, [('csf_odf', 'csffod.mif')])
    ])


    #fod2fixel wmfod.mif wmfixels -fmls_peak_value 0 -fmls_integral 0.10 -afd afd.mif -peak peak.mif -disp disp.mif 
    # OUTPUTS: -afd afd.mif -peak peak.mif -disp disp.mif 
    n_fod2fixel = Node(
        interface= fod2fixelfunc.fod2fixel(
           fmls_peak_value = 0,
           fmls_integral = 0.10,
           out_file = 'wmfixels',
           #afd_file = 'afd.mif',
           peak_file = 'peak.mif',
           disp_file = 'disp.mif'
           
        ),
        name='n_fod2fixel'
    )
    # obtain wm fibre image as input
    wf.connect([
        (n_dwi2fod, n_fod2fixel, [('wm_odf', 'in_file')])
    ])
    # ouputs of fod2fixel
    wf.connect([
        (n_fod2fixel, n_datasink, [('out_file', 'wmfixels')])
    ]) 
    wf.connect([
        (n_fod2fixel, n_datasink, [('afd_file', 'afd.mif')])
    ]) 
    wf.connect([
        (n_fod2fixel, n_datasink, [('peak_file', 'peak.mif')])
    ]) 
    wf.connect([
        (n_fod2fixel, n_datasink, [('disp_file', 'disp.mif')])
    ]) 

    ## Fixel2peaks 
    n_fixel2peaks = Node(
        interface= fixel2peaksfunc.fixel2peaks(
           out_file = 'peaks_wmdirections.mif',
           number = 1
        ),
        name='n_fixel2peaks'
    )
    # obtain directions file in output folder of fod2fixel, as input
    wf.connect([
        (n_fod2fixel, n_fixel2peaks, [('out_file', 'in_file')])
    ])
    # ouputs of fixel2peaks
    wf.connect([
        (n_fixel2peaks, n_datasink, [('out_file', 'peaks_wmdirections.mif')])
    ]) 
   
    #mrmath to find normalised value of peak WM directions
    n_mrmath = Node(
        interface=mrt.MRMath(
            axis = 3,
            operation = 'norm',
            out_file = 'norm_peaks_wmdirections.mif'
        ),
        name='n_mrmath'
    )

    wf.connect([
        (n_fixel2peaks, n_mrmath, [('out_file', 'in_file')])
    ])

    wf.connect([
        (n_mrmath, n_datasink, [('out_file', 'norm_peaks_wmdirections.mif')])
    ]) 

    # mrcalc to divide peak WM direction by normalised value
    n_mrcalc = Node(
        interface=mrcalcfunc.MRCalc(
            operation = 'divide',
            out_file = 'wm_peak_dir.mif'
        ),
        name='n_mrcalc'
    )

    wf.connect([
        (n_fixel2peaks, n_mrcalc, [('out_file', 'in_file1')])
    ])

    wf.connect([
        (n_mrmath, n_mrcalc, [('out_file', 'in_file2')])
    ])

    wf.connect([
        (n_mrcalc, n_datasink, [('out_file', 'WM_peak_dir.mif')])
    ])

    #mrconvert to extract Z component of peak directions
    n_mrconvert2 = Node(
        interface=utils.MRConvert(
            out_file = 'Zpeak_WM_Directions.mif',
            coord = [3, 2]
        ),
        name='n_mrconvert2'
    )

    wf.connect([
        (n_mrcalc, n_mrconvert2, [('out_file', 'in_file')])
    ])

    wf.connect([
        (n_mrconvert2, n_datasink, [('out_file', 'Zpeak_WM_Directions.mif')])
    ]) 

    # mrcalc to find absolute value
    n_mrcalc2 = Node(
        interface=mrcalcfunc.MRCalc(
            operation = 'abs',
            out_file = 'absZpeak_WM_Directions.mif'
        ),
        name='n_mrcalc2'
    )

    wf.connect([
        (n_mrconvert2, n_mrcalc2, [('out_file', 'in_file1')])
    ])

    wf.connect([
        (n_mrcalc2, n_datasink, [('out_file', 'absZpeak_WM_Directions.mif')])
    ]) 

    # mrcalc to get angle by doing inverse cosine
    n_mrcalc3 = Node(
        interface=mrcalcfunc.MRCalc(
            operation = 'acos',
            out_file = 'acosZpeak_WM_Directions.mif'
        ),
        name='n_mrcalc3'
    )

    wf.connect([
        (n_mrcalc2, n_mrcalc3, [('out_file', 'in_file1')])
    ])

    wf.connect([
        (n_mrcalc3, n_datasink, [('out_file', 'acosZpeak_WM_Directions.mif')])
    ]) 
    
    # mrcalc to convert angle to degrees
    n_mrcalc4 = Node(
        interface=mrcalcfunc.MRCalc(
            operation = 'multiply',
            operand = 180,
            out_file = 'Fixel1_Z_angle.mif'
        ),
        name='n_mrcalc4'
    )

    wf.connect([
        (n_mrcalc3, n_mrcalc4, [('out_file', 'in_file1')])
    ])

    wf.connect([
        (n_mrcalc4, n_datasink, [('out_file', 'Fixel1_Z_angle.mif')])
    ]) 

    n_mrcalc5 = Node(
        interface=mrcalcfunc.MRCalc(
            operation = 'divide',
            operand = 3.14159265,
            out_file = 'Fixel1_Z_cos_deg.mif'
        ),
        name='n_mrcalc5'
    )

    wf.connect([
        (n_mrcalc4, n_mrcalc5, [('out_file', 'in_file1')])
    ])

    wf.connect([
        (n_mrcalc5, n_datasink, [('out_file', 'Fixel1_Z_cos_deg.mif')])
    ]) 

    ## B) Tensor-based analysis
    #dwi2tensor
    n_dwi2tensor = Node(
        interface=mrt.FitTensor(
            out_file = 'dti.mif'
        ),
        name='n_dwi2tensor'
    )

    wf.connect([
        (n_dwibiascorrect, n_dwi2tensor, [('out_file', 'in_file')])
    ])

    wf.connect([
        (n_dwi2mask, n_dwi2tensor, [('out_file', 'in_mask')])
    ])

    wf.connect([
        (n_dwi2tensor, n_datasink, [('out_file', 'dt.mif')])
    ]) 

    #tensor2metric 
    n_tensor2metric = Node(
        interface= tensor2metricfunc.tensor2metric(
            modulate = 'none',
            num = 1,
            vector_file = 'eigenvector.mif'
        ),
        name='n_tensor2metric'
    )

    wf.connect([
        (n_dwi2tensor, n_tensor2metric, [('out_file', 'input_file')])
    ])

    wf.connect([
        (n_tensor2metric, n_datasink, [('vector_file', 'eigenvector.mif')])
    ]) 

    #mrconvert to get Z eigenvector
    n_mrconvert3 = Node(
        interface=utils.MRConvert(
            coord = [3, 2],
            out_file = 'eigenvectorZ.mif'
        ),
        name='n_mrconvert3'
    )

    wf.connect([
        (n_tensor2metric, n_mrconvert3, [('vector_file', 'in_file')])
    ])

    wf.connect([
        (n_mrconvert3, n_datasink, [('out_file', 'eigenvectorZ.mif')])
    ]) 

    #ALL SUBSEQUENT STEPS GET ANGLE IN DEGREES
    # mrcalc to find absolute value
    n_mrcalc6 = Node(
        interface=mrcalcfunc.MRCalc(
            operation = 'abs',
            out_file = 'abs_eigenvectorZ.mif'
        ),
        name='n_mrcalc6'
    )

    wf.connect([
        (n_mrconvert3, n_mrcalc6, [('out_file', 'in_file1')])
    ])

    wf.connect([
        (n_mrcalc6, n_datasink, [('out_file', 'abs_eigenvectorZ.mif')])
    ]) 

    # mrcalc to get angle by doing inverse cosine
    n_mrcalc7 = Node(
        interface=mrcalcfunc.MRCalc(
            operation = 'acos',
            out_file = 'acos_eigenvectorZ.mif'
        ),
        name='n_mrcalc7'
    )

    wf.connect([
        (n_mrcalc6, n_mrcalc7, [('out_file', 'in_file1')])
    ])

    wf.connect([
        (n_mrcalc7, n_datasink, [('out_file', 'acos_eigenvectorZ.mif')])
    ]) 

    # mrcalc to convert angle to degrees
    n_mrcalc8 = Node(
        interface=mrcalcfunc.MRCalc(
            operation = 'multiply',
            operand = 180,
            out_file = 'degrees_eigenvectorZ.mif'
        ),
        name='n_mrcalc8'
    )

    wf.connect([
        (n_mrcalc7, n_mrcalc8, [('out_file', 'in_file1')])
    ])

    wf.connect([
        (n_mrcalc8, n_datasink, [('out_file', 'degrees_eigenvectorZ.mif')])
    ]) 

    n_mrcalc9 = Node(
        interface=mrcalcfunc.MRCalc(
            operation = 'divide',
            operand = 3.14159265,
            out_file = 'dti_z_cos_deg.mif'
        ),
        name='n_mrcalc9'
    )

    wf.connect([
        (n_mrcalc8, n_mrcalc9, [('out_file', 'in_file1')])
    ])

    wf.connect([
        (n_mrcalc9, n_datasink, [('out_file', 'dti_z_cos_deg.mif')])
    ]) 

    # Difference image between fixel based and tensor based outputs
    n_mrcalc10 = Node(
        interface=mrcalcfunc.MRCalc(
            operation = 'subtract',
            out_file = 'diff_imag_tensor_minus_fixel.mif'
        ),
        name='n_mrcalc10'
    )

    wf.connect([
        (n_mrcalc9, n_mrcalc10, [('out_file', 'in_file1')])
    ])

    wf.connect([
        (n_mrcalc5, n_mrcalc10, [('out_file', 'in_file2')])
    ])

    wf.connect([
        (n_mrcalc10, n_datasink, [('out_file', 'diff_imag_tensor_minus_fixel.mif')])
    ]) 

#################################################################################3
    return wf


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="DWI processing pipeline",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument(
        '--bids_dir',
        required=True,
        help='bids directory'
    )

    parser.add_argument(
        '--subjects',
        default=None,
        const=None,
        nargs='*',
        help='list of subjects as seen in bids_dir'
    )

    parser.add_argument(
        '--work_dir',
        required=True,
        help='work directory'
    )

    parser.add_argument(
        '--out_dir',
        required=True,
        help='output directory'
    )

    parser.add_argument(
        '--debug',
        dest='debug',
        action='store_true',
        help='debug mode'
    )


    parser.add_argument(
        '--pbs',
        dest='pbs',
        action='store_true',
        help='use PBS graph'
    )

    args = parser.parse_args()

    # environment variables
    os.environ["FSLOUTPUTTYPE"] = "NIFTI_GZ"
    os.environ["PATH"] += os.pathsep + os.path.join(os.path.dirname(os.path.abspath(__file__)), "scripts")

    this_dir = os.path.dirname(os.path.abspath(__file__))
    if "PYTHONPATH" in os.environ: os.environ["PYTHONPATH"] += os.pathsep + this_dir
    else:                          os.environ["PYTHONPATH"]  = this_dir

    if args.debug:
        from nipype import config
        config.enable_debug_mode()
        config.set('execution', 'stop_on_first_crash', 'true')
        config.set('execution', 'remove_unnecessary_outputs', 'false')
        config.set('execution', 'keep_inputs', 'true')
        config.set('logging', 'workflow_level', 'DEBUG')
        config.set('logging', 'interface_level', 'DEBUG')
        config.set('logging', 'utils_level', 'DEBUG')

    if not args.subjects:
        subject_list = [subj for subj in os.listdir(args.bids_dir) if 'sub' in subj]
    else:
        subject_list = args.subjects

    bids_templates = {
        'all_b0_PA': '{subject_id_p}/dwi/all_b0_PA.mif',
        'DWI_all': '{subject_id_p}/dwi/DWI_all.mif',
        'DTI_B0_PA': '{subject_id_p}/dwi/DTI_B0_PA',
    }

    wf = create_DWI_workflow(
        subject_list=subject_list,
        bids_dir=os.path.abspath(args.bids_dir),
        work_dir=os.path.abspath(args.work_dir),
        out_dir=os.path.abspath(args.out_dir),
        bids_templates=bids_templates
    )

    os.makedirs(os.path.abspath(args.work_dir), exist_ok=True)
    os.makedirs(os.path.abspath(args.out_dir), exist_ok=True)

    wf.write_graph(graph2use='flat', format='png', simple_form=False)
    # run workflow

    if args.pbs:
        wf.run(
            plugin='PBSGraph',
            plugin_args={
                'qsub_args': '-A UQ-CAI -q Short -l nodes=1:ppn=1,mem=5GB,vmem=5GB,walltime=00:30:00'
            }
        )
    else:
        wf.run(
            plugin='MultiProc',
            plugin_args={
                'n_procs': int(os.environ["NCPUS"]) if "NCPUS" in os.environ else int(os.cpu_count())
            }
        )