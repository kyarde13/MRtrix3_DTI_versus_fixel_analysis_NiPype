#!/usr/bin/env python3

import os.path
import os
import glob
import nipype.interfaces.mrtrix3 as mrt

from nipype.interfaces.utility import IdentityInterface, Function
from nipype.interfaces.io import SelectFiles, DataSink
from nipype.pipeline.engine import Workflow, Node, MapNode

import argparse

# Import all created functions that implement MRtrix commands
import mrcat_function as mrcatfunc
import preproc_function as preprocfunc
import preprocess as preprocess
import fod2fixel_function as fod2fixelfunc
import fixel2peaks_function as fixel2peaksfunc
import mrcalc_function as mrcalcfunc
import utils as utils
import tensor2metric_function as tensor2metricfunc
import mrthreshold_function as mrthresholdfunc

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

########## IMPLEMENTING MRTRIX COMMANDS FOR IMAGE ANALYSIS #######################################

## 1) Preprocessing of Data
# https://nipype.readthedocs.io/en/latest/api/generated/nipype.interfaces.mrtrix3.preprocess.html
# DWIDenoise to remove Gaussian noise
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
    # output denoised data into 'DWI_all_denoised'
    wf.connect([
        (n_selectfiles, n_datasink, [('all_b0_PA', 'all_b0_PA_unchanged')]),
        (n_denoise, n_datasink, [('out_file', 'DWI_all_denoised')])
    ])

# MRDeGibbs to remove Gibbs ringing artifact
    n_degibbs = Node(
        interface=mrt.MRDeGibbs(
            out_file = 'DWI_all_denoised_degibbs.mif'
        ),
        name='n_degibbs'
    )
    # input denoised data into degibbs function
    wf.connect([
        (n_denoise, n_degibbs, [('out_file', 'in_file')])
    ])
    # output degibbs data into 'DWI_all_denoised_degibbs.mif'
    wf.connect([
        (n_degibbs, n_datasink, [('out_file', 'DWI_all_denoised_degibbs.mif')])
    ])

# DWI Extract to extract b0 volumes from multi-b image data
    n_dwiextract = Node(
        interface=mrt.DWIExtract(
            bzero=True,
            out_file='b0vols.mif'
        ),
        name='n_dwiextract'
    )
    # input degibbs data into dwiextract function
    wf.connect([
        (n_degibbs, n_dwiextract, [('out_file', 'in_file')])
    ])
    # output extracted b0 volume from degibbs data (contains multiple b values)
    wf.connect([
        (n_dwiextract, n_datasink, [('out_file', 'noddi_b0_degibbs')])
    ])

# MRcat to combine b0 volumes from input image and reverse phase encoded data
    n_mrcat = Node(
        interface=mrcatfunc.MRCat(
            #axis=3,
            out_file = 'b0s.mif'
        ),
        name='n_mrcat'
    )
    # input DTI images (all b0 volumes; reverse phase encoded) for concatenating
    wf.connect([
        (n_selectfiles, n_mrcat, [('DTI_B0_PA', 'in_file1')])
    ])
    # input b0 volumes from NODDI data for concatenating
    wf.connect([
        (n_dwiextract, n_mrcat, [('out_file', 'in_file2')])
    ])
    # output the mrcat file into 'noddi_and_PA_b0s.mif'
    wf.connect([
        (n_mrcat, n_datasink, [('out_file', 'noddi_and_PA_b0s.mif')])
    ])

# DWIfslpreproc for image pre-processing using FSL's eddy tool
    n_dwifslpreproc = Node(
        interface=preprocfunc.DWIFslPreProc(
            out_file = 'preprocessedDWIs.mif',
            use_header = True
        ),
        name='n_dwifslpreproc'
    )
    # output of degibbs as input for preprocessing
    wf.connect([
        (n_degibbs, n_dwifslpreproc, [('out_file', 'in_file')])
    ])
    # output of mrcat (extracted b0 volumes) as se_epi input
    wf.connect([
        (n_mrcat, n_dwifslpreproc, [('out_file', 'se_epi_file')])
    ])
    # output of dwifslpreproc into 'preprocessedDWIs.mif'
    wf.connect([
        (n_dwifslpreproc, n_datasink, [('out_file', 'preprocessedDWIs.mif')])
    ])

# DWI bias correct for B1 field inhomogeneity correction
    n_dwibiascorrect = Node(
        interface = preprocess.DWIBiasCorrect(
            use_ants = True
        ),
        name = 'n_dwibiascorrect',
    )
    # input preprocessed data
    wf.connect([
        (n_dwifslpreproc, n_dwibiascorrect, [('out_file', 'in_file')])
    ])
    # output biascorrect data into 'ANTSpreprocessedDWIs.mif'
    wf.connect([
        (n_dwibiascorrect, n_datasink, [('out_file', 'ANTSpreprocessedDWIs.mif')])
    ]) 

# DWI2mask to compute whole brain mask from bias corrected data
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

##################################################################################
## 2) Fixel-based analysis
# DWI2response for etimation of response function for spherical deconvolution
    n_dwi2response = Node(
        interface=mrt.ResponseSD(
            algorithm = 'dhollander',
            wm_file = 'wm_res.txt',
            gm_file = 'gm_res.txt',
            csf_file = 'csf_res.txt'
        ),
        name='n_dwi2response'
    )
    # input bias corrected data for response function estimation
    wf.connect([
        (n_dwibiascorrect, n_dwi2response, [('out_file', 'in_file')])
    ])
    # output WM, GM, CSF response text files
    wf.connect([
        (n_dwi2response, n_datasink, [('wm_file', 'wm_res.txt')])
    ]) 
    wf.connect([
        (n_dwi2response, n_datasink, [('gm_file', 'gm_res.txt')])
    ]) 
    wf.connect([
        (n_dwi2response, n_datasink, [('csf_file', 'csf_res.txt')])
    ]) 

# DWI2fod for fibre orientation distribution estimation (FOD)
    n_dwi2fod = Node(
        interface=mrt.ConstrainedSphericalDeconvolution(
            algorithm = 'msmt_csd',
            wm_odf = 'wmfod.mif',
            gm_odf = 'gmfod.mif',
            csf_odf = 'csffod.mif'
        ),
        name='n_dwi2fod'
    )
    # utilise dwi2fod response files as input
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
    # output WM, GM and CSF FODs for saving 
    wf.connect([
        (n_dwi2fod, n_datasink, [('wm_odf', 'wmfod.mif')])
    ])
    wf.connect([
        (n_dwi2fod, n_datasink, [('gm_odf', 'gmfod.mif')])
    ])
    wf.connect([
        (n_dwi2fod, n_datasink, [('csf_odf', 'csffod.mif')])
    ])

# Mrconvert to extract z component (component w.r.t main field) of WM FOD
    n_mrconvert_fod = Node(
        interface=utils.MRConvert(
            out_file = 'Zwmfod.mif',
            coord = [3, 0]
        ),
        name='n_mrconvert_fod'
    )
    # utilise WM FOD as input
    wf.connect([
        (n_dwi2fod, n_mrconvert_fod, [('wm_odf', 'in_file')])
    ])
    # output z component of WM FOD
    wf.connect([
        (n_mrconvert_fod, n_datasink, [('out_file', 'Zwmfod.mif')])
    ]) 

# MRcat to concatenate all WM, GM, CSF FOD files to see their distributions throughout Brain
    n_mrcat_fod = Node(
        interface=mrcatfunc.MRCat(
            out_file = 'vf.mif'
        ),
        name='n_mrcat_fod'
    )
    # connect Zwmfod, gmfod and csffod as inputs
    wf.connect([
        (n_mrconvert_fod, n_mrcat_fod, [('out_file', 'in_file1')])
    ])
    wf.connect([
        (n_dwi2fod, n_mrcat_fod, [('gm_odf', 'in_file2')])
    ])
    wf.connect([
        (n_dwi2fod, n_mrcat_fod, [('csf_odf', 'in_file3')])
    ])
    # output the mrcat file into file 'vf.mif'
    wf.connect([
        (n_mrcat_fod, n_datasink, [('out_file', 'vf.mif')])
    ]) 

# fod2fixel 
# Perform segmentation of continuous FODs to produce discrete fixels
# OUTPUTS: -afd afd.mif -peak peak.mif -disp disp.mif 
    n_fod2fixel = Node(
        interface= fod2fixelfunc.fod2fixel(
           out_file = 'wmfixels',
           #afd_file = 'afd.mif',
           peak_file = 'peak.mif',
           disp_file = 'disp.mif'
           
        ),
        name='n_fod2fixel'
    )
    # let the peak value parameter be trialed as multiple values
    n_fod2fixel.iterables = ('fmls_peak_value', [0, 0.10, 0.50])
    n_fod2fixel.iterables = ('fmls_integral', [0, 0.10, 0.50])

    # obtain WM fibre image as input
    wf.connect([
        (n_dwi2fod, n_fod2fixel, [('wm_odf', 'in_file')])
    ])
    # ouputs of fod2fixel saved
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

# fixel2peaks to convert data in the fixel directory format into 4D image of 3-vectors
    n_fixel2peaks = Node(
        interface= fixel2peaksfunc.fixel2peaks(
           out_file = 'peaks_wmdirections.mif'
        ),
        name='n_fixel2peaks'
    )
    # look at multiple values for maximum number of fixels in each voxel
    n_fixel2peaks.iterables = ('number', [1, 2, 3])

    # obtain directions file in output folder of fod2fixel, as input
    wf.connect([
        (n_fod2fixel, n_fixel2peaks, [('out_file', 'in_file')])
    ])
    # ouput of fixel2peaks saved in peaks_wmdirections.mif'
    wf.connect([
        (n_fixel2peaks, n_datasink, [('out_file', 'peaks_wmdirections.mif')])
    ]) 
   
# MRmath to find normalised value of peak WM directions
    n_mrmath = Node(
        interface=mrt.MRMath(
            axis = 3,
            operation = 'norm',
            out_file = 'norm_peaks_wmdirections.mif'
        ),
        name='n_mrmath'
    )
    # input peak fixel data
    wf.connect([
        (n_fixel2peaks, n_mrmath, [('out_file', 'in_file')])
    ])
    # output saved into 'norm_peaks_wmdirections.mif'
    wf.connect([
        (n_mrmath, n_datasink, [('out_file', 'norm_peaks_wmdirections.mif')])
    ]) 

# MRcalc to divide peak WM direction by normalised value
    n_mrcalc = Node(
        interface=mrcalcfunc.MRCalc(
            operation = 'divide',
            out_file = 'wm_peak_dir.mif'
        ),
        name='n_mrcalc'
    )
    # fixel2peaks image as input 1
    wf.connect([
        (n_fixel2peaks, n_mrcalc, [('out_file', 'in_file1')])
    ])
    # normalised fixel2peak image as input 2
    wf.connect([
        (n_mrmath, n_mrcalc, [('out_file', 'in_file2')])
    ])
    # save output image as 'WM_peak_dir.mif'
    wf.connect([
        (n_mrcalc, n_datasink, [('out_file', 'WM_peak_dir.mif')])
    ])

# MRconvert to extract Z component of peak directions
    n_mrconvert2 = Node(
        interface=utils.MRConvert(
            out_file = 'Zpeak_WM_Directions.mif',
            coord = [3, 2]
        ),
        name='n_mrconvert2'
    )
    # input normalised peak direction file
    wf.connect([
        (n_mrcalc, n_mrconvert2, [('out_file', 'in_file')])
    ])
    # save ouptut as 'Zpeak_WM_Directions.mif'
    wf.connect([
        (n_mrconvert2, n_datasink, [('out_file', 'Zpeak_WM_Directions.mif')])
    ]) 

# MRcalc to find absolute value of peak fibre directions
    n_mrcalc2 = Node(
        interface=mrcalcfunc.MRCalc(
            operation = 'abs',
            out_file = 'absZpeak_WM_Directions.mif'
        ),
        name='n_mrcalc2'
    )
    # input z peaks image
    wf.connect([
        (n_mrconvert2, n_mrcalc2, [('out_file', 'in_file1')])
    ])
    # save output as 'absZpeak_WM_Directions.mif'
    wf.connect([
        (n_mrcalc2, n_datasink, [('out_file', 'absZpeak_WM_Directions.mif')])
    ]) 

# MRcalc to get angle by doing inverse cosine
    n_mrcalc3 = Node(
        interface=mrcalcfunc.MRCalc(
            operation = 'acos',
            out_file = 'acosZpeak_WM_Directions.mif'
        ),
        name='n_mrcalc3'
    )
    # input normalised z component of peaks image
    wf.connect([
        (n_mrcalc2, n_mrcalc3, [('out_file', 'in_file1')])
    ])
    # save ouput as 'acosZpeak_WM_Directions.mif'
    wf.connect([
        (n_mrcalc3, n_datasink, [('out_file', 'acosZpeak_WM_Directions.mif')])
    ]) 
    
# MRcalc to convert angle of peak fibre (w.r.t z axis), to degrees
    n_mrcalc4 = Node(
        interface=mrcalcfunc.MRCalc(
            operation = 'multiply',
            operand = 180,
            out_file = 'Fixel1_Z_angle.mif'
        ),
        name='n_mrcalc4'
    )
    # input inverse cosine image of peak fibre
    wf.connect([
        (n_mrcalc3, n_mrcalc4, [('out_file', 'in_file1')])
    ])
    # output image as 'Fixel1_Z_angle.mif'
    wf.connect([
        (n_mrcalc4, n_datasink, [('out_file', 'Fixel1_Z_angle.mif')])
    ]) 

# MRcalc to divide by pi to finish converting from radians to degrees
    n_mrcalc5 = Node(
        interface=mrcalcfunc.MRCalc(
            operation = 'divide',
            operand = 3.14159265,
            out_file = 'Fixel1_Z_cos_deg.mif'
        ),
        name='n_mrcalc5'
    )
    # input image multiplied by 180
    wf.connect([
        (n_mrcalc4, n_mrcalc5, [('out_file', 'in_file1')])
    ])
    # save output as 'Fixel1_Z_cos_deg.mif'
    wf.connect([
        (n_mrcalc5, n_datasink, [('out_file', 'Fixel1_Z_cos_deg.mif')])
    ]) 

##################################################################################
## 3) Tensor-based analysis
# dwi2tensor to compute tensor from biascorrected DWI image
    n_dwi2tensor = Node(
        interface=mrt.FitTensor(
            out_file = 'dti.mif'
        ),
        name='n_dwi2tensor'
    )
    # input bias corrected image
    wf.connect([
        (n_dwibiascorrect, n_dwi2tensor, [('out_file', 'in_file')])
    ])
    # utilise mask to only compute tensors for regions of Brain
    wf.connect([
        (n_dwi2mask, n_dwi2tensor, [('out_file', 'in_mask')])
    ])
    # output data into 'dt.mif'
    wf.connect([
        (n_dwi2tensor, n_datasink, [('out_file', 'dt.mif')])
    ]) 

# tensor2metric to convert tensors to generate maps of tensor-derived parameters
    n_tensor2metric = Node(
        interface= tensor2metricfunc.tensor2metric(
            modulate = 'none',
            num = 1,
            vector_file = 'eigenvector.mif'
        ),
        name='n_tensor2metric'
    )
    # input tensor image
    wf.connect([
        (n_dwi2tensor, n_tensor2metric, [('out_file', 'input_file')])
    ])
    # save output eigenvectors of the diffusion tensor
    wf.connect([
        (n_tensor2metric, n_datasink, [('vector_file', 'eigenvector.mif')])
    ]) 

# MRconvert to get eigenvector w.r.t z direction (main field)
    n_mrconvert3 = Node(
        interface=utils.MRConvert(
            coord = [3, 2],
            out_file = 'eigenvectorZ.mif'
        ),
        name='n_mrconvert3'
    )
    # input eigenvector file from tensor2metric
    wf.connect([
        (n_tensor2metric, n_mrconvert3, [('vector_file', 'in_file')])
    ])
    # save output as 'eigenvectorZ.mif'
    wf.connect([
        (n_mrconvert3, n_datasink, [('out_file', 'eigenvectorZ.mif')])
    ]) 

# ALL SUBSEQUENT STEPS GET ANGLE IN DEGREES
# MRcalc to find absolute value of z eigenvector file
    n_mrcalc6 = Node(
        interface=mrcalcfunc.MRCalc(
            operation = 'abs',
            out_file = 'abs_eigenvectorZ.mif'
        ),
        name='n_mrcalc6'
    )
    # z eigenvector image as input
    wf.connect([
        (n_mrconvert3, n_mrcalc6, [('out_file', 'in_file1')])
    ])
    # save output as 'abs_eigenvectorZ.mif'
    wf.connect([
        (n_mrcalc6, n_datasink, [('out_file', 'abs_eigenvectorZ.mif')])
    ]) 

# MRcalc to get angle by doing inverse cosine
    n_mrcalc7 = Node(
        interface=mrcalcfunc.MRCalc(
            operation = 'acos',
            out_file = 'acos_eigenvectorZ.mif'
        ),
        name='n_mrcalc7'
    )
    # input absolute value of z eigenvector image
    wf.connect([
        (n_mrcalc6, n_mrcalc7, [('out_file', 'in_file1')])
    ])
    # save output as 'acos_eigenvectorZ.mif'
    wf.connect([
        (n_mrcalc7, n_datasink, [('out_file', 'acos_eigenvectorZ.mif')])
    ]) 

# MRcalc to convert angle to degrees
    n_mrcalc8 = Node(
        interface=mrcalcfunc.MRCalc(
            operation = 'multiply',
            operand = 180,
            out_file = 'degrees_eigenvectorZ.mif'
        ),
        name='n_mrcalc8'
    )
    # input inverse cosine image of z eigenvector
    wf.connect([
        (n_mrcalc7, n_mrcalc8, [('out_file', 'in_file1')])
    ])
    # save output as 'degrees_eigenvectorZ.mif'
    wf.connect([
        (n_mrcalc8, n_datasink, [('out_file', 'degrees_eigenvectorZ.mif')])
    ]) 

# MRcalc to divide by pi to finish converting from radians to degrees
    n_mrcalc9 = Node(
        interface=mrcalcfunc.MRCalc(
            operation = 'divide',
            operand = 3.14159265,
            out_file = 'dti_z_cos_deg.mif'
        ),
        name='n_mrcalc9'
    )
    # input z eigenvector image multiplied by 180
    wf.connect([
        (n_mrcalc8, n_mrcalc9, [('out_file', 'in_file1')])
    ])
    # save output as 'dti_z_cos_deg.mif'
    wf.connect([
        (n_mrcalc9, n_datasink, [('out_file', 'dti_z_cos_deg.mif')])
    ]) 

# MRcalc to give difference image between fixel based and tensor based outputs
    n_mrcalc10 = Node(
        interface=mrcalcfunc.MRCalc(
            operation = 'subtract',
            out_file = 'diff_imag_tensor_minus_fixel.mif'
        ),
        name='n_mrcalc10'
    )
    # input tensor based image of whole Brain
    wf.connect([
        (n_mrcalc9, n_mrcalc10, [('out_file', 'in_file1')])
    ])
    # input fixel based image of Brain
    wf.connect([
        (n_mrcalc5, n_mrcalc10, [('out_file', 'in_file2')])
    ])
    # output difference image as 'diff_imag_tensor_minus_fixel.mif'
    wf.connect([
        (n_mrcalc10, n_datasink, [('out_file', 'diff_imag_tensor_minus_fixel.mif')])
    ]) 

#####################################################################################
## 4) Tensor based analysis on WM fibres only (NOT WHOLE BRAIN TENSORS)

# MRthreshold to create WM mask from WM FOD (created earlier)
    n_mrthreshold = Node(
        interface=mrthresholdfunc.MRThreshold(
            out_file = 'thresholded_wmfod.mif'
        ),
        name='n_mrthreshold'
    )
    # input WM FOD
    wf.connect([
        (n_dwi2fod, n_mrthreshold, [('wm_odf', 'in_file')])
    ])
    # output thresholded WM FOD
    wf.connect([
        (n_mrthreshold, n_datasink, [('out_file', 'thresholded_wmfod.mif')])
    ]) 

# MRconvert to extract 1st volume of thresholded WM FOD
    n_mrconvert4 = Node(
        interface=utils.MRConvert(
            coord = [3, 0],
            out_file = 'WMmask.mif'
        ),
        name='n_mrconvert4'
    )
    # input thresholded wmfod
    wf.connect([
        (n_mrthreshold, n_mrconvert4, [('out_file', 'in_file')])
    ])
    # save output as 'WMmask.mif'
    wf.connect([
        (n_mrconvert4, n_datasink, [('out_file', 'WMmask.mif')])
    ]) 

# MRcalc to multiple WM mask with dti image to get tensors only of WM regions
    n_mrcalc11 = Node(
        interface=mrcalcfunc.MRCalc(
            operation = 'multiply',
            out_file = 'WM_dt.mif'
        ),
        name='n_mrcalc11'
    )
    # WM mask as input 1
    wf.connect([
        (n_mrconvert4, n_mrcalc11, [('out_file', 'in_file1')])
    ])
    # dti image as input 2
    wf.connect([
        (n_dwi2tensor, n_mrcalc11, [('out_file', 'in_file2')])
    ]) 
    # save output as 'WM_dt.mif'
    wf.connect([
        (n_mrcalc11, n_datasink, [('out_file', 'WM_dt.mif')])
    ]) 

# tensor2metric to convert tensors to generate maps of tensor-derived parameters
    n_tensor2metric2 = Node(
        interface= tensor2metricfunc.tensor2metric(
            modulate = 'none',
            num = 1,
            vector_file = 'WMeigenvector.mif'
        ),
        name='n_tensor2metric2'
    )
    # input tensor image
    wf.connect([
        (n_mrcalc11, n_tensor2metric2, [('out_file', 'input_file')])
    ])
    # save output eigenvectors of the diffusion tensor
    wf.connect([
        (n_tensor2metric2, n_datasink, [('vector_file', 'WMeigenvector.mif')])
    ])

# MRconvert to get eigenvector w.r.t z direction (main field)
    n_mrconvert5 = Node(
        interface=utils.MRConvert(
            coord = [3, 2],
            out_file = 'WMeigenvectorZ.mif'
        ),
        name='n_mrconvert5'
    )
    # input eigenvector file from tensor2metric
    wf.connect([
        (n_tensor2metric2, n_mrconvert5, [('vector_file', 'in_file')])
    ])
    # save output as 'eigenvectorZ.mif'
    wf.connect([
        (n_mrconvert5, n_datasink, [('out_file', 'WMeigenvectorZ.mif')])
    ]) 

# ALL SUBSEQUENT STEPS GET ANGLE IN DEGREES
# MRcalc to find absolute value of z eigenvector file
    n_mrcalc12 = Node(
        interface=mrcalcfunc.MRCalc(
            operation = 'abs',
            out_file = 'WM_abs_eigenvectorZ.mif'
        ),
        name='n_mrcalc12'
    )
    # z eigenvector image as input
    wf.connect([
        (n_mrconvert5, n_mrcalc12, [('out_file', 'in_file1')])
    ])
    # save output as 'WM_abs_eigenvectorZ.mif'
    wf.connect([
        (n_mrcalc12, n_datasink, [('out_file', 'WM_abs_eigenvectorZ.mif')])
    ]) 
    
# MRcalc to get angle by doing inverse cosine
    n_mrcalc13 = Node(
        interface=mrcalcfunc.MRCalc(
            operation = 'acos',
            out_file = 'acos_WMeigenvectorZ.mif'
        ),
        name='n_mrcalc13'
    )
    # input absolute value of z eigenvector image
    wf.connect([
        (n_mrcalc12, n_mrcalc13, [('out_file', 'in_file1')])
    ])
    # save output as 'acos_WMeigenvectorZ.mif'
    wf.connect([
        (n_mrcalc13, n_datasink, [('out_file', 'acos_WMeigenvectorZ.mif')])
    ]) 

# MRcalc to convert angle to degrees
    n_mrcalc14 = Node(
        interface=mrcalcfunc.MRCalc(
            operation = 'multiply',
            operand = 180,
            out_file = 'degrees_WMeigenvectorZ.mif'
        ),
        name='n_mrcalc14'
    )
    # input inverse cosine image of WM z eigenvector
    wf.connect([
        (n_mrcalc13, n_mrcalc14, [('out_file', 'in_file1')])
    ])
    # save output as 'degrees_WMeigenvectorZ.mif'
    wf.connect([
        (n_mrcalc14, n_datasink, [('out_file', 'degrees_WMeigenvectorZ.mif')])
    ]) 

# MRcalc to divide by pi to finish converting from radians to degrees
    n_mrcalc15 = Node(
        interface=mrcalcfunc.MRCalc(
            operation = 'divide',
            operand = 3.14159265,
            out_file = 'WMdti_z_cos_deg.mif'
        ),
        name='n_mrcalc15'
    )
    # input WM z eigenvector image multiplied by 180
    wf.connect([
        (n_mrcalc14, n_mrcalc15, [('out_file', 'in_file1')])
    ])
    # save output as 'WMdti_z_cos_deg.mif'
    wf.connect([
        (n_mrcalc15, n_datasink, [('out_file', 'WMdti_z_cos_deg.mif')])
    ]) 

# MRcalc to give difference image between fixel based and WM tensor based outputs
    n_mrcalc16 = Node(
        interface=mrcalcfunc.MRCalc(
            operation = 'subtract',
            out_file = 'diffImage_WMtensor_minus_fixel.mif'
        ),
        name='n_mrcalc16'
    )
    # input fixel image of Brain
    wf.connect([
        (n_mrcalc15, n_mrcalc16, [('out_file', 'in_file1')])
    ])
    # input tensor image of WM fibres of Brain
    wf.connect([
        (n_mrcalc5, n_mrcalc16, [('out_file', 'in_file2')])
    ])
    # output difference image as 'diff_imag_WMtensor_minus_fixel.mif'
    wf.connect([
        (n_mrcalc16, n_datasink, [('out_file', 'diffImage_WMtensor_minus_fixel.mif')])
    ]) 
######################################################################################
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