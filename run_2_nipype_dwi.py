#!/usr/bin/env python3

import os.path
import os
import glob
from nipype.interfaces.fsl import BET, ImageMaths, ImageStats, MultiImageMaths, CopyGeom, Merge, UnaryMaths
from nipype.interfaces.utility import IdentityInterface, Function
from nipype.interfaces.io import SelectFiles, DataSink
from nipype.pipeline.engine import Workflow, Node, MapNode

import nipype_interface_tgv_DWI as tgv
import nipype_interface_romeo as romeo
import nipype_interface_bestlinreg as bestlinreg
import nipype_interface_applyxfm as applyxfm
import nipype_interface_makehomogeneous as makehomogeneous
import nipype_interface_nonzeroaverage as nonzeroaverage
import nipype_interface_composite as composite

import argparse


def create_DWI_workflow(
    subject_list,
    bids_dir,
    work_dir,
    out_dir,
    atlas_dir,
    bids_templates,
    masking='bet-multiecho',
    homogeneity_filter=True,
    DWI_threads=1
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
        # output: ['mag', 'phs', 'params']
    )
    wf.connect([
        (n_infosource, n_selectfiles, [('subject_id', 'subject_id_p')])
    ])

    # count the number of echoes
    def get_length(in_):
        return len(in_)
    n_num_echoes = Node(
        interface=Function(
            input_names=['in_'],
            output_names=['num_echoes'],
            function=get_length
        ),
        iterfield=['in_'],
        name='get_num_echoes'
    )
    wf.connect([
        (n_selectfiles, n_num_echoes, [('phs', 'in_')])
    ])

    # scale phase data
    mn_stats = MapNode(
        # -R : <min intensity> <max intensity>
        interface=ImageStats(op_string='-R'),
        iterfield=['in_file'],
        name='get_min_max',
        # output: 'out_stat'
    )
    mn_phs_range = MapNode(
        interface=ImageMaths(suffix="_scaled"),
        name='scale_phase',
        iterfield=['in_file']
        # inputs: 'in_file', 'op_string'
        # output: 'out_file'
    )

    def scale_to_pi(min_and_max):
        from math import pi

        min_value = min_and_max[0][0]
        max_value = min_and_max[0][1]
        fsl_cmd = ""

        # set range to [0, max-min]
        fsl_cmd += "-sub %.10f " % min_value
        max_value -= min_value
        min_value -= min_value

        # set range to [0, 2pi]
        fsl_cmd += "-div %.10f " % (max_value / (2*pi))

        # set range to [-pi, pi]
        fsl_cmd += "-sub %.10f" % pi
        return fsl_cmd

    wf.connect([
        (n_selectfiles, mn_stats, [('phs', 'in_file')]),
        (n_selectfiles, mn_phs_range, [('phs', 'in_file')]),
        (mn_stats, mn_phs_range, [(('out_stat', scale_to_pi), 'op_string')])
    ])

    # read echotime and field strengths from json files
    def read_json(in_file):
        import os
        te = 0.001
        b0 = 7
        if os.path.exists(in_file):
            import json
            with open(in_file, 'rt') as fp:
                data = json.load(fp)
                te = data['EchoTime']
                b0 = data['MagneticFieldStrength']
        return te, b0

    mn_params = MapNode(
        interface=Function(
            input_names=['in_file'],
            output_names=['EchoTime', 'MagneticFieldStrength'],
            function=read_json
        ),
        iterfield=['in_file'],
        name='read_json'
    )
    wf.connect([
        (n_selectfiles, mn_params, [('params', 'in_file')])
    ])
    n_params_e01 = Node(
        interface=Function(
            input_names=['in_file'],
            output_names=['EchoTime', 'MagneticFieldStrength'],
            function=read_json
        ),
        name='read_json1'
    )
    wf.connect([
        (n_selectfiles, n_params_e01, [('params1', 'in_file')])
    ])

    def repeat(in_file):
        return in_file


    # brain extraction
    if 'bet' in masking:
        # homogeneity filter
        n_mag = MapNode(
            interface=Function(
                input_names=['in_file'],
                output_names=['out_file'],
                function=repeat
            ),
            iterfield=['in_file'],
            name='repeat_magnitude'
        )
        if homogeneity_filter:
            mn_homogeneity_filter = MapNode(
                interface=makehomogeneous.MakeHomogeneousInterface(),
                iterfield=['in_file'],
                name='make_homogeneous'
                # output : out_file
            )
            wf.connect([
                (n_selectfiles, mn_homogeneity_filter, [('mag', 'in_file')]),
                (mn_homogeneity_filter, n_mag, [('out_file', 'in_file')])
            ])
        else:
            wf.connect([
                (n_selectfiles, n_mag, [('mag', 'in_file')])
            ])

        bet = MapNode(
            interface=BET(frac=0.4, mask=True, robust=True),
            iterfield=['in_file'],
            name='fsl_bet'
            # output: 'mask_file'
        )

        wf.connect([
            (n_mag, bet, [('out_file', 'in_file')])
        ])

        mn_mask = MapNode(
            interface=Function(
                input_names=['in_file'],
                output_names=['mask_file'],
                function=repeat
            ),
            iterfield=['in_file'],
            name='repeat_mask'
        )
        wf.connect([
            (bet, mn_mask, [('mask_file', 'in_file')])
        ])
    elif masking == 'romeo':
        # per-echo ROMEO masks
        n_romeo = MapNode(
            interface=romeo.RomeoInterface(
                weights_threshold=200
            ),
            iterfield=['in_file', 'echo_time'],
            name='romeo_mask'
            # output: 'out_file'
        )
        wf.connect([
            (n_selectfiles, n_romeo, [('phs', 'in_file')]),
            (mn_params, n_romeo, [('EchoTime', 'echo_time')])
        ])

        n_romeo_maths = MapNode(
            interface=ImageMaths(
                suffix='_ero_dil',
                op_string='-ero -dilM'
            ),
            iterfield=['in_file'],
            name='romeo_ero_dil'
            # input  : 'in_file'
            # output : 'out_file'
        )
        wf.connect([
            (n_romeo, n_romeo_maths, [('out_file', 'in_file')])
        ])

        mn_mask = MapNode(
            interface=Function(
                input_names=['in_file'],
                output_names=['mask_file'],
                function=repeat
            ),
            iterfield=['in_file'],
            name='repeat_mask'
        )
        wf.connect([
            (n_romeo_maths, mn_mask, [('out_file', 'in_file')])
        ])

        # first-echo ROMEO mask for composite inclusions
        n_romeo_e01 = Node(
            interface=romeo.RomeoInterface(
                weights_threshold=200
            ),
            name='romeo_e01_mask'
            # output : 'out_file'
        )
        wf.connect([
            (n_selectfiles, n_romeo_e01, [('phs1', 'in_file')]),
            (n_params_e01, n_romeo_e01, [('EchoTime', 'echo_time')])
        ])
        n_romeo_e01_maths = Node(
            interface=ImageMaths(
                suffix='_ero_dil_fillh',
                op_string='-ero -dilM -fillh'
            ),
            name='romeo_e01_ero_dil_fillh'
            # input  : 'in_file'
            # output : 'out_file'
        )
        wf.connect([
            (n_romeo_e01, n_romeo_e01_maths, [('out_file', 'in_file')])
        ])
    
    # DWI processing
    mn_DWI_iterfield = ['phase_file', 'TE', 'b0']
    
    # if using a multi-echo masking method, add mask_file to iterfield
    if masking not in ['bet-firstecho', 'bet-lastecho']: mn_DWI_iterfield.append('mask_file')

    mn_DWI = MapNode(
        interface=tgv.DWIappingInterface(
            iterations=1000,
            alpha=[0.0015, 0.0005],
            erosions=0 if masking == 'romeo' else 5,
            num_threads=DWI_threads,
            out_suffix='_DWI_recon'
        ),
        iterfield=mn_DWI_iterfield,
        name='DWI'
        # output: 'out_file'
    )

    # args for PBS
    mn_DWI.plugin_args = {
        'qsub_args': f'-A UQ-CAI -q Short -l nodes=1:ppn={DWI_threads},mem=20gb,vmem=20gb,walltime=03:00:00',
        'overwrite': True
    }

    wf.connect([
        (mn_params, mn_DWI, [('EchoTime', 'TE')]),
        (mn_params, mn_DWI, [('MagneticFieldStrength', 'b0')]),
        (mn_mask, mn_DWI, [('mask_file', 'mask_file')]),
        (mn_phs_range, mn_DWI, [('out_file', 'phase_file')])
    ])

    # DWI averaging
    n_DWI_average = Node(
        interface=nonzeroaverage.NonzeroAverageInterface(),
        name='DWI_average'
        # input : in_files
        # output : out_file
    )
    wf.connect([
        (mn_DWI, n_DWI_average, [('out_file', 'in_files')])
    ])

    # datasink
    n_datasink = Node(
        interface=DataSink(base_directory=bids_dir, container=out_dir),
        name='datasink'
    )

    wf.connect([
        (n_DWI_average, n_datasink, [('out_file', 'DWI_average')]),
        (mn_DWI, n_datasink, [('out_file', 'DWI_single')]),
        (mn_mask, n_datasink, [('mask_file', 'mask_single')])
    ])

    if masking == 'romeo':
        mn_DWI_filled = MapNode(
            interface=tgv.DWIappingInterface(
                iterations=1000,
                alpha=[0.0015, 0.0005],
                erosions=5,
                num_threads=DWI_threads,
                out_suffix='_DWI_recon_filled'
            ),
            iterfield=['phase_file', 'TE', 'b0'],
            name='DWI_filled'
            # output: 'out_file'
        )

        # args for PBS
        mn_DWI_filled.plugin_args = {
            'qsub_args': f'-A UQ-CAI -q Short -l nodes=1:ppn={DWI_threads},mem=20gb,vmem=20gb,walltime=03:00:00',
            'overwrite': True
        }

        wf.connect([
            (mn_params, mn_DWI_filled, [('EchoTime', 'TE')]),
            (mn_params, mn_DWI_filled, [('MagneticFieldStrength', 'b0')]),
            (n_romeo_e01_maths, mn_DWI_filled, [('out_file', 'mask_file')]),
            (mn_phs_range, mn_DWI_filled, [('out_file', 'phase_file')])
        ])

        # DWI averaging
        n_DWI_filled_average = Node(
            interface=nonzeroaverage.NonzeroAverageInterface(),
            name='DWI_filled_average'
            # input : in_files
            # output : out_file
        )
        wf.connect([
            (mn_DWI_filled, n_DWI_filled_average, [('out_file', 'in_files')])
        ])

        # composite DWI
        n_composite = Node(
            interface=composite.CompositeNiftiInterface(),
            name='DWI_composite'
        )
        wf.connect([
            (n_DWI_average, n_composite, [('out_file', 'in_file1')]),
            (n_DWI_filled_average, n_composite, [('out_file', 'in_file2')])
        ])

        wf.connect([
            (n_romeo_e01_maths, n_datasink, [('out_file', 'mask_filled_single')]),
            (n_DWI_filled_average, n_datasink, [('out_file', 'DWI_filled_average')]),
            (n_composite, n_datasink, [('out_file', 'DWI_composite')]),
            (mn_DWI_filled, n_datasink, [('out_file', 'DWI_filled_single')])
        ])


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