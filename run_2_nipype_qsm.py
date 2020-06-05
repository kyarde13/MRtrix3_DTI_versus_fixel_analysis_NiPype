#!/usr/bin/env python3

import os.path
import os
from nipype.interfaces.fsl import BET, ImageMaths, ImageStats, MultiImageMaths, CopyGeom, Merge, UnaryMaths
from nipype.interfaces.utility import IdentityInterface, Function
from nipype.interfaces.io import SelectFiles, DataSink
from nipype.pipeline.engine import Workflow, Node, MapNode

import nipype_interface_tgv_qsm as tgv
import nipype_interface_romeo as romeo
import nipype_interface_bestlinreg as bestlinreg
import nipype_interface_applyxfm as applyxfm
import nipype_interface_makehomogeneous as makehomogeneous

import argparse


def create_qsm_workflow(
    subject_list,
    bids_dir,
    work_dir,
    out_dir,
    atlas_dir,
    masking='bet',
    bids_templates={
        'mag1': '{subject_id_p}/anat/*gre*E01*magnitude*.nii.gz',
        'mag': '{subject_id_p}/anat/*gre*magnitude*.nii.gz',
        'phs': '{subject_id_p}/anat/*gre*phase*.nii.gz',
        'params': '{subject_id_p}/anat/*gre*phase*.json'
    },
    qsm_threads=1
):

    # create initial workflow
    wf = Workflow(name='qsm', base_dir=work_dir)

    # use infosource to iterate workflow across subject list
    n_infosource = Node(
        interface=IdentityInterface(
            fields=['subject_id']
        ),
        name="infosource"
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
        name='selectfiles'
        # output: ['mag', 'phs', 'params']
    )
    wf.connect([
        (n_infosource, n_selectfiles, [('subject_id', 'subject_id_p')])
    ])

    # scale phase data
    mn_stats = MapNode(
        # -R : <min intensity> <max intensity>
        interface=ImageStats(op_string='-R'),
        iterfield=['in_file'],
        name='stats_node',
        # output: 'out_stat'
    )
    mn_phs_range = MapNode(
        interface=ImageMaths(),
        name='phs_range_node',
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

    def repeat(name):
        return name

    # brain extraction
    if masking == 'bet1echo':
        n_bet = Node(
            interface=BET(frac=0.4, mask=True, robust=True),
            iterfield=['in_file'],
            name='bet'
            # output: 'mask_file'
        )
        wf.connect([
            (n_selectfiles, n_bet, [('mag1', 'in_file')])
        ])

        mn_mask = Node(
            interface=Function(
                input_names=['name'],
                output_names=['mask_file'],
                function=repeat
            ),
            iterfield=['name'],
            name='join'
        )
        wf.connect([
            (n_bet, mn_mask, [('mask_file', 'name')])
        ])
    if masking == 'bet':
        mn_homogeneity_filter = MapNode(
            interface=makehomogeneous.MakeHomogeneousInterface(),
            iterfield=['in_file'],
            name='makehomogeneous'
        )
        wf.connect([
            (n_selectfiles, mn_homogeneity_filter, [('mag', 'in_file')])
        ])

        mn_bet = MapNode(
            interface=BET(frac=0.4, mask=True, robust=True),
            iterfield=['in_file'],
            name='bet'
            # output: 'mask_file'
        )
        wf.connect([
            (mn_homogeneity_filter, mn_bet, [('out_file', 'in_file')])
        ])

        mn_mask = MapNode(
            interface=Function(
                input_names=['name'],
                output_names=['mask_file'],
                function=repeat
            ),
            iterfield=['name'],
            name='join'
        )
        wf.connect([
            (mn_bet, mn_mask, [('mask_file', 'name')])
        ])
    elif masking == 'romeo':
        # ROMEO only operates on stacked .nii files
        n_stacked_magnitude = Node(
            interface=Merge(
                dimension='t',
                output_type='NIFTI'
            ),
            name="stack_magnitude",
            iterfield=['in_files']
            # output: 'merged_file'
        )
        wf.connect([
            (n_selectfiles, n_stacked_magnitude, [('mag', 'in_files')])
        ])
        n_stacked_phase = Node(
            interface=Merge(
                dimension='t',
                output_type='NIFTI'
            ),
            name="stack_phase",
            iterfield=['in_files']
            # output: 'merged_file'
        )
        wf.connect([
            (mn_phs_range, n_stacked_phase, [('out_file', 'in_files')])
        ])

        n_romeo = Node(
            interface=romeo.RomeoInterface(
                weights_threshold=200
            ),
            iterfield=['in_file', 'echo_times'],
            name='romeo'
            # output: 'out_file'
        )
        wf.connect([
            (n_stacked_phase, n_romeo, [('merged_file', 'in_file')]),
            (mn_params, n_romeo, [('EchoTime', 'echo_times')])
        ])
        mn_mask = MapNode(
            interface=Function(
                input_names=['name'],
                output_names=['mask_file'],
                function=repeat
            ),
            iterfield=['name'],
            name='join'
        )
        wf.connect([
            (n_romeo, mn_mask, [('out_file', 'name')])
        ])
    elif masking == 'atlas-based':
        n_selectatlas = Node(
            interface=SelectFiles(
                templates={
                    'template': '*template*',
                    'mask': '*mask*'
                },
                base_directory=atlas_dir
            ),
            name='selectatlas'
            # output: ['template', 'mask']
        )

        mn_bestlinreg = MapNode(
            interface=bestlinreg.NiiBestLinRegInterface(),
            iterfield=['in_fixed', 'in_moving'],
            name='bestlinreg'
            # output: out_transform
        )

        wf.connect([
            (n_selectfiles, mn_bestlinreg, [('mag', 'in_fixed')]),
            (n_selectatlas, mn_bestlinreg, [('template', 'in_moving')])
        ])

        mn_applyxfm = MapNode(
            interface=applyxfm.NiiApplyMincXfmInterface(),
            iterfield=['in_file', 'in_like', 'in_transform'],
            name='applyxfm'
            # output: out_file
        )

        wf.connect([
            (n_selectatlas, mn_applyxfm, [('mask', 'in_file')]),
            (n_selectfiles, mn_applyxfm, [('mag', 'in_like')]),
            (mn_bestlinreg, mn_applyxfm, [('out_transform', 'in_transform')])
        ])

        mn_mask = MapNode(
            interface=Function(
                input_names=['name'],
                output_names=['mask_file'],
                function=repeat
            ),
            iterfield=['name'],
            name='join'
        )
        wf.connect([
            (mn_applyxfm, mn_mask, [('out_file', 'name')])
        ])

    # qsm processing
    mn_qsm_iterfield = ['phase_file', 'TE', 'b0']
    
    # if using a multi-echo masking method, add mask_file to iterfield
    if masking != 'bet1echo': mn_qsm_iterfield.append('mask_file')

    mn_qsm = MapNode(
        interface=tgv.QSMappingInterface(
            iterations=1000,
            alpha=[0.0015, 0.0005],
            erosions=2 if masking == 'romeo' else 5,
            num_threads=qsm_threads
        ),
        iterfield=mn_qsm_iterfield,
        name='qsm_node'
        # output: 'out_file'
    )

    # args for PBS
    mn_qsm.plugin_args = {
        'qsub_args': f'-A UQ-CAI -q Short -l nodes=1:ppn={qsm_threads},mem=20gb,vmem=20gb,walltime=03:00:00',
        'overwrite': True
    }

    wf.connect([
        (mn_params, mn_qsm, [('EchoTime', 'TE')]),
        (mn_params, mn_qsm, [('MagneticFieldStrength', 'b0')]),
        (mn_mask, mn_qsm, [('mask_file', 'mask_file')]),
        (mn_phs_range, mn_qsm, [('out_file', 'phase_file')])
    ])

    # mask processing
    def generate_multiimagemaths_lists(in_files):
        if type(in_files) == str: # fix for single-echo masking
            return in_files, [in_files], '-mul %s '
        in_file = in_files[0]
        operand_files = in_files[1:]
        op_string = '-add %s '
        op_string = len(operand_files) * op_string
        return in_file, operand_files, op_string

    n_generate_add_masks_lists = Node(
        interface=Function(
            input_names=['in_files'],
            output_names=[
                'list_in_file',
                'list_operand_files',
                'list_op_string'
            ],
            function=generate_multiimagemaths_lists
        ),
        name='generate_add_masks_lists_node'
    )

    n_add_masks = Node(
        interface=MultiImageMaths(),
        name="add_masks_node"
        # output: 'out_file'
    )

    wf.connect([
        (mn_mask, n_generate_add_masks_lists, [('mask_file', 'in_files')]),
        (n_generate_add_masks_lists, n_add_masks, [('list_in_file', 'in_file')]),
        (n_generate_add_masks_lists, n_add_masks, [('list_operand_files', 'operand_files')]),
        (n_generate_add_masks_lists, n_add_masks, [('list_op_string', 'op_string')])
    ])

    # qsm post-processing
    n_generate_add_qsms_lists = Node(
        interface=Function(
            input_names=['in_files'],
            output_names=['list_in_file', 'list_operand_files', 'list_op_string'],
            function=generate_multiimagemaths_lists
        ),
        name='generate_add_qsms_lists_node'
        # output: 'out_file'
    )

    n_add_qsms = Node(
        interface=MultiImageMaths(),
        name="add_qsms_node"
        # output: 'out_file'
    )
    wf.connect([
        (mn_qsm, n_generate_add_qsms_lists, [('out_file', 'in_files')]),
        (n_generate_add_qsms_lists, n_add_qsms, [('list_in_file', 'in_file')]),
        (n_generate_add_qsms_lists, n_add_qsms, [('list_operand_files', 'operand_files')]),
        (n_generate_add_qsms_lists, n_add_qsms, [('list_op_string', 'op_string')])
    ])

    # divide qsm by mask
    n_final_qsm = Node(
        interface=ImageMaths(op_string='-div'),
        name="divide_added_qsm_by_added_masks"
        # output: 'out_file'
    )
    wf.connect([
        (n_add_qsms, n_final_qsm, [('out_file', 'in_file')]),
        (n_add_masks, n_final_qsm, [('out_file', 'in_file2')])
    ])

    # datasink
    n_datasink = Node(
        interface=DataSink(base_directory=bids_dir, container=out_dir),
        name='datasink'
    )

    wf.connect([
        (n_add_masks, n_datasink, [('out_file', 'mask_sum')]),
        (n_add_qsms, n_datasink, [('out_file', 'qsm_sum')]),
        (n_final_qsm, n_datasink, [('out_file', 'qsm_final_default')]),
        (mn_qsm, n_datasink, [('out_file', 'qsm_singleEchoes')]),
        (mn_mask, n_datasink, [('mask_file', 'mask_singleEchoes')])
    ])

    return wf


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="QSM processing pipeline",
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
        '--masking',
        default='bet',
        const='bet',
        nargs='?',
        choices=['bet', 'bet1echo', 'romeo', 'atlas-based'],
        help='masking strategy'
    )

    parser.add_argument(
        '--atlas_dir',
        default=os.path.join(
            os.path.dirname(os.path.abspath(__file__)),
            "atlas"
        ),
        const=os.path.join(
            os.path.dirname(os.path.abspath(__file__)),
            "atlas"
        ),
        nargs='?',
        help='atlas directory',
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

    wf = create_qsm_workflow(
        subject_list=subject_list,
        bids_dir=os.path.abspath(args.bids_dir),
        work_dir=os.path.abspath(args.work_dir),
        out_dir=os.path.abspath(args.out_dir),
        masking=args.masking,
        atlas_dir=os.path.abspath(args.atlas_dir),
        qsm_threads=16 if args.pbs else 1
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

    #wf.write_graph(graph2use='flat', format='png', simple_form=False)
    #wf.run(plugin='PBS', plugin_args={'-A UQ-CAI -l nodes=1:ppn=16,mem=5gb,vmem=5gb, walltime=30:00:00'})
