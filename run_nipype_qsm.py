#!/usr/bin/env python3
from os.path import join as opj
import os
from nipype.interfaces.fsl import (BET, FLIRT, ImageMaths)
from nipype.interfaces.utility import IdentityInterface
from nipype.interfaces.io import SelectFiles, DataSink
from nipype.pipeline.engine import Workflow, Node
from nipype import MapNode
os.environ["FSLOUTPUTTYPE"] = "NIFTI"


experiment_dir = '/QRISdata/Q0538/17042_detection_of_concussion/interim'
output_dir = '/QRISdata/Q0538/17042_detection_of_concussion/processed'
working_dir = '/gpfs1/scratch/30days/uqsbollm/17042_detection_of_concussion'

subject_list = ['sub-S008LCBL']

# Infosource - a function free node to iterate over the list of subject names
infosource = Node(IdentityInterface(fields=['subject_id']), name="infosource")
infosource.iterables = [('subject_id', subject_list)]

# SelectFiles - to grab the data (alternative to DataGrabber)
templates = {'mag': '{subject_id}/anat/*gre_M_echo_*.nii.gz',
             'phs': '{subject_id}/anat/*gre_P_echo_*.nii.gz'}
selectfiles = Node(SelectFiles(templates, base_directory=experiment_dir), name='selectfiles')

# Datasink - creates output folder for important outputs
datasink = Node(DataSink(base_directory=experiment_dir, container=output_dir), name='datasink')

# Create a preprocessing workflow
preproc = Workflow(name='preprocessing')
preproc.base_dir = opj(experiment_dir, working_dir)

bet = MapNode(BET(frac=0.4, mask=True, robust=True),
              name='bet', iterfield=['in_file'])

phs_range = MapNode(ImageMaths(op_string='-div 4096 -mul 6.28318530718 -sub 3.14159265359'),
                    name='phs_range', iterfield=['in_file'])

# Connect all components of the preprocessing workflow
preproc.connect([(infosource, selectfiles, [('subject_id', 'subject_id')]),
                 (selectfiles, bet, [('mag', 'in_file')]),
                 (selectfiles, phs_range, [('phs', 'in_file')]),
                 (bet, datasink, [('mask_file', 'preprocessing.@bet')]),
                 (phs_range, datasink, [('out_file', 'preprocessing.@phs_range')]),
                 ])

preproc.run('MultiProc', plugin_args={'n_procs': 9})