# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:
# -*- coding: utf-8 -*-
import os.path as op
from nipype.interfaces.base import CommandLine, traits, TraitedSpec, File, Directory, CommandLineInputSpec, Undefined
from nipype.interfaces.mrtrix3.base import MRTrix3BaseInputSpec, MRTrix3Base

class fixel2peaksInputSpec(MRTrix3BaseInputSpec):
    
    #in_file is pathname to a directions file in the wmfixels folder
    in_file = Directory(argstr="%s/directions.mif", usedefault=True, mandatory=True, position=-2, desc="output folder")
    out_file = File(
        argstr="%s", mandatory=True, position=-1, desc="input dwi image"
    )
    number = traits.Int(
        1, argstr="-number %d", desc="maximum number of fixels in each voxel"
    )
    

class fixel2peaksOutputSpec(TraitedSpec):
    out_file = File(exists=True, desc="output image")

class fixel2peaks(MRTrix3Base):
    """
    Convert data in the fixel directory format into a 4D image of 3-vectors
    Example
    -------
    >>> fixel2peaks.inputs.in_file = 'wmfixels/directions' 
    >>> fixel2peaks.inputs.number = 1
    >>> fixel2peaks.cmdline                            
    'fixel2peaks wmfixels/directions -number 1 peaks_wmdirections.mif
    >>> fixel2peaks.run()                               
    """
    _cmd = "fixel2peaks"
    input_spec = fixel2peaksInputSpec
    output_spec = fixel2peaksOutputSpec

    def _list_outputs(self):
        outputs = self.output_spec().get()
        outputs["out_file"] = op.abspath(self.inputs.out_file)
        return outputs   

    