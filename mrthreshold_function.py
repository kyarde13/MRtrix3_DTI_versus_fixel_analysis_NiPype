# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:
# -*- coding: utf-8 -*-
import os.path as op
from nipype.interfaces.base import CommandLine, traits, TraitedSpec, File, CommandLineInputSpec
from nipype.interfaces.mrtrix3.base import MRTrix3BaseInputSpec, MRTrix3Base

class MRThresholdInputSpec(MRTrix3BaseInputSpec):
    in_file = traits.String(
        exists=True, argstr="%s", mandatory=True, position=-2, desc="input image to be thresholded"
    )
    out_file = File(argstr="%s", mandatory=True, position=-1, desc="output binary mask")

class MRThresholdOutputSpec(TraitedSpec):
    out_file = File(exists=True, desc="output image")

class MRThreshold(MRTrix3Base):
    """
    Create bitwise image by thresholding image intensity
    Example
    -------
    >>> import nipype.interfaces.mrtrix3 as mrt
    >>> mrthreshold = mrt.mrthreshold()
    >>> mrthreshold.inputs.in_file1 = 'a.mif' 
    >>> mrthreshold.cmdline                             # doctest: +ELLIPSIS
    'mrthreshold a.mif mask.mif'
    >>> mrthreshold.run()                               # doctest: +SKIP
    """
    _cmd = "mrthreshold"
    input_spec = MRThresholdInputSpec
    output_spec = MRThresholdOutputSpec

    def _list_outputs(self):
        outputs = self.output_spec().get()
        outputs["out_file"] = op.abspath(self.inputs.out_file)
        return outputs