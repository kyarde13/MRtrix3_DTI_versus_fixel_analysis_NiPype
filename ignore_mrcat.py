# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:
# -*- coding: utf-8 -*-

from nipype.interfaces.base import CommandLine, traits, TraitedSpec, File, CommandLineInputSpec

class MRCatInputSpec(MRTrix3BaseInputSpec):
in_file1 = File(
exists=True, argstr="%s", mandatory=True, position=-2, desc="input image 1"
    )
in_file2 = File(
exists=True, argstr="%s", mandatory=True, position=-3, desc="input image 2"
    )
out_file = File(argstr="%s", mandatory=True, position=-1, desc="output image")
axis = traits.Int(
        0, argstr="-axis %d", desc="specfied axis to perform the operation along"
    )
    )
class MRCatOutputSpec(TraitedSpec):
out_file = File(exists=True, desc="output image")

class MRCat(MRTrix3Base):
"""
    Concatenate several images into one
    along a specified axis
    Example
    -------
    >>> import nipype.interfaces.mrtrix3 as mrt
    >>> mrcat = mrt.mrcat()
    >>> mrcat.inputs.in_file1 = 'a.mif' 
    >>> mrcat.inputs.in_file2 = 'b.mif'
    >>> mrmath.inputs.axis = 3
    >>> mrmath.inputs.out_file = 'concatenated.mif'
    >>> mrcat.cmdline                             # doctest: +ELLIPSIS
    'mrcat a.mif b.mif -axis 3 concatenated.mif'
    >>> mrcat.run()                               # doctest: +SKIP
    """
_cmd = "mrcat"
input1_spec = MRCatInputSpec
input2_spec = MRCatInputSpec
output_spec = MRCatOutputSpec

def _list_outputs(self):
outputs = self.output_spec().get()
outputs["out_file"] = op.abspath(self.inputs.out_file)
return outputs