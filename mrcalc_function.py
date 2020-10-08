# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:
# -*- coding: utf-8 -*-
import os.path as op
from nipype.interfaces.base import CommandLine, traits, TraitedSpec, File, CommandLineInputSpec
from nipype.interfaces.mrtrix3.base import MRTrix3BaseInputSpec, MRTrix3Base

class MRCalcInputSpec(MRTrix3BaseInputSpec):
    in_file1 = traits.String(
        exists=True, argstr="%s", mandatory=True, position=-3, desc="input image 1"
    )
    # Either perform operation with another fiel (in_file2) or a constant number (operand)
    in_file2 = File(
        exists=True, argstr="%s", mandatory=False, position=-4, desc="input image 2"
    )
    operand = traits.Int(
        0, argstr="-axis %d", mandatory=False, position=-4, desc="specfied number to perform operation with"
    )
    out_file = File(argstr="%s", mandatory=True, position=-1, desc="output image")
    operation = traits.Enum(
        "abs",
        "neg",
        "add",
        "subtract",
        "multiply",
        "divide",
        "min",
        "max",
        "acos",
        argstr="-%s",
        position=-2,
        mandatory=True,
        desc="operation to perform",
    )

class MRCalcOutputSpec(TraitedSpec):
    out_file = File(exists=True, desc="output image")

class MRCalc(MRTrix3Base):
    """
    Concatenate several images into one
    along a specified axis
    Example
    -------
    >>> mrcalc.inputs.in_file1 = peaks_wmdirections.mif' 
    >>> mrcalc.inputs.in_file2 = norm_peaks_wmdirections.mif 
    >>> mrcalc.inputs.operation = div
    >>> mrcalc.inputs.out_file = 'wm_peak_dir.mif '
    >>> mrcalc.cmdline                            
    'mrcalc peaks_wmdirections.mif norm_peaks_wmdirections.mif -div wm_peak_dir.mif 
    >>> mrcalc.run()                               
    """
    _cmd = "mrcalc"
    input_spec = MRCalcInputSpec
    output_spec = MRCalcOutputSpec

    def _list_outputs(self):
        outputs = self.output_spec().get()
        outputs["out_file"] = op.abspath(self.inputs.out_file)
        return outputs