# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:
# -*- coding: utf-8 -*-
import os.path as op
from nipype.interfaces.base import CommandLine, traits, TraitedSpec, File, Directory, CommandLineInputSpec, Undefined
from nipype.interfaces.mrtrix3.base import MRTrix3BaseInputSpec, MRTrix3Base

class fod2fixelInputSpec(MRTrix3BaseInputSpec):
    in_file = File(
        exists=True, argstr="%s", position=-5, mandatory=True, desc="input dwi image"
    )
    #out_file is a folder
    out_file = Directory(argstr="%s", usedefault=True, mandatory=True, position=-4, desc="output folder")
    fmls_peak_value = traits.Float(
        argstr="-fmls_peak_value %d", desc="any lobe with a maximal peak amplitude smaller than this threshold will be discarded"
    )
    
    fmls_integral = traits.Float(
        argstr="-fmls_integral %f", desc="any lobe with an integral smaller than this threshold will be discarded"
    )
    afd_file = File(
        usedefault=True, argstr="-afd %s", position=-3, 
        desc="output the total Apparent Fibre Density per fixel (integral of FOD lobe)"
    )
    peak_file = File(
        usedefault=True, argstr="-peak_amp %s", position=-2, 
        desc="output the amplitude of the FOD at the maximal peak per fixel"
    )
    disp_file = File(
        usedefault=True, argstr="-disp %s", position=-1, 
        desc="output a measure of dispersion per fixel as the ratio between FOD lobe integral and maximal peak amplitude"
    )

class fod2fixelOutputSpec(TraitedSpec):
    out_file = File(argstr="%s", desc="output image")
    afd_file = File(argstr="-afd %s", desc="output AFD file")
    peak_file = File(argstr="-peak_amp %s", desc="output peak amplitude file")
    disp_file = File(argstr="-disp %s", desc="output fixel dispersion file")

class fod2fixel(MRTrix3Base):
    """
    Perform segmentation of continuous Fibre Orientation Distributions 
    (FODs) to produce discrete fixels
    Example
    -------
    >>> fod2Fixel.inputs.in_file = 'wmfod.mif' 
    >>> fod2Fixel.inputs.fmls_peak_value = 0
    >>> fod2Fixel.inputs.fmls_integral = 0.1
    >>> fod2fixel.cmdline                            
    'fod2Fixel wmfod.mif -fmls_peak_value 0 -fmls_integral 0.1 -afd afd,mif -peak_amp peak.mif -disp disp.mif
    >>> fod2Fixel.run()                               
    """
    _cmd = "fod2fixel"
    input_spec = fod2fixelInputSpec
    output_spec = fod2fixelOutputSpec

    def _list_outputs(self):
        outputs = self.output_spec().get()
        outputs["out_file"] = op.abspath(self.inputs.out_file)
        if self.inputs.afd_file != Undefined:
            outputs["afd_file "] = op.abspath(self.inputs.afd_file )
        if self.inputs.peak_file != Undefined:
            outputs["peak_file"] = op.abspath(self.inputs.peak_file)
        if self.inputs.disp_file != Undefined:
            outputs["disp_file"] = op.abspath(self.inputs.disp_file)
        return outputs   

    