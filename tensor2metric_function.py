# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:
# -*- coding: utf-8 -*-
import os.path as op
from nipype.interfaces.base import CommandLine, traits, TraitedSpec, File, CommandLineInputSpec
from nipype.interfaces.mrtrix3.base import MRTrix3BaseInputSpec, MRTrix3Base

class tensor2metricInputSpec(MRTrix3BaseInputSpec):
    
    input_file = File(
        exists=True, argstr="%s", mandatory=True, position=-4, desc="input DTI"
    )
    
    vector_file = File(
        argstr="-vector %s", mandatory=False, position=-1, 
        desc="compute the selected eigenvector(s) of the diffusion tensor.")
    
    fa_file = File(
        argstr="-fa %s", mandatory = False, position=-1, 
        desc="compute the fractional anisotropy (FA) of the diffusion tensor."
    )
    adc_file = File(
        argstr="-adc %s", mandatory = False, position=-1, 
        desc="compute the mean apparent diffusion coefficient (ADC) of the diffusion tensor"
    )

    modulate = traits.String(
        argstr="-modulate %s", mandatory = False, position=-3,
        desc="specify how to modulate (deault = FA)"
    )
    num = traits.Int(
        argstr="-num %d", mandatory = False, position=-2, desc="specify the desired eigenvalue/eigenvector(s)"
    )
    

class tensor2metricOutputSpec(TraitedSpec):
    vector_file = File(exists=True, desc="output image")

class tensor2metric(MRTrix3Base):
    """
    Generate maps of tensor-derived parameters
    Example
    -------
    >>> tensor2metric.inputs.in_file = 'dti.mif' 
    >>> tensor2metric.cmdline                            
    tensor2metric dt.mif -fa fa.mif #computes fa map
    tensor2metric dt.mif -adc adc.mif #computes adc map
    tensor2metric dt.mif -modulate none -num 1 -vector eigenvector.mif 
    >>> tensor2metric.run()                               
    """
    _cmd = "tensor2metric"
    input_spec = tensor2metricInputSpec
    output_spec = tensor2metricOutputSpec

    def _list_outputs(self):
        outputs = self.output_spec().get()
        outputs["vector_file"] = op.abspath(self.inputs.vector_file)
        return outputs 