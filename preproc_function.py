# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:
# -*- coding: utf-8 -*-
import os.path as op
from nipype.interfaces.base import CommandLine, traits, TraitedSpec, File, CommandLineInputSpec
from nipype.interfaces.mrtrix3.base import MRTrix3BaseInputSpec, MRTrix3Base

class DWIPreProcInputSpec(MRTrix3BaseInputSpec):
    in_file = File(
        exists=True, argstr="%s", mandatory=True, position=-2, desc="input dwi image"
    )
    out_file = File(argstr="%s", mandatory=True, position=-1, desc="output image")
    use_header = traits.Bool(
        argstr="-rpe_header",
        mandatory=True,
        desc="Phase-encoding information can be found in the image header",
    )
    se_epi_file = File(exists=True, argstr="-se_epi %s", desc="Spin Echo EPI image")
    
class DWIPreProcOutputSpec(TraitedSpec):
    out_file = File(exists=True, desc="output image")

class DWIFslPreProc(MRTrix3Base):
    """
    Perform diffusion image pre-processing using FSL's eddy tool; 
    including inhomogeneity distortion correction 
    using FSL's topup tool if possible
    Example
    -------
    >>> import nipype.interfaces.mrtrix3 as mrt
    >>> dwifslpreproc = mrt.dwifslpreproc()
    >>> dwifslpreproc.inputs.in_file = 'dwi.mif' 
    >>> dwifslpreproc.inputs.se_epi_file = 'b0s.mif'
    >>> dwifslpreproc.inputs.out_file = 'preprocessed.mif'
    >>> dwifslpreproc.cmdline                            
    'dwifslpreproc dwi.mif preprocessed.mif -rpe_header _se_epi b0s.mif'
    >>> dwifslpreproc.run()                               
    """
    _cmd = "dwifslpreproc"
    input_spec = DWIPreProcInputSpec
    output_spec = DWIPreProcOutputSpec

    def _list_outputs(self):
        outputs = self.output_spec().get()
        outputs["out_file"] = op.abspath(self.inputs.out_file)
        return outputs