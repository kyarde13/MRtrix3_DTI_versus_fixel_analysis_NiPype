MRtrix3 was utilised to process the DWI dataset from a single subject. Tensor-based WM fibre representation utilised the main eigenvector to indicate the fibre’s direction, and its orientation with respect to the main magnetic field (z direction).
 
For fixel-based analysis, spherical deconvolution identified multiple fibre bundles within each imaging voxel. Fixels were converted into primary, secondary and tertiary vectors, with x, y and z components. From the z component of the primary fixel, the WM fibre’s orientation with respect to the main magnetic field was obtained. 
 
A Python pipeline was created to perform all MRtrix3 commands in sequence. Some of the open-sourced code available in Nipype was utilised, so that value permutations of various function parameters could be ran.


# 1) install dependencies
setup a miniconda python environment e.g.
```
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
```

logout or source .bashrc
```
source ~/.bashrc
```

install nipype
```
conda install --channel conda-forge nipype
```

# 2) run
- use nipype to run MRtrix pipeline by running run_2_nipype_dwi.py (adjust subject names inside the script)
