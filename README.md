# MRCdwitools

Tools for Diffusion Weighted Imaging (DWI):

## dwifit
- Voxelwise fitting of various DWI models
- Input format is Nifti
- Implementation is with dlib C++ library http://dlib.net
- Uses Nifti file format https://radiopaedia.org/articles/nifti-file-format, 
  while ASCII file format is internally available
- B-value ASCII file has b-values in intergers, separated by space, e. g.:<br>
  0 100 200 300 400 500 1000 2000 3000

### Usage of programs:

#### 1 Convert Nifti to ASCII for fittings
- dwiNifti2ASCII.sh DWI Nifti filename Mask Nifti filename

#### 2 Run one or more Nifti files in series
- run_dwifit_for_all_ASCII.sh model name

#### 3 Construct Nifti from fitted parameter values in ASCII format
- dwiASCII2Nifti.sh DWI Nifti filename Mask Nifti filename B-value ASCII filename


#### The binary for execution can be used as standalone
- dwifit
