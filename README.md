# SyntheticMicroBundle
This repository contains the code for generating synthetic data of beating microbundles based on experimentally-informed Finite Element simulations, as implemented in the paper [`TITLE`](addlink).

![Pipeline](Figures/Pipeline_SyntheticData_RealFrame.png)

## In this Repository
This repository contains the following files:
* [`GenerateSyntheticTextures.py`](GenerateSyntheticTextures.py) -- This code generates synthetic regions of beating microbundles based on displacement results obtained via Finite Element Simulations given a real frame image of a microbundle, a mask of the tissue region of the frame, and FEA displacement results at different points within the tissue domain. 
* [`Add_Noise.py`](Add_Noise.py) -- This code adds a Perlin noise to the generated synthetic microbundle frames given a mask of the texture region (region to add noise), number of octaves, and magnitude ratio relative to the maximum intensity value within the mask region.
* [`fcns.py`](fcns.py) --  This file contains all the user-defined functions needed to run the codes.
## Tutorial
We provide further details about the code structure to ensure that users can make use of it smoothly. We also provide an example folder `Example` showing the required inputs and expected outputs.
### Preparing data for analysis (Input)
The code expects three input folders. The naming of the folders is not important as long as the names are consistent with the folder paths indicated inside the codes. Specifically, we use the following folder structure:

1. `Movies` folder: This folder contains real movies of beating microbundles. Critically, movies should have a `.tif` extension.

2. `Masks` folder: This folder contains binary masks corresponding to each of the real movies. Also, masks should have a `.tif` extension. Synthetic masks outputted by `GenerateSyntheticTextures.py` are also saved to this folder as `Synthetic_Mask_*.tif` (e.g. 'Synthetic_Mask_D1T4_Before.tif')

3. `FEA_Results_*` folder (e.g. 'FEA_Results_Homog_MaxAct0.1_VFA'): This folder contains displacement results extracted at cell centers and saved as a text file for each step of the timeseries Finite Element simulation.

Initially, the `Example` folder should have the following structure:

```bash
|___ Example
|        |___ Movies
|                |___"*.tif"
|        |___ Masks
|                |___"*_Mask.tif"
|        |___ FEA_Results_*
|                |___"disp_all_Step%i.txt"
```

### Running the code (Output)
Once the input data structuring described above is followed, generating synthetic data should be straightforward. Specifically, the following folders and files are outputted:

The output of the code as shown in `Example` folder should be as follows:

```bash
|___ Example
|        |___ Movies
|                |___"*.tif"
|        |___ Masks
|                |___"*_Mask.tif"
|                |___"Synthetic_Mask_*.tif"
|        |___ FEA_Results_$\color{green}{*}$
|                |___"disp_all_Step%i.txt"
|        |___ Frames
|                |___Frames_*
|                        |___"Frame%04.tif"
|        |___ Textures
|                |___"*_Frame%04_TissueTexture.png"
|        |___ *_Frame%04_SyntheticTextures
|                |___ "Warped_Frame%04.tif"
|        |___ Results_*_Frame%04_SyntheticData_SquareMask
|                |___ "Results_G%ix%i_ErrorF_SquareMask.txt"
|        |___ Noisy_*_Frame%04_SyntheticTextures
|                |___ G%ix%i_MagR%f_Oct%f
|                        |___ "Noisy_Frame%04.tif"
|                        |___ "Noisy_Synthetic_Perlin_MagR%f_Oct%f.gif"
```

### Understanding the output files

## Supporting Repositories 
Include link to [`MicroBundleCompute`](https://github.com/elejeune11/MicroBundleCompute)