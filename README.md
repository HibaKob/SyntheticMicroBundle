# SyntheticMicroBundle
This repository contains the code for generating synthetic data of beating microbundles based on experimentally-informed Finite Element Simulations, as implemented in the paper 'TITLE' [Link](addlink).

![Pipeline](Pipeline_SyntheticData_RealFrame.png)

## In this Repository
This repository contains the following files:
* [`GenerateSyntheticTextures.py`](GenerateSyntheticTextures.py) -- This code generates synthetic regions of beating microbundles based on displacement results obtained via Finite Element Simulations given a real frame image of a microbundle, a mask of the tissue region of the frame, and FEA displacement results at different points within the tissue domain. 
* [`Add_Noise.py`](Add_Noise.py) -- This code adds a Perlin noise to the generated synthetic microbundle frames given a mask of the texture region (rgion to add noise), number of octaves, and magnitude ratio relative to the maximum intensity value within the mask region.
* [`Track_Frames.py`](Track_Frames.py) -- This code tracks the displacement of a region in real and synthetic movies of beating microbundles across frames given a movie frame and a mask. 
* [`fcns.py`](fcns.py) --  This file contains all the user-defined functions needed to run the codes.
## Tutorial
We provide further details about the code structure to ensure that users can make use of it smoothly. (Include an example folder?)
### Preparing data for analysis (Input)
The code expects three input folders. The naming of the folders is not important as long as the names are consistent with the folder paths indicated inside the codes. Specifically, we use the following folder structure:

1. `Movies` folder: This folder contains real movies of beating microbundles. Critically, movies should have a `.tif` extension.

2. `Masks` folder: This folder contains binary masks corresponding to each of the real movies. Also, masks should have a `.tif` extension.

3. `FEA_Results_xxxxxx` folder:

### Running the code (Output)

### Understanding the output files

## Supporting Repositories 
Include link to [`MicroBundleCompute`](https://github.com/elejeune11/MicroBundleCompute)