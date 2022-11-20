# SyntheticMicroBundle: Warp Code
This folder contains the files to create synthetic movie frames of beating microbundles based on FEA results. All files were tested in `Python 3.9.13`.
## In this Folder
The  `Warp_code` folder contains the following files:
* [`GenerateSyntheticTextures.py`](GenerateSyntheticTextures.py) -- This code generates synthetic regions of beating microbundles based on displacement results obtained via Finite Element simulations given a real frame image of a microbundle, a mask of the tissue region of the frame, and FEA displacement results at different points within the tissue domain. 
* [`Add_Noise.py`](Add_Noise.py) -- This code adds a Perlin noise to the generated synthetic microbundle frames given a mask of the texture region (region to add noise), number of octaves, and magnitude ratio relative to the maximum intensity value within the mask region.
* [`fcns.py`](fcns.py) --  This file contains all the user-defined functions needed to run the codes.