# SyntheticMicroBundle
This repository contains the code for generating synthetic data of beating microbundles based on experimentally-informed Finite Element simulations, as implemented in the paper [`TITLE`](addlink).

![Pipeline](/Figures/Pipeline_SyntheticData_RealFrame.png)

## In this Repository
This repository contains the following folders:
* [`Warp_code`](Warp_code) -- This folder contains the files needed to create synthetic movie frames of beating microbundles based on FEA data.
* [`FEA_code`](FEA_code) -- This folder contains the FEA code and mesh files used to obtain the FEA data.  
* [`Example`](Example) -- This folder contains a brief tutorial on how to use the provided code files, organize the input files, and understand the output files.
* [`Data`](Data) -- This folder contains a description of the generated synthetic data as well as a link to the full dataset. 

## Supporting Repositories 
This repository is recommended to be used in conjunction with the  [`MicroBundleCompute`](https://github.com/elejeune11/MicroBundleCompute) repository. The latter adds functionalities to generate segmentation masks of the microbundle region of the frames and track this region or a portion of it across the movie frames.