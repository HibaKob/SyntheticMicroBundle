# SyntheticMicroBundle: Example
## Tutorial
In this folder, we share further details about the code structure, the required inputs and expected outputs to ensure that users can make use of it smoothly.
### Preparing data for analysis (Input)
The code expects three input folders. The naming of the folders is not important as long as the names are consistent with the folder paths indicated inside the codes. Specifically, we use the following folder structure:

1. `Movies` folder: This folder contains real movies of beating microbundles. Critically, movies should have a `.tif` extension.

2. `Masks` folder: This folder contains binary masks corresponding to each of the real movies. Also, masks should have a `.tif` extension. Synthetic masks outputted by `GenerateSyntheticTextures.py` are also saved to this folder as `Synthetic_Mask_*.tif` (e.g. "Synthetic_Mask_D1T4_Before.tif")

3. `FEA_Results_*` folder (e.g. "FEA_Results_Homog_MaxAct0.08_VFA_Z_Lite"): This folder contains X, Y, and Z displacement results extracted at mesh cell centers for each step of the timeseries Finite Element (FE) simulation and saved as text files. These files are automatically generated when `FEA_Synthetic_Microbundle.py` file in `FEA_code` folder is run.

    3a. Each `pos_*.txt` file is a 'mxn' array of displacement results at the 'm' mesh cell centers for 'n' steps in the X, Y, and Z directions respectively. For the 'lite' version shared in this folder, we include 25 simulation steps only due to restrictions on file sizes. 

4. `Tissue_Slice_Coordinates.txt` file: This text file contains the contour coordinates (X,Y) of the microbundle used to generate the mesh for the Finite Element simulations.

In addition to the `run_code_*` python files, the `Example` folder should have the following structure initially:

```bash
|___ Example
|        |___ Movies
|                |___"*.tif"
|        |___ Masks
|                |___"*_Mask.tif"
|        |___ FEA_Results_**
|                |___"disp_all_Step%i.txt"
|        |___ "Tissue_Slice_Coordinates.txt"
```

### Running the code (Output)
Once the input data structuring described above is followed, generating synthetic data should be straightforward. Specifically, the following folders and files are outputted according to the structure detailed below:

```bash
|___ Example
|        |___ Movies
|                |___"*.tif"
|        |___ Masks
|                |___"*_Mask.tif"
|                |___"Synthetic_Mask_*.tif"
|        |___ FEA_Results_**
|                |___"disp_all_Step%i.txt"
|        |___ "Tissue_Slice_Coordinates.txt"
|        |___ Frames
|                |___Frames_*   
|                        |___"Frame%04.tif"
|        |___ Textures
|                |___"*_Frame%04_TissueTexture.png"
|        |___ *_Frame%04_SyntheticTextures
|                |___ **_G%ix%i
|                         |___ "Warped_Frame%04.tif"
|        |___ Noisy_*_Frame%04_SyntheticTextures
|                |___ **_G%ix%i_MagR%f_Oct%f
|                        |___ "Noisy_Frame%04.tif"
```

### Understanding the output files
The output of the code mainly includes the synthetic frames of beating microbundles with and without added Perlin noise. The remaining output consists of crucial intermediate steps necessary to obtain the final product. Mainly, outputs are saved in `.tif` or `.png` formats for easier visualization.

* `Synthetic_Mask_*.tif`: A binary mask of the region of interest for generating synthetic data. The size of this mask is equal to that of the original movie frames. Here, the size of the `.tif` image is '512x512'. 

* `Frame%04.tif`: The frames of the movie file saved individually. We note here that we preserve the size of the frames as well as data type ('uint16' in this case). 

* `*_Frame%04_TissueTexture.png`: The microbundle texture isolated from "Frame%04" based on "*_Mask.tif" and saved with transparent background. 

* `Warped_Frame%04.tif`: The region of interest isolated from the microbundle texture image using "Synthetic_Mask_*.tif" and deformed based on displacement results extracted from FE simulations. Each frame is a '512x512' 'uint16' image. 

* `Noisy_Frame%04.tif`: The synthetic frames with Perlin noise added to the textured region. Each frame is a '512x512' 'uint8' image. 