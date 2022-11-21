# SyntheticMicroBundle: Example
## Tutorial
In this folder, we share further details about the code structure, the required inputs and expected outputs to ensure that users can make use of it smoothly.
### Preparing data for analysis (Input)
The code expects three input folders. The naming of the folders is not important as long as the names are consistent with the folder paths indicated inside the codes. Specifically, we use the following folder structure:

1. `Movies` folder: This folder contains real movies of beating microbundles. Critically, movies should have a `.tif` extension.

2. `Masks` folder: This folder contains binary masks corresponding to each of the real movies. Also, masks should have a `.tif` extension. Synthetic masks outputted by `GenerateSyntheticTextures.py` are also saved to this folder as `Synthetic_Mask_*.tif` (e.g. 'Synthetic_Mask_D1T4_Before.tif')

3. `FEA_Results_*` folder (e.g. 'FEA_Results_Homog_MaxAct0.1_VFA'): This folder contains displacement results extracted at cell centers and saved as text files for each step of the timeseries Finite Element simulation. These files are automatically generated when `FEA_Synthetic_Microbundle.py` file in `FEA_code` folder is run.

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

\*this text is surrounded by literal asterisks\*
```bash
|___ Example
|        |___ Movies
|                |___"*.tif"
|        |___ Masks
|                |___"*_Mask.tif"
|                |___"Synthetic_Mask_*.tif"
|        |___ FEA_Results_${**}
|                |___"disp_all_Step%i.txt"
|        |___ Frames
|                |___Frames_*   
|                        |___"Frame%04.tif"
|        |___ Textures
|                |___"*_Frame%04_TissueTexture.png"
|        |___ *_Frame%04_SyntheticTextures
|                |___ **_G%ix%i
                         |___ "Warped_Frame%04.tif"
                         |___ "Synthetic_Frames_G%ix%i.gif"
|        |___ Noisy_*_Frame%04_SyntheticTextures
|                |___ **_G%ix%i_MagR%f_Oct%f
|                        |___ "Noisy_Frame%04.tif"
|                        |___ "Noisy_Synthetic_Perlin_MagR%f_Oct%f.gif"
```

### Understanding the output files