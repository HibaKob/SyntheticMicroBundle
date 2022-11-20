# SyntheticMicroBundle: Example
## Tutorial
In this folder, we share further details about the code structure, the required inputs and expected outputs to ensure that users can make use of it smoothly.
### Preparing data for analysis (Input)
The code expects three input folders. The naming of the folders is not important as long as the names are consistent with the folder paths indicated inside the codes. Specifically, we use the following folder structure:

1. `Movies` folder: This folder contains real movies of beating microbundles. Critically, movies should have a `.tif` extension.

2. `Masks` folder: This folder contains binary masks corresponding to each of the real movies. Also, masks should have a `.tif` extension. Synthetic masks outputted by `GenerateSyntheticTextures.py` are also saved to this folder as `Synthetic_Mask_*.tif` (e.g. 'Synthetic_Mask_D1T4_Before.tif')

3. `FEA_Results_*` folder (e.g. 'FEA_Results_Homog_MaxAct0.1_VFA'): This folder contains displacement results extracted at cell centers and saved as text files for each step of the timeseries Finite Element simulation.

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

The output of the code as shown in this folder should be as follows:

\*this text is surrounded by literal asterisks\*
```bash
|___ Example
|        |___ Movies
|                |___"*.tif"
|        |___ Masks
|                |___"*_Mask.tif"
|                |___"Synthetic_Mask_*.tif"
|        |___ FEA_Results_*
|                |___"disp_all_Step%i.txt"
|        |___ Frames
|                |___Frame_*   
|                        |___"Frame%04.tif"
|        |___ Textures
|                |___"*_Frame%04_TissueTexture.png"
|        |___ *_Frame%04_SyntheticTextures
|                |___ "Warped_Frame%04.tif"
|        |___ Noisy_*_Frame%04_SyntheticTextures
|                |___ G%ix%i_MagR%f_Oct%f
|                        |___ "Noisy_Frame%04.tif"
|                        |___ "Noisy_Synthetic_Perlin_MagR%f_Oct%f.gif"
```

### Understanding the output files