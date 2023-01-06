# SyntheticMicroBundle: FEA Code
This folder contains all the files needed to run Finite Element simulations and extract displacement data used for the generation of synthetic beating microundle data. The results are similar to those in `FEA_Results_Homog_MaxAct0.08_VFA_Z` folder shared in this repository under `Example` folder.

The following files are provided: 

* [`FEA_Synthetic_Microbundle.py`](FEA_Synthetic_Microbundle.py) -- The Finite Element code file run in [`FEniCS 2019.1.0`](https://fenicsproject.org). More details about the implementation of the Finite Element model can be found in paper [`TITLE`](addlink).

* [`MicroTug_RefDim_3D.xdmf`](MicroTug_RefDim_3D.xdmf) & [`MicroTug_RefDim_3D.h5`](`MicroTug_RefDim_3D.h5) -- The microbundle mesh generated in [`Gmsh 4.10.5`](https://gmsh.info) and saved in `.xdmf` format readable by FEniCS.
