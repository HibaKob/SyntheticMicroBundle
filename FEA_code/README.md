# SyntheticMicroBundle: FEA Code
This folder contains all the files needed to run Finite Element simulations and extract displacement data used for the generation of synthetic beating microundle data. The results are similar to those in `FEA_Results_Homog_MaxAct0.08_VFA_Z_Lite` folder shared in this repository under `Example` folder.

The following files are provided: 

* [`FEA_Synthetic_Microbundle.py`](FEA_Synthetic_Microbundle.py) -- The Finite Element code file run in [`FEniCS 2019.1.0`](https://fenicsproject.org). More details about the implementation of the Finite Element model can be found in paper [`MicroBundleCompute: Automated segmentation, tracking, and analysis of subdomain deformation in cardiac microbundles`](https://doi.org/10.1371/journal.pone.0298863). We note here that we run this code in parallel. As an example, the following command with the FEniCS module loaded can be used:
    ```
    mpirun -np 2 python3 FEA_Synthetic_Microbundle.py
    ```

* [`MicroTug_Geom_3D.xdmf`](MicroTug_Geom_3D.xdmf) & [`MicroTug_Geom_3D.h5`](`MicroTug_Geom_3D.h5) -- The microbundle mesh generated in [`Gmsh 4.10.5`](https://gmsh.info) and saved in `.xdmf` format readable by FEniCS.
