# SyntheticMicroBundle: Data

The synthetic microbundle dataset consists of 60 `200x256x256` `.tif` files generated based on textures extracted from real data and warped according to experimentally-informed Finite Element (FE) simulations. We make the entire dataset available [`here`](https://drive.google.com/drive/folders/1nlaBlfXcrup4wdoifaoX3D6US0cQ1Zi8?usp=sharing).

<p align = "center">
<img alt="SyntheticExample" src="/Figures/Synthetic_Frames_G16x16.gif" width="50%" />
</p>

We obtained 15 different microbundle textures from 5 real microbundle movies at 3 different frames: the first frame, a peak frame, and a valley frame. From these textures, we cropped `90x90` square domains centerd at each microbundle domain and then warped each of these 15 square textures based on displacement results of 4 FE simulations: 

1. [`FEA_Results_Homog_MaxAct0.08_VFA_X`](https://drive.google.com/drive/folders/18WdScdh4GIq9YaA5ygX6wx20lC1ZgEqa?usp=sharing) -- Homogeneous activation across the whole microbundle domain with the fiber direction varying linearly in X direction (microbundle length)

2. [`FEA_Results_Homog_MaxAct0.08_VFA_Z`](https://drive.google.com/drive/folders/1FobFgv8s0IjGSRAax_m2sTKDNFiD8kCI?usp=sharing) -- Homogeneous activation across the whole microbundle domain with the fiber direction varying linearly in Z direction (microbundle depth) 

3. [`FEA_Results_Heterog_MaxAct0.08_VFA_X`](https://drive.google.com/drive/folders/1004hc0bogykYXMqm70VJd2l86QVufLD4?usp=sharing) -- Heterogeneous activation where the active microbundle domain has a passive inlcusion in the middle with the fiber direction varying linearly in X direction (microbundle length)

4. [`FEA_Results_Heterog_MaxAct0.08_VFA_Z`](https://drive.google.com/drive/folders/1OhL6X4cpTSrMWCU8bARD9L-BdnRkvcqc?usp=sharing) -- Heterogeneous activation where the active microbundle domain has a passive inlcusion in the middle with the fiber direction varying linearly in Z direction (microbundle depth) 

Each of these 4 folders contains 9 `.txt` result files, 3 for the mesh cell center positions in the X, Y, and Z directions, and 6 for the Green-Lagrange strain results $E_{xx}$, $E_{xy}$, $E_{xz}$, $E_{yy}$, $E_{yz}$, $E_{zz}$. Each of these text files contains a `number of points x number of steps` array. In other words, each column represents the results at the corresponding simulation step. We briefly note here that all FE simulations were run for 200 steps, simulating 8 complete microbundle beats. 

The dataset can be accessed [`here`](https://drive.google.com/drive/folders/1bomqRqcy550tiXtsZ4Om9pYwl1iOco7e?usp=sharing). Each `.tif` file is approximatley 105 MB in size. The table below summarizes the conditions for generating the 60 different synthetic microbundles.

| File Name                | Texture Source | Frame     | FEA data                       |
| -------------------------| -------------- | --------- | -----------------------------  |
| `SyntheticTextures_1.tif`  | Movie 1 | Frame0000 | Homogeneous Activation, Variable Fiber Direction in X |
| `SyntheticTextures_2.tif`  | Movie 1 | Frame0000 | Homogeneous Activation, Variable Fiber Direction in Z |
| `SyntheticTextures_3.tif`  | Movie 1 | Frame0000 | Heterogeneous Activation, Variable Fiber Direction in X |
| `SyntheticTextures_4.tif`  | Movie 1 | Frame0000 | Heterogeneous Activation, Variable Fiber Direction in Z |
| `SyntheticTextures_5.tif`  | Movie 1 | Frame0041 | Homogeneous Activation, Variable Fiber Direction in X |
| `SyntheticTextures_6.tif`  | Movie 1 | Frame0041 | Homogeneous Activation, Variable Fiber Direction in Z |
| `SyntheticTextures_7.tif`  | Movie 1 | Frame0041 | Heterogeneous Activation, Variable Fiber Direction in X |
| `SyntheticTextures_8.tif`  | Movie 1 | Frame0041 | Heterogeneous Activation, Variable Fiber Direction in Z |
| `SyntheticTextures_9.tif`  | Movie 1 | Frame0077 | Homogeneous Activation, Variable Fiber Direction in X |
| `SyntheticTextures_10.tif`  | Movie 1 | Frame0077 | Homogeneous Activation, Variable Fiber Direction in Z |
| `SyntheticTextures_11.tif`  | Movie 1 | Frame0077 | Heterogeneous Activation, Variable Fiber Direction in X |
| `SyntheticTextures_12.tif`  | Movie 1 | Frame0077 | Heterogeneous Activation, Variable Fiber Direction in Z |
| `SyntheticTextures_13.tif`  | Movie 2 | Frame0000 | Homogeneous Activation, Variable Fiber Direction in X |
| `SyntheticTextures_14.tif`  | Movie 2 | Frame0000 | Homogeneous Activation, Variable Fiber Direction in Z |
| `SyntheticTextures_15.tif`  | Movie 2 | Frame0000 | Heterogeneous Activation, Variable Fiber Direction in X |
| `SyntheticTextures_16.tif`  | Movie 2 | Frame0000 | Heterogeneous Activation, Variable Fiber Direction in Z |
| `SyntheticTextures_17.tif`  | Movie 2 | Frame0029 | Homogeneous Activation, Variable Fiber Direction in X |
| `SyntheticTextures_18.tif`  | Movie 2 | Frame0029 | Homogeneous Activation, Variable Fiber Direction in Z |
| `SyntheticTextures_19.tif`  | Movie 2 | Frame0029 | Heterogeneous Activation, Variable Fiber Direction in X |
| `SyntheticTextures_20.tif`  | Movie 2 | Frame0029 | Heterogeneous Activation, Variable Fiber Direction in Z |
| `SyntheticTextures_21.tif`  | Movie 2 | Frame0061 | Homogeneous Activation, Variable Fiber Direction in X |
| `SyntheticTextures_22.tif`  | Movie 2 | Frame0061 | Homogeneous Activation, Variable Fiber Direction in Z |
| `SyntheticTextures_23.tif`  | Movie 2 | Frame0061 | Heterogeneous Activation, Variable Fiber Direction in X |
| `SyntheticTextures_24.tif`  | Movie 2 | Frame0061 | Heterogeneous Activation, Variable Fiber Direction in Z |
| `SyntheticTextures_25.tif`  | Movie 3 | Frame0000 | Homogeneous Activation, Variable Fiber Direction in X |
| `SyntheticTextures_26.tif`  | Movie 3 | Frame0000 | Homogeneous Activation, Variable Fiber Direction in Z |
| `SyntheticTextures_27.tif`  | Movie 3 | Frame0000 | Heterogeneous Activation, Variable Fiber Direction in X |
| `SyntheticTextures_28.tif`  | Movie 3 | Frame0000 | Heterogeneous Activation, Variable Fiber Direction in Z |
| `SyntheticTextures_29.tif`  | Movie 3 | Frame0043 | Homogeneous Activation, Variable Fiber Direction in X |
| `SyntheticTextures_30.tif`  | Movie 3 | Frame0043 | Homogeneous Activation, Variable Fiber Direction in Z |
| `SyntheticTextures_31.tif`  | Movie 3 | Frame0043 | Heterogeneous Activation, Variable Fiber Direction in X |
| `SyntheticTextures_32.tif`  | Movie 3 | Frame0043 | Heterogeneous Activation, Variable Fiber Direction in Z |
| `SyntheticTextures_33.tif`  | Movie 3 | Frame0080 | Homogeneous Activation, Variable Fiber Direction in X |
| `SyntheticTextures_34.tif`  | Movie 3 | Frame0080 | Homogeneous Activation, Variable Fiber Direction in Z |
| `SyntheticTextures_35.tif`  | Movie 3 | Frame0080 | Heterogeneous Activation, Variable Fiber Direction in X |
| `SyntheticTextures_36.tif`  | Movie 3 | Frame0080 | Heterogeneous Activation, Variable Fiber Direction in Z |
| `SyntheticTextures_37.tif`  | Movie 4 | Frame0000 | Homogeneous Activation, Variable Fiber Direction in X |
| `SyntheticTextures_38.tif`  | Movie 4 | Frame0000 | Homogeneous Activation, Variable Fiber Direction in Z |
| `SyntheticTextures_39.tif`  | Movie 4 | Frame0000 | Heterogeneous Activation, Variable Fiber Direction in X |
| `SyntheticTextures_40.tif`  | Movie 4 | Frame0000 | Heterogeneous Activation, Variable Fiber Direction in Z |
| `SyntheticTextures_41.tif`  | Movie 4 | Frame0033 | Homogeneous Activation, Variable Fiber Direction in X |
| `SyntheticTextures_42.tif`  | Movie 4 | Frame0033 | Homogeneous Activation, Variable Fiber Direction in Z |
| `SyntheticTextures_43.tif`  | Movie 4 | Frame0033 | Heterogeneous Activation, Variable Fiber Direction in X |
| `SyntheticTextures_44.tif`  | Movie 4 | Frame0033 | Heterogeneous Activation, Variable Fiber Direction in Z |
| `SyntheticTextures_45.tif`  | Movie 4 | Frame0064 | Homogeneous Activation, Variable Fiber Direction in X |
| `SyntheticTextures_46.tif`  | Movie 4 | Frame0064 | Homogeneous Activation, Variable Fiber Direction in Z |
| `SyntheticTextures_47.tif`  | Movie 4 | Frame0064 | Heterogeneous Activation, Variable Fiber Direction in X |
| `SyntheticTextures_48.tif`  | Movie 4 | Frame0064 | Heterogeneous Activation, Variable Fiber Direction in Z |
| `SyntheticTextures_49.tif`  | Movie 5 | Frame0000 | Homogeneous Activation, Variable Fiber Direction in X |
| `SyntheticTextures_50.tif`  | Movie 5 | Frame0000 | Homogeneous Activation, Variable Fiber Direction in Z |
| `SyntheticTextures_51.tif`  | Movie 5 | Frame0000 | Heterogeneous Activation, Variable Fiber Direction in X |
| `SyntheticTextures_52.tif`  | Movie 5 | Frame0000 | Heterogeneous Activation, Variable Fiber Direction in Z |
| `SyntheticTextures_53.tif`  | Movie 5 | Frame0030 | Homogeneous Activation, Variable Fiber Direction in X |
| `SyntheticTextures_54.tif`  | Movie 5 | Frame0030 | Homogeneous Activation, Variable Fiber Direction in Z |
| `SyntheticTextures_55.tif`  | Movie 5 | Frame0030 | Heterogeneous Activation, Variable Fiber Direction in X |
| `SyntheticTextures_56.tif`  | Movie 5 | Frame0030 | Heterogeneous Activation, Variable Fiber Direction in Z |
| `SyntheticTextures_57.tif`  | Movie 5 | Frame0062 | Homogeneous Activation, Variable Fiber Direction in X |
| `SyntheticTextures_58.tif`  | Movie 5 | Frame0062 | Homogeneous Activation, Variable Fiber Direction in Z |
| `SyntheticTextures_59.tif`  | Movie 5 | Frame0062 | Heterogeneous Activation, Variable Fiber Direction in X |
| `SyntheticTextures_60.tif`  | Movie 5 | Frame0062 | Heterogeneous Activation, Variable Fiber Direction in Z |
