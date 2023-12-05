# EVPSC_DV
A modified version of [EVPSC_CPP](https://github.com/Kecheng96/EVPSC_CPP). In this version, deformation modes containing slip systems and twinning systems are dealed as equitable modes, each of them can be interactive with the other one, instead of dealing them as deformation mode families. Equipped with an Orowan version of slip movement law, and a draft copy of the twinning-de-twinning scheme. A single crystal version is at [SingleCrystalCpp](https://github.com/ShawnWolgu/SingleCrystalCpp) and the detailed description of the deformation laws can be obtained in that repo.

This program is still ongoing so it is still not easy to use, and the code may be renewed sometimes.

## Compilation
To compile this code you can use cmake and refer to the CMakeLists.txt example. Be aware that these libraries are needed:
- Eigen3
- nlohmann-json
- OpenMP
make sure the correct link has been made in your environment.

## Examples
I prepared some examples. The meaning of the model parameters can be readed directly in the code or the docs in the repo [SingleCrystalCpp](https://github.com/ShawnWolgu/SingleCrystalCpp). To run these examples just simply run "EVPSC_CPP" the executible file in your working folder with the input file "EVPSC_CPP.in".

