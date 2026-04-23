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

## Citation

If you use this code or find the methodology helpful for your research, please cite our corresponding paper published in **Journal of the Mechanics and Physics of Solids (JMPS)**:

### Standard Reference
> Kaiwei Wu, Xiaochuan Sun, Bowen Liu, et al. **A multiscale investigation into the electroplastic effects in copper: Experiments and crystal plasticity modeling.** *Journal of the Mechanics and Physics of Solids*, Vol. 212, 106597 (2026).  
> DOI: [10.1016/j.jmps.2026.106597](https://doi.org/10.1016/j.jmps.2026.106597)

### BibTeX
```bibtex
@article{WU2026106597,
  title = {A multiscale investigation into the electroplastic effects in copper: Experiments and crystal plasticity modeling},
  journal = {Journal of the Mechanics and Physics of Solids},
  volume = {212},
  pages = {106597},
  year = {2026},
  issn = {0022-5096},
  doi = {[https://doi.org/10.1016/j.jmps.2026.106597](https://doi.org/10.1016/j.jmps.2026.106597)},
  author = {Kaiwei Wu and Xiaochuan Sun and Bowen Liu and Quan Li and Wei Wen and Qiwei Shi and Zhutian Xu and Linfa Peng and Huamiao Wang and Yin Zhang and Samuel Forest}
}
