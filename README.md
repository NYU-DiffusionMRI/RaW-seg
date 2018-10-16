# Random Walker (RaW) segmentation

The code implements random walker (RaW) segmentations of myelinated axons in a 3d electron microscopy (EM) data, originally developed in [Lee, et al., 2018](https://doi.org/). All functions are implemented in MATLAB and embeded into 3 classes:

* **rawseg:** Random walker segmentation for myelinated axons in a 3d EM data.
* **analyzeseg:** Analyze segmented axons and calculate fiber orientation function (FOD), dispersion angle, axonal diameter, and g-ratio.
* **plotseg** Visualize segmented axons into a 3d box or 3d polyhedrons.

The usage of the classes can be found in demo_\*.m files.

## References
* **Random walker (RaW) segmentation**
  - [Lee, et al., 2018](https://doi.org/)

* **Distortion correction**
  - [Sun, et al., IEEE CVPR, 2010](https://doi.org/10.1109/CVPR.2010.5539939)
  - [Sun, et al., IJCV, 2014](https://doi.org/10.1007/s11263-013-0644-x)
  
* **FOD on a triangulated spherical surface**
  - [Efficient spherical design](http://web.maths.unsw.edu.au/~rsw/Sphere/EffSphDes)
  - [Womersley, et al., arXiv, 2017](https://arxiv.org/abs/1709.01624)

* **[Sherical harmonic transform library](https://www.mathworks.com/matlabcentral/fileexchange/43856-real-complex-spherical-harmonic-transform-gaunt-coefficients-and-rotations)**
  - [Politis, et al., 2016](https://aaltodoc.aalto.fi/handle/123456789/22499)
  
* **[Bingham statistics library](https://github.com/SebastianRiedel/bingham)**
  
* **[Magicwand2](https://www.mathworks.com/matlabcentral/fileexchange/6034-magicwand2)**

## Authors
* [Hong-Hsi Lee](http://www.diffusion-mri.com/people/hong-hsi-lee)
* [Dmitry S Novikov](http://www.diffusion-mri.com/people/dmitry-novikov)
* [Els Fieremans](http://www.diffusion-mri.com/people/els-fieremans)

## License
This project is licensed under the [LICENSE](https://github.com/NYU-DiffusionMRI/RaW-seg/blob/master/LICENSE).
