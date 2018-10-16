# Random Walker (RaW) segmentation

The code implements random walker (RaW) segmentations of myelinated axons in a 3d EM data, originally developed in [Lee, et al., 2018](https://doi.org/). All functions are implemented in MATLAB and embeded into 3 classes:

* **rawseg:** Random walker segmentation for myelinated axons in a 3d EM data.
* **analyzeseg:** Analyze segmented axons and calculate fiber orientation function, dispersion angle, axonal diameter, and g-ratio.
* **plotseg** Visualize segmented axons into a 3d box or 3d polyhedrons.

The usage of the classes can be found in demo_\*.m files.

## References
* **Random walker (RaW) segmentation**
  - [Lee, et al., 2018](https://doi.org/)

* **Optical flow distortion correction**
  - [Sun, et al., 2010](https://doi.org/)
  - [Sun, et al., 2014](https://doi.org/)
  
* **Triangulated spherical surface**
  - [Womersley, et al., 2017] (https://doi.org/)

* **Sherical harmonics toolbox**
  - [Politis, et al., 2016](https://doi.org/)
  
* **Bingham statistics library**
  - https://github.com/SebastianRiedel/bingham
  
* **Magicwand2**
  - https://www.mathworks.com/matlabcentral/fileexchange/6034-magicwand2

## Authors
* [Hong-Hsi Lee](http://www.diffusion-mri.com/people/hong-hsi-lee)
* [Dmitry S Novikov](http://www.diffusion-mri.com/people/dmitry-novikov)
* [Els Fieremans](http://www.diffusion-mri.com/people/els-fieremans)

## License
This project is licensed under the [LICENSE](https://github.com/NYU-DiffusionMRI/monte-carlo-simulation-recipes/blob/example1/LICENSE).
