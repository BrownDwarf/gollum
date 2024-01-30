---
title: 'The `gollum` interface to precomputed synthetic spectral model grids'
tags:
  - Python
  - astronomy
  - spectroscopy
  - stars
  - echelle
authors:
  - name: Sujay Shankar
    orcid: 0000-0002-2290-6810
    affiliation: 1
  - name: Michael A. Gully-Santiago
    orcid: 0000-0002-4020-3457
    affiliation: 1
  - name: Caroline V. Morley
    orcid: 0000-0002-4404-0456
    affiliation: 1
  - name: Jiayi Cao
    orcid: 0000-0002-2466-3816
    affiliation: 1
affiliations:
 - name: The University of Texas at Austin Department of Astronomy, Austin, TX, USA
   index: 1
date: 30 January 2024
bibliography: paper.bib
---

# Summary

The spectra of stars, brown dwarfs, and planets are incredibly complex and information-rich. Complexity can be reduced by distilling spectra down to their fundamental stellar parameters&mdash;effective temperature, surface gravity, and iron abundance, and sometimes others such as alpha element abundance or sedimentation efficiency&mdash;which to a first order control the appearance of stellar and substellar spectra. Synthetic spectra models mimic stellar spectra across a grid of these fundamental parameters. Due to the computational impact of modeling stars, model grids have high resolution in the spectral axis (wavelength coordinates) but are coarsely sampled over fundamental stellar parameters. Comparing these grids to data has been a challenge due to the large number of model grids, their size, and their dimensionality, making the relation between fundamental stellar parameters and spectral appearance somewhat unintuitive at times. This necessitates the development of an intuitive, performant interface allowing astronomers to explore numerous model grids and compare them to data.

`gollum` is a Python package for intuitive analysis and visualization of precomputed synthetic spectra. Its API is designed to have modules dedicated to each model grid it supports, with each module then containing classes for both individual spectra and bulk grid access. The programmatic interface to spectral analysis uses method-chaining to make `gollum` code very readable. The visual interface in the form of interactive dashboards powered by `bokeh` take intuition to the next level with sliders and toggles that allow users to tweak both fundamental stellar parameters and  extrinsic parameters such as radial velocity and rotational broadening, viewing the results in real time [@bokeh2018]. `gollum`'s modularity allows for a wide range of model grids to potentially be supported, and its performance is optimized with libraries such as `numpy`, `scipy`, `astropy`, and `specutils` to allow for quick loading and processing of large amounts of data [@harris2020array; @scipy2020; @astropy13; @astropy18; @specutils].

`gollum` appeals to use cases ranging from entry-level astronomers to seasoned researchers with its combination of intuition and power. It interoperates with `muler`, a similar framework designed for true observed échelle spectra, and thanks to its dashboard, can also be used for by-eye fitting of models to data [@muler]. `gollum`'s programmatic interface is being used to create an extension to the `blase` framework that will allow for inference of fundamental stellar parameters from observed spectra using interpretable machine learning and interpolation techniques, taking inspiration from other frameworks such as `starfish` that also specialize in spectroscopic inference [@blase2022; @czekala15]. `gollum` currently supports PHOENIX, Sonora (Alpha, Bobcat, and Diamondback), and CoolTLUSTY model grids, and support for other model grids is planned [@husser2013; @alpha; @bobcat; @diamondback; @lacy2023].

# Acknowledgements

This research has made use of NASA's Astrophysics Data System Bibliographic Services.  

# References
