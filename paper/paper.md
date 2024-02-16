---
title: '`gollum`: An intuitive programmatic and visual interface for precomputed synthetic spectral model grids'
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
  - name: Kyle Kaplan
    orcid: 0000-0001-6909-3856
    affiliation: 1
  - name: Karina Kimani-Stewart
    orcid: 0000-0002-6825-351X
    affiliation: 2
  - name: Diana Gonzalez-Argúeta
    affiliation: 3
affiliations:
  - name: The University of Texas at Austin | Department of Astronomy, Austin, TX, USA
    index: 1
  - name: Texas Tech University | Department of Physics and Astronomy, Lubbock, TX, USA
    index: 2
  - name: New Jersey Institute of Technology | Department of Physics, Newark, NJ, USA
    index: 3
date: 14 February 2024
bibliography: paper.bib
---

# Summary

The spectra of stars, brown dwarfs, and planets are complex and information-rich. Complexity can be reduced by distilling spectra down to their fundamental stellar parameters: effective temperature, surface gravity, and iron abundance, and sometimes others such as alpha element abundance, carbon-to-oxygen ratio,  or sedimentation efficiency of clouds. To a first order, these control the appearance of stellar and substellar spectra. Synthetic spectral models mimic stellar spectra across a grid of these fundamental parameters. Due to the computational constraints of modeling stellar and substellar spectra, model grids have high resolution in the spectral axis (wavelength coordinates) but are coarsely sampled over fundamental stellar parameters. Comparing these grids to data has been a challenge due to the large number of model grids, their size, their coarse sampling, and their dimensionality, making the relation between fundamental stellar parameters and spectral appearance somewhat unintuitive at times. This motivates the development of an intuitive, performant interface allowing astronomers to explore numerous model grids and compare them to data.

# Statement of need

`gollum` is a Python package for intuitive analysis and visualization of precomputed synthetic spectra. Its API is designed to have modules dedicated to each model grid it supports, with each module then containing classes for both individual spectra and bulk grid access. The programmatic interface to spectral analysis uses method-chaining to make `gollum` code very readable, taking inspiration from frameworks like lightkurve [@lightkurve:2018]. The visual interface in the form of interactive dashboards powered by `bokeh` [@bokeh2018] sports low-latency sliders and toggles. This allows users to tweak both fundamental stellar parameters and extrinsic parameters such as radial velocity and rotational broadening, building their intuition. `gollum`'s modularity allows for a wide range of model grids to potentially be supported, and its performance is optimized with libraries such as `numpy`, `scipy`, `astropy`, and `specutils` to allow for quick loading and processing of large amounts of data [@harris:2020; @virtanen:2020; @astropy:2022; @earl:2023].

`gollum` appeals to use cases ranging from entry-level astronomers to seasoned researchers with its combination of intuition-building and performance. The framework has been demonstrated on high-resolution spectra of brown dwarfs from Keck-NIRSPEC using `bdexda` [@kimani-stewart:2021], has experimental support for starspots, used in `acdc` [@cao:2022], and its programmatic interface has been used in the analysis of IGRINS spectra with `plotspec` [@kaplan2023]. It interoperates with `muler` [@gully-santiago:2022a], a similar framework designed for observed data from échelle spectrographs, and thanks to its dashboard, can also be used for flux calibration and empirical telluric correction. `gollum`'s programmatic interface is being used to create an extension to the `blase` framework [@gully-santiago:2022b; @shankar:2023] that will allow for inference of fundamental stellar parameters from observed spectra using interpretable machine learning and interpolation techniques, taking inspiration from other frameworks such as `starfish` [@czekala:2015] that also specialize in spectroscopic inference. `gollum` currently supports PHOENIX, Sonora (Bobcat, and Diamondback), and CoolTLUSTY model grids, and support for other model grids such as Sonora Elf Owl is planned [@husser:2013; @marley:2021; @morley:2024; @lacy:2023; @mukherjee:2024].

# Acknowledgements

This research has made use of NASA's Astrophysics Data System Bibliographic Services. This material is based upon work supported by the National Aeronautics and Space Administration under Grant Numbers 80NSSC21K0650 for the NNH20ZDA001N-ADAP:D.2 program,
and 80NSSC20K0257 for the XRP program issued through the Science Mission Directorate. S.S. acknowledges the National Science Foundation, which supported the work presented here under grant AST-1908892. 

# References
