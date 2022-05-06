---
title: 'The `gollum` interface to precomputed synthetic spectral model grids'
tags:
  - Python
  - astronomy
  - spectroscopy
  - stars
  - echelle
authors:
  - name: Michael A. Gully-Santiago
    orcid: 0000-0002-4020-3457
    affiliation: 1
  - name: Caroline V. Morley
    orcid: 0000-0002-4404-0456
    affiliation: 1
  - name: More TBD
    affiliation: 2
affiliations:
 - name: The University of Texas at Austin Department of Astronomy, Austin, TX, USA
   index: 1
 - name: Second Affiliation goes here, Somewhere, USA
   index: 2
date: 16 February 2022
bibliography: paper.bib

# Optional fields if submitting to a AAS journal too, see this blog post:
# https://blog.joss.theoj.org/2018/12/a-new-collaboration-with-aas-publishing
#aas-doi: 10.3847/xxxxx <- update this with the DOI from AAS once you know it.
#aas-journal: Astrophysical Journal <- The name of the AAS journal.
---

# Summary

The spectra of stars, brown dwarfs, and planets are amazingly complex and information-rich.  Centuries of effort in astrophysics have distilled that complexity into the properties that control the bulk appearance of stellar and substellar spectra: effective temperature, surface gravity, iron abundance, and sometimes other compositional consituents.  Synthetic spectral models can now predict the appearance of a star's spectrum given merely these few properties.  The computational expense of modeling stars has driven practitioners to precompute the models finely in wavelength coordinates, but coarsely over trios or quartets of these fundamental properties.  Comparing these coarsely sampled grids of precomputed synthetic spectra to data remains a challenge.  The shear data volume of these grids has made it difficult for newcomers to build an intuition for how the spectra change with their inputs.  
The `gollum` framework provides a user-friendly Application Programming Interface (API) to load, manipulate, and visualize these voluminous precomputed synthetic spectral models.  The framework provides a first-of-its-kind interactive dashboard with user-controlled sliders that instantaneously update the fundamental, extrinsic, and compositional properties of the spectra.  `gollum` currently supports the PHOENIX grid of stellar spectra, and the Sonora-Bobcat grid of substellar brown dwarf and free-floating giant exoplanet spectra.  Support for other model grids is planned.  The `specutils`-based Python 3 framework interoperates easily with the astropy ecosystem of tools, including its sibling API for échelle spectra, `muler`.


# Statement of need

Citation to  `starfish` framework [@czekala15]


`gollum` depends on `astropy` [@astropy13; @astropy18], `numpy` [@harris2020array], `specutils`, `scipy` [@scipy2020], and others.

# Interoperation with the `muler` framework

`gollum` can interoperate with `muler` (Gully-Santiago et al. *in prep*)

# Supported model grids

We currently support two precomputed synthetic spectral models: PHOENIX  (XX Cite Husser), and Sonora-Bobcat (Marley et al. XX).

# Dashboard

We have a human in the loop interactive dashboard which allows users to compare data to models. Specifically, the current dashboard is designed to support the latest available Sonora-Bobcat 2021 Models, which allows users to control sliders correlating with intrinsic properties (effective temperature, surface gravity, and metallicity) and extrinsic properties (rotational broadening, radial velocity, and a normalization scalar). 

# Acknowledgements

This research has made use of NASA's Astrophysics Data System Bibliographic Services.  

# References
