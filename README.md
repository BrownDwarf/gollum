# gollum

A microservice for programmatic access to precomputed synthetic spectral model grids in astronomy.

The goal of this repo is to provide a Python Application Programming Interface (API) to several different synthetic spectral models.  `gollum` will be built on the astropy affiliated package [specutils](https://specutils.readthedocs.io/en/stable/), and will be inspired by the API design of [lightkurve](http://docs.lightkurve.org/).  This project is loosely related to the parallel [muler](http://docs.lightkurve.org/) framework that is built on specutils and focuses on data.  This project is all about models.  The code itself and will have some overlap with functionality in [Starfish](https://starfish.readthedocs.io/en/latest/), and this project could one day become a microservice to Starfish, rather than duplicate code.

The package is under active development.  Please contribute by commenting on the Issues, providing Pull Requests, and starring the GitHub repo.
