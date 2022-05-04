# gollum

`v0.2.1`

<a href="https://gollum-astro.readthedocs.io/en/latest/"><img src="https://img.shields.io/badge/Read-the%20docs-blue"></a>
<a href="https://pypi.org/project/gollum/"><img src="https://img.shields.io/badge/pip_install-gollum-yellow"></a>
<a href="https://ui.adsabs.harvard.edu/abs/2013A%26A...553A...6H/abstract"><img src="https://img.shields.io/badge/Works_with-PHOENIX-brightgreen"></a>
<a href="https://zenodo.org/record/1309035#.YL_SQoRKiV4"><img src="https://img.shields.io/badge/Works_with-Sonora_Bobcat-brightgreen"></a>

A microservice for programmatic access to precomputed synthetic spectral model grids in astronomy.

![gollum demo](docs/_static/gollum_resample.png?raw=true "Code Snippet of PHOENIX Synthetic Model Demonstration")

The goal of this repo is to provide a Python Application Programming Interface (API) to several different synthetic spectral models. `gollum` will be built on the astropy affiliated package [specutils](https://specutils.readthedocs.io/en/stable/), and will be inspired by the API design of [lightkurve](http://docs.lightkurve.org/). This project is loosely related to the parallel [muler](http://muler.readthedocs.io/) framework that is built on specutils and focuses on data. This project is all about models. The code itself and will have some overlap with functionality in [Starfish](https://starfish.readthedocs.io/en/latest/), and this project could one day become a microservice to Starfish, rather than duplicate code.

The package is under active development. Please contribute by commenting on the Issues, providing Pull Requests, and starring the GitHub repo.
