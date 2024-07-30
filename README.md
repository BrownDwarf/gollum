# gollum

`v0.4.1`

<a href="https://gollum-astro.readthedocs.io/en/latest/"><img src="https://img.shields.io/badge/Read-the%20docs-blue"></a>
<a href="https://pypi.org/project/gollum/"><img src="https://img.shields.io/badge/pip_install-gollum-yellow"></a>
<a href="https://ui.adsabs.harvard.edu/abs/2013A%26A...553A...6H/abstract"><img src="https://img.shields.io/badge/Works_with-PHOENIX-brightgreen"></a>
<a href="https://zenodo.org/record/1309035#.YL_SQoRKiV4"><img src="https://img.shields.io/badge/Works_with-Sonora_Bobcat-brightgreen"></a>

A microservice for programmatic access to precomputed synthetic spectral model grids in astronomy.

![gollum demo](docs/_static/gollum_resample.png?raw=true "Code Snippet of PHOENIX Synthetic Model Demonstration")

The goal of this repo is to provide a Python Application Programming Interface (API) to several different synthetic spectral models. `gollum` will be built on the astropy affiliated package [specutils](https://specutils.readthedocs.io/en/stable/), and will be inspired by the API design of [lightkurve](http://docs.lightkurve.org/). This project is loosely related to the parallel [muler](http://muler.readthedocs.io/) framework that is built on specutils and focuses on data. This project is all about models. The code itself and will have some overlap with functionality in [Starfish](https://starfish.readthedocs.io/en/latest/), and this project could one day become a microservice to Starfish, rather than duplicate code.

# Dashboard

We have a human-in-the-loop interactive dashboard which allows users to compare data to models. The current version of the dashboard supports the Sonora-Bobcat 2021 Models and the PHOENIX model grid.

![dashboard demo](https://user-images.githubusercontent.com/98151293/167173097-31427d83-f7fc-4146-a520-34e6b97b3b1b.gif)

This dashboard allows users to control sliders correlating with intrinsic properties (effective temperature, surface gravity, and metallicity) and extrinsic properties (rotational broadening, radial velocity, and a normalization scalar). From the selected intrinsic values, the dashboard can find the closest matching model (based on the closest existing point in a jagged 3D array of existing intrinsic values) and display it on screen so that the user can compare it with the real data. The data itself shows up as a blue plot, while the model is red, which will allow users to make by-eye fittings of the models to the data displayed.

There is some latency in the updating of the model's graph when the user moves certain sliders too quickly. This latency comes from the large amount of data points and the effect of the curse of dimensionality when it comes to the search for the nearest grid point based on intrinsic values that the dashboard must do with each update of the sliders. This latency mostly only applies when the user moves the sliders very quickly, however. More gradual movement of the sliders allows for relatively smooth updating of the model spectrum with minimal latency.




The package is under active development. Feel free to contibute by either raising issues here in GitHub or by submitting pull requests. If you have questions or need help, please also use GitHub issues to reach out to the development team and we will do our best to assist you.
