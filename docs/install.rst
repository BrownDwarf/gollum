.. _installation:

**********************************
Installing the development version
**********************************




.. note::

    Conda installation is not yet available.

Pip installation is currently pending, and should eventually work ::

    pip install gollum

Currently we recommend installing the developer version from source to help us with beta testing.


To install `gollum` from source ::

    git clone https://github.com/BrownDwarf/gollum.git
    cd gollum
    conda env create -f environment.yml
    conda activate gollum
    python setup.py develop


You can run the tests in the `tests/` directory to double-check that everything installed correctly::

    py.test -vs



Requirements
============

The project may work with a variety of Python 3 minor versions. We have tested it on Linux Operating Systems.
