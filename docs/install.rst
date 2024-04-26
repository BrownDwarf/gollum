.. _installation:

**********************************
Installing the development version
**********************************




.. note::

    Conda installation is not yet available.

Pip installation is available ::

    pip install gollum

We recommend installing the developer version from source to help us with beta testing and to stay on the bleeding edge with bug fixes and updates. You can create a ready-made `conda` environment for gollum using the code below. ::

    git clone https://github.com/BrownDwarf/gollum.git
    cd gollum
    conda env create -f environment.yml
    conda activate gollum_dev
    pip install --editable .

Alternatively, if you want control over the Python version, you can create a new environment and install the dependencies using the following. You can replace "gollum_dev" and "3.12" with whatever name and compatible Python version you wish. ::
    
    git clone https://github.com/BrownDwarf/gollum.git
    cd gollum
    conda create -n "gollum_dev" python=3.12
    conda activate gollum_dev
    conda install -c conda-forge --file requirements.txt
    pip install --editable .

Eventually you can run the tests in the `tests/` directory to double-check that everything installed correctly.  Currently we are evaluating the best way to specify paths to voluminous libraries to make testing robust and quick across machines.  Stay tuned! ::

    py.test -vs test_precomputed.py  # should work for everyone
    py.test -vs test_phoenix.py # requires local PHOENIX models
    py.test -vs test_sonora.py # requires local Sonora models






Requirements
============

The project appears to work on modern Mac OS, Linux, and Windows operating systems, and officially supports Python 3.8 and above.