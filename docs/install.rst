.. _installation:

**********************************
Installing the development version
**********************************




.. note::

    Conda installation is not yet available.

Pip installation is available ::

    pip install gollum

Currently we recommend installing the developer version from source to help us with beta testing.  You can install gollum into its own isolated conda python environment using the instructions below.  Line 3 will create a conda environment titled *gollum_dev* that contains all the dependencies required to install `gollum`, and even all the affiliated code to make the documentation and run unit tests from scratch.  Very useful!  You can activate this conda environment using the code on line 4 `conda activate gollum_dev`.  You could alternatively activate some other existing conda environment on your computer, but it is not guaranteed to have all of the extra packages.  Our sibling package `muler` offers an identical conda environment to this one, so if you made one for muler already, you can do: `conda activate muler_dev` ::

    git clone https://github.com/BrownDwarf/gollum.git
    cd gollum
    conda env create -f environment.yml
    conda activate gollum_dev # or conda activate muler_dev if you have that already
    python setup.py develop


Eventually you can run the tests in the `tests/` directory to double-check that everything installed correctly.  Currently we are evaluating the best way to specify paths to voluminous libraries to make testing robust and quick across machines.  Stay tuned! ::

    py.test -vs test_precomputed.py  # should work for everyone
    py.test -vs test_phoenix.py # requires local PHOENIX models
    py.test -vs test_sonora.py # requires local Sonora models






Requirements
============

The project appears to work on modern Mac OS, Linux, and Windows operating systems, and has been tested for Python 3.7 and above.  It may work on Python 3.6, and will not work on Python 2.