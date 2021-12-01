.. _modelgrids:

***********************
Downloading model grids
***********************



Sonora
======

Currently we only support the 2018 Sonora Model grids.  In a web browser, navigate to the `Sonora 2018 Zenodo website <https://zenodo.org/record/1309035#.YafR7oDML9A>`_.  Click on the ``Download`` button next to ``spectra.tar     1.2 GB``.  That will download the ``spectra.tar`` file to your ``~/Downloads/`` directory or equivalent default location.  Untar that directory.  The pattern should look something like below, depending on what your operating system defaults are, so copy-ing and past-ing this code is not likely to work--- try each step at a time.  ::

    tar -xzvf spectra.tar
    mv spectra ~/libraries/raw/Sonora
    


The end result should look something like this: ::

    $ ls ~/libraries/raw/Sonora/
    sp_t2100g3160nc_m0.0.gz  sp_t450g316nc_m0.0.gz    sp_t950g31nc_m0.0.gz
    sp_t2100g316nc_m0.0.gz   sp_t450g31nc_m0.0.gz     sp_t950g562nc_m0.0.gz
    sp_t2100g31nc_m0.0.gz    sp_t450g562nc_m0.0.gz    sp_t950g56nc_m0.0.gz


PHOENIX
=======


.. note::

    Downloading all the PHOENIX models can take hours or days! Start the downloading early.

The PHOENIX models total over 100 GB, and generally download at a relatively slow bandwidth.  The files are arranged into sub-directories for metallicity and alpha-element abundance.

You can get all of the PHOENIX models in one-fell-swoop from the command line if you have `wget` ::

    cd ~/Downloads
    wget -r -l 0 ftp://phoenix.astro.physik.uni-goettingen.de/HiResFITS/

If you don't have `wget` on your computer, please help the ``gollum`` grow by filing a GitHub Issue with how you resolved the problem, or what problems you are encountering.

As noted, this process will take a while as each individual file is painstakingly downloaded from a single German computer.  The commandline script, as written, will preserve the directory structure--- that's good! The gollum code demands that the directory structure is preserved.  Once it's all downloaded it should look like this ::


    # From inside this directory: ~/Downloads/phoenix.astro.physik.uni-goettingen.de/
    $ tree 
    .
    └── HiResFITS
        ├── PHOENIX-ACES-AGSS-COND-2011
        │   └── Z-0.0
        │       ├── lte02300-0.00-0.0.PHOENIX-ACES-AGSS-COND-2011-HiRes.fits
        │       ├── lte02300-0.50-0.0.PHOENIX-ACES-AGSS-COND-2011-HiRes.fits
        │       ├── lte02300-1.00-0.0.PHOENIX-ACES-AGSS-COND-2011-HiRes.fits
        │       ├── lte02300-1.50-0.0.PHOENIX-ACES-AGSS-COND-2011-HiRes.fits
                                        [ ... ]
        │       ├── lte04000-2.50-0.0.PHOENIX-ACES-AGSS-COND-2011-HiRes.fits
        │       ├── lte04000-3.00-0.0.PHOENIX-ACES-AGSS-COND-2011-HiRes.fits
        │       └── lte04000-3.50-0.0.PHOENIX-ACES-AGSS-COND-2011-HiRes.fits
        │                                [ ... ]
        │    ├── Z-0.0
        │    ├── Z+0.5
        │    ├── Z-0.5
        │    ├── Z+1.0
        │    ├── Z-1.0
        │    ├── Z-1.5
        │    ├── Z-2.0
        │    ├── Z-3.0
        │    └── Z-4.0
        └── WAVE_PHOENIX-ACES-AGSS-COND-2011.fits


While you `could` leave these directories in say ``~/Downloads/phoenix.astro.physik.uni-goettingen.de/``, I recommend making a more permanent and recognizable home for these models.  In particular `gollum` attempts to search a single default path for models: ``~/libraries/raw/``, where the tilde ``~/`` denotes your home directory ::


    mv ~/Downloads/phoenix.astro.physik.uni-goettingen.de/HiResFITS/PHOENIX-ACES-AGSS-COND-2011 ~/libraries/raw/PHOENIX/
    ls ~/libraries/raw/PHOENIX/
    Z-0.0  Z+0.5  Z-0.5  Z+1.0  Z-1.0  Z-1.5  Z-2.0  Z-3.0  Z-4.0

Finally, you must copy the wavelength file into this directory as well.  Notice that this placement breaks the native directory structure, so you must complete this step in order for ``gollum`` to work. ::

    mv ~/Downloads/phoenix.astro.physik.uni-goettingen.de/HiResFITSWAVE_PHOENIX-ACES-AGSS-COND-2011.fits ~/libraries/raw/PHOENIX/
    ls ~/libraries/raw/PHOENIX/
    WAVE_PHOENIX-ACES-AGSS-COND-2011.fits  Z-0.0  Z+0.5  Z-0.5  Z+1.0  Z-1.0  Z-1.5  Z-2.0  Z-3.0  Z-4.0

