import setuptools


setuptools.setup(
    name="gollum",
    version="0.4.0",
    author="gully",
    author_email="igully@gmail.com",
    description="A Python package for working with precomputed synthetic spectral models",
    long_description="A Python package for working with precomputed synthetic spectral model grids such as PHOENIX and Sonora",
    long_description_content_type="text/markdown",
    url="https://github.com/BrownDwarf/gollum",
    install_requires=[
        "numpy",
        "scipy",
        "astropy",
        "specutils>=1.8",
        "importlib_resources",
        "pandas",
        "matplotlib",
        "tqdm",
        "bokeh",
    ],
    packages=setuptools.find_packages(where="src"),
    package_dir={"": "src"},
    package_data={
        # If any package contains *.csv files, include them:
        "": ["*.csv"]
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.8",  # astropy 5 and up now requires Python 3.8+, sorry!
)
