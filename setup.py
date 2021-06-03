import setuptools


setuptools.setup(
    name="gollum",
    version="0.0.1",
    author="gully",
    author_email="igully@gmail.com",
    description="A Python package for working with precomputed synthetic spectral models",
    long_description="A Python package for working with precomputed synthetic spectral models such as PHOENIX and Sonora-Bobcat",
    long_description_content_type="text/markdown",
    url="https://github.com/BrownDwarf/gollum",
    install_requires=["numpy", "scipy", "astropy", "specutils", "pandas", "matplotlib"],
    packages=setuptools.find_packages(where="src"),
    package_dir={"": "src"},
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.7",
)
