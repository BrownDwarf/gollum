# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
import os
import sys

sys.path.insert(0, os.path.abspath("../src/"))

# -- Project information -----------------------------------------------------

project = "gollum"
copyright = "2021, 2022, gully"
author = "gully"

# The full version, including alpha/beta/rc tags
release = "0.3.0"


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.coverage",
    "sphinx.ext.napoleon",
    "sphinx.ext.todo",
    "sphinx.ext.mathjax",
    "nbsphinx",
    "sphinx.ext.githubpages",
    "sphinx.ext.ifconfig",
    "sphinx_gallery.load_style",
    "numpydoc",
]

# Add any paths that contain templates here, relative to this directory.
templates_path = ["_templates"]

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store", "**/.ipynb_checkpoints"]
nbsphinx_timeout = 60

autosummary_generate = True
html_show_sourcelink = True
numpydoc_show_class_members = False

nbsphinx_thumbnails = {
    "tutorials/best_fit_for_fixed_template": "_static/extrinsic_grid_search.png",
    "tutorials/PHOENIX_tutorial": "_static/gollum_units.png",
    "tutorials/Divide_by_a_blackbody": "_static/blackbody.png",
    "tutorials/gollum_demo_Sonora_and_BDSS": "_static/sliders.png",
    "tutorials/Phoenix_dashboard_demo": "_static/PHOENIX_dashboard_thumbnail.png",
}

# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = "sphinx_material"

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ["_static"]


# Set link name generated in the top bar.
html_title = "Home"

nbsphinx_codecell_lexer = "python"

# Material theme options (see theme.conf for more information)
html_theme_options = {
    "base_url": "https://gollum.readthedocs.io/",
    "nav_title": "gollum docs",
    "nav_links": [
        {"title": "Quickstart", "href": "quickstart", "internal": True},
        {"title": "Installation", "href": "install", "internal": True},
        {"title": "API", "href": "api", "internal": True},
        {"title": "Tutorials", "href": "tutorials/index", "internal": True},
    ],
    "color_primary": "green",
    "color_accent": "yellow",
    "theme_color": "d35400",
    "repo_url": "https://github.com/BrownDwarf/gollum/",
    "repo_name": "gollum",
    "repo_type": "github",
    "master_doc": True,
    "globaltoc_depth": 2,
    "globaltoc_collapse": True,
    "globaltoc_includehidden": True,
    "logo_icon": "&#xE85C",
}

html_sidebars = {
    "**": ["logo-text.html", "globaltoc.html", "localtoc.html", "searchbox.html"]
}


html_use_index = True
html_domain_indices = True

nbsphinx_kernel_name = "python"
nbsphinx_execute = "never"
