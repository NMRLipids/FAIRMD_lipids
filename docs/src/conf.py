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
from datetime import datetime

# Only on Read the Docs
if os.getenv("READTHEDOCS") == "True":
    repo_root = os.environ.get("READTHEDOCS_REPOSITORY_PATH")
    if repo_root:
        # make package importable if autodoc needs it
        if repo_root not in sys.path:
            sys.path.insert(0, repo_root)
        data_path = os.path.join(repo_root, "src", "fairmd", "lipids", "data", "ToyData")
        os.environ.setdefault("NMLDB_DATA_PATH", data_path)

import fairmd.lipids as dbl

# Directory containing this conf.py
here = os.path.dirname(__file__)
# Scripts directory
repo_root = os.path.abspath(os.path.join(here, "..", "..", ".."))

# Add to path:
sys.path.insert(0, repo_root)
year = datetime.now().year

# -- Project information -----------------------------------------------------

project = f"FAIRMD Lipids v{dbl.__version__}"
author = dbl.__author__
copyright = f"""{year}, {author}
    OSI Approved: GNU General Public License v3 (GPLv3)
    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version
   """


# The full version, including alpha/beta/rc tags
release = dbl.__version__
html_context = {
    "copyright_link": dbl.__url__ + "/blob/main/LICENSE.txt",
    "repo_link": dbl.__url__,
}

# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.viewcode",
    "sphinx.ext.napoleon",
    "myst_parser",
]

# Add any paths that contain templates here, relative to this directory.
templates_path = ["_templates"]

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ["**tests**", "*__init__.py"]


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
# html_theme = 'sphinxdoc'
# html_theme = 'classic'
html_theme = "sphinx_rtd_theme"

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
# html_static_path = ["_static"]

source_suffix = {
    ".rst": "restructuredtext",
    ".txt": "markdown",
    ".md": "markdown",
}

autodoc_member_order = "bysource"
