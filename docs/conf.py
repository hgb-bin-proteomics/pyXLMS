# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

import os
import sys

sys.path.insert(0, os.path.abspath("../src/"))

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = "pyXLMS"
copyright = "2024, Micha Johannes Birklbauer"
author = "Micha Johannes Birklbauer"
version = "0.2"
release = "0.2.0"

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.napoleon",
    "sphinx.ext.viewcode",
]

templates_path = ["_templates"]
exclude_patterns = ["build"]
root_doc = "index"
autosummary_generate = True
autodoc_default_options = {"members": True, "inherited-members": True}

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_title = "pyXLMS - A python package to process protein cross-linking data"
html_short_title = "pyXLMS"
html_logo = (
    "https://github.com/hgb-bin-proteomics/MSAnnika/raw/master/logo/icons/icon.png"
)
html_favicon = (
    "https://github.com/hgb-bin-proteomics/MSAnnika/raw/master/logo/icons/favicon.png"
)
html_theme = "pydata_sphinx_theme"
html_theme_options = {
    "logo": {
        "alt_text": "pyXLMS - Home",
        "text": "Python Crosslink Analysis",
        "image_light": "https://github.com/hgb-bin-proteomics/MSAnnika/raw/master/logo/icons/icon.png",
        "image_dark": "https://github.com/hgb-bin-proteomics/MSAnnika/raw/master/logo/icons/icon.png",
    },
    "external_links": [
        {"name": "MS Annika", "url": "https://github.com/hgb-bin-proteomics/MSAnnika"}
    ],
    "icon_links": [
        {
            "name": "GitHub",
            "url": "https://github.com/hgb-bin-proteomics/pyXLMS",
            "icon": "fa-brands fa-github",
            "type": "fontawesome",
        }
    ],
}
