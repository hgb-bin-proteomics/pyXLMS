# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

import os
import sys

sys.path.insert(0, os.path.abspath("../src/"))

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = "pyXLMS-pydantic"
copyright = "2026, Bioinformatics Research Group, FH Oberösterreich Campus Hagenberg"
author = "Micha Johannes Birklbauer"
version = "1.8"
release = "1.8.11"

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.napoleon",
    "sphinx.ext.viewcode",
    "myst_parser",
    "sphinx_copybutton",
]

templates_path = ["_templates"]
exclude_patterns = ["build"]
root_doc = "index"
autosummary_generate = True
autodoc_default_options = {"members": True, "inherited-members": False}
python_maximum_signature_line_length = 88

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_title = "pyXLMS - A python package to process protein cross-linking data"
html_short_title = "pyXLMS-pydantic"
html_logo = "icons/icon.png"
html_favicon = "icons/favicon.png"
html_theme = "pydata_sphinx_theme"
html_show_sourcelink = False
html_theme_options = {
    "logo": {
        "alt_text": "pyXLMS",
        "text": "pyXLMS",
        "image_light": "icons/icon.png",
        "image_dark": "icons/icon.png",
    },
    "header_links_before_dropdown": 4,
    "external_links": [],
    "icon_links": [
        {
            "name": "GitHub",
            "url": "https://github.com/hgb-bin-proteomics/pyXLMS-pydantic",
            "icon": "fa-brands fa-github",
            "type": "fontawesome",
        },
    ],
    "show_toc_level": 2,
    "use_edit_page_button": False,
    "primary_sidebar_end": ["indices.html"],
}
html_context = {
    "github_url": "https://github.com",
    "github_user": "hgb-bin-proteomics",
    "github_repo": "pyXLMS-pydantic",
    "github_version": "master",
    "doc_path": "sphinx",
    "default_mode": "auto",
}
