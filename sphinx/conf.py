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
copyright = "2025, Micha Johannes Birklbauer"
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
    "myst_parser",
]

templates_path = ["_templates"]
exclude_patterns = ["build"]
root_doc = "index"
autosummary_generate = True
autodoc_default_options = {"members": True, "inherited-members": True}
python_maximum_signature_line_length = 88

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
        "alt_text": "pyXLMS",
        "text": "pyXLMS",
        "image_light": "https://github.com/hgb-bin-proteomics/MSAnnika/raw/master/logo/icons/icon.png",
        "image_dark": "https://github.com/hgb-bin-proteomics/MSAnnika/raw/master/logo/icons/icon.png",
    },
    "external_links": [
        {"name": "Documentation", "url": "https://hgb-bin-proteomics.github.io/pyXLMS-docs"},
        {"name": "Contact", "url": "https://github.com/hgb-bin-proteomics/pyXLMS?tab=readme-ov-file#contact"},
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
    "use_edit_page_button": True,
    "primary_sidebar_end": ["indices.html"],
}
html_context = {
    "github_url": "https://github.com",
    "github_user": "hgb-bin-proteomics",
    "github_repo": "pyXLMS",
    "github_version": "master",
    "doc_path": "sphinx",
}
