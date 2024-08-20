# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

import os
import sys
sys.path.insert(0, os.path.abspath('..'))

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'Example'
copyright = '2024, Micha Johannes Birklbauer'
author = 'Micha Johannes Birklbauer'
release = '1.0.0'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.autosummary',
    'sphinx.ext.napoleon',
    'sphinx.ext.viewcode'
]

templates_path = ['_templates']
exclude_patterns = ['build']
root_doc = "index"
autosummary_generate = True
autodoc_default_options = {"members": True, "inherited-members": True}

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'pydata_sphinx_theme'
html_static_path = ['_static']
html_theme_options = {
    "icon_links": [
        {
            "name": "GitHub",
            "url": "https://github.com/michabirklbauer/python_template",
            "icon": "fa-brands fa-github",
            "type": "fontawesome",
        }
    ]
}
