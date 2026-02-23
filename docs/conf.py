# Configuration file for the Sphinx documentation builder.
# https://www.sphinx-doc.org/en/master/usage/configuration.html

import os
import sys

# -- Path setup --------------------------------------------------------------
# Add the project root so Sphinx can find the 'starterDD' package
sys.path.insert(0, os.path.abspath('..'))

# -- Project information -----------------------------------------------------
project = 'starterDD'
copyright = '2026, R. Guasch'
author = 'R. Guasch'
release = '0.1.0'

# -- General configuration ---------------------------------------------------
extensions = [
    'sphinx.ext.autodoc',           # Auto-generate docs from docstrings
    'sphinx.ext.autosummary',       # Generate summary tables
    'sphinx.ext.napoleon',          # Support Google/NumPy-style docstrings
    'sphinx.ext.viewcode',          # Add [source] links to highlighted source code
    'sphinx.ext.intersphinx',       # Link to other projects' documentation
    'sphinx_autodoc_typehints',     # Better type hint rendering
]

# Napoleon settings (for Google/NumPy style docstrings)
napoleon_google_docstrings = True
napoleon_numpy_docstrings = True
napoleon_include_init_with_doc = True
napoleon_include_private_with_doc = False

# Autodoc settings
autodoc_default_options = {
    'members': True,
    'undoc-members': True,
    'show-inheritance': True,
    'member-order': 'bysource',
}
autodoc_typehints = 'description'
autosummary_generate = True

# Mock imports for optional dependencies that may not be installed
autodoc_mock_imports = ['glow']

# Templates path
templates_path = ['_templates']

# Patterns to exclude
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']

# -- Options for HTML output -------------------------------------------------
html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']

# Theme options
html_theme_options = {
    'navigation_depth': 4,
    'titles_only': False,
    'collapse_navigation': False,
}

# Intersphinx mapping (link to external docs)
intersphinx_mapping = {
    'python': ('https://docs.python.org/3', None),
    'numpy': ('https://numpy.org/doc/stable/', None),
}
