# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -------------------------------------------------

project = 'RS-LMTO-ASA'
copyright = '2024, RS-LMTO Development Team'
author = 'RS-LMTO Development Team'

# The full version, including alpha/beta/rc tags
release = '1.0.0'
version = '1.0'

# -- General configuration ------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with sphinx.
extensions = [
    'sphinx.ext.mathjax',
    'sphinx.ext.autosectionlabel',
    'sphinx_rtd_theme',
]

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that should NOT be included
# when building the documentation. You can also use filename patterns like
# `*.rst` to filter out that set of files.
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']

# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and PDF output.  See the Sphinx documentation for
# a list of built-in themes.
html_theme = 'sphinx_rtd_theme'

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are put in a directory named '_static'
# under the build directory.
html_static_path = ['_static']

# Theme options for Read the Docs theme
html_theme_options = {
    'navigation_depth': 4,
    'collapse_navigation': False,
    'sticky_navigation': True,
    'includehidden': True,
    'titles_only': False,
    'prev_next_buttons_location': 'bottom',
    'style_external_links': False,
}

# -- Options for LaTeX output ------------------------------------------------

latex_elements = {
    'papersize': 'letterpaper',
    'pointsize': '10pt',
}

latex_documents = [
    ('index', 'RS-LMTO-ASA.tex', 'RS-LMTO-ASA Documentation', 
     'RS-LMTO Development Team', 'manual'),
]

# -- Options for text output -------------------------------------------------

text_documents = [
    ('index', 'RS-LMTO-ASA.txt', {'title': 'RS-LMTO-ASA Documentation'}),
]

# -- Additional settings -------------------------------------------------

# Use mathjax to render math (MathJax 3 configuration)
mathjax3_config = {
    'tex': {
        'inlineMath': [['$', '$'], ['\\(', '\\)']],
        'displayMath': [['$$', '$$'], ['\\[', '\\]']],
    }
}

# Enable cross-references between rst files
autosectionlabel_prefix_document = True
