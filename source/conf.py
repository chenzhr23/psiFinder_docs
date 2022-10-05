# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'psiFinder'
copyright = '2022, chenzhr23'
author = 'chenzhr23'
release = 'v1.0'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

templates_path = ['_templates']

exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

# html_theme = 'sphinxdoc' #nice
# html_theme = 'pyramid' #nice
# html_theme = 'nature' #nice
# html_theme = 'haiku' #OK
# html_theme = 'classic' #nice
html_theme = 'bizstyle' #nice
# html_theme = 'alabaster' #nice

html_theme_options = {'body_max_width': '90%'}

html_sidebars = {
    '**': [
        'globaltoc.html',
    ]
}

# For md files
source_suffix = {
    '.rst': 'restructuredtext',
    '.txt': 'markdown',
    '.md': 'markdown',
}

