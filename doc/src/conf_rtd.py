import sys
import os
extensions = ['sphinx.ext.autodoc', 'sphinx.ext.imgmath', 'numpydoc',
              'sphinx.ext.intersphinx', 'sphinx.ext.coverage',
              'sphinx.ext.doctest', 'sphinx.ext.autosummary',
              'matplotlib.sphinxext.plot_directive']

templates_path = ['__templates']

source_suffix = '.rst'

#numpydoc_show_class_members = False
numpydoc_class_members_toctree = False

#source_encoding = 'utf-8-sig'

master_doc = 'index'

project = u'MAPP'
copyright = u'2017, Sina Moeini'
author = u'Sina Moeini'

version = u'0.00'
release = u'0.00'

language = None

exclude_patterns = []



pygments_style = 'sphinx'



todo_include_todos = False

# -- HTML output ----------------------------------------------------------
import sphinx_rtd_theme
html_theme = 'sphinx_rtd_theme'
html_theme_path = [sphinx_rtd_theme.get_html_theme_path()]

html_theme_options = {
    'sticky_navigation': False,  # Set to False to disable the sticky nav while scrolling.
    'logo_only': True,  # if we have a html_logo below, this shows /only/ the logo with no title text
    'display_version': False
}

html_logo = "logo.png"

html_static_path = ['__static']

htmlhelp_basename = 'MAPP'

# -- Options for LaTeX output ---------------------------------------------

latex_elements = {

}

latex_documents = [
    (master_doc, 'MAPP.tex', u'MAPP Documentation',
     u'Sina Moeini', 'manual'),
]


# -- Options for manual page output ---------------------------------------

man_pages = [
    (master_doc, 'mapp', u'MAPP Documentation',
     [author], 1)
]



# -- Options for Texinfo output -------------------------------------------

texinfo_documents = [
    (master_doc, 'MAPP', u'MAPP Documentation',
     author, 'MAPP', 'One line description of project.',
     'Miscellaneous'),
]

# -----------------------------------------------------------------------------
# Autosummary
# -----------------------------------------------------------------------------

import glob
autosummary_generate = glob.glob("*.rst")


