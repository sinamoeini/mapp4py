from __future__ import division, absolute_import, print_function

import sys, os, re, mapp

import sphinx
if sphinx.__version__ < "1.0.1":
    raise RuntimeError("Sphinx 1.0.1 or newer required")

needs_sphinx = '1.0'

# -----------------------------------------------------------------------------
# General configuration
# -----------------------------------------------------------------------------


sys.path.insert(0, os.path.abspath('../sphinxext'))

extensions = ['sphinx.ext.autodoc', 'sphinx.ext.pngmath', 'numpydoc',
              'sphinx.ext.intersphinx', 'sphinx.ext.coverage',
              'sphinx.ext.doctest', 'sphinx.ext.autosummary',
              'matplotlib.sphinxext.plot_directive']

templates_path = ['__templates']

source_suffix = '.rst'

project = u'MAPP'
copyright = u'2017, Sina Moeini'
author = u'Sina Moeini'


version = u'beta'
release = u'beta'

master_doc = 'index'
numpydoc_show_class_members = True
numpydoc_class_members_toctree = False

today_fmt = '%B %d, %Y'


default_role = "autolink"

exclude_dirs = []

add_function_parentheses = False

pygments_style = 'sphinx'


# -----------------------------------------------------------------------------
# HTML output
# -----------------------------------------------------------------------------

themedir = os.path.join(os.pardir, 'scipy-sphinx-theme', '_theme')
if not os.path.isdir(themedir):
    raise RuntimeError("Get the scipy-sphinx-theme first, "
                       "via git submodule init && git submodule update")

html_theme = 'scipy'
html_theme_path = [themedir]

html_theme_options = {
    "edit_link": False,
    "sidebar": "left",
    "scipy_org_logo": False,
    "rootlinks": []
}
html_sidebars = {'index': 'indexsidebar.html'}


html_title = "%s v%s Manual" % (project, version)
html_static_path = ['__static']
html_last_updated_fmt = '%b %d, %Y'

html_use_modindex = True
html_copy_source = False
html_domain_indices = False
html_file_suffix = '.html'

html_logo = "logo.png"

htmlhelp_basename = 'mapp'

pngmath_use_preview = True
pngmath_dvipng_args = ['-gamma', '1.5', '-D', '96', '-bg', 'Transparent']


# -----------------------------------------------------------------------------
# Texinfo output
# -----------------------------------------------------------------------------
texinfo_documents = [
    (master_doc, 'MAPP', u'MAPP Documentation',
     author, 'MAPP', 'One line description of project.',
     'Miscellaneous'),
]

# -----------------------------------------------------------------------------
# NumPy extensions
# -----------------------------------------------------------------------------

# If we want to do a phantom import from an XML file for all autodocs
phantom_import_file = 'dump.xml'

# Make numpydoc to generate plots for example sections
numpydoc_use_plots = True

# -----------------------------------------------------------------------------
# Autosummary
# -----------------------------------------------------------------------------

import glob
autosummary_generate = glob.glob("*.rst")

# -----------------------------------------------------------------------------
# Coverage checker
# -----------------------------------------------------------------------------
coverage_ignore_modules = r"""
    """.split()
coverage_ignore_functions = r"""
    test($|_) (some|all)true bitwise_not cumproduct pkgload
    generic\.
    """.split()
coverage_ignore_classes = r"""
    """.split()

coverage_c_path = []
coverage_c_regexes = {}
coverage_ignore_c_items = {}



