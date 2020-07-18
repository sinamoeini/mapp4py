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

project = 'MAPP'
copyright = '2008-2009, The Scipy community'

import numpy
version = re.sub(r'(\d+\.\d+)\.\d+(.*)', r'\1\2','0.0.0')
version = re.sub(r'(\.dev\d+).*?$', r'\1', version)
release = '0.0.0'
print("%s %s" % (version, release))

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

htmlhelp_basename = 'mapp'

pngmath_use_preview = True
pngmath_dvipng_args = ['-gamma', '1.5', '-D', '96', '-bg', 'Transparent']


# -----------------------------------------------------------------------------
# LaTeX output
# -----------------------------------------------------------------------------
_stdauthor = 'Written by the NumPy community'
latex_documents = [
  ('reference/index', 'numpy-ref.tex', 'NumPy Reference',
   _stdauthor, 'manual'),
  ('user/index', 'numpy-user.tex', 'NumPy User Guide',
   _stdauthor, 'manual'),
]

latex_preamble = r'''
\usepackage{amsmath}
\DeclareUnicodeCharacter{00A0}{\nobreakspace}

% In the parameters section, place a newline after the Parameters
% header
\usepackage{expdlist}
\let\latexdescription=\description
\def\description{\latexdescription{}{} \breaklabel}

% Make Examples/etc section headers smaller and more compact
\makeatletter
\titleformat{\paragraph}{\normalsize\py@HeaderFamily}%
            {\py@TitleColor}{0em}{\py@TitleColor}{\py@NormalColor}
\titlespacing*{\paragraph}{0pt}{1ex}{0pt}
\makeatother

% Fix footer/header
\renewcommand{\chaptermark}[1]{\markboth{\MakeUppercase{\thechapter.\ #1}}{}}
\renewcommand{\sectionmark}[1]{\markright{\MakeUppercase{\thesection.\ #1}}}
'''

latex_use_modindex = False


# -----------------------------------------------------------------------------
# Texinfo output
# -----------------------------------------------------------------------------

texinfo_documents = [
  ('mapp', 'MAPP Documentation', _stdauthor, 'MAPP',
   "NumPy: array processing for numbers, strings, records, and objects.",
   'Programming',
   1),
]


# -----------------------------------------------------------------------------
# Intersphinx configuration
# -----------------------------------------------------------------------------
"""
intersphinx_mapping = {
    'python': ('https://docs.python.org/dev', None),
    'scipy': ('https://docs.scipy.org/doc/scipy/reference', None),
    'matplotlib': ('http://matplotlib.org', None)
}
"""

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
autosummary_generate = glob.glob("reference/*.rst")

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


# -----------------------------------------------------------------------------
# Plots
# -----------------------------------------------------------------------------
plot_pre_code = """
import numpy as np
np.random.seed(0)
"""
plot_include_source = True
plot_formats = [('png', 100), 'pdf']

import math
phi = (math.sqrt(5) + 1)/2

plot_rcparams = {
    'font.size': 8,
    'axes.titlesize': 8,
    'axes.labelsize': 8,
    'xtick.labelsize': 8,
    'ytick.labelsize': 8,
    'legend.fontsize': 8,
    'figure.figsize': (3*phi, 3),
    'figure.subplot.bottom': 0.2,
    'figure.subplot.left': 0.2,
    'figure.subplot.right': 0.9,
    'figure.subplot.top': 0.85,
    'figure.subplot.wspace': 0.4,
    'text.usetex': False,
}


# -----------------------------------------------------------------------------
# Source code links
# -----------------------------------------------------------------------------

import inspect
from os.path import relpath, dirname

for name in ['sphinx.ext.linkcode', 'numpydoc.linkcode']:
    try:
        __import__(name)
        extensions.append(name)
        break
    except ImportError:
        pass
else:
    print("NOTE: linkcode extension not found -- no links to source generated")

def linkcode_resolve(domain, info):
    """
    Determine the URL corresponding to Python object
    """
    if domain != 'py':
        return None

    modname = info['module']
    fullname = info['fullname']

    submod = sys.modules.get(modname)
    if submod is None:
        return None

    obj = submod
    for part in fullname.split('.'):
        try:
            obj = getattr(obj, part)
        except:
            return None

    try:
        fn = inspect.getsourcefile(obj)
    except:
        fn = None
    if not fn:
        return None

    try:
        source, lineno = inspect.getsourcelines(obj)
    except:
        lineno = None

    if lineno:
        linespec = "#L%d-L%d" % (lineno, lineno + len(source) - 1)
    else:
        linespec = ""

    fn = relpath(fn, start=dirname(numpy.__file__))

    if 'dev' in numpy.__version__:
        return "http://github.com/numpy/numpy/blob/master/numpy/%s%s" % (
           fn, linespec)
    else:
        return "http://github.com/numpy/numpy/blob/v%s/numpy/%s%s" % (
           numpy.__version__, fn, linespec)
