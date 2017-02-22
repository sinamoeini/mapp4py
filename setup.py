from distutils.core import setup, Extension
import os
from numpy import __path__ as npy_path
npy_lib = npy_path[0] + '/core/lib'
npy_include = npy_path[0] + '/core/include'



cpp_files=[]
cpp_files += ['src/'+ each for each in os.listdir('src') if each.endswith('.cpp')]

#cpp_files =['src/ff_lj.cpp']

module_mapp = Extension('mapp',
                    libraries = ['npymath','mpi','mpi_cxx','util'],
                    library_dirs = [npy_lib,'/usr/local/lib'],
                    include_dirs = [npy_include,'/usr/local/include'],
                    sources = cpp_files,
		    extra_compile_args=['-std=c++11','-Wno-reorder','-O0'])

setup (name ='mapp',
       version = '0.0.0',
       description = 'MIT Atomistic Parallel Package',
       author = 'Sina Moeini',
       author_email = 'sinam@mit.edu',
       url = 'https://github.com/sinamoeini/mapp4py',
       ext_modules = [module_mapp])
