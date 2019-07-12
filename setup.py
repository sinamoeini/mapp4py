from setuptools import setup , find_packages
from setuptools.extension import Extension
import os,sys
import subprocess

mpi_cxx='mpic++'





try:
    from numpy import __path__ as npy_path
except:
    pass;
    raise Exception('mapp4py requires numpy');

def mpi_cxx_params(mpi_cxx):
    result='';
    try:
        if sys.version_info < (3,0):
            result=subprocess.check_output([mpi_cxx,'-show']);
        else:
            result=subprocess.check_output([mpi_cxx,'-show']).decode('utf-8');
    except:
        pass;
        err_msg='could not excute \"'+mpi_cxx+' '+'-show'+'\"';
        raise Exception(err_msg)


    lst=result.split();
    mpi_inc_dir=[];
    mpi_lib_dir=[];
    mpi_lib=[];
    mpi_extra_ld_args=[];
    for i in range(1,lst.__len__()):
        if lst[i].startswith('-I'):
            mpi_inc_dir+=[lst[i].replace('-I','')];
        elif lst[i].startswith('-L'):
            mpi_lib_dir+=[lst[i].replace('-L','')];
        elif lst[i].startswith('-l'):
            mpi_lib+=[lst[i].replace('-l','')];
        elif lst[i].startswith('-Wl,'):
            mpi_extra_ld_args+=[lst[i].replace('-Wl,','')];

    return mpi_inc_dir,mpi_lib_dir,mpi_lib,mpi_extra_ld_args



npy_lib=['npymath']
npy_lib_dir = [npy_path[0]+'/core/lib']
npy_inc_dir = [npy_path[0]+'/core/include']


mpi_inc_dir,mpi_lib_dir,mpi_lib,mpi_extra_ld_args=mpi_cxx_params(mpi_cxx);





libs=mpi_lib+npy_lib
inc_dirs=npy_inc_dir+mpi_inc_dir
lib_dirs=npy_lib_dir+mpi_lib_dir
extra_ld_args=mpi_extra_ld_args

cpp_files=[]
cpp_files+=['src/'+ each for each in os.listdir('src') if each.endswith('.cpp')]

module_mapp = Extension('mapp',
                    include_dirs=inc_dirs,
                    libraries=libs,
                    library_dirs=lib_dirs,
                    sources=cpp_files,
                    extra_compile_args=['-std=c++11','-Wno-unused-but-set-variable','-Wno-reorder','-O3'])




setup (name ='mapp',
       version = '0.0.0.2',
       description = 'MIT Atomistic Parallel Package',
       author = 'Sina Moeini',
       author_email = 'sinam@mit.edu',
       url = 'https://github.com/sinamoeini/mapp4py',
       packages=find_packages(),
       ext_modules = [module_mapp])

