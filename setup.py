from setuptools import setup , find_packages
from setuptools.extension import Extension
import os,sys
import subprocess
from sys import platform
import re

mpi_cxx='mpic++'


def mpi_cxx_params(mpi_cxx):
    result='';
    cmd=mpi_cxx+' -show';
    try:
        if sys.version_info < (3,0):
            result=subprocess.check_output(cmd.split());
        else:
            result=subprocess.check_output(cmd.split()).decode('utf-8');
    except:
        pass;
        err_msg='could not excute \"'+cmd+'\"';
        raise Exception(err_msg)


    lst=result.split();
    mpi_inc_dirs=[];
    mpi_lib_dirs=[];
    mpi_libs=[];
    mpi_ldds=[];
    mpi_flgs=[];
    nargs=lst.__len__();
    i=1;
    while i<nargs:
        if lst[i].startswith('-I'):
            mpi_inc_dirs+=[lst[i].replace('-I','')];
        elif lst[i].startswith('-L'):
            mpi_lib_dirs+=[lst[i].replace('-L','')];
        elif lst[i].startswith('-l'):
            mpi_libs+=[lst[i].replace('-l','')];
        elif lst[i].startswith('-Wl,'):
            mpi_ldds+=[lst[i].replace('-Wl,','')];
        elif lst[i]=='-Xlinker' and i+1<nargs:
            i+=1;
            mpi_ldds+=[lst[i]];
        else:
            mpi_flgs+=[lst[i]];
        i+=1;

    return mpi_inc_dirs,mpi_lib_dirs,mpi_libs,mpi_ldds,mpi_flgs



def get_libmpi_so(mpi_cxx):
    f = open('helloworld.cpp','w');
    
    f.write('#include <stdio.h>\n');
    f.write('#include <cstdlib>\n');
    f.write('#include <iostream>\n');
    f.write('#include <mpi.h>\n');
    f.write('int main(int nargs,char* args[])\n');
    f.write('{\n');
    f.write('    MPI_Init(&nargs,&args);\n');
    f.write('    int rank,size;\n');
    f.write('    MPI_Comm_rank(MPI_COMM_WORLD,&rank);\n');
    f.write('    MPI_Comm_size(MPI_COMM_WORLD,&size);\n');
    f.write('    printf(\"Hello world I am proc %d out of %d\\n\",rank,size);\n');
    f.write('    MPI_Finalize();\n');
    f.write('    return EXIT_SUCCESS;\n');
    f.write('}\n');
    f.close();
    

    
    returncode=0;
    cmd=mpi_cxx+' helloworld.cpp -o helloworld';
    try:
        returncode=subprocess.call(cmd.split());
    except:
        pass;
        raise Exception('could not excute \"'+cmd+'\"');
    if returncode!=0:
        raise Exception('could not excute \"'+cmd+'\"');

    if platform=='darwin':
        cmd='otool -L helloworld';
        ptrn='libmpi.*.dylib'
        regex_ptrn='libmpi.[0-9]*.dylib$'
    else:
        cmd='ldd helloworld';
        ptrn='libmpi.so.*'
        regex_ptrn='libmpi.so.[0-9]*$'


    try:
        if sys.version_info < (3,0):
            result=subprocess.check_output(cmd.split());
        else:
            result=subprocess.check_output(cmd.split()).decode('utf-8');
    except:
        pass;
        raise Exception('could not excute \"'+cmd+'\"')

    def find_ptrn(result,ptrn):
        for line in result.split():
            if re.search(ptrn,line):
                return line;
        return None;

    so_file=find_ptrn(result,regex_ptrn);
    if(so_file==None):
        raise Exception('could not find file \"'+ptrn+'\"');

    cmd='rm -f helloworld.cpp helloworld';
    try:
        returncode=subprocess.call(cmd.split());
    except:
        pass;
        raise Exception('could not excute \"'+cmd+'\"');
    if returncode!=0:
        raise Exception('could not excute \"'+cmd+'\"');
    return so_file.split('/').pop()




def numpy_params():
    try:
        from numpy import __path__ as npy_path
    except:
        pass;
        raise Exception('mapp4py requires numpy');

    npy_libs=['npymath']
    npy_lib_dirs=[npy_path[0]+'/core/lib']
    npy_inc_dirs=[npy_path[0]+'/core/include']
    return npy_inc_dirs,npy_lib_dirs,npy_libs




mpi_inc_dirs,mpi_lib_dirs,mpi_libs,mpi_ldds,mpi_flgs=mpi_cxx_params(mpi_cxx);
npy_inc_dirs,npy_lib_dirs,npy_libs=numpy_params();


libs=mpi_libs+npy_libs
inc_dirs=mpi_inc_dirs+npy_inc_dirs
lib_dirs=mpi_lib_dirs+npy_lib_dirs
ldds=mpi_ldds
flgs=mpi_flgs+['-std=c++11','-Wno-unused-but-set-variable','-Wno-reorder','-O3']

#print(libs)
#print(inc_dirs)
#print(lib_dirs)
#print(ldds)
#print(flgs)

cpp_files=[]
cpp_files+=['src/'+ each for each in os.listdir('src') if each.endswith('.cpp')]

module_mapp = Extension('mapp',
                    include_dirs=inc_dirs,
                    libraries=libs,
                    library_dirs=lib_dirs,
                    sources=cpp_files,
                    extra_compile_args=flgs)




setup (name ='mapp',
       version = '0.0.0.2',
       description = 'MIT Atomistic Parallel Package',
       author = 'Sina Moeini',
       author_email = 'sinam@mit.edu',
       url = 'https://github.com/sinamoeini/mapp4py',
       packages=find_packages(),
       ext_modules = [module_mapp])



