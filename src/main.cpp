#include "api.h"
#include <stdio.h>
#include <iostream>
#include <mpi.h>
#include "MAPP.h"
/*
 note that this file is just running the package without
 importing through python mainly for debugging purposes
 */
#ifdef IS_PY3K
int main(int nargs,char* args[])
{
    MPI_Init(&nargs, &args);
    wchar_t** __args=NULL;
    if(nargs) __args=new wchar_t*[nargs];
    for(int i=0;i<nargs;i++) __args[i]=Py_DecodeLocale(args[i], NULL);
    PyImport_AppendInittab("mapp4py",MAPP_NS::MAPP::init_module);
    Py_Main(nargs,__args);
    delete [] __args;
    int mpi_finalized;
    MPI_Finalized(&mpi_finalized);
    if(!mpi_finalized) MPI_Finalize();
    return EXIT_SUCCESS;
}
#else
namespace MAPP_NS
{ PyMODINIT_FUNC initmapp4py(void);}
int main(int nargs,char* args[])
{
    MPI_Init(&nargs, &args);
    PyImport_AppendInittab("mapp4py",MAPP_NS::initmapp4py);
    Py_Main(nargs,args);
    int mpi_finalized;
    MPI_Finalized(&mpi_finalized);
    if(!mpi_finalized) MPI_Finalize();
    return EXIT_SUCCESS;
}
#endif

