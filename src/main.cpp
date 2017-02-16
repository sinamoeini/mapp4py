#include <stdio.h>
#include <iostream>
#include <mpi.h>
#include <Python/Python.h>
#include "MAPP.h"

int main(int nargs,char* args[])
{
    MPI_Init(&nargs, &args);
    PyImport_AppendInittab("mapp",MAPP_NS::MAPP::init_module);
    
    Py_Main(nargs,args);
    MPI_Finalize();
    return EXIT_SUCCESS;
}
