#include <Python.h>
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#define PY_ARRAY_UNIQUE_SYMBOL ARRAY_API
#include <numpy/arrayobject.h>
#include "MAPP.h"
#include <mpi.h>
#include <dlfcn.h>
#include "comm.h"
#include "example.h"
#include "atoms_styles.h"
#include "md_styles.h"
#include "min_styles.h"
#include "dae_styles.h"
#include "export_styles.h"
#include "newton_gmres.h"
#define GET_FILE(file_name) reinterpret_cast<PyFileObject*>(PySys_GetObject((char*)#file_name))->f_fp
using namespace MAPP_NS;
/*--------------------------------------------*/
PyMethodDef MAPP::methods[]={[0 ... 2]={NULL}};
/*--------------------------------------------*/
void MAPP::setup_methods()
{
    methods[0].ml_name="pause_slave_out";
    methods[0].ml_meth=(PyCFunction)pause_out;
    methods[0].ml_flags=METH_NOARGS;
    methods[0].ml_doc=(char*)R"---(
    pause_slave_out()
    
    Pauses stdout & stderr of non-root processes
    
    Returns
    -------
    None
   
    Notes
    -----
    This function stops python's stdout and stderr of all processesors except the root processor

    
    Examples
    --------
    ::
     
        >>> pause_slave_out()

    )---";

    methods[1].ml_name="resume_slave_out";
    methods[1].ml_meth=(PyCFunction)resume_out;
    methods[1].ml_flags=METH_NOARGS;
    methods[1].ml_doc=(char*)R"---(
    resume_slave_out()
    
    Resumes stdout & stderr of non-root processes
    
    Returns
    -------
    None
   
    Notes
    -----
    This function resumes python's stdout and stderr of all if they have been paused by :py:function:`mapp.pause_slave_out`

    
    Examples
    --------
    ::
     
        >>> resume_slave_out()

    )---";
}
/*--------------------------------------------
 
 --------------------------------------------*/
PyObject* MAPP::pause_out(PyObject*)
{
    MPI_Barrier(MPI_COMM_WORLD);
    if(!glbl_rank) Py_RETURN_NONE;
    if(!__devnull__) __devnull__=fopen(devnull_path,"w");

    GET_FILE(stdout)=__devnull__;
    GET_FILE(stderr)=__devnull__;
    mapp_out=__devnull__;
    mapp_err=__devnull__;
    Py_RETURN_NONE;
}
/*--------------------------------------------
 
 --------------------------------------------*/
PyObject* MAPP::resume_out(PyObject*)
{
    MPI_Barrier(MPI_COMM_WORLD);
    GET_FILE(stdout)=__stdout__;
    GET_FILE(stderr)=__stderr__;
    mapp_out=__stdout__;
    mapp_err=__stderr__;
    if(__devnull__)
    {
        fclose(__devnull__);
        __devnull__=NULL;
    }
    Py_RETURN_NONE;
}
/*--------------------------------------------*/
const char* MAPP::devnull_path;
int MAPP::glbl_rank=0;
FILE* MAPP_NS::MAPP::__devnull__(NULL);
FILE* MAPP_NS::MAPP::__stdout__(NULL);
FILE* MAPP_NS::MAPP::__stderr__(NULL);
FILE* MAPP_NS::MAPP::mapp_out(NULL);
FILE* MAPP_NS::MAPP::mapp_err(NULL);
/*--------------------------------------------*/
PyMODINIT_FUNC initmapp(void)
{return MAPP_NS::MAPP::init_module();}

void MAPP::init_module(void)
{
    int mpi_initialized;
    MPI_Initialized(&mpi_initialized);
    if(!mpi_initialized)
    {
        dlopen("libmpi.so.0",RTLD_GLOBAL | RTLD_LAZY);
        MPI_Init(NULL,NULL);
        Py_AtExit([](){MPI_Finalize();});
    }
    
    PyObject* posixpath=PyImport_ImportModule("posixpath");
    PyObject* devnull_path_op=PyObject_GetAttrString(posixpath,"devnull");
    devnull_path=PyString_AsString(devnull_path_op);
    Py_DECREF(devnull_path_op);
    Py_DECREF(posixpath);
    glbl_rank=Communication::get_rank();
    
    mapp_out=__stdout__=GET_FILE(stdout);
    mapp_err=__stderr__=GET_FILE(stderr);
    
    import_array();
    setup_methods();
    PyObject* module=Py_InitModule3("mapp",methods,"MIT Atomistic Parallel Package");
    if(module==NULL) return;
    
    MAPP_MPI::setup_tp();
    if(PyType_Ready(&MAPP_MPI::TypeObject)<0) return;
    Py_INCREF(&MAPP_MPI::TypeObject);
    PyModule_AddObject(module,"mpi",reinterpret_cast<PyObject*>(&MAPP_MPI::TypeObject));
    
    
    if(LineSearch::setup_tp()<0) return;
    PyModule_AddObject(module,"ls",reinterpret_cast<PyObject*>(&LineSearch::TypeObject));
    
    
    if(LineSearchGoldenSection::setup_tp()<0) return;
    PyModule_AddObject(module,"ls_golden",reinterpret_cast<PyObject*>(&LineSearchGoldenSection::TypeObject));
    
    
    if(LineSearchBrent::setup_tp()<0) return;
    PyModule_AddObject(module,"ls_brent",reinterpret_cast<PyObject*>(&LineSearchBrent::TypeObject));
    
    
    if(LineSearchBackTrack::setup_tp()<0) return;
    PyModule_AddObject(module,"ls_bt",reinterpret_cast<PyObject*>(&LineSearchBackTrack::TypeObject));
    
    
    PyObject* md=MAPP::MD::init_module();
    if(md==NULL) return;
    PyModule_AddObject(module,"md",md);
    
    
    PyObject* dmd=MAPP::DMD::init_module();
    if(dmd==NULL) return;
    PyModule_AddObject(module,"dmd",dmd);
    
    pause_out(NULL);
}
/*--------------------------------------------*/
PyMethodDef MAPP::MD::methods[]={[0 ... 0]={NULL}};
/*--------------------------------------------*/
void MAPP::MD::setup_methods()
{
    /*
    ExamplePython::ml_phonon(methods[0]);
    ExamplePython::ml_phonon_1d(methods[1]);
    ExamplePython::ml_phonon_1dd(methods[2]);*/
}
/*--------------------------------------------
 
 --------------------------------------------*/
PyObject* MAPP::MD::init_module(void)
{
    setup_methods();
    PyObject* module=Py_InitModule3("mapp.md",methods,"Molecular Dynamics (MD) module");
    if(module==NULL) return NULL;
    
    
    if(AtomsMD::setup_tp()<0) return NULL;
    PyModule_AddObject(module,"atoms",reinterpret_cast<PyObject*>(&AtomsMD::TypeObject));
    
    if(MDNVT::setup_tp()<0) return NULL;
    PyModule_AddObject(module,"nvt",reinterpret_cast<PyObject*>(&MDNVT::TypeObject));
    
    
    if(MDNST::setup_tp()<0) return NULL;
    PyModule_AddObject(module,"nst",reinterpret_cast<PyObject*>(&MDNST::TypeObject));
    
    
    if(MDMuVT::setup_tp()<0) return NULL;
    PyModule_AddObject(module,"muvt",reinterpret_cast<PyObject*>(&MDMuVT::TypeObject));
    
    
    if(MinCG::setup_tp()<0) return NULL;
    PyModule_AddObject(module,"min_cg",reinterpret_cast<PyObject*>(&MinCG::TypeObject));
    
    
    if(MinLBFGS::setup_tp()<0) return NULL;
    PyModule_AddObject(module,"min_lbfgs",reinterpret_cast<PyObject*>(&MinLBFGS::TypeObject));
    
    
    if(ExportMD::setup_tp()<0) return NULL;
    PyModule_AddObject(module,"export",reinterpret_cast<PyObject*>(&ExportMD::TypeObject));
    
    
    if(ExportCFGMD::setup_tp()<0) return NULL;
    PyModule_AddObject(module,"export_cfg",reinterpret_cast<PyObject*>(&ExportCFGMD::TypeObject));
    
    
    return module;
}
/*--------------------------------------------*/
PyMethodDef MAPP::DMD::methods[]={[0 ... 3]={NULL}};
/*--------------------------------------------*/
void MAPP::DMD::setup_methods()
{
    
    ExamplePython::ml_mv_c(methods[0]);
    ExamplePython::ml_mv_c2(methods[1]);
    ExamplePython::ml_mv_c3(methods[2]);
    /*
    ExamplePython::ml_alpha(methods[1]);
    ExamplePython::ml_prt(methods[2]);
    ExamplePython::ml_delta_c(methods[3]);
     */
}
/*--------------------------------------------
 
 --------------------------------------------*/
PyObject* MAPP::DMD::init_module(void)
{
    setup_methods();
    PyObject* module=Py_InitModule3("mapp.dmd",methods,"Diffusive Molecular Dynamics (DMD) module");
    if(module==NULL) return NULL;
    
    
    if(AtomsDMD::setup_tp()<0) return NULL;
    PyModule_AddObject(module,"atoms",reinterpret_cast<PyObject*>(&AtomsDMD::TypeObject));
    
    
    if(MinCGDMD::setup_tp()<0) return NULL;
    PyModule_AddObject(module,"min_cg",reinterpret_cast<PyObject*>(&MinCGDMD::TypeObject));
    
    
    if(MinLBFGSDMD::setup_tp()<0) return NULL;
    PyModule_AddObject(module,"min_lbfgs",reinterpret_cast<PyObject*>(&MinLBFGSDMD::TypeObject));
    
    
    if(DAEBDF::setup_tp()<0) return NULL;
    PyModule_AddObject(module,"bdf",reinterpret_cast<PyObject*>(&DAEBDF::TypeObject));
    
    if(NewtonGMRES::setup_tp()<0) return NULL;
    PyModule_AddObject(module,"newton_gmres",reinterpret_cast<PyObject*>(&NewtonGMRES::TypeObject));
    
    if(ExportDMD::setup_tp()<0) return NULL;
    PyModule_AddObject(module,"export",reinterpret_cast<PyObject*>(&ExportDMD::TypeObject));
    
    if(ExportCFGDMD::setup_tp()<0) return NULL;
    PyModule_AddObject(module,"export_cfg",reinterpret_cast<PyObject*>(&ExportCFGDMD::TypeObject));
    
    
    return module;
}
