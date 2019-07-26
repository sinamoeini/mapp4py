#include <Python.h>
#include "py_compat.h"
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#define PY_ARRAY_UNIQUE_SYMBOL ARRAY_API
#include <numpy/arrayobject.h>
#include "MAPP.h"
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <mpi.h>
#ifdef LIBMPI_SONAME
#include <dlfcn.h>
#endif
#include "comm.h"
#include "example.h"
#include "atoms_styles.h"
#include "md_styles.h"
#include "min_styles.h"
#include "dae_styles.h"
#include "export_styles.h"
#ifdef POTFIT
#include "potfit.h"
#endif
#include "import_eam.h"
using namespace MAPP_NS;
/*--------------------------------------------*/
PyMethodDef MAPP::methods[]=EmptyPyMethodDef(6);
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
    This function resumes python's stdout and stderr of all if they have been paused by :py:function:`mapp4py.pause_slave_out`

    
    Examples
    --------
    ::
     
        >>> resume_slave_out()

    )---";
    
    ImportEAM::ml_read_eam(methods[2],methods[3],methods[4]);
    
}
/*--------------------------------------------
 
 --------------------------------------------*/
PyObject* MAPP::pause_out(PyObject*)
{
    MPI_Barrier(MPI_COMM_WORLD);
    if(!glbl_rank) Py_RETURN_NONE;
    if(!__devnull__) __devnull__=fopen(devnull_path,"w");

    GET_FILE_FROM_PY(stdout)=__devnull__;
    GET_FILE_FROM_PY(stderr)=__devnull__;
    mapp_out=__devnull__;
    mapp_err=__devnull__;
    Py_RETURN_NONE;
}
/*--------------------------------------------
 
 --------------------------------------------*/
PyObject* MAPP::resume_out(PyObject*)
{
    MPI_Barrier(MPI_COMM_WORLD);
    GET_FILE_FROM_PY(stdout)=__stdout__;
    GET_FILE_FROM_PY(stderr)=__stderr__;
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
#ifdef MAPP_DEBUG_MODE
FILE* MAPP_NS::MAPP::mapp_debug(NULL);
#endif
#ifdef IS_PY3K
PyModuleDef MAPP::module=EmptyModule;
#endif
MOD_INIT(mapp4py,MAPP_NS::MAPP::init_module())
/*--------------------------------------------*/
PyObject* MAPP::init_module(void)
{
    int mpi_initialized;
    MPI_Initialized(&mpi_initialized);
    if(!mpi_initialized)
    {
#ifdef LIBMPI_SONAME
        int mode = RTLD_NOW | RTLD_GLOBAL;
#ifdef RTLD_NOLOAD
        mode |= RTLD_NOLOAD;
#endif
        dlopen(LIBMPI_SONAME,mode);
#endif
        MPI_Init(NULL,NULL);
    }
    
#ifdef MAPP_DEBUG_MODE
#include <unistd.h>
#include <iostream>
    char debug_file [100];
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    sprintf(debug_file, "Debug/debug-%d",rank);
    mapp_debug=fopen(debug_file,"w");
    fprintf(mapp_debug,"my pid is %d\n",getpid());
#endif
    Py_AtExit([]()
    {
#ifdef MAPP_DEBUG_MODE
        fclose(mapp_debug);
#endif
        int mpi_finalized;
        MPI_Finalized(&mpi_finalized);
        if(!mpi_finalized) MPI_Finalize();
    });
    
    PyObject* posixpath=PyImport_ImportModule("posixpath");
    PyObject* devnull_path_op=PyObject_GetAttrString(posixpath,"devnull");
    devnull_path=__PyString_AsString(devnull_path_op);
    Py_DECREF(devnull_path_op);
    Py_DECREF(posixpath);
    glbl_rank=Communication::get_rank();
    
    mapp_out=__stdout__=GET_FILE_FROM_PY(stdout);
    mapp_err=__stderr__=GET_FILE_FROM_PY(stderr);
    
    import_array2("mapp4py was unable to import numpy.core.multiarray", NULL)
    setup_methods();
#ifdef IS_PY3K
    module.m_name="mapp4py";
    module.m_doc="MIT Atomistic Parallel Package";
    module.m_methods=methods;
    PyObject* module_ob=PyModule_Create(&module);
#else
    PyObject* module_ob=Py_InitModule3("mapp4py",methods,"MIT Atomistic Parallel Package");
#endif

    
    if(module_ob==NULL) return NULL;
    
    MAPP_MPI::setup_tp();
    if(PyType_Ready(&MAPP_MPI::TypeObject)<0) return NULL;
    Py_INCREF(&MAPP_MPI::TypeObject);
    PyModule_AddObject(module_ob,"mpi",reinterpret_cast<PyObject*>(&MAPP_MPI::TypeObject));
    
    
    if(LineSearch::setup_tp()<0) return NULL;
    PyModule_AddObject(module_ob,"ls",reinterpret_cast<PyObject*>(&LineSearch::TypeObject));
    
    
    if(LineSearchGoldenSection::setup_tp()<0) return NULL;
    PyModule_AddObject(module_ob,"ls_golden",reinterpret_cast<PyObject*>(&LineSearchGoldenSection::TypeObject));
    
    
    if(LineSearchBrent::setup_tp()<0) return NULL;
    PyModule_AddObject(module_ob,"ls_brent",reinterpret_cast<PyObject*>(&LineSearchBrent::TypeObject));
    
    
    if(LineSearchBackTrack::setup_tp()<0) return NULL;
    PyModule_AddObject(module_ob,"ls_bt",reinterpret_cast<PyObject*>(&LineSearchBackTrack::TypeObject));

#ifdef POTFIT   
    if(PotFit<ForceFieldEAMPotFitAckOgata,2>::setup_tp()<0) return NULL;
    PyModule_AddObject(module_ob,"potfit",reinterpret_cast<PyObject*>(&PotFit<ForceFieldEAMPotFitAckOgata,2>::TypeObject));
    
    if(PotFit<ForceFieldEAMPotFitAckJP,2>::setup_tp()<0) return NULL;
    PyModule_AddObject(module_ob,"potfit_jp",reinterpret_cast<PyObject*>(&PotFit<ForceFieldEAMPotFitAckJP,2>::TypeObject));
#endif
    
    PyObject* md=MAPP::MD::init_module();
    if(md==NULL) return NULL;
    PyModule_AddObject(module_ob,"md",md);
    
    
    PyObject* dmd=MAPP::DMD::init_module();
    if(dmd==NULL) return NULL;
    PyModule_AddObject(module_ob,"dmd",dmd);
    
    pause_out(NULL);
#ifdef IS_PY3K
    return module_ob;
#else
    return NULL;
#endif
}
/*--------------------------------------------*/
PyMethodDef MAPP::MD::methods[]=EmptyPyMethodDef(1);
/*--------------------------------------------*/
void MAPP::MD::setup_methods()
{
    /*
    ExamplePython::ml_phonon(methods[0]);
    ExamplePython::ml_phonon_1d(methods[1]);
    ExamplePython::ml_phonon_1dd(methods[2]);*/
}
/*--------------------------------------------*/
#ifdef IS_PY3K
PyModuleDef MAPP::MD::module=EmptyModule;
#endif
/*--------------------------------------------*/
PyObject* MAPP::MD::init_module(void)
{
    setup_methods();
#ifdef IS_PY3K
    module.m_name="mapp4py.md";
    module.m_doc="Molecular Dynamics (MD) module";
    module.m_methods=methods;
    PyObject* module_ob=PyModule_Create(&module);
#else
    PyObject* module_ob=Py_InitModule3("mapp4py.md",methods,"Molecular Dynamics (MD) module");
#endif
    if(module_ob==NULL) return NULL;
    
    
    if(AtomsMD::setup_tp()<0) return NULL;
    PyModule_AddObject(module_ob,"atoms",reinterpret_cast<PyObject*>(&AtomsMD::TypeObject));
    
    if(MDNVT::setup_tp()<0) return NULL;
    PyModule_AddObject(module_ob,"nvt",reinterpret_cast<PyObject*>(&MDNVT::TypeObject));
    
    
    if(MDNST::setup_tp()<0) return NULL;
    PyModule_AddObject(module_ob,"nst",reinterpret_cast<PyObject*>(&MDNST::TypeObject));
    
    
    if(MDMuVT::setup_tp()<0) return NULL;
    PyModule_AddObject(module_ob,"muvt",reinterpret_cast<PyObject*>(&MDMuVT::TypeObject));
    
    
    if(MinCGOld::setup_tp()<0) return NULL;
    PyModule_AddObject(module_ob,"min_cg",reinterpret_cast<PyObject*>(&MinCGOld::TypeObject));
    
    if(MinCG::setup_tp()<0) return NULL;
    PyModule_AddObject(module_ob,"min_cg_new",reinterpret_cast<PyObject*>(&MinCG::TypeObject));
    
    
    if(MinLBFGS::setup_tp()<0) return NULL;
    PyModule_AddObject(module_ob,"min_lbfgs",reinterpret_cast<PyObject*>(&MinLBFGS::TypeObject));
    
    
    if(ExportMD::setup_tp()<0) return NULL;
    PyModule_AddObject(module_ob,"export",reinterpret_cast<PyObject*>(&ExportMD::TypeObject));
    
    
    if(ExportCFGMD::setup_tp()<0) return NULL;
    PyModule_AddObject(module_ob,"export_cfg",reinterpret_cast<PyObject*>(&ExportCFGMD::TypeObject));
    
    
    return module_ob;
}
/*--------------------------------------------*/
PyMethodDef MAPP::DMD::methods[]=EmptyPyMethodDef(4);
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
/*--------------------------------------------*/
#ifdef IS_PY3K
PyModuleDef MAPP::DMD::module=EmptyModule;
#endif
/*--------------------------------------------*/
PyObject* MAPP::DMD::init_module(void)
{
    setup_methods();
#ifdef IS_PY3K
    module.m_name="mapp4py.dmd";
    module.m_doc="Diffusive Molecular Dynamics (MD) module";
    module.m_methods=methods;
    PyObject* module_ob=PyModule_Create(&module);
#else
    PyObject* module_ob=Py_InitModule3("mapp4py.dmd",methods,"Diffusive Molecular Dynamics (DMD) module");
#endif
    if(module_ob==NULL) return NULL;
    
    
    if(AtomsDMD::setup_tp()<0) return NULL;
    PyModule_AddObject(module_ob,"atoms",reinterpret_cast<PyObject*>(&AtomsDMD::TypeObject));
    
    
    if(MinCGDMDOld::setup_tp()<0) return NULL;
    PyModule_AddObject(module_ob,"min_cg",reinterpret_cast<PyObject*>(&MinCGDMDOld::TypeObject));
    
    if(MinCGDMD::setup_tp()<0) return NULL;
    PyModule_AddObject(module_ob,"min_cg_new",reinterpret_cast<PyObject*>(&MinCGDMD::TypeObject));
    
    if(MinLBFGSDMDOld::setup_tp()<0) return NULL;
    PyModule_AddObject(module_ob,"min_lbfgs",reinterpret_cast<PyObject*>(&MinLBFGSDMDOld::TypeObject));
    
    
    if(DAEBDF::setup_tp()<0) return NULL;
    PyModule_AddObject(module_ob,"bdf",reinterpret_cast<PyObject*>(&DAEBDF::TypeObject));
    
    if(ExportDMD::setup_tp()<0) return NULL;
    PyModule_AddObject(module_ob,"export",reinterpret_cast<PyObject*>(&ExportDMD::TypeObject));
    
    if(ExportCFGDMD::setup_tp()<0) return NULL;
    PyModule_AddObject(module_ob,"export_cfg",reinterpret_cast<PyObject*>(&ExportCFGDMD::TypeObject));
    
    
    return module_ob;
}
