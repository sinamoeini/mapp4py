#ifndef __MAPP__newton_gmres__
#define __MAPP__newton_gmres__
#include "api.h"
#include "export_styles.h"
namespace MAPP_NS
{
    template<typename> class Vec;
    class NewtonGMRES
    {
    private:
    protected:
        void pre_run_chk(class AtomsDMD*,class ForceFieldDMD*);
        int c_dim;
        bool chng_box;
        bool S_dof[__dim__][__dim__];
        type0 S[__dim__][__dim__];
        type0 a_tol_sqrt_nc_dofs;
        
        class AtomsDMD* atoms;
        class ForceFieldDMD* ff;
        class ExportDMD* xprt;
        class DynamicDMD* dynamic;
    public:
        
        int ntally;
        int max_nsteps;
        type0 a_tol;
        int max_ngmres_iters;
        
        
        NewtonGMRES();
        ~NewtonGMRES();

        void init();
        void fin();
        void run(int);
        
        typedef struct
        {
            PyObject_HEAD
            NewtonGMRES* ngmres;
            ExportDMD::Object* xprt;
        }Object;
        
        
        static PyTypeObject TypeObject;
        static PyObject* __new__(PyTypeObject*,PyObject*, PyObject*);
        static int __init__(PyObject*, PyObject*,PyObject*);
        static PyObject* __alloc__(PyTypeObject*,Py_ssize_t);
        static void __dealloc__(PyObject*);
        
        static PyMethodDef methods[];
        static void setup_tp_methods();
        static void ml_run(PyMethodDef&);
        
        static PyGetSetDef getset[];
        static void setup_tp_getset();
        static void getset_max_ngmres_iters(PyGetSetDef&);
        static void getset_S(PyGetSetDef&);
        static void getset_a_tol(PyGetSetDef&);
        static void getset_ntally(PyGetSetDef&);
        static void getset_export(PyGetSetDef&);
        
        static int setup_tp();
    };
}

#endif
