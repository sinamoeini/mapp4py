#ifndef __MAPP__dea__
#define __MAPP__dea__
#include "api.h"
#include "export_styles.h"
namespace MAPP_NS
{
    template<typename> class Vec;
    class DAE
    {
    private:
    protected:
        int calc_ndofs(class AtomsDMD*);
        type0 a_tol_sqrt_nc_dofs;
        void pre_run_chk(class AtomsDMD*,class ForceFieldDMD*);
        int c_dim;
        int ncs;
        
        
        bool chng_box;
        bool S_dof[__dim__][__dim__];
        type0 S[__dim__][__dim__];
        
        class AtomsDMD* atoms;
        class ForceFieldDMD* ff;
        class ExportDMD* xprt;
        class DynamicDMD* dynamic;
    public:
        int ntally;
        int nreset;
        int max_nsteps;
        type0 a_tol;
        type0 min_dt;
        
        type0 dt;
        type0 t_cur;
        type0 t_fin;
        type0 err;
        
        type0* c;
        type0* c_d;
        
        
        
        DAE();
        virtual ~DAE();
        virtual void init_static();
        virtual void fin_static();
        virtual void run(type0){};
        void min_error();
        
        
        virtual void init();
        virtual void fin();
        
        
        typedef struct
        {
            PyObject_HEAD
            DAE* dae;
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
        static void getset_S(PyGetSetDef&);
        static void getset_a_tol(PyGetSetDef&);
        static void getset_max_nsteps(PyGetSetDef&);
        static void getset_min_dt(PyGetSetDef&);
        static void getset_nreset(PyGetSetDef&);
        static void getset_ntally(PyGetSetDef&);
        static void getset_export(PyGetSetDef&);
        
        static int setup_tp();
    };
}
#endif 
