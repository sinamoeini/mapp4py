#ifndef __MAPP__dmd_implicit__
#define __MAPP__dmd_implicit__
#include "dmd.h"
namespace MAPP_NS
{
    template<typename,class> class GMRES;
    class DMDImplicit:public DMD
    {
    private:
    protected:
        bool nonlin();
        type0 update_c();
        void J_test();
    public:
        int max_niters_nonlin;
        int max_niters_lin;
        
        int nnonlin_acc;
        int nnonlin_rej;
        int ninteg_acc;
        int ninteg_rej;
        int nintpol_acc;
        int nintpol_rej;
        
        type0 err_fac;
        type0 beta;
        
        type0* y_0;
        type0* c_0;
        type0* a;
        type0* del_c;
        type0* F;
        Vec<type0>* F_ptr;
        Vec<type0>* del_c_ptr;
        
        
        
        class GMRES<type0,ForceFieldDMD>* gmres;
        class __GMRES* __gmres;
        
        DMDImplicit();
        virtual ~DMDImplicit();
        void init_static();
        void fin_static();

        
        typedef struct
        {
            PyObject_HEAD
            DMDImplicit* dmd;
        }Object;
        
        
        static PyTypeObject TypeObject;
        static PyObject* __new__(PyTypeObject*,PyObject*, PyObject*);
        static int __init__(PyObject*, PyObject*,PyObject*);
        static PyObject* __alloc__(PyTypeObject*,Py_ssize_t);
        static void __dealloc__(PyObject*);
        
        static PyMethodDef methods[];
        static void setup_tp_methods();
        
        static PyGetSetDef getset[];
        static void setup_tp_getset();
        static void getset_max_niters_nonlin(PyGetSetDef&);
        static void getset_max_niters_lin(PyGetSetDef&);
        
        static void setup_tp();
    };
}
#endif 
