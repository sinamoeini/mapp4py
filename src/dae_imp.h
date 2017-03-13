#ifndef __MAPP__dae_imp__
#define __MAPP__dae_imp__
#include "dae.h"
namespace MAPP_NS
{
    class DAEImplicit:public DAE
    {
    private:
    protected:
        bool newton();
        type0 update_c();
        void J_test();
    public:
        int max_nnewton_iters;
        int max_ngmres_iters;
        
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
        
        
        
        class GMRES* gmres;
        
        DAEImplicit();
        virtual ~DAEImplicit();
        void init_static();
        void fin_static();

        
        typedef struct
        {
            PyObject_HEAD
            DAEImplicit* dae;
            ExportDMD::Object* xprt;
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
        static void getset_max_nnewton_iters(PyGetSetDef&);
        static void getset_max_ngmres_iters(PyGetSetDef&);
        
        static void setup_tp();
    };
}


#endif 
