#ifdef Min_Style
    MinStyle(MinLBFGS,l-bfgs)
#else
#ifndef __MAPP__MinLBFGS__
#define __MAPP__MinLBFGS__
#include "min_cg.h"
namespace MAPP_NS
{
    class MinLBFGS:public MinCG
    {
    private:
    protected:
        int m;
        type0* rho;
        type0* alpha;
        
        VecTens<type0,1>* s;
        VecTens<type0,1>* y;
    public:
        MinLBFGS(AtomsMD*&,ForceFieldMD*&,int);
        ~MinLBFGS();
        void run(int);
        void init();
        void fin();
        
        typedef struct
        {
            PyObject_HEAD
            MinLBFGS* min;
            AtomsMD::Object* atoms;
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
        static void getset_m(PyGetSetDef&);
        
        static void setup_tp();
    };
}

#endif
#endif
