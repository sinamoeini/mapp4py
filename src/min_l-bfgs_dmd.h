#ifndef __MAPP__min_l_bfgs_dmd__
#define __MAPP__min_l_bfgs_dmd__
#include "min_cg_dmd.h"
namespace MAPP_NS
{
    class MinLBFGSDMD:public MinCGDMD
    {
    private:
    protected:
        int m;
        type0* rho;
        type0* alpha;
        
        VecTens<type0,2>* s;
        VecTens<type0,2>* y;
    public:
        MinLBFGSDMD(AtomsDMD*,ForceFieldDMD*,int);
        ~MinLBFGSDMD();
        void run(int);
        void init();
        void fin();
        
        typedef struct
        {
            PyObject_HEAD
            MinLBFGSDMD* min;
            AtomsDMD::Object* atoms;
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
