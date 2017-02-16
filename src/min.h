#ifndef __MAPP__min__
#define __MAPP__min__
#include "global.h"
#include "api.h"
#include "atoms_styles.h"
namespace MAPP_NS
{
    class Min
    {
    private:
    protected:
        static const char* err_msgs[];
        type0 max_dx;
        bool chng_box;
        type0 f_h;
        bool affine;
        type0 e_tol;
        bool H_dof[__dim__][__dim__];
        int ntally;
        
        class Atoms* atoms;
        class ForceField* ff;
    public:
        Min(class Atoms*,class ForceField*);
        virtual ~Min();
        
        typedef struct
        {
            PyObject_HEAD
            Min* min;
            Atoms::Object* atoms;
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
        static void getset_e_tol(PyGetSetDef&);
        static void getset_affine(PyGetSetDef&);
        static void getset_H_dof(PyGetSetDef&);
        static void getset_max_dx(PyGetSetDef&);
        static void getset_ntally(PyGetSetDef&);
        
        static void setup_tp();
    };
}
#endif
