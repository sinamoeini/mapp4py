#ifndef __MAPP__atoms_md__
#define __MAPP__atoms_md__
#include "atoms.h"
namespace MAPP_NS
{
    class AtomsMD: public Atoms
    {
    private:
    protected:
    public:
        Vec<type0>* x_d;
        Vec<elem_type>* elem;
        void create_T(type0,int);
        //void DO(PyFunctionObject*);
        AtomsMD(MPI_Comm&);
        AtomsMD(Communication&);
        ~AtomsMD();
        
        
        typedef struct
        {
            PyObject_HEAD
            class AtomsMD* atoms;
            class ForceFieldMD* ff;
        }Object;
        
        static PyTypeObject TypeObject;
        static PyObject* __new__(PyTypeObject*,PyObject*, PyObject*);
        static int __init__(PyObject*, PyObject*,PyObject*);
        static PyObject* __alloc__(PyTypeObject*,Py_ssize_t);
        static void __dealloc__(PyObject*);
        
        static PyGetSetDef getset[];
        static void setup_tp_getset();
        
        static PyMethodDef methods[];
        static void setup_tp_methods();
        static void ml_create_T(PyMethodDef&);
        static void ml_import_cfg(PyMethodDef&);
        
        static void setup_tp();
    };
}
#endif
