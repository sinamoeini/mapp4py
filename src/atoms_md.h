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
        type0 S_pe[__dim__][__dim__];
        type0 pe;
        
        Vec<type0>* x_d;
        Vec<elem_type>* elem;
        void create_T(type0,int);
        int count_elem(elem_type);
        void DO(PyObject*);
        AtomsMD(MPI_Comm&);
        ~AtomsMD();
        AtomsMD& operator=(const Atoms&);
        void x_d2s_d_dump();
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
        static void getset_S_pe(PyGetSetDef&);
        static void getset_pe(PyGetSetDef&);
        
        static PyMethodDef methods[];
        static void setup_tp_methods();
        static void ml_create_temp(PyMethodDef&);
        static void ml_add_elem(PyMethodDef&);
        static void ml_import_cfg(PyMethodDef&);
        //static void ml_xxx(PyMethodDef&);
        
        static int setup_tp();
    };
}
#endif
