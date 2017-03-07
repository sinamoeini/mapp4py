#ifndef __MAPP__atoms_dmd__
#define __MAPP__atoms_dmd__
#include "atoms.h"
namespace MAPP_NS
{
    class AtomsDMD: public Atoms
    {
    private:
    protected:
    public:
        Vec<type0>* alpha;
        Vec<type0>* c;
        Vec<type0>* c_d;
        Vec<bool>* c_dof;
        Vec<elem_type>* elem;
        void sort_stack(vec**&,int&,vec**&,int&,vec**&,int&,vec**&,int&);
        AtomsDMD(MPI_Comm&,int);
        AtomsDMD(Communication&,int);
        ~AtomsDMD();
        type0 max_alpha;
        int c_dim;
        const int N;
        type0* xi;
        type0* wi;
        
        typedef struct
        {
            PyObject_HEAD
            class AtomsDMD* atoms;
            class ForceFieldDMD* ff;
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
        
        static void setup_tp();
    };
}
#endif 
