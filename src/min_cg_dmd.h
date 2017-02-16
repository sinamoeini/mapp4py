#ifndef __MAPP__min_cg_dmd__
#define __MAPP__min_cg_dmd__
#include "min.h"
#include "min_vec.h"
#include "ls.h"
namespace MAPP_NS
{
    class MinCGDMD:public Min
    {
    private:
    protected:
        type0 max_dalpha;
        
        void prepare_affine_h();
        type0 calc_ndofs();
        type0 ndofs;
        
        VecTens<type0,2> h;
        VecTens<type0,2> x;
        VecTens<type0,2> x0;
        VecTens<type0,2> f;
        VecTens<type0,2> f0;
        
        MPI_Comm& world;
        class AtomsDMD* atoms;
        class ForceFieldDMD* ff;
        class DynamicDMD* dynamic;
        void print_error();
        
        vec* uvecs[2];
    public:
        MinCGDMD(AtomsDMD*,ForceFieldDMD*);
        ~MinCGDMD();
        virtual void run(int);
        void init();
        void fin();
        
        void force_calc();
        
        type0 F(type0);
        type0 dF(type0,type0&);
        void ls_prep(type0&,type0&,type0&);
        void F_reset();
        
        
        static LineSearch* ls;
        
        
        typedef struct
        {
            PyObject_HEAD
            MinCGDMD* min;
            AtomsDMD::Object* atoms;
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
        static void getset_max_dalpha(PyGetSetDef&);
        
        
        static void setup_tp();
        
        
        
    };

}
#endif
