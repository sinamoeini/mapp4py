#ifndef __MAPP__min_cg_dmd__
#define __MAPP__min_cg_dmd__
#include "min_vec.h"
#include "api.h"
#include "ls.h"
#include "atoms_dmd.h"
namespace MAPP_NS
{
    class MinCGDMD
    {
    private:
    protected:
        type0 max_dx;
        type0 max_dalpha;
        type0 MLT[__dim__][__dim__];
        bool chng_box;
        int err;
        type0 f_h;
        type0 curr_energy;
        
        void prepare_affine_h();
        type0 calc_ndofs();
        type0 ndofs;
        
        
        bool affine;
        type0 e_tol;
        
        bool H_dof[__dim__][__dim__];
        
        int ntally;
        
        VecTens<type0,2> h;
        VecTens<type0,2> x;
        VecTens<type0,2> x0;
        VecTens<type0,2> f;
        VecTens<type0,2> f0;
        
        MPI_Comm& world;
        class AtomsDMD*& atoms;
        class ForceFieldDMD*& ff;
        class DynamicDMD* dynamic;
        void print_error();
        
        vec* uvecs[2];
    public:
        MinCGDMD(AtomsDMD*&,ForceFieldDMD*&);
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
        static void getset_e_tol(PyGetSetDef&);
        static void getset_affine(PyGetSetDef&);
        static void getset_H_dof(PyGetSetDef&);
        static void getset_max_dx(PyGetSetDef&);
        static void getset_max_dalpha(PyGetSetDef&);
        static void getset_ntally(PyGetSetDef&);
        
        static void setup_tp();
        
        
        
    };

}
#endif
