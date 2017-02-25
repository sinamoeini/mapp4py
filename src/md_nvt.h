#ifndef __MAPP__md_nvt__
#define __MAPP__md_nvt__
#include "api.h"
#include <mpi.h>
#include <stdio.h>
#include <iostream>
#include <limits.h>
#include "atoms_md.h"
#include "thermostat.h"

namespace MAPP_NS
{
    class MDNVT
    {
    private:
    protected:
        class AtomsMD* atoms;
        class ForceFieldMD* ff;
        
        class DynamicMD* dynamic;
        type0 mvv[__nvoigt__];
        type0 __vec[__nvoigt__];
        type0 __vec_lcl[__nvoigt__];
        
        type0 S_part[__dim__][__dim__];
        
        type0 T_part;
        int ndof_part;
        ThermostatNHC thermo_part;
        
        
        type0 kB;
        type0 T;
        type0 dt;
        type0 dt2;
        
        int ntally;
        
        void update_x_d__x__x_d(type0);
        void update_x();
        void update_x_d(type0);
        void update_x_d();
        void update_x_d_final(type0);
        virtual void change_dt(type0);
        virtual void pre_run_chk(AtomsMD*,ForceFieldMD*);
    public:
        MDNVT(type0,type0);
        virtual ~MDNVT();
        void init();
        virtual void run(int);
        void fin();
        
        typedef struct
        {
            PyObject_HEAD
            MDNVT* md;
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
        static void getset_niters(PyGetSetDef&);
        static void getset_nchains(PyGetSetDef&);
        static void getset_dt(PyGetSetDef&);
        static void getset_t_relax(PyGetSetDef&);
        static void getset_T(PyGetSetDef&);
        static void getset_ntally(PyGetSetDef&);
        
        
        static void setup_tp();
         
    };
}


#endif
