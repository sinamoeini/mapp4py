#ifndef __MAPP__min_cg_dmd__
#define __MAPP__min_cg_dmd__
#include "min.h"
#include "min_vec.h"
#include "ff_dmd.h"
#include "atoms_dmd.h"
#include "dynamic_dmd.h"
#include "MAPP.h"
#include "thermo_dynamics.h"
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
        
        class AtomsDMD* atoms;
        class ForceFieldDMD* ff;
        class DynamicDMD* dynamic;
        class ExportDMD* xprt;
        void print_error();
        
        vec* uvecs[2];
        
    public:
        MinCGDMD();
        MinCGDMD(type0,bool(&)[__dim__][__dim__],bool,type0,type0,class LineSearch*);
        ~MinCGDMD();
        virtual void run(int);
        template<class C>
        void run(C*,int);
        void init();
        void fin();
        
        void force_calc();
        
        type0 F(type0);
        type0 dF(type0,type0&);
        void ls_prep(type0&,type0&,type0&);
        void F_reset();

        
        typedef struct
        {
            PyObject_HEAD
            MinCGDMD* min;
            LineSearch::Object* ls;
            ExportDMD::Object* xprt;
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
        static void getset_export(PyGetSetDef&);
        
        static void setup_tp();
        
        
        
    };

}
using namespace MAPP_NS;
/*--------------------------------------------
 
 --------------------------------------------*/
template<class C>
void MinCGDMD::run(C* ls,int nsteps)
{
    force_calc();
    
    int nevery_xprt=xprt==NULL ? 0:xprt->nevery;
    if(nevery_xprt) xprt->write(0);
    
    type0 S[__dim__][__dim__];
    ThermoDynamics thermo(6,
    "PE",ff->nrgy_strss[0],
    "S[0][0]",S[0][0],
    "S[1][1]",S[1][1],
    "S[2][2]",S[2][2],
    "S[1][2]",S[2][1],
    "S[2][0]",S[2][0],
    "S[0][1]",S[1][0]);
    
    thermo.init();
    Algebra::DyadicV_2_MLT(&ff->nrgy_strss[1],S);
    thermo.print(0);
    
    
    type0 e_prev,e_curr=ff->nrgy_strss[0];
    type0 f0_f0,f_f,f_f0;
    type0 ratio,alpha;
    int err=nsteps==0? MIN_F_MAX_ITER:LS_S;
    h=f;
    f0_f0=f*f;
    int istep=0;
    for(;istep<nsteps && err==LS_S;istep++)
    {
        if(f0_f0==0.0)
        {
            err=LS_F_GRAD0;
            continue;
        }
        
        x0=x;
        f0=f;
        
        e_prev=e_curr;
        
        
        f_h=f*h;
        if(f_h<0.0)
        {
            h=f;
            f_h=f0_f0;
        }
        if(affine) prepare_affine_h();
        err=ls->min(this,e_curr,alpha,1);
        
        if(err!=LS_S)
            continue;
        
        force_calc();
        
        if(e_prev-e_curr<e_tol)
            err=MIN_S_TOLERANCE;
        
        if(istep+1==nsteps)
            err=MIN_F_MAX_ITER;
        
        if((istep+1)%ntally==0)
        {
            Algebra::DyadicV_2_MLT(&ff->nrgy_strss[1],S);
            thermo.print(istep+1);
        }
        
        if(nevery_xprt && (istep+1)%nevery_xprt==0) xprt->write(istep+1);
        
        if(err) continue;
        
        f_f=f*f;
        f_f0=f*f0;
        
        ratio=(f_f-f_f0)/(f0_f0);
        
        h=ratio*h+f;
        f0_f0=f_f;
    }
    
    if(istep%ntally)
    {
        Algebra::DyadicV_2_MLT(&ff->nrgy_strss[1],S);
        thermo.print(istep);
    }
    
    if(nevery_xprt && istep%nevery_xprt) xprt->write(istep);

    thermo.fin();
    fprintf(MAPP::mapp_out,"%s",err_msgs[err]);
}
#endif
