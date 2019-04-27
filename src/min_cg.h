#ifndef __MAPP__min_cg__
#define __MAPP__min_cg__
#include "min.h"
#include "min_vec.h"
#include "ff_md.h"
#include "atoms_md.h"
#include "dynamic_md.h"
#include "MAPP.h"
#include "thermo_dynamics.h"
namespace MAPP_NS
{
    class MinCG:public Min
    {
    private:
    protected:
        
        void prep();
        type0 calc_ndofs();
        type0 ndofs;
        
        VecTens<type0,1> h;
        VecTens<type0,1> x;
        VecTens<type0,1> x0;
        VecTens<type0,1> x_d;
        VecTens<type0,1> f;
        VecTens<type0,1> f0;
        
        type0 S_tmp[__dim__][__dim__];
        
        class AtomsMD* atoms;
        class ForceFieldMD* ff;
        class DynamicMD* dynamic;
        class ExportMD* xprt;
        
    public:
        MinCG(type0,bool(&)[__dim__][__dim__],bool,type0,class LineSearch*);
        ~MinCG();
        virtual void run(int);
        template<class C>
        void run(C*,int);
        virtual void init();
#ifdef POTFIT
        template<class FF,size_t NELEMS>
        friend class PotFit;
        void init(vec*);
#endif
        virtual void fin();
        void ff_test(int,type0,type0,int);
        void force_calc();
        
        type0 F(type0);
        type0 dF(type0,type0&);
        void ls_prep(type0&,type0&,type0&);
        void F_reset();

        
        typedef struct
        {
            PyObject_HEAD
            MinCG* min;
            LineSearch::Object* ls;
            ExportMD::Object* xprt;
        }Object;
        
        static PyTypeObject TypeObject;
        static PyObject* __new__(PyTypeObject*,PyObject*, PyObject*);
        static int __init__(PyObject*, PyObject*,PyObject*);
        static PyObject* __alloc__(PyTypeObject*,Py_ssize_t);
        static void __dealloc__(PyObject*);
        
        static PyMethodDef methods[];
        static void setup_tp_methods();
        static void ml_run(PyMethodDef&);
        static void ml_ff_test(PyMethodDef&);
        
        static PyGetSetDef getset[];
        static void setup_tp_getset();
        static void getset_export(PyGetSetDef&);
        
        static int setup_tp();
        
        
        
    };
}
using namespace MAPP_NS;
/*--------------------------------------------
 
 --------------------------------------------*/
template<class C>
void MinCG::run(C* ls,int nsteps)
{
    ls->reset();
    int step=atoms->step;
    
    force_calc();
    
    int nevery_xprt=xprt==NULL ? 0:xprt->nevery;
    if(nevery_xprt) xprt->write(step);
    
    ThermoDynamics thermo(6,
    "PE",atoms->pe,
    "S[0][0]",atoms->S_pe[0][0],
    "S[1][1]",atoms->S_pe[1][1],
    "S[2][2]",atoms->S_pe[2][2],
    "S[1][2]",atoms->S_pe[2][1],
    "S[2][0]",atoms->S_pe[2][0],
    "S[0][1]",atoms->S_pe[1][0]);
    
    if(ntally) thermo.init();
    if(ntally) thermo.print(step);
    
    type0 e_prev,e_curr=atoms->pe;
    type0 f0_f0,f_f,f_f0;
    type0 ratio,alpha;
    int err=nsteps==0? MIN_F_MAX_ITER:LS_S;
    h=f;
    f0_f0=f*f;
    int istep=0;
    for(;err==LS_S;istep++)
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
        prep();
        err=ls->min(this,e_curr,alpha,1);
        
        if(err!=LS_S)
        {
            // this was a bullshit step so we have to decrease the setup once to compensate
            // the last step was the previous one
            istep--;
            continue;
        }
        
        force_calc();
        
        if(ntally && (istep+1)%ntally==0)
            thermo.print(step+istep+1);
        
        if(nevery_xprt && (istep+1)%nevery_xprt==0) xprt->write(step+istep+1);
        
        //this was a successfull step but the last one
        if(e_prev-e_curr<e_tol) err=MIN_S_TOLERANCE;
        if(istep+1==nsteps) err=MIN_F_MAX_ITER;
        if(err) continue;
        
        f_f=f*f;
        f_f0=f*f0;
        
        ratio=(f_f-f_f0)/(f0_f0);
        
        h=ratio*h+f;
        f0_f0=f_f;
    }
    
    if(ntally && istep%ntally)
        thermo.print(step+istep);
    
    
    if(nevery_xprt && istep%nevery_xprt) xprt->write(step+istep);
    

    if(ntally) thermo.fin();
    
    if(ntally) fprintf(MAPP::mapp_out,"%s",err_msgs[err]);
    
    atoms->step+=istep;
}
#endif
