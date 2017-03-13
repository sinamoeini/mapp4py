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
        MinLBFGSDMD(int,type0,bool(&)[__dim__][__dim__],bool,type0,type0,class LineSearch*);
        ~MinLBFGSDMD();
        void run(int);
        template<class C>
        void run(C*,int);
        void init();
        void fin();
        
        typedef struct
        {
            PyObject_HEAD
            MinLBFGSDMD* min;
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
        
        static PyGetSetDef getset[];
        static void setup_tp_getset();
        static void getset_m(PyGetSetDef&);
        
        static void setup_tp();
    };
}
using namespace MAPP_NS;
/*--------------------------------------------
 
 --------------------------------------------*/
template<class C>
void MinLBFGSDMD::run(C* ls,int nsteps)
{
    int step=atoms->step;
    
    force_calc();
    
    int nevery_xprt=xprt==NULL ? 0:xprt->nevery;
    if(nevery_xprt) xprt->write(step);
    
    type0 S[__dim__][__dim__];
    ThermoDynamics thermo(6,
    "PE",ff->nrgy_strss[0],
    "S[0][0]",S[0][0],
    "S[1][1]",S[1][1],
    "S[2][2]",S[2][2],
    "S[1][2]",S[2][1],
    "S[2][0]",S[2][0],
    "S[0][1]",S[1][0]);
    
    
    if(ntally) thermo.init();
    Algebra::DyadicV_2_MLT(&ff->nrgy_strss[1],S);
    if(ntally) thermo.print(step);
    
    
    type0 e_prev,e_curr=ff->nrgy_strss[0];
        
    type0 alpha_m,gamma;
    type0 inner0,inner1;
    
    
    
    int k=0;
    gamma=1.0;
    int err=nsteps==0? MIN_F_MAX_ITER:LS_S;
    

    int istep=0;
    for(;istep<nsteps && err==LS_S;istep++)
    {
        x0=x;
        h=f0=f;
        
        for(int i=0;i<k;i++)
        {
            alpha[i]=-rho[i]*(s[i]*h);
            h+=alpha[i]*y[i];
        }
        
        h*=gamma;
        
        for(int i=k-1;i>-1;i--)
        h+=(-alpha[i]-rho[i]*(y[i]*h))*s[i];
        
        e_prev=e_curr;
        
        
        f_h=f*h;
        if(f_h<0.0)
        {
            h=f;
            k=0;
            f_h=f*f;
        }
        if(affine) prepare_affine_h();
        err=ls->min(this,e_curr,alpha_m,0);
        
        if(err!=LS_S)
        {
            // this was a bullshit step so we have to decrease the setup once to compensate
            // the last step was the previous one
            istep--;
            continue;
        }
        
        force_calc();
        
        if(ntally && (istep+1)%ntally==0)
        {
            Algebra::DyadicV_2_MLT(&ff->nrgy_strss[1],S);
            thermo.print(step+istep+1);
        }
        
        if(nevery_xprt && (istep+1)%nevery_xprt==0) xprt->write(step+istep+1);
        
        //this was a successfull step but the last one
        if(e_prev-e_curr<e_tol) err=MIN_S_TOLERANCE;
        if(istep+1==nsteps) err=MIN_F_MAX_ITER;
        if(err) continue;
        
        if(m)
        {
            if(k!=m) k++;
            
            s[0].cyclic_shift(k);
            y[0].cyclic_shift(k);
            
            for(int i=m-1;i>0;i--)
                rho[i]=rho[i-1];
            
            s[0]=x-x0;
            y[0]=f0-f;
            
            inner0=s[0]*y[0];
            inner1=y[0]*y[0];
            
            gamma=inner0/inner1;
            rho[0]=1.0/inner0;
        }
        else
            gamma=(x*f0-x*f-x0*f0+x0*f)/(f*f+f0*f0-2.0*(f*f0));
    }
    
    if(ntally && istep%ntally)
    {
        Algebra::DyadicV_2_MLT(&ff->nrgy_strss[1],S);
        thermo.print(step+istep);
    }

    if(nevery_xprt && istep%nevery_xprt) xprt->write(step+istep);
    
    if(ntally) thermo.fin();
    
    fprintf(MAPP::mapp_out,"%s",err_msgs[err]);
    
    atoms->step+=istep;
}
#endif 
