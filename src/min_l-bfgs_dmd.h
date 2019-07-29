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
        void* s_ptr;
        void* y_ptr;
        
    public:
        MinLBFGSDMD(int,type0,bool(&)[__dim__][__dim__],bool,type0,type0,class LineSearch*);
        ~MinLBFGSDMD();
        void run(int);
        void init();
        void fin();
        
        template<bool BC,bool X,bool ALPHA,bool C>
        void  __init();
        template<bool BC,bool X,bool ALPHA,bool C,class LS>
        void  __run(LS*,int);
        template<bool BC,bool X,bool ALPHA,bool C>
        void  __fin();
        
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
        
        static int setup_tp();
    };
}
using namespace MAPP_NS;
/*--------------------------------------------
 
 --------------------------------------------*/
template<bool BC,bool X,bool ALPHA,bool C>
void MinLBFGSDMD::__init()
{
    
    MinDMDHandler<BC,X,ALPHA,C>* __handler_ptr=new (&handler_buff) MinDMDHandler<BC,X,ALPHA,C>(atoms,ff,max_dx,max_dalpha,H_dof);
    __handler_ptr->init();
    typedef typename MinDMDHandler<BC,X,ALPHA,C>::VECTENS1 VECTENS1;
    VECTENS1* __y_ptr=m==0 ? NULL:new VECTENS1[m];
    VECTENS1* __s_ptr=m==0 ? NULL:new VECTENS1[m];
    int c_dim=atoms->c_dim;
    for(int i=0;i<m;i++)
    {
        __y_ptr[i].~VECTENS1();
        new (__y_ptr+i) VECTENS1(atoms,__dim__,X,c_dim,ALPHA,c_dim,C);
        
        __s_ptr[i].~VECTENS1();
        new (__s_ptr+i) VECTENS1(atoms,__dim__,X,c_dim,ALPHA,c_dim,C);
    }
    
    s_ptr=__s_ptr;
    y_ptr=__y_ptr;
    
    static constexpr int N1=MinDMDHandler<BC,X,ALPHA,C>::N1;
    for(int i=0;i<m;i++)
    {
        Algebra::Do<N1>::func([&__s_ptr,&__y_ptr,&i,&__handler_ptr](int j)
        {
            __handler_ptr->dynamic.add_xchng(__y_ptr[i].vecs[j]);
            __handler_ptr->dynamic.add_xchng(__s_ptr[i].vecs[j]);
            
        });
    }
    __handler_ptr->dynamic.init();
    
    rho=m==0 ?NULL:new type0[m];
    alpha=m==0 ?NULL:new type0[m];
    
    if(xprt)
    {
        try
        {
            xprt->atoms=atoms;
            xprt->init();
        }
        catch(std::string& err_msg)
        {
            __fin<BC,X,ALPHA,C>();
            throw err_msg;
        }
    }
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<bool BC,bool X,bool ALPHA,bool C,class LS>
void MinLBFGSDMD::__run(LS* ls,int nsteps)
{
    
    MinDMDHandler<BC,X,ALPHA,C>& handler=*reinterpret_cast<MinDMDHandler<BC,X,ALPHA,C>*>(&handler_buff);
    typedef typename MinDMDHandler<BC,X,ALPHA,C>::VECTENS0 VECTENS0;
    typedef typename MinDMDHandler<BC,X,ALPHA,C>::VECTENS1 VECTENS1;
    VECTENS1& f=handler.f;
    VECTENS1& f0=handler.f0;
    VECTENS1& h=handler.h;
    VECTENS0& x0=handler.x0;
    VECTENS0& x=handler.x;
    VECTENS1* y=reinterpret_cast<VECTENS1*>(y_ptr);
    VECTENS1* s=reinterpret_cast<VECTENS1*>(s_ptr);
    type0& f_h=handler.f_h;
    
    ls->reset();
    int step=atoms->step;
    handler.force_calc();
    
    int nevery_xprt=xprt==NULL ? 0:xprt->nevery;
    if(nevery_xprt) xprt->write(step);
    
    ThermoDynamics thermo(6,
      "FE",atoms->fe,
      "S[0][0]",atoms->S_fe[0][0],
      "S[1][1]",atoms->S_fe[1][1],
      "S[2][2]",atoms->S_fe[2][2],
      "S[1][2]",atoms->S_fe[2][1],
      "S[2][0]",atoms->S_fe[2][0],
      "S[0][1]",atoms->S_fe[1][0]);
    
    if(ntally) thermo.init();
    if(ntally) thermo.print(step);
    
    
    type0 e_prev,e_curr=atoms->pe;
    type0 alpha_m,gamma;
    type0 inner0,inner1,scl0;
    int k=0;
    gamma=1.0;
    int err=nsteps==0? MIN_F_MAX_ITER:LS_S;
    int istep=0;
    for(;err==LS_S;istep++)
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
        {
            scl0=-alpha[i]-rho[i]*(y[i]*h);
            h+=scl0*s[i];
        }
        
        e_prev=e_curr;
        
        
        f_h=f*h;
        if(f_h<0.0)
        {
            h=f;
            k=0;
            f_h=f*f;
        }
        handler.prep();
        err=ls->min(&handler,e_curr,alpha_m,0);
        if(err!=LS_S)
        {
            // this was a bullshit step so we have to decrease the setup once to compensate
            // the last step was the previous one
            istep--;
            continue;
        }
        
        handler.force_calc();
        
        if(ntally && (istep+1)%ntally==0)
            thermo.print(step+istep+1);
        
        if(nevery_xprt && (istep+1)%nevery_xprt==0) xprt->write(step+istep+1);
        
        //this was a successfull step but the last one
        if(e_prev-e_curr<e_tol) err=MIN_S_TOLERANCE;
        if(istep+1==nsteps) err=MIN_F_MAX_ITER;
        if(err) continue;
        
        if(m)
        {
            if(k!=m) k++;
            
            cyclic_shift(s,k);
            cyclic_shift(y,k);
            
            for(int i=m-1;i>0;i--)
                rho[i]=rho[i-1];
            
            s[0]=alpha_m*h;
            y[0]=f0-f;
            
            inner0=s[0]*y[0];
            inner1=y[0]*y[0];
            
            gamma=inner0/inner1;
            rho[0]=1.0/inner0;
        }
        else
            gamma=alpha_m*(h*f0-h*f)/(f*f+f0*f0-2.0*(f*f0));
    }
    
    if(ntally && istep%ntally)
        thermo.print(step+istep);
    
    if(nevery_xprt && istep%nevery_xprt) xprt->write(step+istep);
    
    if(ntally) thermo.fin();
    
    if(ntally) fprintf(MAPP::mapp_out,"%s",err_msgs[err]);
    
    atoms->step+=istep;
    
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<bool BC,bool X,bool ALPHA,bool C>
void MinLBFGSDMD::__fin()
{
    delete [] alpha;
    alpha=NULL;
    delete [] rho;
    rho=NULL;
    
    if(xprt)
    {
        xprt->fin();
        xprt->atoms=NULL;
    }
    MinDMDHandler<BC,X,ALPHA,C>* __handler_ptr=reinterpret_cast<MinDMDHandler<BC,X,ALPHA,C>*>(&handler_buff);
    __handler_ptr->dynamic.fin();
    typedef typename MinDMDHandler<BC,X,ALPHA,C>::VECTENS1 VECTENS1;
    VECTENS1* __y_ptr=reinterpret_cast<VECTENS1*>(y_ptr);
    VECTENS1* __s_ptr=reinterpret_cast<VECTENS1*>(s_ptr);
    delete [] __s_ptr;
    delete [] __y_ptr;
    
    s_ptr=NULL;
    y_ptr=NULL;
    
    __handler_ptr->fin();
    __handler_ptr->~MinDMDHandler();    
}

#endif
