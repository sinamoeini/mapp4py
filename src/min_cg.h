#ifndef __MAPP__min_cg__
#define __MAPP__min_cg__
#include "min.h"
#include "new_min_vec.h"
#include "ff_md.h"
#include "atoms_md.h"
#include "dynamic_md.h"
#include "MAPP.h"
#include "thermo_dynamics.h"
namespace MAPP_NS
{
    
    
    
    namespace MinMDHelper
    {
        template<bool...Bs0>
        class CondB
        {
        public:
            template<class F,class LS>
            static void run(F& f,int n,LS* ls)
            {
                return f.template __run<Bs0...>(ls,n);
            }
            
            template<class F,class LS, class... Bs>
            static void run(F& f,int n,LS* ls,bool b0,Bs...bs)
            {
                if(b0==true)
                    return CondB<Bs0...,true>::run(f,n,ls,bs...);
                else
                    return CondB<Bs0...,false>::run(f,n,ls,bs...);
            }
            
            
            template<class F>
            static void init(F& f)
            {
                return f.template __init<Bs0...>();
            }
            
            template<class F, class... Bs>
            static void init(F& f,bool b0,Bs...bs)
            {
                if(b0==true)
                    return CondB<Bs0...,true>::init(f,bs...);
                else
                    return CondB<Bs0...,false>::init(f,bs...);
            }
            
            template<class F>
            static void fin(F& f)
            {
                return f.template __fin<Bs0...>();
            }
            
            template<class F, class... Bs>
            static void fin(F& f,bool b0,Bs...bs)
            {
                if(b0==true)
                    return CondB<Bs0...,true>::fin(f,bs...);
                else
                    return CondB<Bs0...,false>::fin(f,bs...);
            }
            
        };
        
        
        template<class LS0,class... LSs>
        class CondLS
        {
        public:
            template<class F,class LS00, class... Bs>
            static void run(F& f,int n,LS00* ls,Bs...bs)
            {
                if(dynamic_cast<LS0*>(ls)!=NULL)
                    return CondB<>::run(f,n,dynamic_cast<LS0*>(ls),bs...);
                return CondLS<LSs...>::run(f,n,ls,bs...);
                
            }
        };
        
        template<class LS0>
        class CondLS<LS0>
        {
        public:
            template<class F,class LS00, class... Bs>
            static void run(F& f,int n,LS00* ls,Bs...bs)
            {
                return CondB<>::run(f,n,dynamic_cast<LS0*>(ls),bs...);
            }
        };
        
        
        
        
        
        
        
        
        template<bool X>
        class GetMaxAlphaX
        {
        public:
            template<class F>
            static void func(F& f,type0&,Vec<type0>**)
            {}
        };
        
        
        template<>
        class GetMaxAlphaX<true>
        {
        public:
            template<class F>
            static void func(F& f,type0& __max_a_lcl,Vec<type0>** __vs)
            {
                f.max_alpha_lcl_x(__max_a_lcl,(*__vs)->begin());
            }
        };

        template<bool BC,bool X>
        class Prep
        {
        public:
            template<class F>
            static void func(F&)
            {}
        };
        
        template<>
        class Prep<true,true>
        {
        public:
            template<class F>
            static void func(F& f)
            {
                f.prep_x_d();
            }
            
        };
        
        template<>
        class Prep<true,false>
        {
        public:
            template<class F>
            static void func(F& f)
            {
                f.prep_x_d_affine();
            }
            
        };
        
        template<bool BC>
        class PostForceCalcBC
        {
        public:
            template<class F>
            static void func(F&)
            {}
        };
        
        template<>
        class PostForceCalcBC<true>
        {
        public:
            template<class F>
            static void func(F& f)
            {
                f.post_force_calc_box();
            }
        };
    }
    
    template<bool BC,bool X>
    class MinMDHandler
    {
    protected:
    public:
        static constexpr int N0=VTN<BC||X>::N;
        static constexpr int N1=VTN<X>::N;
        typedef VT<BC,VTN<BC||X>::N> VECTENS0;
        typedef VT<BC,VTN<X>::N> VECTENS1;
        
        class ForceFieldMD* ff;
        class AtomsMD* atoms;
        type0 max_dx;
        bool H_dof[__dim__][__dim__];
        
        
        type0 f_h;
        class NewDynamicMD<BC,BC||X>* dynamic;
        VECTENS0 x,x0,x_d;
        VECTENS1 h,f,f0;
        
        MinMDHandler(AtomsMD*,ForceFieldMD*,type0,bool (&)[__dim__][__dim__]);
        MinMDHandler();
        
        void add_extra_vec(VECTENS1& __v)
        {
            new (&__v) VECTENS1(atoms,__dim__,X);
            Algebra::Do<N1>::func([this,&__v](int i){dynamic->add_xchng(__v.vecs[i]);});
        }
        
        void rm_extra_vec(VECTENS1& __v){__v.~VECTENS1();}
        void init();
        void fin();
        
        
        void max_alpha_lcl_x(type0&,type0*);
        void ls_prep(type0&,type0&,type0&);
        type0 F(type0);
        void F_reset();
        type0 force_calc();
        void prep();
        void prep_x_d();
        void prep_x_d_affine();
        void post_force_calc_box();
    };
}

using namespace MAPP_NS;
/*--------------------------------------------
 
 --------------------------------------------*/
template<bool BC,bool X>
MinMDHandler<BC,X>::MinMDHandler(AtomsMD* __atoms,
ForceFieldMD* __ff,type0 __max_dx,
bool (&__H_dof)[__dim__][__dim__]):
ff(__ff),
atoms(__atoms),
max_dx(__max_dx),
dynamic(NULL)
{
    Algebra::V_eq<__dim__*__dim__>(__H_dof[0],H_dof[0]);
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<bool BC,bool X>
MinMDHandler<BC,X>::MinMDHandler():
atoms(NULL),
ff(NULL),
max_dx(0.0),
dynamic(NULL)
{
    Algebra::set<__dim__*__dim__>(H_dof[0],false);
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<bool BC,bool X>
void MinMDHandler<BC,X>::init()
{
    
    f.~VECTENS1();
    new (&f) VECTENS1(atoms,ff->F_H,ff->f,X);
    
    f0.~VECTENS1();
    new (&f0) VECTENS1(atoms,__dim__,X);
    h.~VECTENS1();
    new (&h) VECTENS1(atoms,__dim__,X);
    
    x.~VECTENS0();
    new (&x) VECTENS0(atoms,atoms->H,atoms->x,BC||X);
    x0.~VECTENS0();
    new (&x0) VECTENS0(atoms,__dim__,BC||X);
    x_d.~VECTENS0();
    new (&x_d) VECTENS0(atoms,h,__dim__,X||BC,X,X&&!BC);
    
    dynamic=new NewDynamicMD<BC,BC||X>(atoms,ff,{},{atoms->x_dof},{});
    
    Algebra::Do<N1>::func([this](int i)
    {
      dynamic->add_xchng(h.vecs[i]);
      dynamic->add_xchng(f0.vecs[i]);
    });
    
    Algebra::Do<N0>::func([this](int i)
    {dynamic->add_xchng(x0.vecs[i]);});
    
    if(BC) dynamic->add_xchng(x_d.vecs[0]);
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<bool BC,bool X>
void MinMDHandler<BC,X>::fin()
{
    delete dynamic;
    dynamic=NULL;
    
    x_d.~VECTENS0();
    x0.~VECTENS0();
    x.~VECTENS0();
    
    h.~VECTENS1();
    f0.~VECTENS1();
    
    
    f.~VECTENS1();
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<bool BC,bool X>
void MinMDHandler<BC,X>::max_alpha_lcl_x
(type0& max_a_lcl,type0* x_d_xvec)
{
    type0 max_x_d_lcl=0.0;
    const int natms_lcl=atoms->natms_lcl;
    int n=natms_lcl*__dim__;
    for(int i=0;i<n;i++)
        max_x_d_lcl=MAX(max_x_d_lcl,fabs(x_d_xvec[i]));
    
    max_a_lcl=MIN(fabs(max_dx/max_x_d_lcl),max_a_lcl);
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<bool BC,bool X>
void MinMDHandler<BC,X>::ls_prep(type0& dfa,type0& h_norm,type0& max_a)
{
    h_norm=h*h;
    
    dfa=-f_h;
    
    if(h_norm==0.0)
    {
        max_a=0.0;
        dfa=0.0;
        return;
    }
    
    if(dfa>=0.0)
    {
        max_a=0.0;
        dfa=1.0;
        return;
    }
    
    h_norm=sqrt(h_norm);
    
    type0 max_a_lcl=std::numeric_limits<type0>::infinity();
    MinMDHelper::GetMaxAlphaX<BC||X>::func(*this,max_a_lcl,x_d.vecs);
    MPI_Allreduce(&max_a_lcl,&max_a,1,Vec<type0>::MPI_T,MPI_MIN,atoms->world);
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<bool BC,bool X>
type0 MinMDHandler<BC,X>::F(type0 __alpha)
{
    
    x=x0+__alpha*x_d;
    dynamic->update();
    return ff->value();
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<bool BC,bool X>
void MinMDHandler<BC,X>::F_reset()
{
    x=x0;
    dynamic->update();
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<bool BC,bool X>
type0 MinMDHandler<BC,X>::force_calc()
{
    type0 ans=*ff->derivative();
    MinMDHelper::PostForceCalcBC<BC>::func(*this);
    return ans;
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<bool BC,bool X>
void MinMDHandler<BC,X>::prep()
{
    MinMDHelper::Prep<BC,X>::func(*this);
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<bool BC,bool X>
void MinMDHandler<BC,X>::prep_x_d()
{
    const int natms_lcl=atoms->natms_lcl;
    type0* xvec=x0.vecs[0]->begin();
    type0* hvec=h.vecs[0]->begin();
    type0* x_dvec=x_d.vecs[0]->begin();
    Algebra::MLT_mul_MLT(atoms->H,h.A,x_d.A);
    memcpy(x_dvec,hvec,natms_lcl*__dim__*sizeof(type0));
    for(int iatm=0;iatm<natms_lcl;iatm++,xvec+=__dim__,x_dvec+=__dim__)
        Algebra::V_mul_MLT_add_in(xvec,h.A,x_dvec);
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<bool BC,bool X>
void MinMDHandler<BC,X>::prep_x_d_affine()
{
    const int natms_lcl=atoms->natms_lcl;
    type0* xvec=x0.vecs[0]->begin();
    type0* x_dvec=x_d.vecs[0]->begin();
    Algebra::MLT_mul_MLT(atoms->H,h.A,x_d.A);
    for(int iatm=0;iatm<natms_lcl;iatm++,xvec+=__dim__,x_dvec+=__dim__)
        Algebra::V_mul_MLT(xvec,h.A,x_dvec);
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<bool BC,bool X>
void MinMDHandler<BC,X>::post_force_calc_box()
{
    Algebra::DoLT<__dim__>::func([this](int i,int j)
    {f.A[i][j]=H_dof[i][j] ? f.A[i][j]:0.0;});
}
/*------------------------------------------------------------------------------------------------*/
namespace MAPP_NS
{
    class MinCG:public Min
    {
    private:
    protected:
        bool B_DOF,X_DOF;
        void* handler_ptr;

        class AtomsMD* atoms;
        class ForceFieldMD* ff;
        class ExportMD* xprt;
        void pre_run_chk(AtomsMD*,ForceFieldMD*);
    public:
        MinCG(type0,bool(&)[__dim__][__dim__],bool,type0,class LineSearch*);
        ~MinCG();
        virtual void run(int);
        virtual void init();
        virtual void fin();

        template<bool BC,bool X>
        void  __init();
        template<bool BC,bool X,class LS>
        void  __run(LS*,int);
        template<bool BC,bool X>
        void  __fin();
        
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
        
        static PyGetSetDef getset[];
        static void setup_tp_getset();
        static void getset_x_dof(PyGetSetDef&);
        static void getset_export(PyGetSetDef&);
        
        static int setup_tp();
        
        
        
    };

}
/*--------------------------------------------
 
 --------------------------------------------*/
template<bool BC,bool X>
void  MinCG::__init()
{
    
    MinMDHandler<BC,X>* __handler_ptr=new MinMDHandler<BC,X>(atoms,ff,max_dx,H_dof);
    handler_ptr=__handler_ptr;
    
    __handler_ptr->init();
    __handler_ptr->dynamic->init();
    
    
    
    if(xprt)
    {
        try
        {
            xprt->atoms=atoms;
            xprt->init();
        }
        catch(std::string& err_msg)
        {
            __fin<BC,X>();
            throw err_msg;
        }
    }
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<bool BC,bool X,class LS>
void  MinCG::__run(LS* ls,int nsteps)
{
    
    MinMDHandler<BC,X> &handler=*reinterpret_cast<MinMDHandler<BC,X>*>(handler_ptr);

    typename MinMDHandler<BC,X>::VECTENS1& f=handler.f;
    typename MinMDHandler<BC,X>::VECTENS1& f0=handler.f0;
    typename MinMDHandler<BC,X>::VECTENS1& h=handler.h;
    typename MinMDHandler<BC,X>::VECTENS0& x0=handler.x0;
    typename MinMDHandler<BC,X>::VECTENS0& x=handler.x;
    type0& f_h=handler.f_h;
    
    ls->reset();
    int step=atoms->step;
    
    handler.force_calc();
    
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
        handler.prep();
        err=ls->min(&handler,e_curr,alpha,1);
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
/*--------------------------------------------
 
 --------------------------------------------*/
template<bool BC,bool X>
void  MinCG::__fin()
{
    if(xprt)
    {
        xprt->fin();
        xprt->atoms=NULL;
    }
    MinMDHandler<BC,X>* __handler_ptr=reinterpret_cast<MinMDHandler<BC,X>*>(handler_ptr);
    __handler_ptr->dynamic->fin();
    __handler_ptr->fin();
    delete __handler_ptr;
    handler_ptr=NULL;
}
#endif
