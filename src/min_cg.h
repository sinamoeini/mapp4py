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
        NewDynamicMD<BC,BC||X> dynamic;
        VECTENS0 x,x0,x_d;
        VECTENS1 h,f,f0;
        
        MinMDHandler(AtomsMD*,ForceFieldMD*,type0,bool (&)[__dim__][__dim__]);
        ~MinMDHandler(){};
        void add_extra_vec(VECTENS1& __v)
        {
            new (&__v) VECTENS1(atoms,__dim__,X);
            Algebra::Do<N1>::func([this,&__v](int i){dynamic.add_xchng(__v.vecs[i]);});
        }
        
        template<class V,class...Vs>
        void add_xchng(V*& v,Vs*&... vs)
        {
            dynamic.add_xchng(v);
            add_xchng(vs...);
        }
        
        
        void add_xchng(){}
        
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
dynamic(atoms,ff,{},{atoms->x_dof},{})
{
    Algebra::V_eq<__dim__*__dim__>(__H_dof[0],H_dof[0]);
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
    
    
    
    Algebra::Do<N1>::func([this](int i)
    {
      dynamic.add_xchng(h.vecs[i]);
      dynamic.add_xchng(f0.vecs[i]);
    });
    
    Algebra::Do<N0>::func([this](int i)
    {dynamic.add_xchng(x0.vecs[i]);});
    
    if(BC) dynamic.add_xchng(x_d.vecs[0]);
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<bool BC,bool X>
void MinMDHandler<BC,X>::fin()
{
   
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
    dynamic.update();
    return ff->value();
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<bool BC,bool X>
void MinMDHandler<BC,X>::F_reset()
{
    x=x0;
    dynamic.update();
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
    template<bool>
    friend class MinHelper::ThermoHandler;
    template<bool>
    friend class MinHelper::ExportHandler;
    private:
    protected:
        static constexpr int MaxS=MinHelper::MaxSizeAlign<MinMDHandler,2>::MaxS;
        static constexpr int MaxA=MinHelper::MaxSizeAlign<MinMDHandler,2>::MaxA;
        typedef std::aligned_storage<MaxS,MaxA>::type MemT;
        MemT handler_buff;

        class AtomsMD* atoms;
        class ForceFieldMD* ff;
        class ExportMD* xprt;
        void pre_run_chk(AtomsMD*,ForceFieldMD*);
        void create_thermo();
        void destroy_thermo();
    public:
        MinCG(type0,bool(&)[__dim__][__dim__],bool,type0,class LineSearch*);
        ~MinCG();
        template<class...Vs>
        void init(Vs*&...);
        virtual void run(int);
        virtual void init();
        virtual void fin();


        template<bool BC,bool X,class...Vs>
        void  __init(Vs*&...);
        template<bool BC,bool X,bool OUT,bool XOUT,class LS>
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
template<class...Vs>
void MinCG::init(Vs*&... vs)
{
    chng_box=false;
    Algebra::DoLT<__dim__>::func([this](int i,int j)
    {
        if(H_dof[i][j]) chng_box=true;
    });
    
    try
    {
        MinHelper::CondB<>::init(*this,chng_box,!affine,vs...);
    }
    catch(std::string& err_msg)
    {
        throw err_msg;
    }
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<bool BC,bool X,class...Vs>
void  MinCG::__init(Vs*&... vs)
{
    
    MinMDHandler<BC,X>* __handler_ptr=new (&handler_buff) MinMDHandler<BC,X>(atoms,ff,max_dx,H_dof);
    __handler_ptr->init();
    __handler_ptr->dynamic.init();
    
    
    
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
    
    __handler_ptr->add_xchng(vs...);
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<bool BC,bool X,bool OUT,bool XOUT,class LS>
void  MinCG::__run(LS* ls,int nsteps)
{
    
    MinMDHandler<BC,X> &handler=*reinterpret_cast<MinMDHandler<BC,X>*>(&handler_buff);
    typedef typename MinMDHandler<BC,X>::VECTENS0 VECTENS0;
    typedef typename MinMDHandler<BC,X>::VECTENS1 VECTENS1;
    VECTENS1& f=handler.f;
    VECTENS1& f0=handler.f0;
    VECTENS1& h=handler.h;
    VECTENS0& x0=handler.x0;
    VECTENS0& x=handler.x;
    type0& f_h=handler.f_h;
    
    ls->reset();
    int step=atoms->step;
    handler.force_calc();
    
    MinHelper::ThermoHandler<OUT>::init(*this,step);
    MinHelper::ExportHandler<XOUT>::init(*this,step);
    
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
        
        MinHelper::ThermoHandler<OUT>::print(*this,step,istep);
        MinHelper::ExportHandler<XOUT>::print(*this,step,istep);
        
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
    
    MinHelper::ThermoHandler<OUT>::fin(*this,step,istep,MAPP::mapp_out,err_msgs[err]);
    MinHelper::ExportHandler<XOUT>::fin(*this,step,istep);
    
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
    MinMDHandler<BC,X>* __handler_ptr=reinterpret_cast<MinMDHandler<BC,X>*>(&handler_buff);
    __handler_ptr->dynamic.fin();
    __handler_ptr->fin();
    __handler_ptr->~MinMDHandler();

}
#endif
