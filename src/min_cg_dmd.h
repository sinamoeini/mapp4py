#ifndef __MAPP__min_cg_dmd__
#define __MAPP__min_cg_dmd__
#include "min.h"
#include "new_min_vec.h"
#include "ff_dmd.h"
#include "atoms_dmd.h"
#include "dynamic_dmd.h"
#include "MAPP.h"
#include "thermo_dynamics.h"
namespace MAPP_NS
{    
    namespace MinDMDHelper
    {
        template<bool C>
        class GetMaxAlphaC
        {
        public:
            template<class F>
            static void func(F&,type0&,Vec<type0>**)
            {
            }
        };
        
        template<>
        class GetMaxAlphaC<true>
        {
        public:
            template<class F>
            static void func(F& f,type0& __max_a_lcl,Vec<type0>** __vs)
            {
                f.max_alpha_lcl_c(__max_a_lcl,(*__vs)->begin());
            }
        };
        
        
        template<bool ALPHA,bool C>
        class GetMaxAlphaAlpha
        {
        public:
            template<class F>
            static void func(F& f,type0& __max_a_lcl,Vec<type0>** __vs)
            {
                GetMaxAlphaC<C>::func(f,__max_a_lcl,__vs);
            }
        };
        
        template<bool C>
        class GetMaxAlphaAlpha<true,C>
        {
        public:
            template<class F>
            static void func(F& f,type0& __max_a_lcl,Vec<type0>** __vs)
            {
                f.max_alpha_lcl_alpha(__max_a_lcl,(*__vs)->begin());
                GetMaxAlphaC<C>::func(f,__max_a_lcl,__vs+1);
            }
        };

        template<bool X,bool ALPHA,bool C>
        class GetMaxAlphaX
        {
        public:
            template<class F>
            static void func(F& f,type0& __max_a_lcl,Vec<type0>** __vs)
            {
                GetMaxAlphaAlpha<ALPHA,C>::func(f,__max_a_lcl,__vs);
            }
        };
        
        template<bool ALPHA,bool C>
        class GetMaxAlphaX<true,ALPHA,C>
        {
        public:
            template<class F>
            static void func(F& f,type0& __max_a_lcl,Vec<type0>** __vs)
            {
                f.max_alpha_lcl_x(__max_a_lcl,(*__vs)->begin());
                GetMaxAlphaAlpha<ALPHA,C>::func(f,__max_a_lcl,__vs+1);
            }
        };
        
        
        template<bool C>
        class Update
        {
        public:
            template<class __F>
            static void F(__F& f,type0& __alpha)
            {
                f.x=f.x0+__alpha*f.x_d;
                f.dynamic.update();
            }
            template<class F>
            static void F_reset(F& f)
            {
                f.x=f.x0;
                f.dynamic.update();
            }
        };
        
        template<>
        class Update<true>
        {
        public:
            template<class __F>
            static void F(__F& f,type0& __alpha)
            {
                f.x=f.x0+__alpha*f.x_d;
                f.aprop_c();
                f.dynamic.update(f.atoms->c);
            }
            template<class F>
            static void F_reset(F& f)
            {
                f.x=f.x0;
                f.dynamic.update(f.atoms->c);
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
            static void func(F& f)
            {
            }
            
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

    
    
    
    
    

    
    
  
    
    

    
    
    
    template<bool BC,bool X,bool ALPHA,bool C>
    class MinDMDHandler
    {
    protected:
    public:
        static constexpr int N0=VTN<BC||X,ALPHA,C>::N;
        static constexpr int N1=VTN<X,ALPHA,C>::N;
        typedef VT<BC,VTN<BC||X,ALPHA,C>::N> VECTENS0;
        typedef VT<BC,VTN<X,ALPHA,C>::N> VECTENS1;
        
        class ForceFieldDMD* ff;
        class AtomsDMD* atoms;
        type0 max_dx,max_dalpha;
        bool H_dof[__dim__][__dim__];
        
        
        type0 f_h;
        int c_dim;
        NewDynamicDMD<BC,BC||X,ALPHA> dynamic;
        VECTENS0 x,x0,x_d;
        VECTENS1 h,f,f0;
        
        
        MinDMDHandler(AtomsDMD*,ForceFieldDMD*,type0,type0,bool (&)[__dim__][__dim__]);
        ~MinDMDHandler(){}
        void add_extra_vec(VECTENS1& __v)
        {
            new (&__v) VECTENS1(atoms,__dim__,X,c_dim,ALPHA,c_dim,C);
            Algebra::Do<N1>::func([this,&__v](int i){dynamic.add_xchng(__v.vecs[i]);});
        }
        
        void rm_extra_vec(VECTENS1& __v){__v.~VECTENS1();}
        void init();
        void fin();


        
        
        void max_alpha_lcl_x(type0&,type0*);
        void max_alpha_lcl_alpha(type0&,type0*);
        void max_alpha_lcl_c(type0&,type0*);
        void ls_prep(type0&,type0&,type0&);
        type0 F(type0);
        void F_reset();
        type0 force_calc();
        void prep();
        void prep_x_d();
        void prep_x_d_affine();
        void post_force_calc_box();
        void aprop_c();
        

    };
}
using namespace MAPP_NS;
/*--------------------------------------------
 
 --------------------------------------------*/
template<bool BC,bool X,bool ALPHA,bool C>
MinDMDHandler<BC,X,ALPHA,C>::MinDMDHandler(
AtomsDMD* __atoms,
ForceFieldDMD* __ff,
type0 __max_dx,
type0 __max_dalpha,
bool (&__H_dof)[__dim__][__dim__]
):
ff(__ff),
atoms(__atoms),
max_dx(__max_dx),
max_dalpha(__max_dalpha),
c_dim(__atoms->c_dim),
dynamic(__atoms,__ff,{},{__atoms->x_dof,__atoms->alpha_dof,__atoms->c_dof},{})
{
    Algebra::V_eq<__dim__*__dim__>(__H_dof[0],H_dof[0]);
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<bool BC,bool X,bool ALPHA,bool C>
void MinDMDHandler<BC,X,ALPHA,C>::init()
{
    
    f.~VECTENS1();
    new (&f) VECTENS1(atoms,ff->F_H,ff->f,X,ff->f_alpha,ALPHA,ff->f_c,C);
    
    f0.~VECTENS1();
    new (&f0) VECTENS1(atoms,__dim__,X,c_dim,ALPHA,c_dim,C);
    h.~VECTENS1();
    new (&h) VECTENS1(atoms,__dim__,X,c_dim,ALPHA,c_dim,C);
    
    x.~VECTENS0();
    new (&x) VECTENS0(atoms,atoms->H,atoms->x,BC||X,atoms->alpha,ALPHA,atoms->c,C);
    x0.~VECTENS0();
    new (&x0) VECTENS0(atoms,__dim__,BC||X,c_dim,ALPHA,c_dim,C);
    x_d.~VECTENS0();
    new (&x_d) VECTENS0(atoms,h,__dim__,X||BC,X,X&&!BC,
                        c_dim,ALPHA,ALPHA,ALPHA,
                        c_dim,C,C,C);
    
    
    Algebra::Do<N1>::func([this](int i)
    {
        dynamic.add_xchng(h.vecs[i]);
        dynamic.add_xchng(f0.vecs[i]);
    });
    
    Algebra::Do<N0>::func([this](int i)
    {
        dynamic.add_xchng(x0.vecs[i]);
    });
    if(BC) dynamic.add_xchng(x_d.vecs[0]);
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<bool BC,bool X,bool ALPHA,bool C>
void MinDMDHandler<BC,X,ALPHA,C>::fin()
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
template<bool BC,bool X,bool ALPHA,bool C>
void MinDMDHandler<BC,X,ALPHA,C>::max_alpha_lcl_x
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
template<bool BC,bool X,bool ALPHA,bool C>
void MinDMDHandler<BC,X,ALPHA,C>::max_alpha_lcl_alpha(type0& max_a_lcl,type0* x_d_alphavec)
{
    
    int n=atoms->natms_lcl*c_dim;
    type0 max_x_d_alpha_lcl=0.0;
    type0* alphavec=atoms->alpha->begin();
    for(int i=0;i<n;i++)
    {
        if(x_d_alphavec[i]<0.0)
            max_a_lcl=MIN(-alphavec[i]/x_d_alphavec[i],max_a_lcl);
        max_x_d_alpha_lcl=MAX(max_x_d_alpha_lcl,fabs(x_d_alphavec[i]));
    }
    max_a_lcl=MIN(max_a_lcl,fabs(max_dalpha/max_x_d_alpha_lcl));
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<bool BC,bool X,bool ALPHA,bool C>
void MinDMDHandler<BC,X,ALPHA,C>::max_alpha_lcl_c(type0& max_a_lcl,type0* x_d_cvec)
{
    type0 max_c_ratio=max_a_lcl,dc_i=0.0,dh_i=0.0;
    type0* cvec=atoms->c->begin();
    type0 x_d_cv,cv;
    type0 lo=1.0+2.0*std::numeric_limits<type0>::epsilon();
    type0 volatile cj;
    int natms_lcl=atoms->natms_lcl;
    for(int i=0;i<natms_lcl;i++,x_d_cvec+=c_dim,cvec+=c_dim)
    {
        x_d_cv=0.0;
        cv=1.0;
        for(int j=0;j<c_dim;j++)
        {
            //printf("%.16lf %.16lf\n",cvec[j],x_d_cvec[j]);
            if(cvec[j]<0.0) continue;
            x_d_cv-=x_d_cvec[j];
            cv-=cvec[j];
            cj=1.0+cvec[j];
            if(x_d_cvec[j]<0.0 && cj+(max_c_ratio*x_d_cvec[j])<lo)
            {
                
                max_c_ratio=(cj-lo)/(-x_d_cvec[j]);
                dc_i=cj-lo;
                dh_i=-x_d_cvec[j];
            }
        }
        
        cj=1.0+cv;
        if(x_d_cv<0.0 && cj+(max_c_ratio*x_d_cv)<lo)
        {
            max_c_ratio=(cj-lo)/(-x_d_cv);
            dc_i=cj-lo;
            dh_i=-x_d_cv;
        }
        
        bool chk=false;
        while(!chk)
        {
            chk=true;
            cv=1.0;
            for(int j=0;j<c_dim && chk;j++)
            {
                if(cvec[j]<0.0) continue;
                cj=1.0+cvec[j];
                cj+=(max_c_ratio*x_d_cvec[j]);
                if(cj<lo && x_d_cvec[j]<0.0) chk=false;
                --cj;
                cv-=cj;
            }
            cj=1.0+cv;
            if(cj<lo && x_d_cv<0.0) chk=false;
            if(chk) continue;
            dc_i-=std::numeric_limits<type0>::epsilon();
            max_c_ratio=dc_i/dh_i;
        }
        
    }
    
    
    max_a_lcl=MIN(max_a_lcl,max_c_ratio);
    
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<bool BC,bool X,bool ALPHA,bool C>
void MinDMDHandler<BC,X,ALPHA,C>::ls_prep(type0& dfa,type0& h_norm,type0& max_a)
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
    MinDMDHelper::GetMaxAlphaX<BC||X,ALPHA,C>::func(*this,max_a_lcl,x_d.vecs);
    MPI_Allreduce(&max_a_lcl,&max_a,1,Vec<type0>::MPI_T,MPI_MIN,atoms->world);
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<bool BC,bool X,bool ALPHA,bool C>
type0 MinDMDHandler<BC,X,ALPHA,C>::F(type0 __alpha)
{
    MinDMDHelper::Update<C>::F(*this,__alpha);
    return ff->new_value<C>();
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<bool BC,bool X,bool ALPHA,bool C>
void MinDMDHandler<BC,X,ALPHA,C>::F_reset()
{
    return MinDMDHelper::Update<C>::F_reset(*this);
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<bool BC,bool X,bool ALPHA,bool C>
type0 MinDMDHandler<BC,X,ALPHA,C>::force_calc()
{
    type0 ans=*ff->new_derivative<C>();
    MinDMDHelper::PostForceCalcBC<BC>::func(*this);
    return ans;
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<bool BC,bool X,bool ALPHA,bool C>
void MinDMDHandler<BC,X,ALPHA,C>::prep()
{
    MinDMDHelper::Prep<BC,X>::func(*this);
    
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<bool BC,bool X,bool ALPHA,bool C>
void MinDMDHandler<BC,X,ALPHA,C>::prep_x_d()
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
template<bool BC,bool X,bool ALPHA,bool C>
void MinDMDHandler<BC,X,ALPHA,C>::prep_x_d_affine()
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
template<bool BC,bool X,bool ALPHA,bool C>
void MinDMDHandler<BC,X,ALPHA,C>::post_force_calc_box()
{
    type0 (&S_fe)[__dim__][__dim__]=atoms->S_fe;
    type0 neg_v=-atoms->vol;
    Algebra::DoLT<__dim__>::func([this,&S_fe,&neg_v](int i,int j)
    {f.A[i][j]=H_dof[i][j] ? S_fe[i][j]*neg_v:0.0;});
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<bool BC,bool X,bool ALPHA,bool C>
void MinDMDHandler<BC,X,ALPHA,C>::aprop_c()
{
    const int n=atoms->natms_lcl*c_dim;
    type0* c_vec=atoms->c->begin();
    type0 volatile cj;
    for(int i=0;i<n;i++)
        if(c_vec[i]>=0.0)
        {
            cj=c_vec[i];
            --++cj;
            c_vec[i]=cj;
        }
}
/*------------------------------------------------------------------------------------------------*/
namespace MAPP_NS
{
    class MinCGDMD:public Min
    {
    template<bool>
    friend class MinHelper::ThermoHandler;
    template<bool>
    friend class MinHelper::ExportHandler;
    friend class DAE;
    private:
    protected:
        static constexpr int MaxS=MinHelper::MaxSizeAlign<MinDMDHandler,4>::MaxS;
        static constexpr int MaxA=MinHelper::MaxSizeAlign<MinDMDHandler,4>::MaxA;
        typedef std::aligned_storage<MaxS,MaxA>::type MemT;
        bool ALPHA_DOF,C_DOF;
        type0 max_dalpha;
        MemT handler_buff;

        class AtomsDMD* atoms;
        class ForceFieldDMD* ff;
        class ExportDMD* xprt;
        void pre_run_chk(AtomsDMD*,ForceFieldDMD*);
        void create_thermo();
        void destroy_thermo();
    public:
        MinCGDMD(type0,bool(&)[__dim__][__dim__],bool,type0,type0,class LineSearch*);
        ~MinCGDMD();
        virtual void run(int);
        virtual void init();
        virtual void fin();
        
        template<bool BC,bool X,bool ALPHA,bool C>
        void  __init();
        template<bool BC,bool X,bool ALPHA,bool C,bool OUT,bool XOUT,class LS>
        void  __run(LS*,int);
        template<bool BC,bool X,bool ALPHA,bool C>
        void  __fin();
        
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
        static void getset_x_dof(PyGetSetDef&);
        static void getset_alpha_dof(PyGetSetDef&);
        static void getset_c_dof(PyGetSetDef&);
        static void getset_export(PyGetSetDef&);
        
        static int setup_tp();
        
        
        
    };

}
using namespace MAPP_NS;
/*--------------------------------------------
 
 --------------------------------------------*/
template<bool BC,bool X,bool ALPHA,bool C>
void MinCGDMD::__init()
{
    
    MinDMDHandler<BC,X,ALPHA,C>* __handler_ptr=new (&handler_buff) MinDMDHandler<BC,X,ALPHA,C>(atoms,ff,max_dx,max_dalpha,H_dof);
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
            __fin<BC,X,ALPHA,C>();
            throw err_msg;
        }
    }
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<bool BC,bool X,bool ALPHA,bool C,bool OUT,bool XOUT,class LS>
void MinCGDMD::__run(LS* ls,int nsteps)
{
    
    MinDMDHandler<BC,X,ALPHA,C>& handler=*reinterpret_cast<MinDMDHandler<BC,X,ALPHA,C>*>(&handler_buff);
    typedef typename MinDMDHandler<BC,X,ALPHA,C>::VECTENS0 VECTENS0;
    typedef typename MinDMDHandler<BC,X,ALPHA,C>::VECTENS1 VECTENS1;
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
    
    type0 e_prev,e_curr=atoms->fe;
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
template<bool BC,bool X,bool ALPHA,bool C>
void MinCGDMD::__fin()
{
    if(xprt)
    {
        xprt->fin();
        xprt->atoms=NULL;
    }
    MinDMDHandler<BC,X,ALPHA,C>* __handler_ptr=reinterpret_cast<MinDMDHandler<BC,X,ALPHA,C>*>(&handler_buff);
    __handler_ptr->dynamic.fin();
    __handler_ptr->fin();
    __handler_ptr->~MinDMDHandler();
}
#endif
