#ifndef __MAPP__min_cg3_dmd__
#define __MAPP__min_cg3_dmd__
#include "min.h"
#include "new_min_vec.h"
#include "ff_dmd.h"
#include "atoms_dmd.h"
#include "dynamic_dmd.h"
#include "MAPP.h"
#include "thermo_dynamics.h"
namespace MAPP_NS
{
    template<bool C>
    class MinDMDHandlerGetMaxAlphaC
    {
    public:
        template<class F>
        static void func(F& f,type0& __max_a_lcl,Vec<type0>** __vs)
        {
            f.max_alpha_lcl_c(__max_a_lcl,(*__vs)->begin());
        }
    };
    template<>
    class MinDMDHandlerGetMaxAlphaC<false>
    {
    public:
        template<class F>
        static void func(F&,type0&,Vec<type0>**)
        {
        }
    };
    
    template<bool ALPHA,bool C>
    class MinDMDHandlerGetMaxAlphaAlpha
    {
    public:
        template<class F>
        static void func(F& f,type0& __max_a_lcl,Vec<type0>** __vs)
        {
            f.max_alpha_lcl_alpha(__max_a_lcl,(*__vs)->begin());
            MinDMDHandlerGetMaxAlphaC<C>::func(f,__max_a_lcl,__vs+1);
        }
    };
    template<bool C>
    class MinDMDHandlerGetMaxAlphaAlpha<false,C>
    {
    public:
        template<class F>
        static void func(F& f,type0& __max_a_lcl,Vec<type0>** __vs)
        {
            MinDMDHandlerGetMaxAlphaC<C>::func(f,__max_a_lcl,__vs);
        }
    };
    
    
    template<bool X,bool ALPHA,bool C>
    class MinDMDHandlerGetMaxAlphaX
    {
    public:
        template<class F>
        static void func(F& f,type0& __max_a_lcl,Vec<type0>** __vs)
        {
            f.max_alpha_lcl_x(__max_a_lcl,(*__vs)->begin());
            MinDMDHandlerGetMaxAlphaAlpha<ALPHA,C>::func(f,__max_a_lcl,__vs+1);
        }
    };
    
    
    template<bool ALPHA,bool C>
    class MinDMDHandlerGetMaxAlphaX<false,ALPHA,C>
    {
    public:
        template<class F>
        static void func(F& f,type0& __max_a_lcl,Vec<type0>** __vs)
        {
            MinDMDHandlerGetMaxAlphaAlpha<ALPHA,C>::func(f,__max_a_lcl,__vs);
        }
    };
    
    
    template<bool C>
    class MinDMDHandlerUpdate
    {
    public:
        template<class __F>
        static void F(__F& f,type0& __alpha)
        {
            f.x=f.x0+__alpha*f.x_d;
            f.aprop_c();
            f.dynamic->update(f.atoms->c);
        }
        template<class F>
        static void F_reset(F& f)
        {
            f.x=f.x0;
            f.dynamic->update(f.atoms->c);
        }
    };
    
    template<>
    class MinDMDHandlerUpdate<false>
    {
    public:
        template<class __F>
        static void F(__F& f,type0& __alpha)
        {
            f.x=f.x0+__alpha*f.x_d;
            f.dynamic->update();
        }
        template<class F>
        static void F_reset(F& f)
        {
            f.x=f.x0;
            f.dynamic->update();
        }
    };
    
    
    template<bool BC,bool X>
    class MinDMDHandlerPrep
    {
    public:
        template<class F>
        static void func(F& f)
        {}
        
    };
    
    template<>
    class MinDMDHandlerPrep<true,true>
    {
    public:
        template<class F>
        static void func(F& f)
        {
            f.prep_x_d();
        }
        
    };
    
    template<>
    class MinDMDHandlerPrep<true,false>
    {
    public:
        template<class F>
        static void func(F& f)
        {
            f.prep_x_d_affine();
        }
        
    };
    
    
    template<bool BC>
    class MinDMDHandlerPostForceCalcBC
    {
    public:
        template<class F>
        static void func(F& f)
        {
            f.post_force_calc_box();
        }
        
    };
    
    template<>
    class MinDMDHandlerPostForceCalcBC<false>
    {
    public:
        template<class F>
        static void func(F& f)
        {
        }
        
    };
    
    template<bool BC,bool X,bool ALPHA,bool C>
    class MinDMDHandler
    {
    protected:
    public:
        void prep_x_d()
        {
            const int natms_lcl=atoms->natms_lcl;
            type0* xvec=x0.vecs[0]->begin();
            type0* hvec=h.vecs[0]->begin();
            type0* x_dvec=x_d.vecs[0]->begin();
            Algebra::MLT_mul_MLT(atoms->H,h.A,x_d.A);
            
            for(int iatm=0;iatm<natms_lcl;iatm++,xvec+=__dim__,x_dvec+=__dim__,hvec+=__dim__)
            {
                Algebra::V_eq<__dim__>(hvec,x_dvec);
                Algebra::V_mul_MLT_add_in(xvec,h.A,x_dvec);
            }
        }
        void prep_x_d_affine()
        {
            const int natms_lcl=atoms->natms_lcl;
            type0* xvec=x0.vecs[0]->begin();
            type0* x_dvec=x_d.vecs[0]->begin();
            Algebra::MLT_mul_MLT(atoms->H,h.A,x_d.A);
            for(int iatm=0;iatm<natms_lcl;iatm++,xvec+=__dim__,x_dvec+=__dim__)
                Algebra::V_mul_MLT(xvec,h.A,x_dvec);
        }
        
        
        void post_force_calc_box()
        {
            type0 (&S_fe)[__dim__][__dim__]=atoms->S_fe;
            type0 neg_v=-atoms->vol;
            Algebra::DoLT<__dim__>::func([this,&S_fe,&neg_v](int i,int j)
            {f.A[i][j]=H_dof[i][j] ? S_fe[i][j]*neg_v:0.0;});
        }
        
        void aprop_c()
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
        class NewDynamicDMD<BC,BC||X,ALPHA>* dynamic;
        VECTENS0 x,x0,x_d;
        VECTENS1 h,f,f0;
        
        MinDMDHandler(
        AtomsDMD* __atoms,
        ForceFieldDMD* __ff,
        type0 __max_dx,
        type0 __max_dalpha,
        bool (&__H_dof)[__dim__][__dim__]
        ):
        atoms(__atoms),
        ff(__ff),
        max_dx(__max_dx),
        max_dalpha(__max_dalpha),
        c_dim(__atoms->c_dim),
        dynamic(NULL)
        {
            Algebra::V_eq<__dim__*__dim__>(__H_dof[0],H_dof[0]);
        }
        
        MinDMDHandler():
        atoms(NULL),
        ff(NULL),
        max_dx(0.0),
        max_dalpha(0.0),
        c_dim(0.0),
        dynamic(NULL)
        {
            Algebra::set<__dim__*__dim__>(H_dof[0],false);
        }
        
        void add_extra_vec(VECTENS1& __v)
        {
            new (&__v) VECTENS1(atoms,__dim__,X,c_dim,ALPHA,c_dim,C);
            Algebra::Do<N1>::func([this,&__v](int i){dynamic->add_xchng(__v.vecs[i]);});
        }
        
        void rm_extra_vec(VECTENS1& __v){__v.~VECTENS1();}
        void init()
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
            
            dynamic=new NewDynamicDMD<BC,BC||X,ALPHA>(atoms,ff,{},{atoms->x_dof,atoms->alpha_dof,atoms->c_dof},{});

            Algebra::Do<N1>::func([this](int i)
            {
                dynamic->add_xchng(h.vecs[i]);
                dynamic->add_xchng(f0.vecs[i]);
            });
            
            Algebra::Do<N0>::func([this](int i)
            {
                dynamic->add_xchng(x0.vecs[i]);
            });
            if(BC) dynamic->add_xchng(x_d.vecs[0]);
            
        }
        void fin()
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
        

        
        
        void max_alpha_lcl_x(type0& max_a_lcl,type0* x_d_xvec)
        {
            type0 max_x_d_lcl=0.0;
            const int natms_lcl=atoms->natms_lcl;
            int n=natms_lcl*__dim__;
            for(int i=0;i<n;i++)
                max_x_d_lcl=MAX(max_x_d_lcl,fabs(x_d_xvec[i]));
            
            max_a_lcl=MIN(fabs(max_dx/max_x_d_lcl),max_a_lcl);
        }
        
        void max_alpha_lcl_alpha(type0& max_a_lcl,type0* x_d_alphavec)
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
        
        void max_alpha_lcl_c(type0& max_a_lcl,type0* x_d_cvec)
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
        void ls_prep(type0& dfa,type0& h_norm,type0& max_a)
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
            MinDMDHandlerGetMaxAlphaX<BC||X,ALPHA,C>::func(*this,max_a_lcl,x_d.vecs);
            MPI_Allreduce(&max_a_lcl,&max_a,1,Vec<type0>::MPI_T,MPI_MIN,atoms->world);
        }
        type0 F(type0 __alpha)
        {
            MinDMDHandlerUpdate<C>::F(*this,__alpha);
            return ff->new_value<C>();
        }
        void F_reset()
        {
            return MinDMDHandlerUpdate<C>::F_reset(*this);
        }
        
        type0 force_calc()
        {
            type0 ans=*ff->new_derivative<C>();
            MinDMDHandlerPostForceCalcBC<BC>::func(*this);
            return ans;
        }
        void prep()
        {
            MinDMDHandlerPrep<BC,X>::func(*this);
            
        }
        

        

    };
    
    
    
    
    
}



    
    
namespace MAPP_NS
{
    class MinCG3DMD:public Min
    {
    private:
    protected:
        bool B_DOF,X_DOF,ALPHA_DOF,C_DOF;
        type0 max_dalpha;
        

        class AtomsDMD* atoms;
        class ForceFieldDMD* ff;
        class ExportDMD* xprt;
        void pre_run_chk(AtomsDMD*,ForceFieldDMD*);
        
    public:
        MinCG3DMD(type0,bool(&)[__dim__][__dim__],bool,type0,type0,class LineSearch*);
        ~MinCG3DMD();
        virtual void run(int);
        virtual void init();
        virtual void fin();

        
        template<bool BC,bool X,bool ALPHA,bool C,class LS>
        void  run(LS* ls,int nsteps)
        {
            ls->reset();
            MinDMDHandler<BC,X,ALPHA,C> handler(atoms,ff,max_dx,max_dalpha,H_dof);
            handler.init();
            handler.dynamic->init();
            
            
            
            if(xprt)
            {
                try
                {
                    xprt->atoms=atoms;
                    xprt->init();
                }
                catch(std::string& err_msg)
                {
                    xprt->fin();
                    xprt->atoms=NULL;
                    handler.dynamic->fin();
                    handler.fin();
                    throw err_msg;
                }
            }
            
            
            
            
            
            
            
            
            
            
            typename MinDMDHandler<BC,X,ALPHA,C>::VECTENS1& f=handler.f;
            typename MinDMDHandler<BC,X,ALPHA,C>::VECTENS1& f0=handler.f0;
            typename MinDMDHandler<BC,X,ALPHA,C>::VECTENS1& h=handler.h;
            typename MinDMDHandler<BC,X,ALPHA,C>::VECTENS0& x0=handler.x0;
            typename MinDMDHandler<BC,X,ALPHA,C>::VECTENS0& x=handler.x;
            type0& f_h=handler.f_h;
            
            
            int step=atoms->step;
            
            type0 e_curr=handler.force_calc();
            
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
            
            
            type0 e_prev;
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
            
            
            
            
            
            
            
            
            handler.dynamic->fin();
            handler.fin();
        }
        

        
        typedef struct
        {
            PyObject_HEAD
            MinCG3DMD* min;
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

#endif
