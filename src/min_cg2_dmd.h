#ifndef __MAPP__min_cg2_dmd__
#define __MAPP__min_cg2_dmd__
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
    class MinDMD2HandlerGetMaxAlphaMu
    {
    public:
        template<class F>
        static void func(F& f,type0& __max_a_lcl,Vec<type0>** __vs)
        {
            f.max_alpha_lcl_c(__max_a_lcl,(*__vs)->begin());
        }
    };
    template<>
    class MinDMD2HandlerGetMaxAlphaMu<false>
    {
    public:
        template<class F>
        static void func(F&,type0&,Vec<type0>**)
        {
        }
    };
    
    template<bool ALPHA,bool C>
    class MinDMD2HandlerGetMaxAlphaAlpha
    {
    public:
        template<class F>
        static void func(F& f,type0& __max_a_lcl,Vec<type0>** __vs)
        {
            f.max_alpha_lcl_alpha(__max_a_lcl,(*__vs)->begin());
            MinDMD2HandlerGetMaxAlphaMu<C>::func(f,__max_a_lcl,__vs+1);
        }
    };
    template<bool C>
    class MinDMD2HandlerGetMaxAlphaAlpha<false,C>
    {
    public:
        template<class F>
        static void func(F& f,type0& __max_a_lcl,Vec<type0>** __vs)
        {
            MinDMD2HandlerGetMaxAlphaMu<C>::func(f,__max_a_lcl,__vs);
        }
    };
    
    
    template<bool X,bool ALPHA,bool C>
    class MinDMD2HandlerGetMaxAlphaX
    {
    public:
        template<class F>
        static void func(F& f,type0& __max_a_lcl,Vec<type0>** __vs)
        {
            f.max_alpha_lcl_x(__max_a_lcl,(*__vs)->begin());
            MinDMD2HandlerGetMaxAlphaAlpha<ALPHA,C>::func(f,__max_a_lcl,__vs+1);
        }
    };
    
    
    template<bool ALPHA,bool C>
    class MinDMD2HandlerGetMaxAlphaX<false,ALPHA,C>
    {
    public:
        template<class F>
        static void func(F& f,type0& __max_a_lcl,Vec<type0>** __vs)
        {
            MinDMD2HandlerGetMaxAlphaAlpha<ALPHA,C>::func(f,__max_a_lcl,__vs);
        }
    };
    
    
    template<bool C>
    class MinDMD2HandlerUpdate
    {
    public:
        template<class __F>
        static void F(__F& f,type0& __alpha)
        {
            f.x=f.x0+__alpha*f.x_d;
            f.eval_c();
            f.dynamic->update(f.atoms->c);
        }
        template<class F>
        static void F_reset(F& f)
        {
            f.x=f.x0;
            f.eval_c();
            f.dynamic->update(f.atoms->c);
        }
    };
    
    template<>
    class MinDMD2HandlerUpdate<false>
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
    class MinDMD2HandlerPrep
    {
    public:
        template<class F>
        static void func(F& f)
        {}
        
    };
    
    template<>
    class MinDMD2HandlerPrep<true,true>
    {
    public:
        template<class F>
        static void func(F& f)
        {
            f.prep_x_d();
        }
        
    };
    
    template<>
    class MinDMD2HandlerPrep<true,false>
    {
    public:
        template<class F>
        static void func(F& f)
        {
            f.prep_x_d_affine();
        }
        
    };
    
    
    template<bool BC>
    class MinDMD2HandlerPostForceCalcBC
    {
    public:
        template<class F>
        static void func(F& f)
        {
            f.post_force_calc_box();
        }
        
    };
    
    template<>
    class MinDMD2HandlerPostForceCalcBC<false>
    {
    public:
        template<class F>
        static void func(F& f)
        {
        }
        
    };
    
    template<bool BC,bool X,bool ALPHA,bool MU>
    class MinDMD2Handler
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
            memcpy(x_dvec,hvec,natms_lcl*__dim__*sizeof(type0));
            for(int iatm=0;iatm<natms_lcl;iatm++,xvec+=__dim__,x_dvec+=__dim__)
                Algebra::V_mul_MLT_add_in(xvec,h.A,x_dvec);
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
        
        
        static constexpr int N0=VTN<BC||X,ALPHA,MU>::N;
        static constexpr int N1=VTN<X,ALPHA,MU>::N;
        typedef VT<BC,VTN<BC||X,ALPHA,MU>::N> VECTENS0;
        typedef VT<BC,VTN<X,ALPHA,MU>::N> VECTENS1;
        
        
        class ForceFieldDMD* ff;
        class AtomsDMD* atoms;
        type0 max_dx,max_dalpha;
        bool H_dof[__dim__][__dim__];
        
        Vec<type0>* beta_mu;
        type0* beta_mu_0;
        
        type0 f_h;
        int c_dim;
        class NewDynamicDMD<BC,BC||X,ALPHA>* dynamic;
        VECTENS0 x,x0,x_d;
        VECTENS1 h,f,f0;
        
        MinDMD2Handler(
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
        beta_mu(NULL),
        beta_mu_0(NULL),
        c_dim(__atoms->c_dim),
        dynamic(NULL)
        {
            Algebra::V_eq<__dim__*__dim__>(__H_dof[0],H_dof[0]);
        }
        
        MinDMD2Handler():
        atoms(NULL),
        ff(NULL),
        max_dx(0.0),
        max_dalpha(0.0),
        c_dim(0.0),
        dynamic(NULL),
        beta_mu(NULL),
        beta_mu_0(NULL)
        {
            Algebra::set<__dim__*__dim__>(H_dof[0],false);
        }
        
        void add_extra_vec(VECTENS1& __v)
        {
            new (&__v) VECTENS1(atoms,__dim__,X,c_dim,ALPHA,c_dim,MU);
            Algebra::Do<N1>::func([this,&__v](int i){dynamic->add_xchng(__v.vecs[i]);});
        }
        
        void rm_extra_vec(VECTENS1& __v){__v.~VECTENS1();}
        void init()
        {
            if(MU || ALPHA) beta_mu=new Vec<type0>(atoms,c_dim);
            
            
            
            f.~VECTENS1();
            new (&f) VECTENS1(atoms,ff->F_H,ff->f,X,ff->f_alpha,ALPHA,ff->f_c,MU);
            
            f0.~VECTENS1();
            new (&f0) VECTENS1(atoms,__dim__,X,c_dim,ALPHA,c_dim,MU);
            h.~VECTENS1();
            new (&h) VECTENS1(atoms,__dim__,X,c_dim,ALPHA,c_dim,MU);
            
            x.~VECTENS0();
            new (&x) VECTENS0(atoms,atoms->H,atoms->x,BC||X,atoms->alpha,ALPHA,beta_mu,MU);
            x0.~VECTENS0();
            new (&x0) VECTENS0(atoms,__dim__,BC||X,c_dim,ALPHA,c_dim,MU);
            x_d.~VECTENS0();
            new (&x_d) VECTENS0(atoms,h,__dim__,X||BC,X,X&&!BC,
                                c_dim,ALPHA,ALPHA,ALPHA,
                                c_dim,MU,MU,MU);
            
            dynamic=new NewDynamicDMD<BC,BC||X,ALPHA>(atoms,ff,{},{atoms->x_dof,atoms->alpha_dof,atoms->c_dof},{});
            if(MU || ALPHA) dynamic->add_updt(beta_mu);
            
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
            if(MU || ALPHA) eval_beta_mu();
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
            
            delete beta_mu;
            beta_mu=NULL;
            delete [] beta_mu_0;
            beta_mu_0=NULL;
        }
        void eval_beta_mu_0()
        {
            size_t nelems=atoms->elements.nelems;
            type0* masses=atoms->elements.masses;
            type0 kBT=(atoms->temp*atoms->kB);
            type0 h=atoms->hP;
            if(nelems) beta_mu_0=new type0[nelems];
            for(size_t i=0;i<nelems;i++)
                beta_mu_0[i]=3.0*log(h/(M_PI*sqrt(masses[i]*kBT)));
        }
        void eval_beta_mu()
        {
            type0* __beta_mu=beta_mu->begin();
            type0* __alpha=atoms->alpha->begin();
            type0* __c=atoms->c->begin();
            elem_type* __elem=atoms->elem->begin();
            int __natms_lcl=atoms->natms_lcl;
            type0 cv;
            for(int i=0;i<__natms_lcl;i++,__beta_mu+=c_dim,__alpha+=c_dim,__c+=c_dim,__elem+=c_dim)
            {
                cv=1.0;
                for(int j=0;j<c_dim;j++)
                {
                    if(__c[j]<0.0) continue;
                    cv-=__c[j];
                }
                
                for(int j=0;j<c_dim;j++)
                {
                    if(__c[j]<0.0)
                    {
                        __beta_mu[j]=0.0;
                        continue;
                    }
                    //__beta_mu[j]=beta_mu_0[__elem[j]]+log(__c[j]/cv)-3.0*log(__alpha[j]);
                    __beta_mu[j]=log(__c[j]/cv)-3.0*log(__alpha[j]);
                }
            }
        }

        void eval_c()
        {
            type0* __beta_mu=beta_mu->begin();
            type0* __alpha=atoms->alpha->begin();
            type0* __c=atoms->c->begin();
            elem_type* __elem=atoms->elem->begin();
            int __natms_lcl=atoms->natms_lcl;
            type0 q,c0,c1,c_dim_0,l;
            //type0 volatile cj;
            type0 eps=std::numeric_limits<type0>::epsilon();
            for(int i=0;i<__natms_lcl;i++,__beta_mu+=c_dim,__alpha+=c_dim,__c+=c_dim,__elem+=c_dim)
            {
                
                q=1.0;
                c_dim_0=0.0;
                for(int j=0;j<c_dim;j++)
                {
                    if(__c[j]<0.0) continue;
                    //__c[j]=exp(__beta_mu[j]-beta_mu_0[__elem[j]])*Algebra::pow<3>(__alpha[j]);
                    __c[j]=exp(__beta_mu[j])*Algebra::pow<3>(__alpha[j]);
                    
                    q+=__c[j];
                    ++c_dim_0;
                }
                
                l=1.0;
                c1=MIN(MAX(1.0+l/q,1.0+eps),2.0-c_dim_0*eps);
                c0=c1;
                --c_dim_0;
                
                for(int j=0;j<c_dim-1;j++)
                {
                    if(__c[j]<0.0) continue;
                    l+=__c[j];
                    c1=MIN(MAX(1.0+l/q,c0+eps),2.0-c_dim_0*eps);
                    __c[j]=c1-c0;
                    
                    c0=c1;
                    --c_dim_0;

                }
                
                if(__c[c_dim-1]>=0.0)
                {
                    __c[c_dim-1]=2.0-c0;
                }
                
                
                /*
                q=0.0;
                type0 n_dim=1.0;
                for(int j=0;j<c_dim;j++)
                {
                    if(__c[j]<0.0) continue;
                    //__c[j]=exp(__beta_mu[j]-beta_mu_0[__elem[j]])*Algebra::pow<3>(__alpha[j]);
                    n_dim++;
                    q+=__beta_mu[j]+3.0*log(__alpha[j]);;
                }
                q/=n_dim;
                
                
                type0 z=exp(-q);
                for(int j=0;j<c_dim;j++)
                {
                    if(__c[j]<0.0) continue;
                    __c[j]=exp(__beta_mu[j]+3.0*log(__alpha[j])-q);
                    z+=__c[j];
                }
                
                for(int j=0;j<c_dim;j++)
                {
                    if(__c[j]<0.0) continue;
                    cj=__c[j]/z;
                    --++cj;
                    __c[j]=cj;
                }
                
                if(std::isnan(__c[0]) || std::isnan(__c[1]))
                {
                    
                }
                */
                
            }
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
                    max_a_lcl=MIN((std::numeric_limits<type0>::epsilon()-alphavec[i])/x_d_alphavec[i],max_a_lcl);
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
            MinDMD2HandlerGetMaxAlphaX<BC||X,ALPHA,MU>::func(*this,max_a_lcl,x_d.vecs);
            MPI_Allreduce(&max_a_lcl,&max_a,1,Vec<type0>::MPI_T,MPI_MIN,atoms->world);
        }
        type0 F(type0 __alpha)
        {
            MinDMD2HandlerUpdate<ALPHA || MU>::F(*this,__alpha);
            return ff->value_gp();
        }
        void F_reset()
        {
            return MinDMD2HandlerUpdate<ALPHA || MU>::F_reset(*this);
        }
        
        type0 force_calc()
        {
            type0 ans=*ff->derivative_gp();
            MinDMD2HandlerPostForceCalcBC<BC>::func(*this);
            return ans;
        }
        void prep()
        {
            MinDMD2HandlerPrep<BC,X>::func(*this);
            
        }
        

        

    };
    
    
    
    
    
}



    
    
namespace MAPP_NS
{
    class MinCG2DMD:public Min
    {
    private:
    protected:
        bool B_DOF,X_DOF,ALPHA_DOF,MU_DOF;
        type0 max_dalpha,max_dbetamu;
        

        class AtomsDMD* atoms;
        class ForceFieldDMD* ff;
        class ExportDMD* xprt;
        void pre_run_chk(AtomsDMD*,ForceFieldDMD*);
        
    public:
        MinCG2DMD(type0,bool(&)[__dim__][__dim__],bool,type0,type0,type0,class LineSearch*);
        ~MinCG2DMD();
        virtual void run(int);
        virtual void init();
        virtual void fin();
        
        template<bool BC,bool X,bool ALPHA,bool C,class LS>
        void  run(LS* ls,int nsteps)
        {
            ls->reset();
            MinDMD2Handler<BC,X,ALPHA,C> handler(atoms,ff,max_dx,max_dalpha,H_dof);
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
            
            
            
            
            
            
            
            
            
            
            typename MinDMD2Handler<BC,X,ALPHA,C>::VECTENS1& f=handler.f;
            typename MinDMD2Handler<BC,X,ALPHA,C>::VECTENS1& f0=handler.f0;
            typename MinDMD2Handler<BC,X,ALPHA,C>::VECTENS1& h=handler.h;
            typename MinDMD2Handler<BC,X,ALPHA,C>::VECTENS0& x0=handler.x0;
            typename MinDMD2Handler<BC,X,ALPHA,C>::VECTENS0& x=handler.x;
            type0& f_h=handler.f_h;
            
            
            int step=atoms->step;
            
            type0 e_curr=handler.force_calc();
            
            int nevery_xprt=xprt==NULL ? 0:xprt->nevery;
            if(nevery_xprt) xprt->write(step);
            
            ThermoDynamics thermo(6,
              "GP",atoms->gp,
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
            MinCG2DMD* min;
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
        static void getset_max_dbetamu(PyGetSetDef&);
        static void getset_x_dof(PyGetSetDef&);
        static void getset_alpha_dof(PyGetSetDef&);
        static void getset_mu_dof(PyGetSetDef&);
        static void getset_export(PyGetSetDef&);
        
        static int setup_tp();
        
        
        
    };

}

#endif
