#ifndef __MAPP__min_cg3_dmd__
#define __MAPP__min_cg3_dmd__
#include "min.h"
#include "min_vec.h"
#include "ff_dmd.h"
#include "atoms_dmd.h"
#include "dynamic_dmd.h"
#include "MAPP.h"
#include "thermo_dynamics.h"
namespace MAPP_NS
{
    template<bool BC,bool X,bool ALPHA,bool C>
    class MinDMDMan_t
    {
    public:
        typedef VecTens<type0,3> VECTENS0;
        typedef VecTens<type0,3> VECTENS1;
    };
    template<>
    class MinDMDMan_t<true,true,true,true>
    {
    public:
        typedef VecTens<type0,3> VECTENS0;
        typedef VecTens<type0,3> VECTENS1;
    };
    template<>
    class MinDMDMan_t<false,true,true,true>
    {
    public:
        typedef VecTens<type0,3> VECTENS0;
        typedef VecTens<type0,3> VECTENS1;
    };
    template<>
    class MinDMDMan_t<true,false,true,true>
    {
    public:
        typedef VecTens<type0,3> VECTENS0;
        typedef VecTens<type0,2> VECTENS1;
    };
    template<>
    class MinDMDMan_t<false,false,true,true>
    {
    public:
        typedef VecTens<type0,2> VECTENS0;
        typedef VecTens<type0,2> VECTENS1;
    };
    
    template<bool BC,bool X,bool ALPHA,bool C>
    class MinDMDMan
    {
    public:
        typedef typename MinDMDMan_t<BC,X,ALPHA,C>::VECTENS0 VECTENS0;
        typedef typename MinDMDMan_t<BC,X,ALPHA,C>::VECTENS1 VECTENS1;
        class ForceFieldDMD* ff;
        class AtomsDMD* atoms;
        class DynamicDMD* dynamic;
        int c_dim;
        type0 max_dx,max_dalpha;
        bool H_dof[__dim__][__dim__];
        VECTENS0 x,x0,x_d;
        VECTENS1 h,f,f0;
        void add_extra_vec(VECTENS1&){}
        void rm_extra_vec(VECTENS1& v){v.~VECTENS1();}
        void init(){}
        void fin()
        {
            
            dynamic->fin();
            delete dynamic;
            dynamic=NULL;
            
            
            x_d.~VECTENS0();
            x0.~VECTENS0();
            x.~VECTENS0();
            
            h.~VECTENS1();
            f0.~VECTENS1();
            f.~VECTENS1();
        }
        void force_calc(){}
        void prep(){}
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
        
        void post_force_calc_c()
        {
            int natms_lcl=atoms->natms_lcl;
            type0* c=atoms->c->begin();
            type0* mu=ff->mu->begin();
            type0 cv;
            type0 logcv;
            type0 kBT=atoms->temp*atoms->kB;
            for(int i=0;i<natms_lcl;i++,c+=c_dim,mu+=c_dim)
            {
                cv=1.0;
                for(int j=0;j<c_dim;j++)
                {
                    if(c[j]>=0.0)
                        cv-=c[j];
                }
                logcv=log(cv);
                for(int j=0;j<c_dim;j++)
                {
                    if(c[j]>=0.0)
                        mu[j]=-mu[j]-kBT*(log(c[j])-logcv);
                }
            }
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
        type0 F(type0);
        
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
            type0 lo=1.0+std::numeric_limits<type0>::epsilon();
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
                
                //printf("%.16lf %.16lf\n",cv,x_d_cv);
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
                        //printf("%.52lf\n",cj);
                        if(cj<lo) chk=false;
                        --cj;
                        cv-=cj;
                    }
                    cj=1.0+cv;
                    if(cj<lo) chk=false;
                    if(chk) continue;
                    dc_i-=std::numeric_limits<type0>::epsilon();
                    max_c_ratio=dc_i/dh_i;
                }
            }
            
            
            max_a_lcl=MIN(max_a_lcl,max_c_ratio);
            
        }
        
    };
    
    
    template<> void MinDMDMan<true,true,true,true>::init();
    template<> void MinDMDMan<false,true,true,true>::init();
    template<> void MinDMDMan<true,false,true,true>::init();
    template<> void MinDMDMan<false,false,true,true>::init();
    
    template<> void MinDMDMan<true,true,true,true>::add_extra_vec(VECTENS1&);
    template<> void MinDMDMan<false,true,true,true>::add_extra_vec(VECTENS1&);
    template<> void MinDMDMan<true,false,true,true>::add_extra_vec(VECTENS1&);
    template<> void MinDMDMan<false,false,true,true>::add_extra_vec(VECTENS1&);
    
    
    template<> void MinDMDMan<true,true,true,true>::prep();
    template<> void MinDMDMan<false,true,true,true>::prep();
    template<> void MinDMDMan<true,false,true,true>::prep();
    template<> void MinDMDMan<false,false,true,true>::prep();
    
    template<> void MinDMDMan<true,true,true,true>::force_calc();
    template<> void MinDMDMan<false,true,true,true>::force_calc();
    template<> void MinDMDMan<true,false,true,true>::force_calc();
    template<> void MinDMDMan<false,false,true,true>::force_calc();
    
    
    template<> type0 MinDMDMan<true,true,true,true>::F(type0);
    template<> type0 MinDMDMan<false,true,true,true>::F(type0);
    template<> type0 MinDMDMan<true,false,true,true>::F(type0);
    template<> type0 MinDMDMan<false,false,true,true>::F(type0);
}



    
    
namespace MAPP_NS
{
    class MinCG3DMD:public Min
    {
    private:
    protected:
        type0 max_dalpha;
        
        void prep();
        type0 calc_ndofs();
        type0 ndofs;
        
        VecTens<type0,3> h;
        VecTens<type0,3> x;
        VecTens<type0,3> x0;
        VecTens<type0,3> x_d;
        VecTens<type0,3> f;
        VecTens<type0,3> f0;
        
        class AtomsDMD* atoms;
        class ForceFieldDMD* ff;
        class DynamicDMD* dynamic;
        class ExportDMD* xprt;
        void print_error();

#ifdef NEW_UPDATE
#else
        vec* uvecs[3];
#endif
        void pre_run_chk(AtomsDMD*,ForceFieldDMD*);
        
    public:
        MinCG3DMD(type0,bool(&)[__dim__][__dim__],bool,type0,type0,class LineSearch*);
        ~MinCG3DMD();
        virtual void run(int);
        template<class C>
        void run(C*,int);
        virtual void init();
        virtual void fin();
        
        void force_calc();
        
        type0 F(type0);
        type0 dF(type0,type0&);
        void ls_prep(type0&,type0&,type0&);
        void F_reset();
        void refine(int,int);
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
        static void ml_refine(PyMethodDef&);
        
        static PyGetSetDef getset[];
        static void setup_tp_getset();
        static void getset_max_dalpha(PyGetSetDef&);
        static void getset_export(PyGetSetDef&);
        
        static int setup_tp();
        
        
        
    };

}
using namespace MAPP_NS;
/*--------------------------------------------
 
 --------------------------------------------*/
template<class C>
void MinCG3DMD::run(C* ls,int nsteps)
{
    int step=atoms->step;
    
    force_calc();
    
    int nevery_xprt=xprt==NULL ? 0:xprt->nevery;
    if(nevery_xprt) xprt->write(step);
    
    ThermoDynamics thermo(6,
    "PE",atoms->fe,
    "S[0][0]",atoms->S_fe[0][0],
    "S[1][1]",atoms->S_fe[1][1],
    "S[2][2]",atoms->S_fe[2][2],
    "S[1][2]",atoms->S_fe[2][1],
    "S[2][0]",atoms->S_fe[2][0],
    "S[0][1]",atoms->S_fe[1][0]);
    
    if(ntally) thermo.init();
    if(ntally) thermo.print(step);
    
    
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
