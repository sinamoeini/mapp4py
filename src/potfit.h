#ifndef __MAPP__potfit__
#define __MAPP__potfit__
#include "atoms_md.h"
#include "min_cg_fit.h"
#include "ff_eam_fit.h"
#include "thermo_dynamics.h"
namespace MAPP_NS
{
    class PotFit
    {
    private:
        class HydrideConst* hyd_const;
        int ntrial;
        static const char* err_msgs[];
        void sort_AR_ij(type0*&,type0*&,bool*&,size_t);
        
        
        int quadratic(type0 b,type0 c,type0* ans)
        {
            type0 delta_sq=b*b-4.0*c;
            if(delta_sq<=0.0) return 0;
            
            type0 delta=sqrt(delta_sq);
            ans[0]=0.5*(-b+delta);
            ans[1]=0.5*(-b-delta);
            return 2;
        }
        
        int cubic(type0 b,type0 c,type0 d,type0* ans)
        {
            type0 p=c-b*b/3.0;
            type0 q=2.0*b*b*b/27.0-b*c/3.0+d;
            
            if(p==0.0) return pow(q,1.0/3.0);
            if(q==0.0) return 0.0;
            
            type0 t=sqrt(fabs(p)/3.0);
            type0 g=1.5*q/(p*t);
            if(p>0.0)
            {
                ans[0]=-2.0*t*sinh(asinh(g)/3.0)-b/3.0;
                return 1;
            }
            
            
            if(4.0*p*p*p+27.0*q*q<0.0)
            {
                type0 theta=acos(g)/3.0;
                ans[0]=2.0*t*cos(theta)-b/3.0;
                ans[1]=2.0*t*cos(theta+2.0*M_PI/3.0)-b/3.0;
                ans[2]=2.0*t*cos(theta-2.0*M_PI/3.0)-b/3.0;
                return 3;
            }
            
            if(q>0.0)
            {
                ans[0]=-2.0*t*cosh(acosh(-g)/3.0)-b/3.0;
                return 1;
            }
            
            ans[0]=2.0*t*cosh(acosh(g)/3.0)-b/3.0;
            return 1;
        }
        
        type0 cubic(type0 b,type0 c,type0 d)
        {
            type0 p=c-b*b/3.0;
            type0 q=2.0*b*b*b/27.0-b*c/3.0+d;
            
            if(p==0.0) return pow(q,1.0/3.0);
            if(q==0.0) return 0.0;
            
            type0 t=sqrt(fabs(p)/3.0);
            type0 g=1.5*q/(p*t);
            if(p>0.0)
                return -2.0*t*sinh(asinh(g)/3.0)-b/3.0;
            
            
            if(4.0*p*p*p+27.0*q*q<0.0)
                return 2.0*t*cos(acos(g)/3.0)-b/3.0;
            
            if(q>0.0)
                return -2.0*t*cosh(acosh(-g)/3.0)-b/3.0;
            
            return 2.0*t*cosh(acosh(g)/3.0)-b/3.0;
        }
        
        int quartic(type0 b,type0 c,type0 d,type0 e,type0* ans)
        {
            
            type0 p=c-0.375*b*b;
            type0 q=0.125*b*b*b-0.5*b*c+d;
            if(q==0)
            {
                type0 m=0.25*p*p+0.01171875*b*b*b*b-e+0.25*b*d-0.0625*b*b*c;
                if(m<0.0) return 0;
                int nroots=0;
                
                if(m-0.5*p>=0.0)
                {
                    type0 absy=sqrt(m-0.5*p);
                    ans[nroots++]=absy-0.25*b;
                    ans[nroots++]=-absy-0.25*b;
                }
                if(m+0.5*p>=0.0)
                {
                    type0 absy=sqrt(m+0.5*p);
                    ans[nroots++]=absy-0.25*b;
                    ans[nroots++]=-absy-0.25*b;
                }
                
                return nroots;
            }
            type0 m=cubic(p,0.25*p*p+0.01171875*b*b*b*b-e+0.25*b*d-0.0625*b*b*c,-0.125*q*q);
            
            if(m<0.0) return 0;
            type0 sqrt_2m=sqrt(2.0*m);
            int nroots=0;
            if(-m-p+q/sqrt_2m>=0.0)
            {
                type0 delta=sqrt(2.0*(-m-p+q/sqrt_2m));
                ans[nroots++]=0.5*(-sqrt_2m+delta)-0.25*b;
                ans[nroots++]=0.5*(-sqrt_2m-delta)-0.25*b;
            }
            
            if(-m-p-q/sqrt_2m>=0.0)
            {
                type0 delta=sqrt(2.0*(-m-p-q/sqrt_2m));
                ans[nroots++]=0.5*(sqrt_2m+delta)-0.25*b;
                ans[nroots++]=0.5*(sqrt_2m-delta)-0.25*b;
            }
            
            return nroots;
        }
        void slpine_alpha_max(const type0& a,const type0& b,const type0& c,const type0& d,
                              const type0& da,const type0& db,const type0& dc,const type0& dd,
                              const type0& x_lo,const type0& x_hi,type0& alpha_max)
        {
            type0 f0=((a*x_lo+b)*x_lo+c)*x_lo+d;
            type0 df0=((da*x_lo+db)*x_lo+dc)*x_lo+dd;
            type0 f1=((a*x_hi+b)*x_hi+c)*x_hi+d;
            type0 df1=((da*x_hi+db)*x_hi+dc)*x_hi+dd;
            if(df0<0.0)
                alpha_max=MIN(alpha_max,-f0/df0);
            if(df1<0.0)
                alpha_max=MIN(alpha_max,-f1/df1);
            
            
            type0 A=a*db-b*da;
            type0 B=2.0*(a*dc-c*da);
            type0 C=3.0*(a*dd-d*da) + (b*dc-c*db);
            type0 D=2.0*(b*dd-d*db);
            type0 E=c*dd-d*dc;
            type0 ans[4];
            int n;
            if(A!=0.0)
            {
                B/=A;
                C/=A;
                D/=A;
                E/=A;
                n=quartic(B,C,D,E,ans);
                
                
            }
            else if(B!=0.0)
            {
                C/=B;
                D/=B;
                E/=B;
                n=cubic(C,D,E,ans);
            }
            else if(C!=0.0)
            {
                D/=C;
                E/=C;
                n=quadratic(D,E,ans);
            }
            else if(D!=0.0)
            {
                n=1;
                ans[0]=-E/D;
            }
            else
            {
                if(a*da<0.0)
                    alpha_max=MIN(alpha_max,-a/da);
                else
                    alpha_max=0.0;
                return;
            }
            
            
            type0 df,dg;
            for(int i=0;i<n;i++)
            {
                if(x_lo<ans[i] && ans[i]<x_hi)
                {
                    dg=(3.0*da*ans[i]+2.0*db)*ans[i]+dc;
                    df=(3.0*a*ans[i]+2.0*b)*ans[i]+c;
                    
                    if(dg!=0.0 && -df/dg >=0.0)
                        alpha_max=MIN(alpha_max,-df/dg);
                }
            }
            
            
            
            
        }
        void spline(type0*& A,type0*& dA,type0*& R,size_t& sz,type0& alpha_max)
        {
            type0 a=0.0,b=0.0,c=0.0,d=0.0;
            type0 da=0.0,db=0.0,dc=0.0,dd=0.0;
            type0 r_lo=0.0,r_hi=0.0;
            for(size_t i=0;i<sz;i++)
            {
                
                
                r_lo=R[0]-R[i];
                a+=A[i];
                b-=A[i]*3.0*r_lo;
                c+=A[i]*3.0*r_lo*r_lo;
                d-=A[i]*r_lo*r_lo*r_lo;
                
                da+=dA[i];
                db-=dA[i]*3.0*r_lo;
                dc+=dA[i]*3.0*r_lo*r_lo;
                dd-=dA[i]*r_lo*r_lo*r_lo;
                
                if(dA[i]!=0.0)
                {
                    size_t j=i+1;
                    while(j<sz && dA[j]==0.0)
                        ++j;
                    if(j==sz)
                        r_hi=R[0];
                    else
                        r_hi=R[0]-R[j];
                    
                    slpine_alpha_max(a,b,c,d,da,db,dc,dd,r_lo,r_hi,alpha_max);
                }
                
                
                
            }
        }
        
        void spline2(type0*& A,type0*& dA,type0*& R,size_t& sz,type0& alpha_max)
        {
            
            for(size_t i=0;i<sz;i++)
            {
                if(dA[i]<0.0)
                {
                    type0 max_a= -A[i]/dA[i];
                    while(A[i]+max_a*dA[i]<0.0)
                        max_a=nextafter(max_a,-1.0);
                    alpha_max=MIN(max_a,alpha_max);
                }
                
                
                
                
                
            }
        }
        
        
        
        type0 dot(type0*& a,type0*& b)
        {
            type0 ans=0.0;
            for(int i=0;i<nvars;i++) ans+=a[i]*b[i];
            return ans;
        }
        type0 val(type0 (&cur)[1+__nvoigt__],type0 (&target)[1+__nvoigt__],type0 (&coef)[1+__nvoigt__])
        {
            type0 ans=0.0;
            Algebra::Do<1+__nvoigt__>::func([&cur,&target,&coef,&ans](int i){ans+=coef[i]*(cur[i]-target[i])*(cur[i]-target[i]);});
            return ans;
        }
        type0 dval(type0 (&comp_deriv)[1+__nvoigt__],type0 (&cur)[1+__nvoigt__],type0 (&target)[1+__nvoigt__],type0 (&coef)[1+__nvoigt__])
        {
            type0 ans=0.0;
            Algebra::Do<1+__nvoigt__>::func([&cur,&target,&coef,&ans,&comp_deriv](int i)
            {ans+=2.0*coef[i]*(cur[i]-target[i])*comp_deriv[i];});
            return ans;
        }
        
        void calc_f(type0 (&cur)[1+__nvoigt__])
        {
            
            for(int i=0;i<nvars;i++) f_lcl[i]=0.0;
            
            size_t nelems=2;
            for(size_t ielem=0;ielem<nelems;ielem++)
                for(size_t jelem=0;jelem<nelems;jelem++)
                    for(size_t i=0;i<rho_sz[ielem][jelem];i++)
                        if(rho_A_dof[ielem][jelem][i])
                            rho_A_f[ielem][jelem][i]-=dval(drho_A[ielem][jelem][i],cur,target,coef);
            
            type0 sum;
            for(size_t ielem=0;ielem<nelems;ielem++)
                for(size_t i=0;i<rho_sz[ielem][0];i++)
                    if(rho_A_dof[ielem][0][i])
                    {
                        sum=0.0;
                        for(size_t jelem=0;jelem<nelems;jelem++)
                            sum+=rho_A_f[ielem][jelem][i];
                        for(size_t jelem=0;jelem<nelems;jelem++)
                            rho_A_f[ielem][jelem][i]=sum;
                    }
            
            
            for(size_t ielem=0;ielem<nelems;ielem++)
                for(size_t jelem=0;jelem<ielem+1;jelem++)
                    for(size_t i=0;i<phi_sz[ielem][jelem];i++)
                        if(phi_A_dof[ielem][jelem][i])
                        {
                            phi_A_f[ielem][jelem][i]-=dval(dphi_A[ielem][jelem][i],cur,target,coef);
                            phi_A_f[jelem][ielem][i]=phi_A_f[ielem][jelem][i];
                        }
            
            
            for(size_t ielem=0;ielem<nelems;ielem++)
                for(size_t i=0;i<3;i++)
                    if(F_A_dof[ielem][i])
                        F_A_f[ielem][i]-=dval(dF_A[ielem][i],cur,target,coef);
            
        }
        
        
        
    protected:
    public:
        PotFit(type0(***&&)[2],size_t**&&,
               type0(***&&)[2],size_t**&&,
               type0(*&&)[3],
               
               size_t*&&,size_t,size_t*&&,size_t
               ,type0*&&,type0**&&,type0*&&,type0**&&,type0**&&,type0**&&,type0***&&,
               
               
               std::string*&&,std::string*&,int*&,type0(*&)[1+__nvoigt__],type0(*&)[1+__nvoigt__],PyObject**,size_t,MPI_Comm&);
        
        ~PotFit();
        
        void alloc(type0(***&)[2],type0(***&)[2],type0(*&)[3]);
        void dealloc();
        
        int nconfigs;
        type0* errs;
        int* roots;
        std::string* names_str;
        const char** names;
        MPI_Comm world;
        MPI_Comm* my_world;
        int my_rank,my_conf,my_lcl_rank;
        AtomsMD* atoms;
        MinCGFit* min;
        LineSearchBrent* min_ls;
        ForceFieldEAMFit* ff;
        
        VecTens<type0,1> X0;
        VecTens<type0,1> Xorig;

        
        type0 find_max_alpha_rho();
        type0 find_max_alpha_phi();
        type0 find_max_alpha_F();
       
        
        type0 f_h,tot_err;
        
        type0 iter();
        void min_cg(int);
        void store_x0();
        void restore_x0();
        void full_reset();
        
        
        int max_ntrials;
        type0 max_rhobF;
        type0 max_AF;
        type0 max_alphaF;
        type0 max_Arho;
        type0 max_Aphi;
        
        

        
        int nvars;
        type0* f_lcl;
        type0* f;
        type0* f0;
        type0* h;
        type0* h0;
        type0* x0;
        type0* x;
        bool* dof;
        type0 (* dA)[1+__nvoigt__];
        
        type0 *** rho_R;
        type0 *** rho_A;
        type0 *** rho_A0;
        type0 *** rho_A_h;
        type0 *** rho_A_f;
        bool*** rho_A_dof;
        type0 (*** drho_A)[1+__nvoigt__];
        size_t** rho_sz;
        
        type0 *** phi_R;
        type0 *** phi_A;
        type0 *** phi_A0;
        type0 *** phi_A_h;
        type0 *** phi_A_f;
        bool*** phi_A_dof;
        type0 (*** dphi_A)[1+__nvoigt__];
        size_t** phi_sz;
        
        type0** F_A;
        type0** F_A0;
        type0** F_A_h;
        type0 **F_A_f;
        bool **F_A_dof;
        type0 (**dF_A)[1+__nvoigt__];
        
        
        
        
        
        type0 target[1+__nvoigt__];
        type0 coef[1+__nvoigt__];
        type0 curr_diff[1+__nvoigt__];
        
        type0 S_target[__dim__][__dim__];
        type0 pe_target;
        type0 S_coef[__dim__][__dim__];
        type0 pe_coef;
        
        ThermoDynamics get_thermo();
        type0 F(type0);
        void ls_prep(type0&,type0&,type0&);
        void F_reset();
        size_t get_rFeH(type0*&,int*&);
        
        
        typedef struct
        {
            PyObject_HEAD
            class PotFit* potfit;
        }Object;
        
        
        
        static PyTypeObject TypeObject;
        static PyObject* __new__(PyTypeObject*,PyObject*, PyObject*);
        static int __init__(PyObject*, PyObject*,PyObject*);
        static PyObject* __alloc__(PyTypeObject*,Py_ssize_t);
        static void __dealloc__(PyObject*);
        
        
        static PyGetSetDef getset[];
        static void setup_tp_getset();
        static void getset_rho_A(PyGetSetDef&);
        static void getset_rho_R(PyGetSetDef&);
        static void getset_rho_dof(PyGetSetDef&);
        static void getset_phi_A(PyGetSetDef&);
        static void getset_phi_R(PyGetSetDef&);
        static void getset_phi_dof(PyGetSetDef&);
        static void getset_F_A(PyGetSetDef&);
        static void getset_F_dof(PyGetSetDef&);
        static void getset_max_ntrials(PyGetSetDef&);
        static void getset_max_rhobF(PyGetSetDef&);
        static void getset_max_AF(PyGetSetDef&);
        static void getset_max_alphaF(PyGetSetDef&);
        static void getset_max_Arho(PyGetSetDef&);
        static void getset_max_Aphi(PyGetSetDef&);
        //static void getset_RFeH(PyGetSetDef&);
        
        
        static PyMethodDef methods[];
        static void setup_tp_methods();
        static void ml_reset(PyMethodDef&);
        static void ml_min(PyMethodDef&);
        
        
        static int setup_tp();
        
        
        

        


    };
    
    
}


#endif
