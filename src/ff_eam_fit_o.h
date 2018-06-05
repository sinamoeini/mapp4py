#ifndef __MAPP__ff_eam_fit_o__
#define __MAPP__ff_eam_fit_o__
#include "ff_md.h"
#include "xmath.h"
#define nrho_H 2
#ifdef FeH_SPLINE
#define nphi_FeH 18
#else
#define nphi_FeH 4
#endif
#define nphi_HH 4
#define nF_H 5
#define rho_H_offset 0
#define phi_FeH_offset rho_H_offset+nrho_H
#define phi_HH_offset phi_FeH_offset+nphi_FeH
#define F_H_offset phi_HH_offset+nphi_HH
#define nHvars nrho_H+nphi_FeH+nphi_HH+nF_H
//A_rho_H=v;
//A_phi_FeH=v+2;
//A_phi_HH=v+5;
//A_F_H=v+9;


namespace MAPP_NS
{
    class ForceFieldEAMFitO: public ForceFieldMD
    {
    private:
        Vec<type0>* rho_ptr;
        //Vec<type0>* S_ptr;
    protected:
        
        void force_calc();
        
        
        void energy_calc();
        void pre_xchng_energy(GCMC*);
        type0 xchng_energy(GCMC*);
        void post_xchng_energy(GCMC*);
    public:
        ForceFieldEAMFitO(class AtomsMD*,type0(&)[nrho_H],type0(&)[nphi_FeH],type0(&)[nphi_HH],type0(&)[nF_H]);
        ~ForceFieldEAMFitO();
        void init();
        void fin();
        void init_xchng();
        void fin_xchng();
        
        static type0 rc_rho_H;
        static type0 rc_rho_Fe;
        static type0 rc_phi_FeFe;
        static type0 rc_phi_FeH;
        static type0 rc_phi_HH;
        
        
        type0 v[nHvars];
        type0* A_rho_H;
        type0* A_phi_FeH;
        type0* A_phi_HH;
        type0* A_F_H;
        
        type0 dv[nHvars];
        type0* dA_rho_H;
        type0* dA_phi_FeH;
        type0* dA_phi_HH;
        type0* dA_F_H;
        
        
        type0 dv_lcl[nHvars];
        type0* dA_rho_H_lcl;
        type0* dA_phi_FeH_lcl;
        type0* dA_phi_HH_lcl;
        type0* dA_F_H_lcl;
        
        bool dof[nHvars];
        
        
        static type0 A_F_Fe[3];
        static type0 AR_rho_Fe[3][2];
        static type0 AR_phi_FeFe[13][2];
#ifdef FeH_SPLINE
        static type0 R_phi_FeH[nphi_FeH];
#endif
        static type0 B_phi_FeFe[4];
        static type0 r0_FeFe,r1_FeFe;
        
        
        
        
        
        type0 E_HH(type0);
        type0 dE_HH(type0);
        type0 E_HH_(type0);
        type0 dE_HH_(type0);
        type0 calc_F(elem_type,type0);
        type0 calc_dF(elem_type,type0);
        type0 calc_ddF(elem_type,type0);
        type0 calc_rho(elem_type,type0);
        type0 calc_drho(elem_type,type0);
        type0 calc_phi(elem_type,elem_type,type0);
        type0 calc_dphi(elem_type,elem_type,type0);
        
        
        void calc_Dphi_FeH(type0,bool);
        void calc_Dphi_HH(type0,bool);
        void calc_DF_H(type0);
        void calc_DF_Fe(type0);
        void calc_Drho_H(type0,type0);
        void calc_Drho_Fe(type0,type0);
        
        void gradient();
        type0 mean_rho_H();
        size_t get_rFeH(type0*&,int*&);
        
        
        
        
        
        
        
        
        template<const int N>
        type0 spline(type0 (&AR)[N][2],type0 r)
        {
            int i=N-1;
            type0 ans=0.0;
            while(i>-1 && r<AR[i][1])
            {
                ans+=AR[i][0]*Algebra::pow<3>(AR[i][1]-r);
                --i;
            }
            return ans;
            
        }
        template<const int N>
        type0 dspline(type0 (&AR)[N][2],type0 r)
        {
            int i=N-1;
            type0 ans=0.0;
            while(i>-1 && r<AR[i][1])
            {
                ans-=3.0*AR[i][0]*Algebra::pow<2>(AR[i][1]-r);
                --i;
            }
            return ans;
            
        }
        
#ifdef FeH_SPLINE
        template<const int N>
        type0 spline(type0*& A,type0 (&R)[N],type0 r)
        {
            int i=N-1;
            type0 ans=0.0;
            while(i>-1 && r<R[i])
            {
                ans+=A[i]*Algebra::pow<3>(R[i]-r);
                --i;
            }
            return ans;
            
        }
        template<const int N>
        type0 dspline(type0*& A,type0 (&R)[N],type0 r)
        {
            int i=N-1;
            type0 ans=0.0;
            while(i>-1 && r<R[i])
            {
                ans-=3.0*A[i]*Algebra::pow<2>(R[i]-r);
                --i;
            }
            return ans;
        }
        
        template<const int N>
        void Dspline(type0& coef,type0*& D,type0 (&R)[N],type0 r)
        {
            int i=N-1;
            while(i>-1 && r<R[i])
            {
                D[i]+=coef*Algebra::pow<3>(R[i]-r);
                --i;
            }
            
        }
#endif
        
        
        
        
        /*
         1769.6687067528123, 4963.443749027828, 2727.5091948962486, 274.2110421849655
         28.571747468384945, 8.41348676233098, 3.597361579691342, 1.8000200905082515
         */
        type0 xi(type0 r)
        {
            return (
            1769.6687067528123*exp(-28.571747468384945*r)+
            4963.443749027828*exp(-8.41348676233098*r)+
            2727.5091948962486*exp(-3.597361579691342*r)+
            274.2110421849655*exp(-1.8000200905082515*r)
                    )/r;
        }
        type0 dxi(type0 r)
        {
            return -(
                    1769.6687067528123*(28.571747468384945*r+1.0)*exp(-28.571747468384945*r)+
                    4963.443749027828*(8.41348676233098*r+1.0)*exp(-8.41348676233098*r)+
                    2727.5091948962486*(3.597361579691342*r+1.0)*exp(-3.597361579691342*r)+
                    274.2110421849655*(1.8000200905082515*r+1.0)*exp(-1.8000200905082515*r)
                    )/(r*r);
        }
        
        
        
        
        static void ml_new(PyMethodDef&);
        
        
    };
    
    
}


#endif
