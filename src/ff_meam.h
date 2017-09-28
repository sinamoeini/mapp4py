#ifndef __MAPP__ff_meam__
#define __MAPP__ff_meam__
#include "ff_md.h"
struct _object;
typedef _object PyObject;
namespace MAPP_NS
{
    class ForceFieldMEAM: public ForceFieldMD
    {
    private:        
    protected:
        void force_calc();
        void energy_calc();
        void pre_xchng_energy(GCMC*);
        type0 xchng_energy(GCMC*);
        void post_xchng_energy(GCMC*);
    public:
        ForceFieldMEAM(class AtomsMD*);
        ~ForceFieldMEAM();
        
        static void ml_new(PyMethodDef&);

        void init();
        void fin();
        void init_xchng();
        void fin_xchng();

    };
}




/*--------------------------------------------
 
 --------------------------------------------*/
#include "xmath.h"
namespace MAPP_NS
{
    namespace MEAMMath
    {
        
        template<const int d>
        inline void d_fij(
        const type0,const type0,
        const type0,const type0,
        const type0* RESTRICT,const type0* RESTRICT,
        const type0* RESTRICT,type0* RESTRICT,
        type0&,type0&){};
        

        
        template<>
        inline void d_fij<1>(
        const type0 prefac_i,const type0 prefac_j,
        const type0 fpair_prefac_i,const type0 fpair_prefac_j,
        const type0* RESTRICT A_i,const type0* RESTRICT A_j,
        const type0* RESTRICT R_ij,type0* RESTRICT F_ij,
        type0& spair,type0& fpair)
        {
            F_ij[0]+=prefac_i*A_i[0]-prefac_j*A_j[0];
            F_ij[1]+=prefac_i*A_i[1]-prefac_j*A_j[1];
            F_ij[2]+=prefac_i*A_i[2]-prefac_j*A_j[2];
            
            type0 m_i=prefac_i*(A_i[0]*R_ij[0]+A_i[1]*R_ij[1]+A_i[2]*R_ij[2]);
            type0 m_j=-prefac_j*(A_j[0]*R_ij[0]+A_j[1]*R_ij[1]+A_j[2]*R_ij[2]);
            fpair+=fpair_prefac_i*m_i+fpair_prefac_j*m_j;
            spair+=m_i+m_j;
            
            
        }
        
        template<>
        inline void d_fij<2>(
        const type0 prefac_i,const type0 prefac_j,
        const type0 fpair_prefac_i,const type0 fpair_prefac_j,
        const type0* RESTRICT A_i,const type0* RESTRICT A_j,
        const type0* RESTRICT R_ij,type0* RESTRICT F_ij,
        type0& spair,type0& fpair)
        {
            type0 m_i=0.0,m_j=0.0;
            type0 y_ij,y_ji;
            y_ij=prefac_i*(A_i[0]*R_ij[0]+A_i[1]*R_ij[1]+A_i[2]*R_ij[2]);
            y_ji=-prefac_j*(A_j[0]*R_ij[0]+A_j[1]*R_ij[1]+A_j[2]*R_ij[2]);
            F_ij[0]+=2.0*(y_ij-y_ji);
            m_i+=y_ij*R_ij[0];
            m_j-=y_ji*R_ij[0];
            
            
            y_ij=prefac_i*(A_i[1]*R_ij[0]+A_i[3]*R_ij[1]+A_i[4]*R_ij[2]);
            y_ji=-prefac_j*(A_j[1]*R_ij[0]+A_j[3]*R_ij[1]+A_j[4]*R_ij[2]);
            F_ij[1]+=2.0*(y_ij-y_ji);
            m_i+=y_ij*R_ij[1];
            m_j-=y_ji*R_ij[1];
            
            
            y_ij=prefac_i*(A_i[2]*R_ij[0]+A_i[4]*R_ij[1]+A_i[5]*R_ij[2]);
            y_ji=-prefac_j*(A_j[2]*R_ij[0]+A_j[4]*R_ij[1]+A_j[5]*R_ij[2]);
            F_ij[2]+=2.0*(y_ij-y_ji);
            m_i+=y_ij*R_ij[2];
            m_j-=y_ji*R_ij[2];
            
            
            fpair+=fpair_prefac_i*m_i+fpair_prefac_j*m_j;
            spair+=m_i+m_j;
        }
        
        template<>
        inline void d_fij<3>(
        const type0 prefac_i,const type0 prefac_j,
        const type0 fpair_prefac_i,const type0 fpair_prefac_j,
        const type0* RESTRICT A_i,const type0* RESTRICT A_j,
        const type0* RESTRICT R_ij,type0* RESTRICT F_ij,
        type0& spair,type0& fpair)
        {
            type0 m_i=0.0,m_j=0.0;
            type0 y_ij,y_ji;
            /*
             A[0]*r[0]*r[0]+
             A[1]*r[0]*r[1]*2.0+
             A[2]*r[0]*r[2]*2.0+
             A[3]*r[1]*r[1]+
             A[4]*r[1]*r[2]*2.0+
             A[5]*r[2]*r[2];
             */
            y_ij=prefac_i*(
            A_i[0]*R_ij[0]*R_ij[0]+
            A_i[1]*R_ij[0]*R_ij[1]*2.0+
            A_i[2]*R_ij[0]*R_ij[2]*2.0+
            A_i[3]*R_ij[1]*R_ij[1]+
            A_i[4]*R_ij[1]*R_ij[2]*2.0+
            A_i[5]*R_ij[2]*R_ij[2]
            );
            
            y_ji=prefac_j*(
            A_j[0]*R_ij[0]*R_ij[0]+
            A_j[1]*R_ij[0]*R_ij[1]*2.0+
            A_j[2]*R_ij[0]*R_ij[2]*2.0+
            A_j[3]*R_ij[1]*R_ij[1]+
            A_j[4]*R_ij[1]*R_ij[2]*2.0+
            A_j[5]*R_ij[2]*R_ij[2]
            );
            F_ij[0]+=3.0*(y_ij-y_ji);
            m_i+=y_ij*R_ij[0];
            m_j-=y_ji*R_ij[0];
            
            
            
            /*
             A[1]*r[0]*r[0]+
             A[3]*r[0]*r[1]*2.0+
             A[4]*r[0]*r[2]*2.0+
             A[6]*r[1]*r[1]+
             A[7]*r[1]*r[2]*2.0+
             A[8]*r[2]*r[2];
             */
            
            y_ij=prefac_i*(
            A_i[1]*R_ij[0]*R_ij[0]+
            A_i[3]*R_ij[0]*R_ij[1]*2.0+
            A_i[4]*R_ij[0]*R_ij[2]*2.0+
            A_i[6]*R_ij[1]*R_ij[1]+
            A_i[7]*R_ij[1]*R_ij[2]*2.0+
            A_i[8]*R_ij[2]*R_ij[2]
            );
            y_ji=prefac_j*(
            A_j[1]*R_ij[0]*R_ij[0]+
            A_j[3]*R_ij[0]*R_ij[1]*2.0+
            A_j[4]*R_ij[0]*R_ij[2]*2.0+
            A_j[6]*R_ij[1]*R_ij[1]+
            A_j[7]*R_ij[1]*R_ij[2]*2.0+
            A_j[8]*R_ij[2]*R_ij[2]
            );
            F_ij[1]+=3.0*(y_ij-y_ji);
            m_i+=y_ij*R_ij[1];
            m_j-=y_ji*R_ij[1];
            
            /*
             A[2]*r[0]*r[0]+
             A[4]*r[0]*r[1]*2.0+
             A[5]*r[0]*r[2]*2.0+
             A[7]*r[1]*r[1]+
             A[8]*r[1]*r[2]*2.0+
             A[9]*r[2]*r[2];
             */
            y_ij=prefac_i*(
            A_i[2]*R_ij[0]*R_ij[0]+
            A_i[4]*R_ij[0]*R_ij[1]*2.0+
            A_i[5]*R_ij[0]*R_ij[2]*2.0+
            A_i[6]*R_ij[1]*R_ij[1]+
            A_i[8]*R_ij[1]*R_ij[2]*2.0+
            A_i[9]*R_ij[2]*R_ij[2]
            );
            y_ji=prefac_j*(
            A_j[2]*R_ij[0]*R_ij[0]+
            A_j[4]*R_ij[0]*R_ij[1]*2.0+
            A_j[5]*R_ij[0]*R_ij[2]*2.0+
            A_j[7]*R_ij[1]*R_ij[1]+
            A_j[8]*R_ij[1]*R_ij[2]*2.0+
            A_j[9]*R_ij[2]*R_ij[2]
            );
            F_ij[2]+=3.0*(y_ij-y_ji);
            m_i+=y_ij*R_ij[2];
            m_j-=y_ji*R_ij[2];
            
            
            fpair+=fpair_prefac_i*m_i+fpair_prefac_j*m_j;
            spair+=m_i+m_j;
        }
        
        template<>
        inline void d_fij<0>(
        const type0 prefac_i,const type0 prefac_j,
        const type0 fpair_prefac_i,const type0 fpair_prefac_j,
        const type0* RESTRICT A_i,const type0* RESTRICT A_j,
        const type0* RESTRICT R_ij,type0* RESTRICT F_ij,
        type0& spair,type0& fpair)
        {
            type0 m_i=prefac_i**A_i;
            type0 m_j=prefac_j**A_j;
            fpair+=fpair_prefac_i*m_i+fpair_prefac_j*m_j;
            spair+=m_i+m_j;
            
            
        }
        
        
        
        enum{RHO_OFFSET=0,Y1_OFFSET=4,Y2_OFFSET=7,W2_OFFSET=13,Y3_OFFSET=14,W3_OFFSET=24};
        
        inline void force_part(
        type0* RESTRICT beta_ov_r0_i,type0* RESTRICT beta_ov_r0_j,
        type0* RESTRICT t0_i,type0* RESTRICT t0_j,
        type0* RESTRICT s_i,type0* RESTRICT s_j,
        type0 prefac_i,type0 prefac_j,
        type0 s_prefac_i,type0 s_prefac_j,
        type0 rho_prefac_i,type0 rho_prefac_j,
        type0* RESTRICT rhoa_i,type0* RESTRICT rhoa_j,
        type0* RESTRICT t_i,type0* RESTRICT t_j,
        type0* RESTRICT rho_i,type0* RESTRICT rho_j,
        type0* RESTRICT R_ij,type0* RESTRICT F_ij,
        type0 phi_ij,type0 dphi_ij,
        type0 r_ij,type0 S_ij,type0 dS_ij_ij)
        {
            type0 dF_ij[3]={[0 ... 2]=0.0};
            type0 fpair=0.0;
            type0 spair=0.0;
            type0 r_ij_inv=1.0/r_ij;
            type0 r_ij_inv_sq=r_ij_inv*r_ij_inv;
            type0 r_ij_inv_cb=r_ij_inv*r_ij_inv_sq;
            
            d_fij<1>(rhoa_j[1]*prefac_i*t_i[0]*2.0*r_ij_inv,
                     rhoa_i[1]*prefac_j*t_j[0]*2.0*r_ij_inv,
                     
                     beta_ov_r0_j[1]+r_ij_inv,
                     beta_ov_r0_i[1]+r_ij_inv,
                     
                     rho_i+Y1_OFFSET,
                     rho_j+Y1_OFFSET,
                     
                     R_ij,dF_ij,spair,fpair);
            
            d_fij<2>(rhoa_j[2]*prefac_i*t_i[1]*2.0*r_ij_inv_sq,
                     rhoa_i[2]*prefac_j*t_j[1]*2.0*r_ij_inv_sq,
                     
                     beta_ov_r0_j[2]+2.0*r_ij_inv,
                     beta_ov_r0_i[2]+2.0*r_ij_inv,
                     
                     rho_i+Y2_OFFSET,
                     rho_j+Y2_OFFSET,
                     
                     R_ij,dF_ij,spair,fpair);
            
            d_fij<3>(rhoa_j[3]*prefac_i*t_i[2]*2.0*r_ij_inv_cb,
                     rhoa_i[3]*prefac_j*t_j[2]*2.0*r_ij_inv_cb,
                     
                     beta_ov_r0_j[3]+3.0*r_ij_inv,
                     beta_ov_r0_i[3]+3.0*r_ij_inv,
                     
                     rho_i+Y3_OFFSET,
                     rho_j+Y3_OFFSET,
                     
                     R_ij,dF_ij,spair,fpair);
            
            
            d_fij<0>(-rhoa_j[2]*prefac_i*t_i[1]*2.0/3.0,
                     -rhoa_i[2]*prefac_j*t_j[1]*2.0/3.0,
                     
                     beta_ov_r0_j[2],
                     beta_ov_r0_i[2],
                     
                     rho_i+W2_OFFSET,
                     rho_j+W2_OFFSET,
                     
                     R_ij,dF_ij,spair,fpair);
            d_fij<1>(-rhoa_j[3]*prefac_i*t_i[2]*0.4*r_ij_inv,
                     -rhoa_i[3]*prefac_j*t_j[2]*0.4*r_ij_inv,
                     
                     beta_ov_r0_j[3]+r_ij_inv,
                     beta_ov_r0_i[3]+r_ij_inv,
                     
                     rho_i+W3_OFFSET,
                     rho_j+W3_OFFSET,
                     
                     R_ij,dF_ij,spair,fpair);
            

            fpair*=-S_ij;
            
            
            
            
            
            
            
            
            type0 mm_i=rhoa_j[0]*(rho_prefac_i+
            (t0_j[0]-t_i[0])*(prefac_i*rho_i[1]*rho_i[1]-s_prefac_i*s_i[0])+
            (t0_j[1]-t_i[1])*(prefac_i*rho_i[2]*rho_i[2]-s_prefac_i*s_i[1])+
            (t0_j[2]-t_i[2])*(prefac_i*rho_i[3]*rho_i[3]-s_prefac_i*s_i[2]));
        
            type0 mm_j=rhoa_i[0]*(rho_prefac_j+
            (t0_i[0]-t_j[0])*(prefac_j*rho_j[1]*rho_j[1]-s_prefac_j*s_j[0])+
            (t0_i[1]-t_j[1])*(prefac_j*rho_j[2]*rho_j[2]-s_prefac_j*s_j[1])+
            (t0_i[2]-t_j[2])*(prefac_j*rho_j[3]*rho_j[3]-s_prefac_j*s_j[2]));
            
            
            
            
            
            spair+=mm_i+mm_j;
            fpair-=S_ij*(beta_ov_r0_j[0]*mm_i+beta_ov_r0_i[0]*mm_j);
            
            spair+=phi_ij;
            fpair+=dphi_ij*S_ij;
            
            fpair+=dS_ij_ij*spair;
            
            dF_ij[0]*=S_ij;
            dF_ij[0]+=fpair*r_ij_inv*R_ij[0];
            dF_ij[1]*=S_ij;
            dF_ij[1]+=fpair*r_ij_inv*R_ij[1];
            dF_ij[2]*=S_ij;
            dF_ij[2]+=fpair*r_ij_inv*R_ij[2];
            
            
            
        }
        
    }
}
#endif
