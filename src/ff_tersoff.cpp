/*--------------------------------------------
 Created by Sina on 07/15/13.
 Copyright (c) 2013 MIT. All rights reserved.
 --------------------------------------------*/
#include "ff_tersoff.h"
#include "neighbor_md.h"
#include "atoms_md.h"
#include "elements.h"
#include "memory.h"
#include "pgcmc.h"
#include "api.h"
using namespace MAPP_NS;
/*--------------------------------------------
 constructor
 --------------------------------------------*/
ForceFieldTersoff::
ForceFieldTersoff(AtomsMD* atoms):
ForceFieldMD(atoms)
{
    gcmc_n_cutoff=1;
    gcmc_n_vars=1;
    gcmc_tag_enabled=false;
    neighbor->pair_wise=false;
    
}
/*--------------------------------------------
 destructor
 --------------------------------------------*/
ForceFieldTersoff::~ForceFieldTersoff()
{
}
/*--------------------------------------------
 initiate before a run
 --------------------------------------------*/
void ForceFieldTersoff::init()
{
    pre_init();
}
/*--------------------------------------------
 after a run
 --------------------------------------------*/
void ForceFieldTersoff::fin()
{
    post_fin();
}
/*--------------------------------------------
 init xchng
 --------------------------------------------*/
void ForceFieldTersoff::init_xchng()
{

}
/*--------------------------------------------
 fin xchng
 --------------------------------------------*/
void ForceFieldTersoff::fin_xchng()
{

}
/*--------------------------------------------
 pre xchng energy
 --------------------------------------------*/
void ForceFieldTersoff::pre_xchng_energy(GCMC* gcmc)
{
    
}
/*--------------------------------------------
 xchng energy
 --------------------------------------------*/
type0 ForceFieldTersoff::xchng_energy(GCMC* gcmc)
{
    return 0.0;
}
/*--------------------------------------------
 post xchng energy
 --------------------------------------------*/
void ForceFieldTersoff::post_xchng_energy(GCMC*)
{
}
/*--------------------------------------------
 force and energy calculation
 --------------------------------------------*/
#include "xmath.h"
void ForceFieldTersoff::__force_calc()
{
    /*
    const type0* x=atoms->x->begin();
    //type0* fvec=f->begin();
    elem_type* evec=atoms->elem->begin();
    
    elem_type ielem,jelem,kelem;
    type0 rsq_ij;
    
    
    
    int** neighbor_list=neighbor->neighbor_list;
    int* neighbor_list_size=neighbor->neighbor_list_size;
    
    type0 x_i[__dim__];
    type0 dx_ij[__dim__];
    const int natms_lcl=atoms->natms_lcl;
    for(int iatm=0;iatm<natms_lcl;iatm++)
    {
        
        ielem=evec[iatm];
        Algebra::V_eq<__dim__>(x+iatm*__dim__,x_i);
        //type0 f_i[__dim__]{DESIG(__dim__,0.0)};
        const int list_size=neighbor_list_size[iatm];
        for(int j=0,jatm;j<list_size;j++)
        {
            jatm=neighbor_list[iatm][j];
            if(jatm<iatm) continue;
            jelem=evec[jatm];
            rsq_ij=Algebra::DX_RSQ(x_i,x+jatm*__dim__,dx_ij);
        }
    }
    
    
    
    
    
    
    
    
    
    
    
    
    
    type0 r_ij,r_ik;
    type0 r_ij_hat[__dim__];
    type0 r_ik_hat[__dim__];
    
    type0** lambda_2_i;
    type0* lambda_2_ij;
    type0 cos_theta_ijk;
    type0 zeta_ij;
    
    for(int iatm=0;iatm<natms_lcl;iatm++)
    {
        
        ielem=evec[iatm];
        lambda_2_i=lambda_2[ielem];
        Algebra::V_eq<__dim__>(x+iatm*__dim__,x_i);
        //type0 f_i[__dim__]{DESIG(__dim__,0.0)};
        const int list_size=neighbor_list_size[iatm];
        for(int j=0,jatm;j<list_size;j++)
        {
            jatm=neighbor_list[iatm][j];
            jelem=evec[jatm];
            r_ij=Algebra::DX_HAT_R(x_i,x+jatm*__dim__,r_ij_hat);
            lambda_2_ij=lambda_2_i[jelem];
            
            
            
            zeta_ij=0.0;
            for(int k=0,katm;k<list_size;k++)
            {
                if(k==j) continue;
                katm=neighbor_list[iatm][k];
                kelem=evec[katm];
                r_ik=Algebra::DX_HAT_R(x_i,x+katm*__dim__,r_ik_hat);
                
                cos_theta_ijk=Algebra::V_mul_V<__dim__>(r_ij_hat,r_ik_hat);
                //zeta_ij+=zeta_ij_k(r_ij,r_ik,cos_theta_ijk);
                
            }
            
            
            
        }
    }
    
    
    
    
    */
    
}
/*--------------------------------------------

 --------------------------------------------*/
type0 ForceFieldTersoff::zeta_ij_k(
const type0& rc_inner,const type0& rc_outter,
const type0& lambda_2,const type0& m,
const type0& gamma,const type0& c_sq_inv,const type0&d_sq,const type0& cos_theta0,
const type0& r_ij,const type0& r_ik,const type0& cos_theta_ijk,
const type0*& r_ij_hat,const type0*& r_ik_hat)
{
    
    /*
    type0 r0=lambda_2*(r_ij-r_ik);
    type0 r1=pow(r0,m);
    //type0 r2=(r1/r0)*lambda_2*m;
    r1=exp(r1);
    
    type0 t=cos_theta_ijk-cos_theta0;
    t*=t;
    t+=d_sq;
    type0 g=(c_sq_inv+1.0/t);
   // type0 d_log_g=-2.0*cos_theta_ijk/(t*t*g);
    g*=gamma;
    
    type0 fc=1.0;
    type0 d_log_fc=0.0;
    
    if(r_ik>rc_inner)
    {
        t=M_PI/(rc_outter-rc_inner);
        type0 s=(r_ik-0.5*(rc_outter+rc_inner))*t;
        fc=0.5-0.5*sin(s);
        d_log_fc=-0.5*t*cos(s)/fc;
    }
    
    
    
    return fc*r1*g;
     */
    return 0.0;
}
/*--------------------------------------------
 only energy calculation this is useful for
 minimization/linesearch methods that do not
 use derivatives of energy
 --------------------------------------------*/
void ForceFieldTersoff::__energy_calc()
{
}
/*--------------------------------------------
 python constructor
 --------------------------------------------*/
void ForceFieldTersoff::ml_new(PyMethodDef& tp_methods)
{

}
