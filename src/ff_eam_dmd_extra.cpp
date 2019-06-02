#include "ff_eam_dmd.h"
#include <stdlib.h>
#include "neighbor_dmd.h"
#include "elements.h"
#include "memory.h"
#include "xmath.h"
#include "atoms_dmd.h"
#include "dynamic_dmd.h"
#include <limits>
#define PI_IN_SQ 0.564189583547756286948079451561
/*--------------------------------------------
 force calculation for grand potential
 --------------------------------------------*/
void ForceFieldEAMDMD::__force_calc_gp()
{
    if(max_pairs0<neighbor->no_pairs)
    {
        delete [] rho_phi;
        delete [] drho_phi_dr;
        delete [] drho_phi_dalpha;
        
        max_pairs0=neighbor->no_pairs;
        size_t no_0=max_pairs0*3;
        Memory::alloc(rho_phi,no_0);
        Memory::alloc(drho_phi_dr,no_0);
        Memory::alloc(drho_phi_dalpha,no_0);
    }
    
    for(size_t i=0;i<max_pairs0*3;i++) rho_phi[i]=drho_phi_dr[i]=drho_phi_dalpha[i]=0.0;
    
    
    
    
    
    type0 r,r_inv;
    size_t m;
    type0* coef;
    type0 tmp0,tmp1;
    type0 fpair,apair;
    
    type0 p,cv_i;
    
    type0 dx_ij[__dim__];
    
    type0 const* c=atoms->c->begin();
    
    
    
    type0* dE=dE_ptr->begin();
    type0* __f_c=f_c->begin();
    type0* rho=E_ptr->begin();
    const int n=atoms->natms_lcl*c_dim;
    for(int i=0;i<n;i++) rho[i]=0.0;
    
    elem_type const* elem_vec=atoms->elem->begin();
    
    type0 alpha_ij,rsq;
    elem_type elem_i,elem_j;
    
    
    
    const type0* x=atoms->x->begin();
    const type0* alpha=atoms->alpha->begin();
    int** neighbor_list=neighbor->neighbor_list;
    int* neighbor_list_size=neighbor->neighbor_list_size;
    size_t istart=0;
    for(int i=0;i<n;i++)
    {
        type0 c_i=c[i];
        if(c_i<0.0) continue;
        elem_i=elem_vec[i];
        const int neigh_sz=neighbor_list_size[i];
        for(int j,__j=0;__j<neigh_sz;__j++,istart+=3)
        {
            j=neighbor_list[i][__j];
            elem_j=elem_vec[j];
            rsq=Algebra::RSQ<__dim__>(x+(i/c_dim)*__dim__,x+(j/c_dim)*__dim__);
            r=sqrt(rsq);
            type0 alpha_ij_sq=alpha[i]*alpha[i]+alpha[j]*alpha[j];
            alpha_ij=sqrt(alpha_ij_sq);
            if(r-alpha_ij*xi[N-1]>=cut[elem_i][elem_j]) continue;
            type0 upper=(r+cut[elem_i][elem_j])/alpha_ij;
            type0 lower=(r-cut[elem_i][elem_j])/alpha_ij;
            type0 __r,p,tmp0,xi2;
            type0* coef;
            
            
            type0 H[3][5]{DESIG2(3,5,0.0)};
            for(int l=0;l<N;l++)
            {
                if(xi[l]<=lower && xi[l]>=upper) continue;
                
                xi2=xi[l]*xi[l];
                
                __r=r-xi[l]*alpha_ij;
                p=fabs(__r)*dr_inv;
                m=static_cast<size_t>(p);
                m=MIN(m,nr-2);
                p-=m;
                p=MIN(p,1.0);
                
                coef=r_phi_arr[elem_i][elem_j][m];
                tmp0=((coef[3]*p+coef[2])*p+coef[1])*p+coef[0];
                if(__r<0.0) tmp0*=-1.0;
                tmp0*=wi_0[l];
                H[0][0]+=tmp0;
                H[0][1]+=tmp0*xi[l];
                H[0][2]+=tmp0*xi2;
                
                
                coef=r_rho_arr[elem_i][elem_j][m];
                tmp0=((coef[3]*p+coef[2])*p+coef[1])*p+coef[0];
                if(__r<0.0) tmp0*=-1.0;
                tmp0*=wi_0[l];
                H[1][0]+=tmp0;
                H[1][1]+=tmp0*xi[l];
                H[1][2]+=tmp0*xi2;
                
                coef=r_rho_arr[elem_j][elem_i][m];
                tmp0=((coef[3]*p+coef[2])*p+coef[1])*p+coef[0];
                if(__r<0.0) tmp0*=-1.0;
                tmp0*=wi_0[l];
                H[2][0]+=tmp0;
                H[2][1]+=tmp0*xi[l];
                H[2][2]+=tmp0*xi2;
            }
            
            r_inv=1.0/r;
            type0 r2_inv=r_inv*r_inv;
            type0 alpha_inv=1.0/alpha_ij;
            type0 alpha2_inv=alpha_inv*alpha_inv;
            type0 alpha4_inv=alpha2_inv*alpha2_inv;
            type0 r2alpha2_inv=r2_inv*alpha2_inv;
            
            
            Algebra::Do<3>::func([&istart,&H,&r_inv,&r2_inv,&alpha2_inv,&r2alpha2_inv,&alpha_ij_sq,&alpha_inv,&rsq,&alpha4_inv,this]
            (int i)
            {
                rho_phi[istart+i]=PI_IN_SQ*r_inv*H[i][0];
                drho_phi_dr[istart+i]=-r2_inv*(rho_phi[istart+i]+2.0*PI_IN_SQ*alpha_inv*H[i][1]);
                drho_phi_dalpha[istart+i]=-alpha2_inv*(rho_phi[istart+i]-2.0*PI_IN_SQ*r_inv*H[i][2]);
            });
            
            
            rho[i]+=c[j]*rho_phi[istart+2];
            
            if(j<n)
            {
                rho[j]+=c_i*rho_phi[istart+1];
                __vec_lcl[0]+=c_i*c[j]*rho_phi[istart];
            }
            else
                __vec_lcl[0]+=0.5*c_i*c[j]*rho_phi[istart];
        }
        
        
        
        p=rho[i]*drho_inv;
        m=static_cast<size_t>(p);
        m=MIN(m,nrho-2);
        p-=m;
        p=MIN(p,1.0);
        coef=F_arr[elem_i][m];
        
        tmp0=(((coef[4]*p+coef[3])*p+coef[2])*p+coef[1])*p+coef[0];
        tmp1=(((4.0*coef[4]*p+3.0*coef[3])*p+2.0*coef[2])*p+coef[1])*drho_inv;
        
        if(rho[i]>rho_max) tmp0+=tmp1*(rho[i]-rho_max);
        
        rho[i]=tmp0;
        dE[i]=tmp1;
        __f_c[i]=-tmp0;
        
        __vec_lcl[0]+=c_i*tmp0;
    }
    

    
    update(dE_ptr);

    type0* __f_x=f->begin();
    type0* __f_alpha=f_alpha->begin();
    type0 f_i[__dim__]={DESIG(__dim__,0.0)};
    type0 x_i[__dim__];
    istart=0;
    for(int i=0;i<n;i++)
    {
        if(i%c_dim==0)
        {
            Algebra::zero<__dim__>(f_i);
            Algebra::V_eq<__dim__>(x+(i/c_dim)*__dim__,x_i);
        }
        
        type0 c_i=c[i];
        type0 alpha_i=alpha[i];
        type0 dE_i=dE[i];
        type0 f_alpha_i=0.0;
        type0 __f_c_i=0.0;
        const int neigh_sz=neighbor_list_size[i];
        for(int j,__j=0;__j<neigh_sz;__j++,istart+=3)
        {
            j=neighbor_list[i][__j];
            Algebra::DX<__dim__>(x_i,x+(j/c_dim)*__dim__,dx_ij);
            
            fpair=-(drho_phi_dr[istart+2]*dE_i+drho_phi_dr[istart+1]*dE[j]+drho_phi_dr[istart])*c_i*c[j];
            apair=-(drho_phi_dalpha[istart+2]*dE_i+drho_phi_dalpha[istart+1]*dE[j]+drho_phi_dalpha[istart])*c_i*c[j];
            __f_c_i-=c[j]*(rho_phi[istart]+rho_phi[istart+1]*dE[j]);
            if(j<n) __f_c[j]-=c_i*(rho_phi[istart]+rho_phi[istart+2]*dE_i);
            
            Algebra::V_add_x_mul_V<__dim__>(fpair,dx_ij,f_i);
            f_alpha_i+=alpha_i*apair;
            if(j<n)
            {
                Algebra::V_add_x_mul_V<__dim__>(-fpair,dx_ij,__f_x+(j/c_dim)*__dim__);
                __f_alpha[j]+=alpha[j]*apair;
            }
            else
                fpair*=0.5;
            Algebra::DyadicV<__dim__>(-fpair,dx_ij,&__vec_lcl[1]);
        }
        
        __f_alpha[i]+=f_alpha_i;
        __f_c[i]+=__f_c_i;
        
        if((i+1)%c_dim==0)
            Algebra::V_add<__dim__>(f_i,__f_x+__dim__*(i/c_dim));
    }
    
    __f_alpha=f_alpha->begin();
    type0* __alpha=atoms->alpha->begin();
    __f_c=f_c->begin();
    elem_type* __elem=atoms->elem->begin();
    c=atoms->c->begin();
    const int __natms=atoms->natms_lcl;
    type0 ave_f_c_i,log_cv_i;
    // so far the energy that we have calculated is included in both GP and FE
    type0 dbeta_gp_lcl=0.0,dbeta_fe_lcl=0.0,dfe_lcl=0.0,tst_lcl=0.0;
    for(int i=0;i<__natms;i++,c+=c_dim,__f_c+=c_dim,__f_alpha+=c_dim,__alpha+=c_dim,__elem+=c_dim)
    {
        cv_i=1.0;
        for(int ic=0;ic<c_dim;ic++)
        {
            if(c[ic]<0.0) continue;
            cv_i-=c[ic];
            dbeta_fe_lcl+=calc_ent(c[ic])-3.0*c[ic]*log(__alpha[ic]);
            dfe_lcl+=c[ic]*mu_0[__elem[ic]];
        }
        ave_f_c_i=0.0;
        log_cv_i=log(cv_i);
        
        for(int ic=0;ic<c_dim;ic++)
        {
            if(c[ic]<0.0) continue;
            __f_c[ic]+=kBT*(1.5+1.0/cv_i);
            //__f_c[ic]+=kBT*(1.5);
            ave_f_c_i+=c[ic]*__f_c[ic];
        }
        
        for(int ic=0;ic<c_dim;ic++)
        {
            if(c[ic]<0.0) continue;
            __f_c[ic]=c[ic]*(__f_c[ic]-ave_f_c_i);
            if(__alpha[ic]>0.0)__f_alpha[ic]+=3.0*__f_c[ic]/__alpha[ic];
            __f_c[ic]/=kBT;
            tst_lcl+=(beta*mu_0[__elem[ic]]+3.0+log(c[ic])-log_cv_i-3.0*log(__alpha[ic]))*__f_c[ic];
        }
        
        dbeta_gp_lcl+=log_cv_i+1.5*(cv_i-1.0);
        //dbeta_gp_lcl+=1.5*(cv_i-1.0);
        dbeta_fe_lcl+=calc_ent(cv_i);
    }
    
    //__vec_lcl[2+__nvoigt__]=__vec_lcl[0];
    __vec_lcl[1+__nvoigt__]=kBT*dbeta_fe_lcl+dfe_lcl+__vec_lcl[0];
    __vec_lcl[0]+=kBT*dbeta_gp_lcl;
    __vec_lcl[2+__nvoigt__]=__vec_lcl[0]+tst_lcl;
    
    //printf("ff: %.16lf %.16lf %.16lf\n",alpha[0],alpha[1],__vec_lcl[0]);
    
}
/*--------------------------------------------
 energy calculation for grand potential
 --------------------------------------------*/
void ForceFieldEAMDMD::__energy_calc_gp()
{
    type0* rho=ddE_ptr->begin();
    const int n=atoms->natms_lcl*c_dim;
    for(int i=0;i<n;i++) rho[i]=0.0;
    
    elem_type const* elem_vec=atoms->elem->begin();
    
    type0 alpha_ij,rsq,r,r_inv,p,tmp0,tmp1;
    elem_type elem_i,elem_j;
    size_t m;
    
    const type0* x=atoms->x->begin();
    const type0* alpha=atoms->alpha->begin();
    int** neighbor_list=neighbor->neighbor_list;
    int* neighbor_list_size=neighbor->neighbor_list_size;
    type0 const* c=atoms->c->begin();
    type0 const* coef;
    
    size_t istart=0;
    for(int i=0;i<n;i++)
    {
        elem_i=elem_vec[i];
        const int neigh_sz=neighbor_list_size[i];
        for(int j,__j=0;__j<neigh_sz;__j++,istart+=3)
        {
            j=neighbor_list[i][__j];
            elem_j=elem_vec[j];
            rsq=Algebra::RSQ<__dim__>(x+(i/c_dim)*__dim__,x+(j/c_dim)*__dim__);
            r=sqrt(rsq);
            alpha_ij=sqrt(alpha[i]*alpha[i]+alpha[j]*alpha[j]);
            if(r-alpha_ij*xi[N-1]>=cut[elem_i][elem_j]) continue;
            
            type0 upper=(r+cut[elem_i][elem_j])/alpha_ij;
            type0 lower=(r-cut[elem_i][elem_j])/alpha_ij;
            type0 __r,p,tmp0,tmp1;
            type0* coef;
            
            r_inv=1.0/r;
            
            type0 __arr[3]{DESIG(3,0.0)};
            for(int l=0;l<N;l++)
            {
                if(xi[l]<=lower && xi[l]>=upper) continue;
                
                __r=r-xi[l]*alpha_ij;
                p=fabs(__r)*dr_inv;
                m=static_cast<size_t>(p);
                m=MIN(m,nr-2);
                p-=m;
                p=MIN(p,1.0);
                
                coef=r_phi_arr[elem_i][elem_j][m];
                tmp0=((coef[3]*p+coef[2])*p+coef[1])*p+coef[0];
                tmp1=((3.0*coef[3]*p+2.0*coef[2])*p+coef[1])*dr_inv;
                if(__r<0.0) tmp0*=-1.0;
                __arr[0]+=wi_0[l]*tmp0;
                
                
                coef=r_rho_arr[elem_i][elem_j][m];
                tmp0=((coef[3]*p+coef[2])*p+coef[1])*p+coef[0];
                tmp1=((3.0*coef[3]*p+2.0*coef[2])*p+coef[1])*dr_inv;
                if(__r<0.0) tmp0*=-1.0;
                __arr[1]+=wi_0[l]*tmp0;
                
                coef=r_rho_arr[elem_j][elem_i][m];
                tmp0=((coef[3]*p+coef[2])*p+coef[1])*p+coef[0];
                tmp1=((3.0*coef[3]*p+2.0*coef[2])*p+coef[1])*dr_inv;
                if(__r<0.0) tmp0*=-1.0;
                __arr[2]+=wi_0[l]*tmp0;
            }
            
            tmp0=PI_IN_SQ*r_inv;
            
            Algebra::Do<3>::func([&__arr,&tmp0,this](int i)
                                 {
                                     __arr[i]*=tmp0;
                                 });
            
            rho[i]+=c[j]*__arr[2];
            
            if(j<n)
            {
                rho[j]+=c[i]*__arr[1];
                __vec_lcl[0]+=c[i]*c[j]*__arr[0];
            }
            else
                __vec_lcl[0]+=0.5*c[i]*c[j]*__arr[0];
        }
        
        
        if(c[i]<0.0) continue;
        tmp0=rho[i];
        p=tmp0*drho_inv;
        m=static_cast<size_t>(p);
        m=MIN(m,nrho-2);
        p-=m;
        p=MIN(p,1.0);
        coef=F_arr[elem_i][m];
        tmp1=(((coef[4]*p+coef[3])*p+coef[2])*p+coef[1])*p+coef[0];
        if(rho[i]>rho_max)
            tmp1+=(((4.0*coef[4]*p+3.0*coef[3])*p+2.0*coef[2])*p+coef[1])*drho_inv*(tmp0-rho_max);
        
        if(c[i]!=0.0)
            __vec_lcl[0]+=c[i]*tmp1;
    }
    
    const int __natms=atoms->natms_lcl;
    type0 cv_i;
    type0 dbeta_gp_lcl=0.0;
    for(int i=0;i<__natms;i++)
    {
        cv_i=1.0;
        for(int ic=0;ic<c_dim;ic++)
            if(c[i*c_dim+ic]>0.0)
                cv_i-=c[i*c_dim+ic];
        //printf("%.16lf\n",cv_i);
        dbeta_gp_lcl+=log(cv_i)+1.5*(cv_i-1.0);
        //dbeta_gp_lcl+=1.5*(cv_i-1.0);
    }
    __vec_lcl[0]+=kBT*dbeta_gp_lcl;
    //printf("%.16lf %.16lf %.16lf\n",alpha[0],alpha[1],__vec_lcl[0]);
}
