#include "ff_eam_dmd.h"
#include <stdlib.h>
#include "neighbor_dmd.h"
#include "elements.h"
#include "memory.h"
#include "xmath.h"
#include "atoms_dmd.h"
#include "MAPP.h"
#include "dynamic_dmd.h"
#include "import_eam.h"
#include <limits>
#define PI_IN_SQ 0.564189583547756286948079451561
using namespace MAPP_NS;
/*--------------------------------------------
 constructor
 --------------------------------------------*/
ForceFieldEAMDMD::
ForceFieldEAMDMD(AtomsDMD* atoms,
type0 __dr,type0 __drho,size_t __nr,size_t __nrho,
type0(***&& __r_phi_arr)[4],type0(***&& __rho_arr)[4],type0(**&& __F_arr)[5],
type0**&& __cut,type0*&& __r_crd,type0*&& __zeta):
ForceFieldDMD(atoms),
nr(__nr),
nrho(__nrho),
dr(__dr),
drho(__drho),
F_arr(__F_arr),
r_phi_arr(__r_phi_arr),
r_rho_arr(__rho_arr),
dE_ptr(NULL),
ddE_ptr(NULL),
cv_ptr(NULL),
vec0(NULL),
vec1(NULL),
vec2(NULL),
vec3(NULL),
rho_phi(NULL),
drho_phi_dr(NULL),
drho_phi_dalpha(NULL),
ddrho_phi_drdr(NULL),
ddrho_phi_drdalpha(NULL),
ddrho_phi_dalphadalpha(NULL),
max_pairs0(0),
max_pairs1(0),
N(atoms->N),
M_IJ(NULL)
{
    
    __r_phi_arr=NULL;
    __rho_arr=NULL;
    __F_arr=NULL;
    
    dr_inv=1.0/dr;
    drho_inv=1.0/drho;
    rho_max=static_cast<type0>(nrho)*drho;

    

    
    Memory::alloc(xi,N);
    Memory::alloc(wi_0,N);
    Memory::alloc(wi_1,N);
    Memory::alloc(zeta,nelems);
    Memory::alloc(c_1,nelems);
    
    memcpy(xi,atoms->xi,N*sizeof(type0));
    memcpy(wi_0,atoms->wi,N*sizeof(type0));
    for(int i=0;i<N;i++)
        wi_1[i]=wi_0[i]*xi[i];


    for(elem_type elem_i=0;elem_i<nelems;elem_i++)
    {
        r_crd[elem_i]=__r_crd[elem_i];
        zeta[elem_i]=__zeta[elem_i];
    }
    Memory::dealloc(__r_crd);
    Memory::dealloc(__zeta);
    
    for(elem_type elem_i=0;elem_i<nelems;elem_i++)
    {
        for(size_t i=0;i<nr;i++)
            r_rho_arr[elem_i][0][i][0]*=static_cast<type0>(i)*dr;
        ImportEAM::interpolate(nr,dr,r_rho_arr[elem_i][0]);
        for(elem_type elem_j=1;elem_j<nelems;elem_j++)
        {
            if(r_rho_arr[elem_i][elem_j]!=r_rho_arr[elem_i][elem_j-1])
            {
                for(size_t i=0;i<nr;i++)
                    r_rho_arr[elem_i][elem_j][i][0]*=static_cast<type0>(i)*dr;
                ImportEAM::interpolate(nr,dr,r_rho_arr[elem_i][elem_j]);
            }
        }
    }
    

    type0 tmp;
    for(size_t i=0;i<nelems;i++)
        for(size_t j=0;j<i+1;j++)
        {
            cut[i][j]=cut[j][i]=__cut[i][j];
            tmp=__cut[i][j];
            cut_sq[i][j]=cut_sq[j][i]=tmp*tmp;
        }
    Memory::dealloc(__cut);
}
/*--------------------------------------------
 destructor
 --------------------------------------------*/
ForceFieldEAMDMD::~ForceFieldEAMDMD()
{
    Memory::dealloc(F_arr);
    Memory::dealloc(r_rho_arr);
    Memory::dealloc(r_phi_arr);

    Memory::dealloc(zeta);
    Memory::dealloc(c_1);
    Memory::dealloc(wi_1);
    Memory::dealloc(wi_0);
    Memory::dealloc(xi);

}
/*--------------------------------------------
 force calculation
 --------------------------------------------*/
void ForceFieldEAMDMD::__force_calc()
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
    type0* mu_vec=f_c->begin();
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
            
            
            Algebra::Do<3>::func([&istart,&H,&r_inv,&r2_inv,&alpha2_inv,&alpha_inv,this]
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
        mu_vec[i]=tmp0;

        __vec_lcl[0]+=c_i*(tmp0+mu_0[elem_i]-3.0*kBT*log(alpha[i]));
        __vec_lcl[1+__nvoigt__]+=calc_ent(c_i);
    }

    const int __natms=atoms->natms_lcl;
    
    for(int i=0;i<__natms;i++)
    {
        cv_i=1.0;
        for(int ic=0;ic<c_dim;ic++)
            if(c[i*c_dim+ic]>0.0)
                cv_i-=c[i*c_dim+ic];
        __vec_lcl[1+__nvoigt__]+=calc_ent(cv_i);
    }
    
    __vec_lcl[0]+=kBT*__vec_lcl[1+__nvoigt__];
    
    update(dE_ptr);

    type0* fvec=f->begin();
    type0* f_alphavec=f_alpha->begin();
    type0 f_i[__dim__]={DESIG(__dim__,0.0)};
    type0 x_i[__dim__];
    Algebra::zero<__dim__>(x_i);
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
        type0 mu_i=0.0;
        const int neigh_sz=neighbor_list_size[i];
        for(int j,__j=0;__j<neigh_sz;__j++,istart+=3)
        {
            j=neighbor_list[i][__j];
            Algebra::DX<__dim__>(x_i,x+(j/c_dim)*__dim__,dx_ij);
            
            fpair=-(drho_phi_dr[istart+2]*dE_i+drho_phi_dr[istart+1]*dE[j]+drho_phi_dr[istart])*c_i*c[j];
            apair=-(drho_phi_dalpha[istart+2]*dE_i+drho_phi_dalpha[istart+1]*dE[j]+drho_phi_dalpha[istart])*c_i*c[j];
            mu_i+=c[j]*(rho_phi[istart]+rho_phi[istart+1]*dE[j]);
            if(j<n) mu_vec[j]+=c_i*(rho_phi[istart]+rho_phi[istart+2]*dE_i);
            
            Algebra::V_add_x_mul_V<__dim__>(fpair,dx_ij,f_i);
            f_alpha_i+=alpha_i*apair;
            if(j<n)
            {
                Algebra::V_add_x_mul_V<__dim__>(-fpair,dx_ij,fvec+(j/c_dim)*__dim__);
                f_alphavec[j]+=alpha[j]*apair;
            }
            else
                fpair*=0.5;
            Algebra::DyadicV<__dim__>(-fpair,dx_ij,&__vec_lcl[1]);
        }
        
        f_alpha_i+=3.0*kBT*c_i/alpha_i;
        f_alphavec[i]+=f_alpha_i;
        mu_vec[i]+=mu_i;
        
        if((i+1)%c_dim==0)
            Algebra::V_add<__dim__>(f_i,fvec+__dim__*(i/c_dim));
    }
    
    type0 norm_sq_lcl=0.0;
    for(int i=0;i<n;i++)
        norm_sq_lcl+=f_alphavec[i]*f_alphavec[i];
    const int natms_lcl=atoms->natms_lcl;
    for(int i=0;i<natms_lcl;i++,fvec+=__dim__)
        norm_sq_lcl+=Algebra::RSQ<__dim__>(fvec,fvec);
}
/*--------------------------------------------
 force calculation
 --------------------------------------------*/
void ForceFieldEAMDMD::__prepJ_n_res(Vec<type0>* fvec_ptr,Vec<type0>* f_alphavec_ptr)
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
    
    if(max_pairs1<neighbor->no_pairs)
    {
        delete [] ddrho_phi_drdr;
        delete [] ddrho_phi_dalphadalpha;
        delete [] ddrho_phi_drdalpha;
        max_pairs1=neighbor->no_pairs;
        size_t no_0=max_pairs1*3;
        Memory::alloc(ddrho_phi_drdr,no_0);
        Memory::alloc(ddrho_phi_drdalpha,no_0);
        Memory::alloc(ddrho_phi_dalphadalpha,no_0);
    }
    
    for(size_t i=0;i<max_pairs0*3;i++) rho_phi[i]=drho_phi_dr[i]=drho_phi_dalpha[i]=0.0;
    for(size_t i=0;i<max_pairs1*3;i++) ddrho_phi_drdr[i]=ddrho_phi_dalphadalpha[i]=ddrho_phi_drdalpha[i]=0.0;
    
    type0 r,r_inv;
    size_t m;
    type0* coef;
    
    type0 p;
    
    const int n=atoms->natms_lcl*c_dim;
    const type0* c=atoms->c->begin();
    const type0* x=atoms->x->begin();
    const type0* alpha=atoms->alpha->begin();
    type0* E=E_ptr->begin();
    type0* dE=dE_ptr->begin();
    type0* ddE=ddE_ptr->begin();
    type0* rho=E;
    for(int i=0;i<n;i++) rho[i]=0.0;
    
    elem_type const* elem_vec=atoms->elem->begin();
    
    type0 alpha_ij,alpha_ij_sq,rsq,__rho;
    elem_type elem_i,elem_j;
    
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
            alpha_ij_sq=alpha[i]*alpha[i]+alpha[j]*alpha[j];
            alpha_ij=sqrt(alpha_ij_sq);
            if(r-alpha_ij*xi[N-1]>=cut[elem_i][elem_j]) continue;
            type0 upper=(r+cut[elem_i][elem_j])/alpha_ij;
            type0 lower=(r-cut[elem_i][elem_j])/alpha_ij;
            type0 __r,p,tmp0,xi2,xi3,xi4;
            type0* coef;
            
            
            
            type0 H[3][5]{DESIG2(3,5,0.0)};
            for(int l=0;l<N;l++)
            {
                if(xi[l]<=lower && xi[l]>=upper) continue;
                
                xi2=xi[l]*xi[l];
                xi3=xi[l]*xi2;
                xi4=xi2*xi2;
                
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
                H[0][3]+=tmp0*xi3;
                H[0][4]+=tmp0*xi4;
                
                
                coef=r_rho_arr[elem_i][elem_j][m];
                tmp0=((coef[3]*p+coef[2])*p+coef[1])*p+coef[0];
                if(__r<0.0) tmp0*=-1.0;
                tmp0*=wi_0[l];
                H[1][0]+=tmp0;
                H[1][1]+=tmp0*xi[l];
                H[1][2]+=tmp0*xi2;
                H[1][3]+=tmp0*xi3;
                H[1][4]+=tmp0*xi4;
                
                coef=r_rho_arr[elem_j][elem_i][m];
                tmp0=((coef[3]*p+coef[2])*p+coef[1])*p+coef[0];
                if(__r<0.0) tmp0*=-1.0;
                tmp0*=wi_0[l];
                H[2][0]+=tmp0;
                H[2][1]+=tmp0*xi[l];
                H[2][2]+=tmp0*xi2;
                H[2][3]+=tmp0*xi3;
                H[2][4]+=tmp0*xi4;
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
                ddrho_phi_drdr[istart+i]=-r2_inv*(3.0*drho_phi_dr[istart+i]-2.0*drho_phi_dalpha[istart+i]);
                ddrho_phi_drdalpha[istart+i]=-r2alpha2_inv*(3.0*rho_phi[istart+i]+3.0*rsq*drho_phi_dr[istart+i]+alpha_ij_sq*drho_phi_dalpha[istart+i]+4.0*PI_IN_SQ*alpha_inv*H[i][3]);
                ddrho_phi_dalphadalpha[istart+i]=-alpha4_inv*(3.0*rho_phi[istart+i]+6.0*alpha_ij_sq*drho_phi_dalpha[istart+i]-4.0*PI_IN_SQ*r_inv*H[i][4]);
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
        
        __rho=rho[i];
        p=__rho*drho_inv;
        m=static_cast<size_t>(p);
        m=MIN(m,nrho-2);
        p-=m;
        p=MIN(p,1.0);
        coef=F_arr[elem_i][m];

        
        E[i]=(((coef[4]*p+coef[3])*p+coef[2])*p+coef[1])*p+coef[0];
        dE[i]=(((4.0*coef[4]*p+3.0*coef[3])*p+2.0*coef[2])*p+coef[1])*drho_inv;
        
        
        if(__rho>rho_max)
        {
            E[i]+=dE[i]*(__rho-rho_max);
            ddE[i]=0.0;
        }
        else
        {
            ddE[i]=(((12.0*coef[4]*p+6.0*coef[3])*p+2.0*coef[2]))*drho_inv*drho_inv;
        }
        
        __vec_lcl[0]+=c_i*(E[i]+mu_0[elem_i]-3.0*kBT*log(alpha[i]));
        __vec_lcl[1+__nvoigt__]+=calc_ent(c_i);
    }
    
    const int __natms=atoms->natms_lcl;
    type0 cv_i;
    for(int i=0;i<__natms;i++)
    {
        cv_i=1.0;
        for(int ic=0;ic<c_dim;ic++)
            if(c[i*c_dim+ic]>0.0)
                cv_i-=c[i*c_dim+ic];
        __vec_lcl[1+__nvoigt__]+=calc_ent(cv_i);
    }
    __vec_lcl[0]+=kBT*__vec_lcl[1+__nvoigt__];
    
    update(dE_ptr);

    
    type0* f_coef=f_c->begin();
    if(c_dim!=1)
    {
        type0 c_tot;
        for(int i=0;i<__natms;i++)
        {
            c_tot=0.0;
            for(int ic=0;ic<c_dim;ic++)
                if(c[i*c_dim+ic]>0.0)
                    c_tot+=c[i*c_dim+ic];
            if(c_tot==0.0)
            {
                for(int ic=0;ic<c_dim;ic++)
                {
                    if(c[i*c_dim+ic]==0.0)
                        f_coef[i*c_dim+ic]=1.0;
                    else
                        f_coef[i*c_dim+ic]=0.0;
                }
            }
            else
            {
                for(int ic=0;ic<c_dim;ic++)
                {
                    if(c[i*c_dim+ic]>=0.0)
                        f_coef[i*c_dim+ic]=c[i*c_dim+ic]/c_tot;
                    else
                        f_coef[i*c_dim+ic]=0.0;
                }
            }
        }
        
    }
    else
    {
        for(int i=0;i<n;i++)
            if(c[i]>=0.0)
                f_coef[i]=1.0;
    }
    
    
    const int natms_lcl=atoms->natms_lcl;
    dE=dE_ptr->begin();
    type0* fvec=fvec_ptr->begin();
    type0* f_alphavec=f_alphavec_ptr->begin();
    for(int i=0;i<n;i++) f_alphavec[i]=0.0;
    for(int i=0;i<natms_lcl*__dim__;i++) fvec[i]=0.0;
    
    
    type0 f_i[__dim__]={DESIG(__dim__,0.0)};
    type0 x_i[__dim__];
    Algebra::zero<__dim__>(x_i);
    type0 dx_ij[__dim__];
    type0 fpair,apair;
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
        const int neigh_sz=neighbor_list_size[i];
        for(int j,__j=0;__j<neigh_sz;__j++,istart+=3)
        {
            j=neighbor_list[i][__j];
            Algebra::DX<__dim__>(x_i,x+(j/c_dim)*__dim__,dx_ij);
            
            fpair=(drho_phi_dr[istart+2]*dE_i+drho_phi_dr[istart+1]*dE[j]+drho_phi_dr[istart]);
            apair=(drho_phi_dalpha[istart+2]*dE_i+drho_phi_dalpha[istart+1]*dE[j]+drho_phi_dalpha[istart]);
            
            Algebra::V_add_x_mul_V<__dim__>(fpair*c[j]*f_coef[i],dx_ij,f_i);
            f_alpha_i+=alpha_i*apair*c[j];
            
            if(j<n)
            {
                Algebra::V_add_x_mul_V<__dim__>(-fpair*c_i*f_coef[j],dx_ij,fvec+(j/c_dim)*__dim__);
                f_alphavec[j]+=alpha[j]*apair*c_i;
            }
            else
                fpair*=0.5;
            
            Algebra::DyadicV<__dim__>(fpair*c_i*c[j],dx_ij,&__vec_lcl[1]);
            
        }
        
        f_alpha_i-=3.0*kBT/alpha_i;
        f_alphavec[i]+=f_alpha_i;
        
        if((i+1)%c_dim==0)
            Algebra::V_add<__dim__>(f_i,fvec+__dim__*(i/c_dim));
    }
}
/*--------------------------------------------
 
 --------------------------------------------*/
void ForceFieldEAMDMD::__J(Vec<type0>* dx_ptr,Vec<type0>* dalpha_ptr,Vec<type0>* Adx_ptr,Vec<type0>* Adalpha_ptr)
{
    type0* deltadE=E_ptr->begin();
    const int n=atoms->natms_lcl*c_dim;
    const int nn=atoms->natms_lcl*__dim__;
    
    
    
    type0 const* c=atoms->c->begin();
    const type0* x=atoms->x->begin();
    const type0* alpha=atoms->alpha->begin();
    type0* ddE=ddE_ptr->begin();
    type0* dE=dE_ptr->begin();
    int** neighbor_list=neighbor->neighbor_list;
    int* neighbor_list_size=neighbor->neighbor_list_size;
    
    
    type0* dx=dx_ptr->begin();
    type0* dalpha=dalpha_ptr->begin();
    
    
    type0* Adx=Adx_ptr->begin();
    type0* Adalpha=Adalpha_ptr->begin();
    
    type0* f_coef=f_c->begin();
    
    type0 dx_ij[__dim__];
    type0 ddx_ij[__dim__];
    type0 x_i[__dim__];
    Algebra::zero<__dim__>(x_i);
    type0 dx_i[__dim__];
    Algebra::zero<__dim__>(dx_i);
    type0 Adx_i[__dim__];
    Algebra::zero<__dim__>(Adx_i);
    
    type0 pair_alpha_0,pair_alpha_1,df_pair,f_pair,a,b,alpha_i,dalpha_i,Adalpha_i;
    for(int i=0;i<n;i++)
        Adalpha[i]=deltadE[i]=0.0;
    for(int i=0;i<nn;i++)
        Adx[i]=0.0;
    
    
    size_t istart=0;
    for(int i=0;i<n;i++)
    {
        if(i%c_dim==0)
        {
            Algebra::V_eq<__dim__>(x+(i/c_dim)*__dim__,x_i);
            Algebra::V_eq<__dim__>(dx+(i/c_dim)*__dim__,dx_i);
            Algebra::zero<__dim__>(Adx_i);
        }
        
        alpha_i=alpha[i];
        dalpha_i=dalpha[i];
        Adalpha_i=0.0;
        
        type0 c_i=c[i];
        if(c_i<0.0) continue;
        
        const int neigh_sz=neighbor_list_size[i];
        for(int j,__j=0;__j<neigh_sz;__j++,istart+=3)
        {
            j=neighbor_list[i][__j];
            
            
            Algebra::DX<__dim__>(x_i,x+(j/c_dim)*__dim__,dx_ij);
            Algebra::DX<__dim__>(dx_i,dx+(j/c_dim)*__dim__,ddx_ij);
        
            a=Algebra::V_mul_V<__dim__>(dx_ij,ddx_ij);
            b=alpha_i*dalpha_i+alpha[j]*dalpha[j];

            pair_alpha_0=-(
            (ddrho_phi_drdalpha[istart]+ddrho_phi_drdalpha[istart+2]*dE[i]+ddrho_phi_drdalpha[istart+1]*dE[j])*a
            +(ddrho_phi_dalphadalpha[istart]+ddrho_phi_dalphadalpha[istart+2]*dE[i]+ddrho_phi_dalphadalpha[istart+1]*dE[j])*b);
            
            pair_alpha_1=-(drho_phi_dalpha[istart]+drho_phi_dalpha[istart+2]*dE[i]+drho_phi_dalpha[istart+1]*dE[j]);
            
            
            
            df_pair=-(
            (ddrho_phi_drdalpha[istart]+ddrho_phi_drdalpha[istart+2]*dE[i]+ddrho_phi_drdalpha[istart+1]*dE[j])*b
            +(ddrho_phi_drdr[istart]+ddrho_phi_drdr[istart+2]*dE[i]+ddrho_phi_drdr[istart+1]*dE[j])*a);
            
            f_pair=-(drho_phi_dr[istart]+drho_phi_dr[istart+2]*dE[i]+drho_phi_dr[istart+1]*dE[j]);
            
            
            
            
            deltadE[i]+=(a*drho_phi_dr[istart+2]+b*drho_phi_dalpha[istart+2])*ddE[i]*c[j];
            Algebra::V_add_x_mul_V<__dim__>(df_pair*c[j]*f_coef[i],dx_ij,Adx_i);
            Algebra::V_add_x_mul_V<__dim__>(f_pair*c[j]*f_coef[i],ddx_ij,Adx_i);
            Adalpha_i+=c[j]*(pair_alpha_0*alpha_i+pair_alpha_1*dalpha_i);
            
            if(j<n)
            {
                deltadE[j]+=(a*drho_phi_dr[istart+1]+b*drho_phi_dalpha[istart+1])*ddE[j]*c[i];
                Algebra::V_add_x_mul_V<__dim__>(-df_pair*c[i]*f_coef[j],dx_ij,Adx+(j/c_dim)*__dim__);
                Algebra::V_add_x_mul_V<__dim__>(-f_pair*c[i]*f_coef[j],ddx_ij,Adx+(j/c_dim)*__dim__);
                Adalpha[j]+=c[i]*(pair_alpha_0*alpha[j]+pair_alpha_1*dalpha[j]);
            }
            else
            {
                f_pair*=0.5;
                df_pair*=0.5;
            }
            
            Algebra::DyadicV<__dim__>(-df_pair*c_i*c[j],dx_ij,__vec_lcl);
            Algebra::DyadicV<__dim__>(-f_pair*c_i*c[j],ddx_ij,dx_ij,__vec_lcl);
            
        }
        
        
        Adalpha[i]+=Adalpha_i;
        if((i+1)%c_dim==0)
            Algebra::V_add<__dim__>(Adx_i,Adx+__dim__*(i/c_dim));
    }
    
    update(E_ptr);

    deltadE=E_ptr->begin();
    istart=0;
    
    for(int i=0;i<n;i++)
    {
        if(i%c_dim==0)
        {
            Algebra::zero<__dim__>(Adx_i);
            Algebra::V_eq<__dim__>(x+(i/c_dim)*__dim__,x_i);
        }
        
        alpha_i=alpha[i];
        dalpha_i=dalpha[i];
        Adalpha_i=0.0;
        
        type0 c_i=c[i];
        if(c_i<0.0) continue;
        
        const int neigh_sz=neighbor_list_size[i];
        for(int j,__j=0;__j<neigh_sz;__j++,istart+=3)
        {
            j=neighbor_list[i][__j];
            Algebra::DX<__dim__>(x_i,x+(j/c_dim)*__dim__,dx_ij);
            
            df_pair=-(drho_phi_dr[istart+2]*deltadE[i]+drho_phi_dr[istart+1]*deltadE[j]);
            pair_alpha_0=-(drho_phi_dalpha[istart+2]*deltadE[i]+drho_phi_dalpha[istart+1]*deltadE[j]);
            
            Algebra::V_add_x_mul_V<__dim__>(df_pair*c[j]*f_coef[i],dx_ij,Adx_i);
            Adalpha_i+=alpha_i*pair_alpha_0*c[j];
            
            if(j<n)
            {
                Algebra::V_add_x_mul_V<__dim__>(-df_pair*c_i*f_coef[j],dx_ij,Adx+(j/c_dim)*__dim__);
                Adalpha[j]+=alpha[j]*pair_alpha_0*c_i;
            }
            else
                df_pair*=0.5;
            
            Algebra::DyadicV<__dim__>(-df_pair*c_i*c[j],dx_ij,__vec_lcl);

        }
        
        
        Adalpha[i]+=Adalpha_i-3.0*kBT*dalpha_i/(alpha_i*alpha_i);
        if((i+1)%c_dim==0)
            Algebra::V_add<__dim__>(Adx_i,Adx+__dim__*(i/c_dim));
    }
}
/*--------------------------------------------
 energy calculation
 --------------------------------------------*/
void ForceFieldEAMDMD::__energy_calc()
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
            type0 __r,p,tmp0;
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
            
            Algebra::Do<3>::func([&__arr,&tmp0](int i)
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
            __vec_lcl[0]+=c[i]*(tmp1+mu_0[elem_i]-3.0*kBT*log(alpha[i]));
        __vec_lcl[1+__nvoigt__]+=calc_ent(c[i]);
    }
    const int __natms=atoms->natms_lcl;
    type0 cv_i;
    for(int i=0;i<__natms;i++)
    {
        cv_i=1.0;
        for(int ic=0;ic<c_dim;ic++)
            if(c[i*c_dim+ic]>0.0)
                cv_i-=c[i*c_dim+ic];
        __vec_lcl[1+__nvoigt__]+=calc_ent(cv_i);
    }
    __vec_lcl[0]+=kBT*__vec_lcl[1+__nvoigt__];
}
/*--------------------------------------------
 init
 --------------------------------------------*/
void ForceFieldEAMDMD::init()
{
    pre_init();
    set_temp();
    
    
    cv_ptr=new Vec<type0>(atoms,1);
    E_ptr=new Vec<type0>(atoms,c_dim);
    dE_ptr=new Vec<type0>(atoms,c_dim);
    ddE_ptr=new Vec<type0>(atoms,c_dim);
}
/*--------------------------------------------
 fin
 --------------------------------------------*/
void ForceFieldEAMDMD::fin()
{
    
    Memory::dealloc(rho_phi);
    Memory::dealloc(drho_phi_dr);
    Memory::dealloc(drho_phi_dalpha);
    
    
    Memory::dealloc(ddrho_phi_drdr);
    Memory::dealloc(ddrho_phi_dalphadalpha);
    Memory::dealloc(ddrho_phi_drdalpha);
    max_pairs0=max_pairs1=0;
    
    delete ddE_ptr;
    delete dE_ptr;
    delete E_ptr;
    delete cv_ptr;
    
    dE_ptr=ddE_ptr=cv_ptr=vec0=vec1=vec2=vec3=NULL;
    post_fin();
}
/*--------------------------------------------
 create the sparse matrices
 --------------------------------------------*/
void ForceFieldEAMDMD::init_static()
{
    Memory::alloc(M_IJ,4*neighbor->no_pairs_2nd);
    vec0=new Vec<type0>(atoms,c_dim);
    vec1=new Vec<type0>(atoms,1);
    vec2=new Vec<type0>(atoms,c_dim);
    vec3=new Vec<type0>(atoms,c_dim);
    fcoef_ptr=new Vec<type0>(atoms,1);
    
    const int __natms_lcl=atoms->natms_lcl;
    type0* f_coef=fcoef_ptr->begin();
    type0* c=atoms->c->begin();
    if(c_dim!=1)
    {
        type0 c_tot;
        for(int i=0;i<__natms_lcl;i++)
        {
            c_tot=0.0;
            for(int ic=0;ic<c_dim;ic++)
                if(c[i*c_dim+ic]>0.0)
                    c_tot+=c[i*c_dim+ic];
            if(c_tot==0.0)
            {
                for(int ic=0;ic<c_dim;ic++)
                {
                    if(c[i*c_dim+ic]==0.0)
                        f_coef[i*c_dim+ic]=1.0;
                    else
                        f_coef[i*c_dim+ic]=0.0;
                }
            }
            else
            {
                for(int ic=0;ic<c_dim;ic++)
                {
                    if(c[i*c_dim+ic]>=0.0)
                        f_coef[i*c_dim+ic]=c[i*c_dim+ic]/c_tot;
                    else
                        f_coef[i*c_dim+ic]=0.0;
                }
            }
        }
        
    }
    else
    {
        const int n=__natms_lcl*c_dim;
        for(int i=0;i<n;i++)
            if(c[i]>=0.0)
                f_coef[i]=1.0;
    }
}
/*--------------------------------------------
 create the sparse matrices
 --------------------------------------------*/
void ForceFieldEAMDMD::fin_static()
{
    delete fcoef_ptr;
    delete vec3;
    delete vec2;
    delete vec1;
    delete vec0;
    Memory::dealloc(M_IJ);
}
/*--------------------------------------------
 set the temperature in the simulation
 --------------------------------------------*/
void ForceFieldEAMDMD::set_temp()
{
    for(size_t i=0;i<nelems;i++)
        c_1[i]=sqrt(0.5*kBT/atoms->elements.masses[i])/M_PI;
}
/*--------------------------------------------
 calculate cv
 calculate mu
 
 external vecs
 atoms->alpha;
 atoms->c;
 atoms->elem;

 internal vecs
 cv_ptr;
 E_ptr;
 dE_ptr;
 mu_ptr;
 --------------------------------------------*/
void ForceFieldEAMDMD::calc_mu()
{
    /*--------------------------------
     2 things are hapenning here:
     0. calculate c_v for local atoms
     1. calculate the non interacting 
     part of energy
     --------------------------------*/
    //external vecs
    const type0* alpha=atoms->alpha->begin();
    type0 const* c=atoms->c->begin();
    elem_type const* elem_vec=atoms->elem->begin();
    
    //internal vecs
    type0* cv=cv_ptr->begin();
    
    int natms_lcl=atoms->natms_lcl;
    type0 en=0.0;
    for(int i=0;i<natms_lcl;i++)
    {
        type0 cv_i=1.0;
        for(int j=0;j<c_dim;j++)
            if(c[i*c_dim+j]>0.0)
            {
                __vec_lcl[1+__nvoigt__]+=calc_ent(c[i*c_dim+j]);
                en+=c[i*c_dim+j]*
                (mu_0[elem_vec[i*c_dim+j]]-3.0*kBT*log(alpha[i*c_dim+j]));
                cv_i-=c[i*c_dim+j];
            }
        cv[i]=cv_i;
        __vec_lcl[1+__nvoigt__]+=calc_ent(cv_i);
    }
    
    __vec_lcl[0]+=en;
    __vec_lcl[0]+=kBT*__vec_lcl[1+__nvoigt__];
    
    /*--------------------------------
     we will calculate the rest of c_v
     here
     --------------------------------*/
    int nall=natms_lcl+atoms->natms_ph;
    for(int i=natms_lcl;i<nall;i++)
    {
        type0 cv_i=1.0;
        for(int j=0;j<c_dim;j++)
            if(c[i*c_dim+j]>0.0)
                cv_i-=c[i*c_dim+j];
        cv[i]=cv_i;
    }
    
    
    
    /*--------------------------------
     calculate the electron densities
     --------------------------------*/
    //internal vecs
    type0* rho=E_ptr->begin();
    type0* mu_vec=f_c->begin();
    
    int n=atoms->natms_lcl*c_dim;
    for(int i=0;i<n;i++) mu_vec[i]=rho[i]=0.0;
    
    size_t m;
    type0 p,__E,__dE,__rho;
    type0* coef;
    
    //internal vecs
    type0* dE=dE_ptr->begin();
    type0* ddE=ddE_ptr->begin();
    int** neighbor_list=neighbor->neighbor_list;
    int* neighbor_list_size=neighbor->neighbor_list_size;
    size_t istart=0;
    for(int i=0;i<n;i++)
    {
        if(c[i]<0.0) continue;
        const int neigh_sz=neighbor_list_size[i];
        for(int j,__j=0;__j<neigh_sz;__j++,istart+=3)
        {
            j=neighbor_list[i][__j];
            rho[i]+=c[j]*rho_phi[istart+2];

            if(j<n)
                rho[j]+=c[i]*rho_phi[istart+1];
        }
        
        __rho=rho[i];
        p=__rho*drho_inv;
        m=static_cast<size_t>(p);
        m=MIN(m,nrho-2);
        p-=m;
        p=MIN(p,1.0);
        coef=F_arr[elem_vec[i]][m];
        
        __E=(((coef[4]*p+coef[3])*p+coef[2])*p+coef[1])*p+coef[0];
        __dE=(((4.0*coef[4]*p+3.0*coef[3])*p+2.0*coef[2])*p+coef[1])*drho_inv;
        ddE[i]=(((12.0*coef[4]*p+6.0*coef[3])*p+2.0*coef[2]))*drho_inv*drho_inv;
        if(__rho>rho_max)
            __E+=__dE*(__rho-rho_max);
        
        __vec_lcl[0]+=c[i]*__E;
        rho[i]=__E;
        mu_vec[i]+=__E;
        dE[i]=__dE;
    }
    
    update(dE_ptr);

    
    type0* fcoef=fcoef_ptr->begin();
    type0* falpha=f_alpha->begin();
    type0* fvec=f->begin();
    memset(falpha,0,sizeof(type0)*n);
    memset(fvec,0,sizeof(type0)*__dim__*atoms->natms_lcl);
    /*
    type0 pair_0,pair_1;
    
    istart=0;
    for(int i=0;i<n;i++)
    {
        const int neigh_sz=neighbor_list_size[i];
        for(int j,__j=0;__j<neigh_sz;__j++,istart+=3)
        {
            j=neighbor_list[i][__j];
            mu[i]+=c[j]*(rho_phi[istart]+rho_phi[istart+1]*dE[j]);
            pair_0=(drho_phi_dr[istart]+drho_phi_dr[istart+1]*dE[j]+drho_phi_dr[istart+2]*dE[i]);
            pair_1=(drho_phi_dalpha[istart]+drho_phi_dalpha[istart+1]*dE[j]+drho_phi_dalpha[istart+2]*dE[i]);
            
            
            if(j<n)
            {
                mu[j]+=c[i]*(rho_phi[istart]+rho_phi[istart+2]*dE[i]);
                __vec_lcl[0]+=c[i]*c[j]*rho_phi[istart];
            }
            else
                __vec_lcl[0]+=0.5*c[i]*c[j]*rho_phi[istart];
        }
    }
    dynamic->update(mu_ptr);
    */
    
    
    
    
    
    type0* x=atoms->x->begin();
    type0 f_i[__dim__]={DESIG(__dim__,0.0)};
    type0 x_i[__dim__];
    Algebra::zero<__dim__>(x_i);
    type0 dx_ij[__dim__];
    type0 fpair,apair;
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
        const int neigh_sz=neighbor_list_size[i];
        for(int j,__j=0;__j<neigh_sz;__j++,istart+=3)
        {
            j=neighbor_list[i][__j];
            mu_vec[i]+=c[j]*(rho_phi[istart]+rho_phi[istart+1]*dE[j]);
            Algebra::DX<__dim__>(x_i,x+(j/c_dim)*__dim__,dx_ij);
            
            fpair=(drho_phi_dr[istart+2]*dE_i+drho_phi_dr[istart+1]*dE[j]+drho_phi_dr[istart]);
            apair=(drho_phi_dalpha[istart+2]*dE_i+drho_phi_dalpha[istart+1]*dE[j]+drho_phi_dalpha[istart]);
            
            Algebra::V_add_x_mul_V<__dim__>(fpair*c[j]*fcoef[i],dx_ij,f_i);
            f_alpha_i+=alpha_i*apair*c[j];
            
            if(j<n)
            {
                mu_vec[j]+=c[i]*(rho_phi[istart]+rho_phi[istart+2]*dE[i]);
                Algebra::V_add_x_mul_V<__dim__>(-fpair*c_i*fcoef[j],dx_ij,fvec+(j/c_dim)*__dim__);
                falpha[j]+=alpha[j]*apair*c_i;
                __vec_lcl[0]+=c[i]*c[j]*rho_phi[istart];
            }
            else
                __vec_lcl[0]+=0.5*c[i]*c[j]*rho_phi[istart];
        }
        
        f_alpha_i-=3.0*kBT/alpha_i;
        falpha[i]+=f_alpha_i;
        
        if((i+1)%c_dim==0)
            Algebra::V_add<__dim__>(f_i,fvec+__dim__*(i/c_dim));
        
    }

    update(f_c);

    
}
/*--------------------------------------------
 force calculation
 --------------------------------------------*/
void ForceFieldEAMDMD::__force_calc_static()
{
    int n=atoms->natms_lcl*c_dim;
    type0 const* E=E_ptr->begin();
    type0 const* c=atoms->c->begin();
    type0 const* cv=cv_ptr->begin();
    type0 const* alpha=atoms->alpha->begin();
    elem_type const* elem=atoms->elem->begin();
    for(int i=0;i<n;i++)
    {
        if(i%c_dim==0)
            __vec_lcl[1+__nvoigt__]+=calc_ent(cv[i]);
        if(c[i]>=0.0)
        {
            __vec_lcl[0]+=c[i]*(E[i]+mu_0[elem[i]]-3.0*kBT*log(alpha[i]));
            __vec_lcl[1+__nvoigt__]+=calc_ent(c[i]);
        }
    }
    
    __vec_lcl[0]+=kBT*__vec_lcl[1+__nvoigt__];
    
    type0 const* dE=dE_ptr->begin();
    type0 const* x=atoms->x->begin();
    

    
    type0 x_i[__dim__];
    Algebra::zero<__dim__>(x_i);
    type0 dx_ij[__dim__];
    
    type0 fpair;
    int** neighbor_list=neighbor->neighbor_list;
    int* neighbor_list_size=neighbor->neighbor_list_size;
    size_t istart=0;
    for(int i=0;i<n;i++)
    {
        if(i%c_dim==0) Algebra::V_eq<__dim__>(x+(i/c_dim)*__dim__,x_i);
        type0 c_i=c[i];
        type0 dE_i=dE[i];
        const int neigh_sz=neighbor_list_size[i];
        
        for(int j,__j=0;__j<neigh_sz;__j++,istart+=3)
        {
            j=neighbor_list[i][__j];
            Algebra::DX<__dim__>(x_i,x+(j/c_dim)*__dim__,dx_ij);
            fpair=-(drho_phi_dr[istart+2]*dE_i+drho_phi_dr[istart+1]*dE[j]+drho_phi_dr[istart])*c_i*c[j];
            if(j<n)
                __vec_lcl[0]+=rho_phi[istart]*c_i*c[j];
            else
            {
                fpair*=0.5;
                __vec_lcl[0]+=0.5*rho_phi[istart]*c_i*c[j];
            }
            
            Algebra::DyadicV<__dim__>(-fpair,dx_ij,&__vec_lcl[1]);
        }
    }
}
/*--------------------------------------------
 create the sparse matrices
 
 external vecs
 atoms->x;
 atoms->alpha;
 atoms->c;
 atoms->elem;
 atoms->c_d;
 
 
 internal vecs
 type0* mu=mu_ptr->begin();
 type0* cv=cv_ptr->begin();
 
 F=c-scalar*c_d(c)+a
 --------------------------------------------*/
void ForceFieldEAMDMD::__c_d_calc()
{
    calc_mu();
    
    //external vecs
    const type0* x=atoms->x->begin();
    const type0* alpha=atoms->alpha->begin();
    const type0* c=atoms->c->begin();
    const elem_type* elem_vec=atoms->elem->begin();
    type0* __c_d=c_d->begin();
    
    
    //internal vecs
    type0* cv=cv_ptr->begin();
    type0* mu_vec=f_c->begin();
    
    
    type0 rsq,alpha_i,alpha_j,alpha_sq_i,alpha_sq_j,gamma_i,gamma_j,mu_i,mu_ji;
    type0 Q_i,gamma_ij_inv,theta_i,theta_j,d_i,d_j,coef0;
    type0 dc_ij,c_1_ielem,zeta_ielem;
    size_t istart=0;
    int** neighbor_list=neighbor->neighbor_list_2nd;
    int* neighbor_list_size=neighbor->neighbor_list_size_2nd;
    const int n=atoms->natms_lcl*c_dim;
    
    type0 c_d_norm_lcl=0.0;
    for(int i=0;i<n;i++)
    {
        if(c[i]<0.0) continue;
        alpha_i=alpha[i];
        alpha_sq_i=alpha_i*alpha_i;
        mu_i=mu_vec[i];
        c_1_ielem=c_1[elem_vec[i]];
        zeta_ielem=zeta[elem_vec[i]];
        const int neigh_sz=neighbor_list_size[i];
        for(int j,__j=0;__j<neigh_sz;__j++,istart+=4)
        {
            j=neighbor_list[i][__j];
            if(!dof_c_empty && !(atoms->c_dof->begin()[i] && atoms->c_dof->begin()[j]))
            {
                M_IJ[istart]=0.0;
                M_IJ[istart+1]=0.0;
                M_IJ[istart+2]=0.0;
                M_IJ[istart+3]=0.0;
                continue;
            }
            alpha_j=alpha[j];
            alpha_sq_j=alpha_j*alpha_j;
            rsq=Algebra::RSQ<__dim__>(x+(i/c_dim)*__dim__,x+(j/c_dim)*__dim__);
            gamma_i=rsq/(zeta_ielem*alpha_sq_i);
            gamma_j=rsq/(zeta_ielem*alpha_sq_j);
            mu_ji=beta*(mu_vec[j]-mu_i);
            calc_Q(gamma_i,gamma_j,mu_ji,Q_i,gamma_ij_inv,theta_i);
            
            
            coef0=-c_1_ielem*rsq*gamma_ij_inv*exp(-Q_i+0.5*mu_ji);
            
            /*
            exp_mod_Q_i=exp(-Q_i);
            exp_mod_Q_j=exp(-Q_i+mu_ji);

            d_i=d_j=-c_1_ielem*rsq*gamma_ij_inv;
            d_i/=alpha_sq_i*alpha_i;
            d_j/=alpha_sq_j*alpha_j;
            d_i*=exp(-Q_i);
            d_j*=exp(-Q_i+mu_ji);;
             */
            theta_j=theta_i+1.0;
            
            dc_ij=coef0*(c[i]*cv[j/c_dim]*exp(-0.5*mu_ji)/(alpha_sq_i*alpha_i)-c[j]*cv[i/c_dim]*exp(0.5*mu_ji)/(alpha_sq_j*alpha_j));
            
            //M_IJ[istart]=beta*(c[i]*cv[j/c_dim]*d_i*theta_i-c[j]*cv[i/c_dim]*d_j*theta_j);
            
            //M_IJ[istart+1]=d_i;
            //M_IJ[istart+2]=d_j;
            
            d_i=exp(-0.5*mu_ji)/(alpha_sq_i*alpha_i);
            d_j=exp(0.5*mu_ji)/(alpha_sq_j*alpha_j);
            
            M_IJ[istart]=beta*(theta_i*c[i]*cv[j/c_dim]*d_i-theta_j*c[j]*cv[i/c_dim]*d_j);
            M_IJ[istart+1]=d_i;
            M_IJ[istart+2]=d_j;
            M_IJ[istart+3]=coef0;
            
            
            __c_d[i]+=dc_ij;

            if(j<n)
                __c_d[j]-=dc_ij;
        }
        c_d_norm_lcl=MAX(c_d_norm_lcl,fabs(__c_d[i]));
    }
    
    MPI_Allreduce(&c_d_norm_lcl,&c_d_norm,1,Vec<type0>::MPI_T,MPI_MAX,atoms->world);
}
/*--------------------------------------------
 create the sparse matrices
 
 external vecs
 atoms->c;
 
 internal vecs
 E_ptr;
 dE_ptr;
 vec0;
 vec1;
 cv_ptr;
 
 input vecs
 x_ptr;
 Ax_ptr;
 
 F_i=c_i-beta*(d c_i)/dt
 J_ij=d F_i/(d c_j)
 Ax_i=-J_ij*x_j;
 --------------------------------------------*/
void ForceFieldEAMDMD::__J(Vec<type0>* x_ptr,Vec<type0>* Ax_ptr)
{
    update(x_ptr);

    //external vecs
    type0* c=atoms->c->begin();
    
    //internal vecs
    type0* ddE=ddE_ptr->begin();
    type0* dE=dE_ptr->begin();
    type0* b0=vec0->begin();
    
    //input vecs
    type0* x=x_ptr->begin();
    type0* Ax=Ax_ptr->begin();
    
    int** neighbor_list=neighbor->neighbor_list;
    int* neighbor_list_size=neighbor->neighbor_list_size;
    size_t istart=0;
    const int n=atoms->natms_lcl*c_dim;
    for(int i=0;i<n;i++) b0[i]=Ax[i]=0.0;
    for(int i=0;i<n;i++)
    {
        const int neigh_sz=neighbor_list_size[i];
        type0 ci_dEEi=c[i]*ddE[i];
        for(int j,__j=0;__j<neigh_sz;__j++,istart+=3)
        {
            j=neighbor_list[i][__j];
            Ax[i]+=ci_dEEi*rho_phi[istart+2]*x[j];
            if(j<n)
                Ax[j]+=c[j]*ddE[j]*rho_phi[istart+1]*x[i];
        }
    }

    update(Ax_ptr);

    
    type0 tmp0;
    istart=0;
    for(int i=0;i<n;i++)
    {
        const int neigh_sz=neighbor_list_size[i];
        for(int j,__j=0;__j<neigh_sz;__j++,istart+=3)
        {
            j=neighbor_list[i][__j];
            tmp0=rho_phi[istart]+dE[i]*rho_phi[istart+2]+dE[j]*rho_phi[istart+1];
            b0[i]+=rho_phi[istart+1]*Ax[j]+tmp0*x[j];
            if(j<n)
                b0[j]+=rho_phi[istart+2]*Ax[i]+tmp0*x[i];
        }
    }
    
    update(vec0);

    
    //internal vecs
    type0* xv=vec1->begin();
    type0* cv=cv_ptr->begin();
    
    const int nall=atoms->natms_lcl+atoms->natms_ph;
    for(int i=0;i<nall;i++)
    {
        xv[i]=0.0;
        for(int j=0;j<c_dim;j++)
            if(c[i*c_dim+j]>=0.0)
                xv[i]+=x[i*c_dim+j];
    }

    for(int i=0;i<n;i++)
        Ax[i]=0.0;
    
    
    neighbor_list=neighbor->neighbor_list_2nd;
    neighbor_list_size=neighbor->neighbor_list_size_2nd;
    istart=0;
    for(int i=0;i<n;i++)
    {
        const int neigh_sz=neighbor_list_size[i];
        for(int j,__j=0;__j<neigh_sz;__j++,istart+=4)
        {
            j=neighbor_list[i][__j];
            
            tmp0=M_IJ[istart]*(b0[j]-b0[i]);
            tmp0+=M_IJ[istart+1]*(cv[j/c_dim]*x[i]-c[i]*xv[j/c_dim]);
            tmp0-=M_IJ[istart+2]*(cv[i/c_dim]*x[j]-c[j]*xv[i/c_dim]);
            tmp0*=M_IJ[istart+3];
            Ax[i]+=tmp0;
            
            if(j<n) Ax[j]-=tmp0;
        }
    }
    
}
/*--------------------------------------------
 gamma_i = x_ij*x_ij/(alpha_i*alpha_i)
 gamma_j = x_ij*x_ij/(alpha_j*alpha_j)
 mu_ji = beta*(mu_j-mu_i)
 
 Q=beta*U
 
 alpha_Q_sq=alpha_U_sq/(x_ij*x_ij)
 --------------------------------------------*/
inline void ForceFieldEAMDMD::
calc_Q(type0& gamma_i,type0& gamma_j,type0& mu_ji
,type0& Q,type0& alpha_Q_sq)
{
    type0 z0;
    type0 a=5.0*(gamma_i-gamma_j-6.0*mu_ji);
    type0 xi=-2.0*(gamma_i+gamma_j);
    type0 d=sqrt((xi-a)*(xi-a)-8.0*a*gamma_i);
    z0=-4.0*gamma_i/((xi-a)-d);
    
    /*
    if(z0<0.0 || z0>1.0)
    {
        Error::abort("could not find z0");
    }
    */
    
    Q=z0*z0*(1.0-z0)*(1.0-z0)*((1.0-z0)*gamma_i+z0*gamma_j);
    Q+=z0*z0*z0*(6.0*z0*z0-15.0*z0+10.0)*mu_ji;
    alpha_Q_sq=1.0/((1.0-z0)*gamma_i+z0*gamma_j);
}
/*--------------------------------------------
 gamma_i = x_ij*x_ij/(alpha_i*alpha_i)
 gamma_j = x_ij*x_ij/(alpha_j*alpha_j)
 mu_ji = beta*(mu_j-mu_i)
 
 Q=beta*U
 
 dz0: derivative of z0 w.r.t. mu_ji
 dQ:  derivative of Q  w.r.t. mu_ji
 
 alpha_Q_sq=alpha_U_sq/(x_ij*x_ij)
 dalpha_Q_sq/dmu_ji
 --------------------------------------------*/
inline void ForceFieldEAMDMD::
calc_Q(type0& gamma_i,type0& gamma_j,type0& mu_ji
,type0& Q,type0& alpha_Q_sq,type0& theta)
{
    type0 z0,dz0,dQ;
    type0 a=5.0*(gamma_i-gamma_j-6.0*mu_ji);
    type0 xi=-2.0*(gamma_i+gamma_j);
    type0 delta_sqrt=sqrt((xi-a)*(xi-a)-8.0*a*gamma_i);
    z0=-4.0*gamma_i/((xi-a)-delta_sqrt);
    dz0=30.0*z0*(1.0-z0)/delta_sqrt;
    /*
     here 
    
     
     */
    dQ=z0*z0*z0*(6.0*z0*z0-15.0*z0+10.0);
    Q=z0*z0*(1.0-z0)*(1.0-z0)*((1.0-z0)*gamma_i+z0*gamma_j);
    Q+=dQ*mu_ji;
    
    alpha_Q_sq=1.0/((1.0-z0)*gamma_i+z0*gamma_j);
    theta=alpha_Q_sq*(gamma_i-gamma_j)*dz0-dQ;

}
/*--------------------------------------------
 python constructor
 --------------------------------------------*/
void ForceFieldEAMDMD::ml_new(PyMethodDef& method_0,PyMethodDef& method_1,PyMethodDef& method_2)
{
    method_0.ml_flags=METH_VARARGS | METH_KEYWORDS;
    method_0.ml_name="ff_eam_funcfl";
    method_0.ml_meth=(PyCFunction)(PyCFunctionWithKeywords)(
    [](PyObject* self,PyObject* args,PyObject* kwds)->PyObject*
    {
        AtomsDMD::Object* __self=reinterpret_cast<AtomsDMD::Object*>(self);
        FuncAPI<std::string*,type0*,type0*,std::string*> f("ff_eam_funcfl",{"funcfl_files","r_crd","C","elems"});
        f.noptionals=1;
        f.logics<1>()[0]=VLogics("gt",0.0);
        f.logics<2>()[0]=VLogics("gt",0.0);
        
        const std::string* names=__self->atoms->elements.names;
        size_t nelems=__self->atoms->elements.nelems;
        if(f(args,kwds)) return NULL;
        if(f.remap<3,0,1,2>("elements present in system",names,nelems)) return NULL;
        
        size_t nr,nrho;
        type0 dr,drho;
        type0** r_c;
        type0(** F)[5]=NULL;
        type0(*** r_phi)[4]=NULL;
        type0(*** rho)[4]=NULL;
        try
        {
            ImportEAM::funcfl(nelems,f.val<0>(),dr,drho,nr,nrho,r_phi,rho,F,r_c);
        }
        catch(char* err_msg)
        {
            PyErr_SetString(PyExc_TypeError,err_msg);
            delete [] err_msg;
            return NULL;
        }
        
        for(size_t i=0;i<nelems;i++)
            if(f.v<1>()[i]>r_c[i][i])
            {
                Memory::dealloc(r_c);
                Memory::dealloc(F);
                Memory::dealloc(r_phi);
                Memory::dealloc(rho);
                PyErr_Format(PyExc_TypeError,"r_crd[%zu] should be less than r_c[%zu][%zu]",i,i,i);
                return NULL;
            }
        
        delete __self->ff;
        __self->ff=new ForceFieldEAMDMD(__self->atoms,dr,drho,nr,nrho,std::move(r_phi),std::move(rho),std::move(F),std::move(r_c),f.mov<1>(),f.mov<2>());
        Py_RETURN_NONE;
    });
    method_0.ml_doc=(char*)R"---(
    ff_eam_funcfl(funcfl_files,r_crd,C,elems=None)
   
    Tabulated EAM force field given by FuncFL file/s
    
    Assigns EAM force field to system
    
    Parameters
    ----------
    funcfl_files : string[nelems]
        list of relative paths to DYNAMO files with FuncFL format
    r_crd : double[nelems]
        cutoff radius for mass exchange
    C : double[nelems]
        dimensionless factor for each element, see :ref:`here <master-ref>`
    elems : string[nelems]
        mapping elements
    
    Returns
    -------
    None
   
    Notes
    -----
    This is tabulated form of Embedded Atom Method (EAM) potential
    
    
    Examples
    --------
    Ni
    
    ::
     
        >>> from mapp import dmd
        >>> sim=dmd.cfg(5,"configs/Ni-DMD.cfg")
        >>> sim.ff_eam_funcfl("potentials/niu3.eam")
    
    

    )---";
    
    method_1.ml_flags=METH_VARARGS | METH_KEYWORDS;
    method_1.ml_name="ff_eam_setfl";
    method_1.ml_meth=(PyCFunction)(PyCFunctionWithKeywords)(
    [](PyObject* self,PyObject* args,PyObject* kwds)->PyObject*
    {
        AtomsDMD::Object* __self=reinterpret_cast<AtomsDMD::Object*>(self);
        FuncAPI<std::string,type0*,type0*,std::string*> f("ff_eam_setfl",{"setfl_file","r_crd","C","elems"});
        f.noptionals=1;
        f.logics<1>()[0]=VLogics("gt",0.0);
        f.logics<2>()[0]=VLogics("gt",0.0);
        
        const std::string* names=__self->atoms->elements.names;
        size_t nelems=__self->atoms->elements.nelems;
        if(f(args,kwds)) return NULL;
        if(f.remap<3,2,1>("elements present in system",names,nelems)) return NULL;
        
        
        size_t nr,nrho;
        type0 dr,drho;
        type0** r_c;
        type0(** F)[5]=NULL;
        type0(*** r_phi)[4]=NULL;
        type0(*** rho)[4]=NULL;
        try
        {
            ImportEAM::setfl(nelems,__self->atoms->elements.names,f.val<0>(),dr,drho,nr,nrho,r_phi,rho,F,r_c);
        }
        catch(char* err_msg)
        {
            PyErr_SetString(PyExc_TypeError,err_msg);
            delete [] err_msg;
            return NULL;
        }
        
        for(size_t i=0;i<nelems;i++)
            if(f.v<1>()[i]>r_c[i][i])
            {
                Memory::dealloc(r_c);
                Memory::dealloc(F);
                Memory::dealloc(r_phi);
                Memory::dealloc(rho);
                PyErr_Format(PyExc_TypeError,"r_crd[%zu] should be less than r_c[%zu][%zu]",i,i,i);
                return NULL;
            }
        
        delete __self->ff;
        __self->ff=new ForceFieldEAMDMD(__self->atoms,dr,drho,nr,nrho,std::move(r_phi),std::move(rho),std::move(F),std::move(r_c),f.mov<1>(),f.mov<2>());
        Py_RETURN_NONE;
    });
    method_1.ml_doc=(char*)R"---(
    ff_eam_setfl(setfl_file,r_crd,C,elems=None)
   
    Tabulated EAM force field given by a single SetFL file
    
    Assigns EAM force field to system
    
    Parameters
    ----------
    setfl_file : string
        relative path to DYNAMO file with SetFL format
    r_crd : double[nelems]
        cutoff radius for mass exchange
    C : double[nelems]
        dimensionless factor for each element, see :ref:`here <master-ref>`
    elems : string[nelems]
        mapping elements
    
    Returns
    -------
    None
   
    Notes
    -----
    This is tabulated form of Embedded Atom Method (EAM) potential
    
    
    Examples
    --------
    Cu
    
    ::
     
        >>> from mapp import dmd
        >>> sim=dmd.cfg(5,"configs/Cu-DMD.cfg")
        >>> sim.ff_eam_setfl("potentials/Cu_mishin.eam.alloy")
    
    

    )---";
    
    method_2.ml_flags=METH_VARARGS | METH_KEYWORDS;
    method_2.ml_name="ff_eam_fs";
    method_2.ml_meth=(PyCFunction)(PyCFunctionWithKeywords)(
    [](PyObject* self,PyObject* args,PyObject* kwds)->PyObject*
    {
        AtomsDMD::Object* __self=reinterpret_cast<AtomsDMD::Object*>(self);
        FuncAPI<std::string,type0*,type0*,std::string*> f("ff_eam_fs",{"fs_file","r_crd","C","elems"});
        f.noptionals=1;
        f.logics<1>()[0]=VLogics("gt",0.0);
        f.logics<2>()[0]=VLogics("gt",0.0);
        
        const std::string* names=__self->atoms->elements.names;
        size_t nelems=__self->atoms->elements.nelems;
        if(f(args,kwds)) return NULL;
        if(f.remap<3,2,1>("elements present in system",names,nelems)) return NULL;
        
        size_t nr,nrho;
        type0 dr,drho;
        type0** r_c;
        type0(** F)[5]=NULL;
        type0(*** r_phi)[4]=NULL;
        type0(*** rho)[4]=NULL;
        try
        {
            ImportEAM::fs(nelems,__self->atoms->elements.names,f.val<0>(),dr,drho,nr,nrho,r_phi,rho,F,r_c);
        }
        catch(char* err_msg)
        {
            PyErr_SetString(PyExc_TypeError,err_msg);
            delete [] err_msg;
            return NULL;
        }
        
        for(size_t i=0;i<nelems;i++)
            if(f.v<1>()[i]>r_c[i][i])
            {
                Memory::dealloc(r_c);
                Memory::dealloc(F);
                Memory::dealloc(r_phi);
                Memory::dealloc(rho);
                PyErr_Format(PyExc_TypeError,"r_crd[%zu] should be less than r_c[%zu][%zu]",i,i,i);
                return NULL;
            }
        
        delete __self->ff;
        __self->ff=new ForceFieldEAMDMD(__self->atoms,dr,drho,nr,nrho,std::move(r_phi),std::move(rho),std::move(F),std::move(r_c),f.mov<1>(),f.mov<2>());
        Py_RETURN_NONE;
    });
    method_2.ml_doc=(char*)R"---(
    ff_eam_fs(fs_file,r_crd,C,elems=None)
   
    Tabulated Finnis-Sinclair EAM
    
    Assigns Finnis-Sinclair EAM force field to system. For explanation of the parameter see the Notes section.
    
    Parameters
    ----------
    fs_file : string
        relative path to DYNAMO file with fs format
    r_crd : double[nelems]
        cutoff radius for mass exchange
    C : double[nelems]
        dimensionless factor for each element, see :ref:`here <master-ref>`
    elems : string[nelems]
        mapping elements
    
    Returns
    -------
    None
   
    Notes
    -----
    This is tabulated form of Finnis-Sinclair Embedded Atom Method (EAM) potential
    
    Consider the general form of EAM potential:
    
    
    .. math:: U=\frac{1}{2}\sum_{i}\sum_{j\neq i}\phi_{\gamma \delta}{(x_{ij}) }+\sum_i E_\gamma \left(\sum_{j\neq i} \rho_{\delta\gamma}(x_{ij}) \right),
    
    
    From here on, greek superscipts/subscripts are used to refer to elements present the system; :math:`\gamma`, and :math:`\delta` denote type of atom :math:`i` and atom :math:`j`, repectively. :math:`E`, :math:`\rho`, and :math:`\phi` are embedding, electron density, and pair functions, respectively. Also :math:`x_{ij}` refers to distance between :math:`i` and :math:`j`. Now the multi component formulation of DMD free energy would be:
    
    .. math::
       F=&\frac{1}{2}\sum_{i,\gamma,j\neq i,\delta}c_i^\gamma c_j^\delta \omega_{\gamma \delta}\left(x_{ij}\right)+\sum_{i,\gamma} c_i^\gamma E_\gamma \left(\sum_{j\neq i,\delta} c_j^{\delta}\psi_{\delta\gamma}\left(x_{ij}\right)\right)-3k_BT\sum_{i,\gamma} c_i^\gamma \log\left(\sqrt{\pi e}\alpha_i^\gamma/\Lambda_\gamma\right)\\
       &+k_BT\sum_{i,\gamma} c_i^\gamma\log c_i^\gamma+k_BT\sum_{i}c_i^v\log (c_i^v),
    
    where
    
    .. math:: \Lambda_\gamma=\frac{h}{\sqrt{2\pi m_\gamma k_BT}}, \quad c_i^v=1-\sum_{\gamma}c_i^\gamma,
    
    .. math:: \omega_{\gamma\delta}(x_{ij})= \frac{1}{\left(\alpha^{\gamma\delta}_{ij}\sqrt{\pi}\right)^{3}}\int d^3\mathbf{x} \exp{\biggl[-\left(\frac{\mathbf{x}_{ij} -\mathbf{x}}{\alpha^{\gamma\delta}_{ij}}\right)^2\biggr]}\phi_{\gamma\delta}(|\mathbf{x}|)
    
    .. math:: \psi_{\gamma\delta}(x_{ij})=\frac{1}{\left(\alpha^{\gamma\delta}_{ij}\sqrt{\pi}\right)^{3}}\int d^3\mathbf{x} \exp{\biggl[-\left(\frac{\mathbf{x}_{ij} -\mathbf{x}}{\alpha^{\gamma\delta}_{ij}}\right)^2\biggr]}\rho_{\gamma\delta}(|\mathbf{x}|)
    
    .. math:: \alpha^{\gamma\delta}_{ij}=\sqrt{{\alpha_i^{\gamma}}^2+{\alpha_j^{\delta}}^2},\quad \mathbf{x}_{ji}=\mathbf{x}_j-\mathbf{x}_i
    
    The numerical evaluation of these integrals is discussed :ref:`here <integ-ref>`.
    
    
    
    Recalling that
    
    .. math:: \langle U \rangle &= \frac{\partial }{\partial \beta}\left(\beta F \right)
    
    the average potential energy is
    
    
    .. math:: \langle U \rangle = \frac{1}{2}\sum_{i,\gamma,j\neq i,\delta}c_i^\gamma c_j^\delta \omega_{\gamma \delta}\left(x_{ij}\right)+\sum_{i,\gamma} c_i^\gamma E_\gamma \left(\sum_{j\neq i,\delta} c_j^{\delta}\psi_{\delta\gamma}\left(x_{ij}\right)\right) -\frac{3}{2} k_B T
    
    
    Examples
    --------
    Iron Hydrogrn mixture
    
    >>> from mapp import dmd
    >>> sim=dmd.cfg(5,"configs/FeH-DMD.cfg")
    >>> sim.ff_eam_fs("potentials/FeH.eam.fs")
    )---";
}


