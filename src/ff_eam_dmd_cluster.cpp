#include "ff_eam_dmd_cluster.h"
#include <stdlib.h>
#include "neighbor_dmd_sc.h"
#include "elements.h"
#include "memory.h"
#include "xmath.h"
#include "atoms_dmd.h"
#include "MAPP.h"
#include "dynamic_dmd.h"
#include "import_eam.h"
#include <limits>
#define PI_IN_SQ 0.564189583547756286948079451561
#define TOL 1.0e-9

#define CONST 1.0
using namespace MAPP_NS;
/*--------------------------------------------
 constructor
 --------------------------------------------*/
ForceFieldEAMDMDCLUSTER::
ForceFieldEAMDMDCLUSTER(AtomsDMD* atoms,
type0 __dr,type0 __drho,size_t __nr,size_t __nrho,
type0(***&& __r_phi_arr)[4],type0(***&& __rho_arr)[4],type0(**&& __F_arr)[5],
type0**&& __cut,type0*&& __r_crd,type0*&& __zeta):
ForceFieldDMD(atoms),
dr(__dr),
drho(__drho),
nr(__nr),
nrho(__nrho),
r_phi_arr(__r_phi_arr),
r_rho_arr(__rho_arr),
F_arr(__F_arr),
rho_phi(NULL),
drho_phi_dr(NULL),
drho_phi_dalpha(NULL),
B_pair(NULL),
is_pair(NULL),
max_pairs(0),
vec0(NULL),
vec1(NULL),
vec2(NULL),
vec3(NULL),
mu_ptr(NULL),
dE_ptr(NULL),
rho_ptr(NULL),
d_ptr(NULL),
b_ptr(NULL),
theta_ptr(NULL),
N(atoms->N)
{
    
    __r_phi_arr=NULL;
    __rho_arr=NULL;
    __F_arr=NULL;
    
    dr_inv=1.0/dr;
    drho_inv=1.0/drho;
    rho_max=static_cast<type0>(nrho)*drho;

    kbT=beta=-1.0;
    

    
    Memory::alloc(xi,N);
    Memory::alloc(wi_0,N);
    Memory::alloc(wi_1,N);
    Memory::alloc(zeta,nelems);
    Memory::alloc(c_0,nelems);
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
    
    delete neighbor;
    neighbor=new NeighborDMDSC(atoms,cut_sk,rsq_crd);
}
/*--------------------------------------------
 destructor
 --------------------------------------------*/
ForceFieldEAMDMDCLUSTER::~ForceFieldEAMDMDCLUSTER()
{
    Memory::dealloc(F_arr);
    Memory::dealloc(r_rho_arr);
    Memory::dealloc(r_phi_arr);

    Memory::dealloc(zeta);
    Memory::dealloc(c_1);
    Memory::dealloc(c_0);
    Memory::dealloc(wi_1);
    Memory::dealloc(wi_0);
    Memory::dealloc(xi);

}
/*--------------------------------------------
 force calculation
 --------------------------------------------*/
void ForceFieldEAMDMDCLUSTER::calc_pair(type0 r,type0* alpha_i,type0* alpha_j,type0*& __rho_phi,type0*&,type0*&)
{}
/*--------------------------------------------
 force calculation
 --------------------------------------------*/
void ForceFieldEAMDMDCLUSTER::force_calc()
{
    if(max_pairs<neighbor->no_pairs)
    {
        delete [] rho_phi;
        delete [] drho_phi_dr;
        delete [] drho_phi_dalpha;
        delete [] B_pair;
        delete [] is_pair;
        
        max_pairs=neighbor->no_pairs;
        size_t no_0=3*max_pairs*c_dim*c_dim;
        Memory::alloc(rho_phi,no_0);
        Memory::alloc(drho_phi_dr,no_0);
        Memory::alloc(drho_phi_dalpha,no_0);
        Memory::alloc(B_pair,c_dim*c_dim*max_pairs);
        Memory::alloc(is_pair,max_pairs);
    }
    
    for(int i=0;i<max_pairs*c_dim*c_dim*3;i++) rho_phi[i]=drho_phi_dr[i]=drho_phi_dalpha[i]=0.0;
    for(int i=0;i<max_pairs*c_dim*c_dim;i++) B_pair[i]=0.0;
    for(int i=0;i<max_pairs;i++) is_pair[i]=false;
    
    
    
    
    
    type0 r,r_inv;
    size_t m;
    type0* coef;
    type0 tmp0,tmp1;
    type0 p,cv_i;
    
    type0 const* c=atoms->c->begin();
    
    
    
    type0* dE=dE_ptr->begin();
    type0* rho=rho_ptr->begin();
    for(int i=0;i<atoms->natms_lcl*c_dim;i++) rho[i]=0.0;
    
    type0 alpha_ij_sq,alpha_ij,rsq;
    
    
    
    const type0* x=atoms->x->begin();
    const type0* alpha=atoms->alpha->begin();
    int** neighbor_list=neighbor->neighbor_list;
    int* neighbor_list_size=neighbor->neighbor_list_size;
    int natms_lcl=atoms->natms_lcl;
    type0 ent_lcl=0.0,vib_lcl=0.0;
    size_t istart=0,iistart;
    for(int i=0;i<natms_lcl;i++)
    {
        const int neigh_sz=neighbor_list_size[i];
        for(int j,__j=0;__j<neigh_sz;__j++)
        {
            j=neighbor_list[i][__j];
            rsq=Algebra::RSQ<__dim__>(x+i*__dim__,x+j*__dim__);
            r=sqrt(rsq);

            for(int ic=0;ic<c_dim;ic++)
                for(int jc=0;jc<c_dim;jc++,istart+=3)
                {
                    alpha_ij_sq=alpha[i*c_dim+ic]*alpha[i*c_dim+ic]+alpha[j*c_dim+jc]*alpha[j*c_dim+jc];
                    alpha_ij=sqrt(alpha_ij_sq);
                    if(r-alpha_ij*xi[N-1]>=cut[ic][jc]) continue;
                    
                    is_pair[istart/(3*c_dim*c_dim)]=true;
                    
                    type0 upper=(r+cut[ic][jc])/alpha_ij;
                    type0 lower=(r-cut[ic][jc])/alpha_ij;
                    type0 __r,p,tmp0,xi2;
                    type0* coef;
                    
                    
                    type0 H[3][3]{[0 ... 2]={[0 ... 2]=0.0}};
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
                        
                        coef=r_phi_arr[ic][jc][m];
                        tmp0=((coef[3]*p+coef[2])*p+coef[1])*p+coef[0];
                        if(__r<0.0) tmp0*=-1.0;
                        tmp0*=wi_0[l];
                        H[0][0]+=tmp0;
                        H[0][1]+=tmp0*xi[l];
                        H[0][2]+=tmp0*xi2;
                        
                        
                        coef=r_rho_arr[ic][jc][m];
                        tmp0=((coef[3]*p+coef[2])*p+coef[1])*p+coef[0];
                        if(__r<0.0) tmp0*=-1.0;
                        tmp0*=wi_0[l];
                        H[1][0]+=tmp0;
                        H[1][1]+=tmp0*xi[l];
                        H[1][2]+=tmp0*xi2;
                        
                        coef=r_rho_arr[jc][ic][m];
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
                    
                    
                    rho[i*c_dim+ic]+=c[j*c_dim+jc]*rho_phi[istart+2];
                    
                    if(j<natms_lcl)
                    {
                        rho[j*c_dim+jc]+=c[i*c_dim+ic]*rho_phi[istart+1];
                        __vec_lcl[0]+=c[i*c_dim+ic]*c[j*c_dim+jc]*rho_phi[istart];
                    }
                    else
                        __vec_lcl[0]+=0.5*c[i*c_dim+ic]*c[j*c_dim+jc]*rho_phi[istart];
                }
        }
        
        
        cv_i=1.0;
        for(int ic=0;ic<c_dim;ic++)
        {
            p=rho[i*c_dim+ic]*drho_inv;
            m=static_cast<size_t>(p);
            m=MIN(m,nrho-2);
            p-=m;
            p=MIN(p,1.0);
            coef=F_arr[ic][m];
            
            tmp0=(((coef[4]*p+coef[3])*p+coef[2])*p+coef[1])*p+coef[0];
            tmp1=(((4.0*coef[4]*p+3.0*coef[3])*p+2.0*coef[2])*p+coef[1])*drho_inv;
            
            if(rho[i*c_dim+ic]>rho_max) tmp0+=tmp1*(rho[i*c_dim+ic]-rho_max);
            
            dE[i*c_dim+ic]=tmp1;
            __vec_lcl[0]+=c[i*c_dim+ic]*tmp0;
            
            
            vib_lcl+=c[i*c_dim+ic]*(c_0[ic]-3.0*kbT*log(alpha[i*c_dim+ic]));
            ent_lcl+=calc_ent(c[i*c_dim+ic]);
            cv_i-=c[i*c_dim+ic];
        }
        
        ent_lcl+=calc_ent(cv_i);
    }
    
    dynamic->update(dE_ptr);
    
    ___sc_loop();
    __vec_lcl[0]+=vib_lcl;
    __vec_lcl[1+__nvoigt__]+=ent_lcl;
    __vec_lcl[0]+=kbT*ent_lcl;
    
    
    
    type0* fvec=f->begin();
    type0* f_alphavec=f_alpha->begin();
    type0 dx_ij[__dim__];
    type0 fpair,apair;
    istart=iistart=0;
    for(int i=0;i<natms_lcl;i++)
    {
        
        const int neigh_sz=neighbor_list_size[i];
        for(int j,__j=0;__j<neigh_sz;__j++)
        {
            j=neighbor_list[i][__j];
            Algebra::DX<__dim__>(x+i*__dim__,x+j*__dim__,dx_ij);
            fpair=0.0;
            
            if(j<natms_lcl)
            {
                for(int ic=0;ic<c_dim;ic++)
                    for(int jc=0;jc<c_dim;jc++,istart+=3,iistart++)
                    {
                        fpair-=(drho_phi_dr[istart+2]*dE[i*c_dim+ic]+drho_phi_dr[istart+1]*dE[j*c_dim+jc]+drho_phi_dr[istart])*B_pair[iistart];
                        apair=-(drho_phi_dalpha[istart+2]*dE[i*c_dim+ic]+drho_phi_dalpha[istart+1]*dE[j*c_dim+jc]+drho_phi_dalpha[istart])*B_pair[iistart];
                        f_alphavec[i*c_dim+ic]+=apair*alpha[i*c_dim+ic];
                        f_alphavec[j*c_dim+jc]+=apair*alpha[j*c_dim+jc];
                        
                    }
                
                Algebra::V_add_x_mul_V<__dim__>(fpair,dx_ij,fvec+i*__dim__);
                Algebra::V_add_x_mul_V<__dim__>(-fpair,dx_ij,fvec+j*__dim__);
            }
            else
            {
                for(int ic=0;ic<c_dim;ic++)
                    for(int jc=0;jc<c_dim;jc++,istart+=3,iistart++)
                    {
                        fpair-=(drho_phi_dr[istart+2]*dE[i*c_dim+ic]+drho_phi_dr[istart+1]*dE[j*c_dim+jc]+drho_phi_dr[istart])*B_pair[iistart];
                        apair=-(drho_phi_dalpha[istart+2]*dE[i*c_dim+ic]+drho_phi_dalpha[istart+1]*dE[j*c_dim+jc]+drho_phi_dalpha[istart])*B_pair[iistart];
                        f_alphavec[i*c_dim+ic]+=apair*alpha[i*c_dim+ic];
                    }
                Algebra::V_add_x_mul_V<__dim__>(fpair,dx_ij,fvec+i*__dim__);
                fpair*=0.5;
            }
            
            Algebra::DyadicV(-fpair,dx_ij,&__vec_lcl[1]);
        }
        
        
        for(int ic=0;ic<c_dim;ic++)
            f_alphavec[i*c_dim+ic]+=3.0*kbT*c[i*c_dim+ic]/alpha[i*c_dim+ic];
        
    }


}
/*--------------------------------------------
 energy calculation
 --------------------------------------------*/
void ForceFieldEAMDMDCLUSTER::energy_calc()
{
    if(max_pairs<neighbor->no_pairs)
    {
        delete [] rho_phi;
        delete [] drho_phi_dr;
        delete [] drho_phi_dalpha;
        delete [] B_pair;
        delete [] is_pair;
        
        max_pairs=neighbor->no_pairs;
        size_t no_0=3*max_pairs*c_dim*c_dim;
        Memory::alloc(rho_phi,no_0);
        Memory::alloc(drho_phi_dr,no_0);
        Memory::alloc(drho_phi_dalpha,no_0);
        Memory::alloc(B_pair,c_dim*c_dim*max_pairs);
        Memory::alloc(is_pair,max_pairs);
    }
    
    for(int i=0;i<max_pairs*c_dim*c_dim*3;i++) rho_phi[i]=drho_phi_dr[i]=drho_phi_dalpha[i]=0.0;
    for(int i=0;i<max_pairs*c_dim*c_dim;i++) B_pair[i]=0.0;
    for(int i=0;i<max_pairs;i++) is_pair[i]=false;
    
    
    
    
    
    type0 r,r_inv;
    size_t m;
    type0* coef;
    type0 tmp0,tmp1;
    type0 p,cv_i;
    
    type0 const* c=atoms->c->begin();
    
    
    
    type0* dE=dE_ptr->begin();
    type0* rho=rho_ptr->begin();
    int natms_lcl=atoms->natms_lcl;
    for(int i=0;i<natms_lcl*c_dim;i++) rho[i]=0.0;
    
    type0 alpha_ij_sq,alpha_ij,rsq;
    
    
    
    const type0* x=atoms->x->begin();
    const type0* alpha=atoms->alpha->begin();
    int** neighbor_list=neighbor->neighbor_list;
    int* neighbor_list_size=neighbor->neighbor_list_size;
    
    type0 ent_lcl=0.0,vib_lcl=0.0;
    size_t istart=0;
    for(int i=0;i<natms_lcl;i++)
    {
        const int neigh_sz=neighbor_list_size[i];
        for(int j,__j=0;__j<neigh_sz;__j++)
        {
            j=neighbor_list[i][__j];
            rsq=Algebra::RSQ<__dim__>(x+i*__dim__,x+j*__dim__);
            r=sqrt(rsq);

            for(int ic=0;ic<c_dim;ic++)
                for(int jc=0;jc<c_dim;jc++,istart+=3)
                {
                    alpha_ij_sq=alpha[i*c_dim+ic]*alpha[i*c_dim+ic]+alpha[j*c_dim+jc]*alpha[j*c_dim+jc];
                    alpha_ij=sqrt(alpha_ij_sq);
                    if(r-alpha_ij*xi[N-1]>=cut[ic][jc]) continue;
                    
                    is_pair[istart/(3*c_dim*c_dim)]=true;
                    type0 upper=(r+cut[ic][jc])/alpha_ij;
                    type0 lower=(r-cut[ic][jc])/alpha_ij;
                    type0 __r,p,tmp0,xi2;
                    type0* coef;
                    
                    
                    type0 H[3][1]{[0 ... 2]={[0 ... 0]=0.0}};
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
                        
                        coef=r_phi_arr[ic][jc][m];
                        tmp0=((coef[3]*p+coef[2])*p+coef[1])*p+coef[0];
                        if(__r<0.0) tmp0*=-1.0;
                        tmp0*=wi_0[l];
                        H[0][0]+=tmp0;

                        
                        
                        coef=r_rho_arr[ic][jc][m];
                        tmp0=((coef[3]*p+coef[2])*p+coef[1])*p+coef[0];
                        if(__r<0.0) tmp0*=-1.0;
                        tmp0*=wi_0[l];
                        H[1][0]+=tmp0;

                        
                        coef=r_rho_arr[jc][ic][m];
                        tmp0=((coef[3]*p+coef[2])*p+coef[1])*p+coef[0];
                        if(__r<0.0) tmp0*=-1.0;
                        tmp0*=wi_0[l];
                        H[2][0]+=tmp0;

                    }
                    
                    r_inv=1.0/r;

                    
                    
                    Algebra::Do<3>::func([&istart,&H,&r_inv,this]
                    (int i)
                    {
                        rho_phi[istart+i]=PI_IN_SQ*r_inv*H[i][0];
                    });
                    
                    
                    rho[i*c_dim+ic]+=c[j*c_dim+jc]*rho_phi[istart+2];
                    
                    if(j<natms_lcl)
                    {
                        rho[j*c_dim+jc]+=c[i*c_dim+ic]*rho_phi[istart+1];
                        __vec_lcl[0]+=c[i*c_dim+ic]*c[j*c_dim+jc]*rho_phi[istart];
                    }
                    else
                        __vec_lcl[0]+=0.5*c[i*c_dim+ic]*c[j*c_dim+jc]*rho_phi[istart];
                }
        }
        
        
        cv_i=1.0;
        for(int ic=0;ic<c_dim;ic++)
        {
            p=rho[i*c_dim+ic]*drho_inv;
            m=static_cast<size_t>(p);
            m=MIN(m,nrho-2);
            p-=m;
            p=MIN(p,1.0);
            coef=F_arr[ic][m];
            
            tmp0=(((coef[4]*p+coef[3])*p+coef[2])*p+coef[1])*p+coef[0];
            tmp1=(((4.0*coef[4]*p+3.0*coef[3])*p+2.0*coef[2])*p+coef[1])*drho_inv;
            
            if(rho[i*c_dim+ic]>rho_max) tmp0+=tmp1*(rho[i*c_dim+ic]-rho_max);
            
            dE[i*c_dim+ic]=tmp1;
            __vec_lcl[0]+=c[i*c_dim+ic]*tmp0;
            
            
            vib_lcl+=c[i*c_dim+ic]*(c_0[ic]-3.0*kbT*log(alpha[i*c_dim+ic]));
            ent_lcl+=calc_ent(c[i*c_dim+ic]);
            cv_i-=c[i*c_dim+ic];
        }
        
        ent_lcl+=calc_ent(cv_i);
    }
    
    dynamic->update(dE_ptr);
    
    ___sc_loop();
    __vec_lcl[0]+=vib_lcl;
    __vec_lcl[1+__nvoigt__]+=ent_lcl;
    __vec_lcl[0]+=kbT*ent_lcl;
}
/*--------------------------------------------
 
 --------------------------------------------*/
void ForceFieldEAMDMDCLUSTER::sc_loop()
{
    
    
    type0* dE=dE_ptr->begin();
    type0* rho=rho_ptr->begin();
    
    type0* coef;
    size_t m,iistart,istart;
    type0 prev_en,curr_en;
    MPI_Allreduce(__vec_lcl,&prev_en,1,Vec<type0>::MPI_T,MPI_SUM,world);
    type0 p,en_diff=1.0,tmp0,tmp1;
    int natms_lcl=atoms->natms_lcl;
    int** neighbor_list=neighbor->neighbor_list;
    int* neighbor_list_size=neighbor->neighbor_list_size;
    type0 const* c=atoms->c->begin();
    type0* b=b_ptr->begin();
    type0 pe_lcl=0.0,ent_corr_lcl=0.0;
    type0 curr_en_lcl=0.0;
    while(en_diff>TOL)
    {
        for(int i=0;i<natms_lcl*c_dim;i++) rho[i]=0.0;
        for(int i=0;i<natms_lcl;i++) b[i]=0.0;
        
        
        istart=0;
        iistart=0;
        pe_lcl=ent_corr_lcl=0.0;
        for(int i=0;i<natms_lcl;i++)
        {
            const int neigh_sz=neighbor_list_size[i];
            for(int j,__j=0;__j<neigh_sz;__j++)
            {
                j=neighbor_list[i][__j];
                
                if(! is_pair[iistart/4])
                {
                    
                    B_pair[iistart]=c[i*c_dim]*c[j*c_dim];
                    B_pair[iistart+1]=c[i*c_dim]*c[j*c_dim+1];
                    B_pair[iistart+2]=c[i*c_dim+1]*c[j*c_dim];
                    B_pair[iistart+3]=c[i*c_dim+1]*c[j*c_dim+1];
                    
                    
                    if(j<natms_lcl)
                    {
                        for(int ic=0;ic<c_dim;ic++)
                        for(int jc=0;jc<c_dim;jc++,istart+=3,iistart++)
                        {
                            pe_lcl+=c[i*c_dim+ic]*c[j*c_dim+jc]*rho_phi[istart];
                            rho[i*c_dim+ic]+=c[j*c_dim+jc]*rho_phi[istart+2];
                            rho[j*c_dim+jc]+=c[i*c_dim+ic]*rho_phi[istart+1];
                            
                        }
                    }
                    else{
                        for(int ic=0;ic<c_dim;ic++)
                        for(int jc=0;jc<c_dim;jc++,istart+=3,iistart++)
                        {
                            pe_lcl+=0.5*c[i*c_dim+ic]*c[j*c_dim+jc]*rho_phi[istart];
                            rho[i*c_dim+ic]+=c[j*c_dim+jc]*rho_phi[istart+2];
                        }
                    }
                    
                    continue;
                }
                
                type0 m=0.0;
                for(int ic=0;ic<c_dim;ic++)
                    for(int jc=0;jc<c_dim;jc++,istart+=3)
                    {
                        if(ic==jc)
                            m+=-beta*(rho_phi[istart+2]*dE[i*c_dim+ic]+rho_phi[istart+1]*dE[j*c_dim+jc]+rho_phi[istart]);
                        else
                            m-=-beta*(rho_phi[istart+2]*dE[i*c_dim+ic]+rho_phi[istart+1]*dE[j*c_dim+jc]+rho_phi[istart]);
                    }
                istart-=3*c_dim*c_dim;
                
                
                type0 u;
                
                if(m!=0.0)
                {
                    
                    type0 __a,__b,__c,__delta;
                    
                    if((c[i*c_dim]==0.0 && c[j*c_dim]==0.0) || (c[i*c_dim]==1.0 && c[j*c_dim]==1.0))
                    {
                        u=c[i*c_dim];
                    }
                    else
                    {
                        __a=0.0;
                        
                        if(m<0.0)
                        {
                            __a=exp(m)-1.0;
                            __b=-1.0-(c[j*c_dim]+c[i*c_dim])*__a;
                            __c=c[j*c_dim]*c[i*c_dim]*(__a+1.0);
                        }
                        else
                        {
                            __a=1.0-exp(-m);
                            __b=-(c[j*c_dim]+c[i*c_dim]-1.0)*__a-1.0;
                            __c=c[j*c_dim]*c[i*c_dim];
                        }
                        
                        type0 __x0=0.0,__x1=0.0;
                        if(__a==0.0)
                        {
                            
                            
                            u=-__c/__b;
                        }
                        else
                        {
                            __delta=sqrt(__b*__b-4.0*__a*__c);
                            
                            __x1=0.5*(-__b+__delta)/__a;
                            __x0=0.5*(-__b-__delta)/__a;
                            if(__x0<=c[i*c_dim] && __x0<=c[j*c_dim] && __x0>=0.0)
                                u=__x0;
                            else
                                u=__x1;
                        }
                        
                    }

                    
                    B_pair[iistart]=u;
                    B_pair[iistart+1]=c[i*c_dim]-u;
                    B_pair[iistart+2]=c[j*c_dim]-u;
                    B_pair[iistart+3]=1.0-c[i*c_dim]-c[j*c_dim]+u;
                    
                    
                    if(j<natms_lcl)
                    {
                        ent_corr_lcl+=calc_ent(B_pair[iistart]);
                        b[i]+=u;
                        b[j]+=u;
                        for(int ic=0;ic<c_dim;ic++)
                            for(int jc=0;jc<c_dim;jc++,istart+=3,iistart++)
                            {
                                pe_lcl+=B_pair[iistart]*rho_phi[istart];
                                if(c[i*c_dim+ic]>0.0)
                                    rho[i*c_dim+ic]+=B_pair[iistart]*rho_phi[istart+2]/c[i*c_dim+ic];
                                if(c[j*c_dim+jc]>0.0)
                                    rho[j*c_dim+jc]+=B_pair[iistart]*rho_phi[istart+1]/c[j*c_dim+jc];
                                
                                
                            }
                    }
                    else
                    {
                        b[i]+=u;
                        ent_corr_lcl+=0.5*calc_ent(B_pair[iistart]);
                        for(int ic=0;ic<c_dim;ic++)
                            for(int jc=0;jc<c_dim;jc++,istart+=3,iistart++)
                            {
                                pe_lcl+=0.5*B_pair[iistart]*rho_phi[istart];
                                if(c[i*c_dim+ic]>0.0)
                                    rho[i*c_dim+ic]+=B_pair[iistart]*rho_phi[istart+2]/c[i*c_dim+ic];

                            }
                    }
                    
                }
                else
                {
                    B_pair[iistart]=c[i*c_dim]*c[j*c_dim];
                    B_pair[iistart+1]=c[i*c_dim]*c[j*c_dim+1];
                    B_pair[iistart+2]=c[i*c_dim+1]*c[j*c_dim];
                    B_pair[iistart+3]=c[i*c_dim+1]*c[j*c_dim+1];
                    
                    
                    if(j<natms_lcl)
                    {
                        b[i]+=c[i*c_dim]*c[j*c_dim];
                        b[j]+=c[i*c_dim]*c[j*c_dim];
                        for(int ic=0;ic<c_dim;ic++)
                            for(int jc=0;jc<c_dim;jc++,istart+=3,iistart++)
                            {
                                pe_lcl+=c[i*c_dim+ic]*c[j*c_dim+jc]*rho_phi[istart];
                                rho[i*c_dim+ic]+=c[j*c_dim+jc]*rho_phi[istart+2];
                                rho[j*c_dim+jc]+=c[i*c_dim+ic]*rho_phi[istart+1];
                                
                            }
                    }
                    else
                    {
                        b[i]+=c[i*c_dim]*c[j*c_dim];
                        for(int ic=0;ic<c_dim;ic++)
                            for(int jc=0;jc<c_dim;jc++,istart+=3,iistart++)
                            {
                                pe_lcl+=0.5*c[i*c_dim+ic]*c[j*c_dim+jc]*rho_phi[istart];
                                rho[i*c_dim+ic]+=c[j*c_dim+jc]*rho_phi[istart+2];
                            }
                    }
                }

                
            }
            
            for(int ic=0;ic<c_dim;ic++)
            {
                
                p=rho[i*c_dim+ic]*drho_inv;
                m=static_cast<size_t>(p);
                m=MIN(m,nrho-2);
                p-=m;
                p=MIN(p,1.0);
                coef=F_arr[ic][m];
                
                tmp0=(((coef[4]*p+coef[3])*p+coef[2])*p+coef[1])*p+coef[0];
                tmp1=(((4.0*coef[4]*p+3.0*coef[3])*p+2.0*coef[2])*p+coef[1])*drho_inv;
                
                if(rho[i*c_dim+ic]>rho_max) tmp0+=tmp1*(rho[i*c_dim+ic]-rho_max);
                
                dE[i*c_dim+ic]=tmp1;
                pe_lcl+=c[i*c_dim+ic]*tmp0;
                
            }
        }
        
        
        
        
        curr_en_lcl=pe_lcl+ent_corr_lcl*kbT;
        MPI_Allreduce(&curr_en_lcl,&curr_en,1,Vec<type0>::MPI_T,MPI_SUM,world);
        
        en_diff=fabs(curr_en-prev_en);
        prev_en=curr_en;
        dynamic->update(dE_ptr);
    }
    
    
    dynamic->update(b_ptr);
    type0* theta=theta_ptr->begin();
    type0* d=d_ptr->begin();
    for(int i=0;i<natms_lcl+atoms->natms_ph;i++) theta[i]=c[2*i];
    en_diff=1.0;
    type0 curr_fe_lcl=0.0,curr_fe=0.0;

    for(int i=0;i<natms_lcl;i++) d[i]=0.0;
    iistart=0;
    for(int i=0;i<natms_lcl;i++)
    {
        const int neigh_sz=neighbor_list_size[i];
        for(int j,__j=0;__j<neigh_sz;__j++,iistart++)
        {
            j=neighbor_list[i][__j];
            if(!  is_pair[iistart]) continue;
            d[i]+=theta[j];
            if(j<natms_lcl)
            d[j]+=theta[i];
        }
        if(theta[i]!=0.0 || theta[i]!=1.0)
        {
            curr_fe_lcl+=(1.0-c[i*c_dim])*log(1.0-theta[i])+(c[i*c_dim]-b[i])*log(theta[i]);
        }
        curr_fe_lcl+=0.5*(theta[i]*d[i]-b[i]);
    }
    
    MPI_Allreduce(&curr_fe_lcl,&curr_fe,1,Vec<type0>::MPI_T,MPI_SUM,world);
    dynamic->update(d_ptr);
    type0 prev_fe=curr_fe;

    
    while(en_diff>TOL)
    {
        
        
    
        for(int i=0;i<natms_lcl;i++)
        {
            if(d[i]!=0.0)
            {
                
                type0 delta_sq=(d[i]+b[i]-1.0)*(d[i]+b[i]-1.0)-4.0*(b[i]-c[i*c_dim])*d[i];
                
                type0 delta=sqrt(delta_sq);
                type0 theta_i0=0.5*(d[i]+b[i]-1.0-delta)/d[i];
                type0 theta_i1=0.5*(d[i]+b[i]-1.0+delta)/d[i];
                bool theta0_acc=false;
                bool theta1_acc=false;
                if(0.0<=theta_i0 && theta_i0<=1.0)
                    theta0_acc=true;
                if(0.0<=theta_i1 && theta_i1<=1.0)
                    theta1_acc=true;
                if(theta0_acc && theta1_acc)
                {
                    type0 d0=fabs(c[i*c_dim]-theta_i0);
                    type0 d1=fabs(c[i*c_dim]-theta_i1);
                    
                    if(d0<d1)
                        theta[i]=theta_i0;
                    else theta[i]=theta_i1;
                }
                else
                {
                    if(theta0_acc)
                        theta[i]=theta_i0;
                    else
                        theta[i]=theta_i1;
                }
                
                //printf("theta 1 is %e %e\n",theta_i0,theta_i1);
                
            }
            else
                theta[i]=(c[i*c_dim]-b[i])/(1.0-b[i]);
            
            
            
        }
        
        
        for(int i=0;i<natms_lcl;i++) d[i]=0.0;
        curr_fe_lcl=0.0;
        iistart=0;
        for(int i=0;i<natms_lcl;i++)
        {
            const int neigh_sz=neighbor_list_size[i];
            for(int j,__j=0;__j<neigh_sz;__j++,iistart++)
            {
                j=neighbor_list[i][__j];
                if(!  is_pair[iistart]) continue;
                d[i]+=theta[j];
                if(j<natms_lcl)
                    d[j]+=theta[i];
            }
            
            if(theta[i]!=0.0 && theta[i]!=1.0)
            {
                curr_fe_lcl+=(1.0-c[i*c_dim])*log(1.0-theta[i])+(c[i*c_dim]-b[i])*log(theta[i]);
                //printf("kkkkunt\n");
                
                
            }
            curr_fe_lcl+=0.5*(theta[i]*d[i]-b[i]);
            //printf("theta[%d] %e %e %e| %e %e\n",i,theta[i],d[i],b[i],curr_fe_lcl,0.5*(theta[i]*d[i]-b[i]));
        }
        
        dynamic->update(d_ptr);
        MPI_Allreduce(&curr_fe_lcl,&curr_fe,1,Vec<type0>::MPI_T,MPI_SUM,world);
        
        en_diff=fabs(curr_fe-prev_fe);
        prev_fe=curr_fe;
    }
    
    
    
    for(int i=0;i<natms_lcl;i++)
        printf("theta[%d] %e %e %e | %e\n",i,theta[i],theta[i]*d[i],b[i],theta[i]+b[i]*(1.0-theta[i])-theta[i]*(1.0-theta[i])*d[i]);

    __vec_lcl[1+__nvoigt__]=ent_corr_lcl+curr_fe_lcl;
    
    
    
    
    __vec_lcl[0]=pe_lcl+kbT*(ent_corr_lcl+curr_fe_lcl);
    
    
    
    
    
}
/*--------------------------------------------
 
 --------------------------------------------*/
void ForceFieldEAMDMDCLUSTER::_sc_loop()
{
    
    
    type0* dE=dE_ptr->begin();
    type0* rho=rho_ptr->begin();
    
    type0* coef;
    size_t m,iistart,istart;
    type0 prev_en,curr_en;
    MPI_Allreduce(__vec_lcl,&prev_en,1,Vec<type0>::MPI_T,MPI_SUM,world);
    type0 p,en_diff=1.0,tmp0,tmp1;
    int natms_lcl=atoms->natms_lcl;
    int** neighbor_list=neighbor->neighbor_list;
    int* neighbor_list_size=neighbor->neighbor_list_size;
    type0 const* c=atoms->c->begin();
    type0 pe_lcl=0.0,ent_corr_lcl=0.0;
    type0 curr_en_lcl=0.0;
    while(en_diff>TOL)
    {
        for(int i=0;i<natms_lcl*c_dim;i++) rho[i]=0.0;
        
        
        istart=0;
        iistart=0;
        pe_lcl=ent_corr_lcl=0.0;
        for(int i=0;i<natms_lcl;i++)
        {
            const int neigh_sz=neighbor_list_size[i];
            for(int j,__j=0;__j<neigh_sz;__j++)
            {
                j=neighbor_list[i][__j];
                
                if(! is_pair[iistart/4])
                {
                    
                    B_pair[iistart]=c[i*c_dim]*c[j*c_dim];
                    B_pair[iistart+1]=c[i*c_dim]*c[j*c_dim+1];
                    B_pair[iistart+2]=c[i*c_dim+1]*c[j*c_dim];
                    B_pair[iistart+3]=c[i*c_dim+1]*c[j*c_dim+1];
                    
                    
                    if(j<natms_lcl)
                    {
                        for(int ic=0;ic<c_dim;ic++)
                        for(int jc=0;jc<c_dim;jc++,istart+=3,iistart++)
                        {
                            pe_lcl+=c[i*c_dim+ic]*c[j*c_dim+jc]*rho_phi[istart];
                            rho[i*c_dim+ic]+=c[j*c_dim+jc]*rho_phi[istart+2];
                            rho[j*c_dim+jc]+=c[i*c_dim+ic]*rho_phi[istart+1];
                            
                        }
                    }
                    else{
                        for(int ic=0;ic<c_dim;ic++)
                        for(int jc=0;jc<c_dim;jc++,istart+=3,iistart++)
                        {
                            pe_lcl+=0.5*c[i*c_dim+ic]*c[j*c_dim+jc]*rho_phi[istart];
                            rho[i*c_dim+ic]+=c[j*c_dim+jc]*rho_phi[istart+2];
                        }
                    }
                    
                    continue;
                }
                
                type0 m=0.0;
                for(int ic=0;ic<c_dim;ic++)
                    for(int jc=0;jc<c_dim;jc++,istart+=3)
                    {
                        if(ic==jc)
                            m+=-beta*(rho_phi[istart+2]*dE[i*c_dim+ic]+rho_phi[istart+1]*dE[j*c_dim+jc]+rho_phi[istart]);
                        else
                            m-=-beta*(rho_phi[istart+2]*dE[i*c_dim+ic]+rho_phi[istart+1]*dE[j*c_dim+jc]+rho_phi[istart]);
                    }
                istart-=3*c_dim*c_dim;
                
                
                type0 u;
                
                if(m!=0.0)
                {
                    
                    type0 __a,__b,__c,__delta;
                    
                    if((c[i*c_dim]==0.0 && c[j*c_dim]==0.0) || (c[i*c_dim]==1.0 && c[j*c_dim]==1.0))
                    {
                        u=c[i*c_dim];
                    }
                    else
                    {
                        __a=0.0;
                        
                        if(m<0.0)
                        {
                            __a=exp(m)-1.0;
                            __b=-1.0-(c[j*c_dim]+c[i*c_dim])*__a;
                            __c=c[j*c_dim]*c[i*c_dim]*(__a+1.0);
                        }
                        else
                        {
                            __a=1.0-exp(-m);
                            __b=-(c[j*c_dim]+c[i*c_dim]-1.0)*__a-1.0;
                            __c=c[j*c_dim]*c[i*c_dim];
                        }
                        
                        type0 __x0=0.0,__x1=0.0;
                        if(__a==0.0)
                        {
                            
                            
                            u=-__c/__b;
                        }
                        else
                        {
                            __delta=sqrt(__b*__b-4.0*__a*__c);
                            
                            __x1=0.5*(-__b+__delta)/__a;
                            __x0=0.5*(-__b-__delta)/__a;
                            if(__x0<=c[i*c_dim] && __x0<=c[j*c_dim] && __x0>=0.0)
                                u=__x0;
                            else
                                u=__x1;
                        }
                        
                    }

                    
                    B_pair[iistart]=u;
                    B_pair[iistart+1]=c[i*c_dim]-u;
                    B_pair[iistart+2]=c[j*c_dim]-u;
                    B_pair[iistart+3]=1.0-c[i*c_dim]-c[j*c_dim]+u;
                    
                    type0 eps_ij=-1.0;
                    if(c[i*c_dim]*c[i*c_dim]!=0.0)
                        eps_ij=u/(c[i*c_dim]*c[j*c_dim])-1.0;
                    //printf("%e\n",eps_ij);
                    if(j<natms_lcl)
                    {
                        ent_corr_lcl+=(calc_ent(1.0+eps_ij)-eps_ij)*(c[i*c_dim]*c[j*c_dim]);

                        for(int ic=0;ic<c_dim;ic++)
                            for(int jc=0;jc<c_dim;jc++,istart+=3,iistart++)
                            {
                                pe_lcl+=B_pair[iistart]*rho_phi[istart];
                                if(c[i*c_dim+ic]>0.0)
                                    rho[i*c_dim+ic]+=B_pair[iistart]*rho_phi[istart+2]/c[i*c_dim+ic];
                                if(c[j*c_dim+jc]>0.0)
                                    rho[j*c_dim+jc]+=B_pair[iistart]*rho_phi[istart+1]/c[j*c_dim+jc];
                                
                                
                            }
                    }
                    else
                    {
                        ent_corr_lcl+=0.5*(calc_ent(1.0+eps_ij)-eps_ij)*(c[i*c_dim]*c[j*c_dim]);
                        
                        
                        for(int ic=0;ic<c_dim;ic++)
                            for(int jc=0;jc<c_dim;jc++,istart+=3,iistart++)
                            {
                                pe_lcl+=0.5*B_pair[iistart]*rho_phi[istart];
                                if(c[i*c_dim+ic]>0.0)
                                    rho[i*c_dim+ic]+=B_pair[iistart]*rho_phi[istart+2]/c[i*c_dim+ic];

                            }
                    }
                    
                }
                else
                {
                    B_pair[iistart]=c[i*c_dim]*c[j*c_dim];
                    B_pair[iistart+1]=c[i*c_dim]*c[j*c_dim+1];
                    B_pair[iistart+2]=c[i*c_dim+1]*c[j*c_dim];
                    B_pair[iistart+3]=c[i*c_dim+1]*c[j*c_dim+1];
                    
                    
                    if(j<natms_lcl)
                    {
                        for(int ic=0;ic<c_dim;ic++)
                            for(int jc=0;jc<c_dim;jc++,istart+=3,iistart++)
                            {
                                pe_lcl+=c[i*c_dim+ic]*c[j*c_dim+jc]*rho_phi[istart];
                                rho[i*c_dim+ic]+=c[j*c_dim+jc]*rho_phi[istart+2];
                                rho[j*c_dim+jc]+=c[i*c_dim+ic]*rho_phi[istart+1];
                                
                            }
                    }
                    else
                    {
                        for(int ic=0;ic<c_dim;ic++)
                            for(int jc=0;jc<c_dim;jc++,istart+=3,iistart++)
                            {
                                pe_lcl+=0.5*c[i*c_dim+ic]*c[j*c_dim+jc]*rho_phi[istart];
                                rho[i*c_dim+ic]+=c[j*c_dim+jc]*rho_phi[istart+2];
                            }
                    }
                }

                
            }
            
            for(int ic=0;ic<c_dim;ic++)
            {
                
                p=rho[i*c_dim+ic]*drho_inv;
                m=static_cast<size_t>(p);
                m=MIN(m,nrho-2);
                p-=m;
                p=MIN(p,1.0);
                coef=F_arr[ic][m];
                
                tmp0=(((coef[4]*p+coef[3])*p+coef[2])*p+coef[1])*p+coef[0];
                tmp1=(((4.0*coef[4]*p+3.0*coef[3])*p+2.0*coef[2])*p+coef[1])*drho_inv;
                
                if(rho[i*c_dim+ic]>rho_max) tmp0+=tmp1*(rho[i*c_dim+ic]-rho_max);
                
                dE[i*c_dim+ic]=tmp1;
                pe_lcl+=c[i*c_dim+ic]*tmp0;
                
            }
        }
        
        
        
        
        curr_en_lcl=pe_lcl+ent_corr_lcl*kbT;
        MPI_Allreduce(&curr_en_lcl,&curr_en,1,Vec<type0>::MPI_T,MPI_SUM,world);
        
        en_diff=fabs(curr_en-prev_en);
        prev_en=curr_en;
        dynamic->update(dE_ptr);
    }
    
    
    __vec_lcl[1+__nvoigt__]=ent_corr_lcl;

    __vec_lcl[0]=pe_lcl+kbT*ent_corr_lcl;
    
}
/*--------------------------------------------
 
 --------------------------------------------*/
void ForceFieldEAMDMDCLUSTER::__sc_loop()
{
    
    
    type0* dE=dE_ptr->begin();
    type0* rho=rho_ptr->begin();
    
    type0* coef;
    size_t m,iistart,istart;
    type0 ent_corr_lcl=0.0,prev_en,curr_en;
    MPI_Allreduce(__vec_lcl,&prev_en,1,Vec<type0>::MPI_T,MPI_SUM,world);
    type0 p,en_diff=1.0,tmp0,tmp1;
    int natms_lcl=atoms->natms_lcl;
    int** neighbor_list=neighbor->neighbor_list;
    int* neighbor_list_size=neighbor->neighbor_list_size;
    type0 const* c=atoms->c->begin();
    type0 tt=0.0;
    while(en_diff>TOL)
    {
        for(int i=0;i<natms_lcl*c_dim;i++) rho[i]=0.0;
        
        __vec_lcl[0]=__vec_lcl[1+__nvoigt__]=0.0;
        istart=0;
        iistart=0;
        ent_corr_lcl=0.0;
        for(int i=0;i<natms_lcl;i++)
        {
            const int neigh_sz=neighbor_list_size[i];
            for(int j,__j=0;__j<neigh_sz;__j++)
            {
                j=neighbor_list[i][__j];
                
                type0 oo=c[i*c_dim]*c[j*c_dim]*(1.0-c[i*c_dim])*(1.0-c[j*c_dim]);
                type0 delta_min=MAX(-c[i*c_dim]*c[j*c_dim],-(1.0-c[i*c_dim])*(1.0-c[j*c_dim]));
                type0 delta_max=MIN(c[i*c_dim]*(1.0-c[j*c_dim]),(1.0-c[i*c_dim])*c[j*c_dim]);
                
                type0 m=0.0;
                for(int ic=0;ic<c_dim;ic++)
                    for(int jc=0;jc<c_dim;jc++,istart+=3)
                    {
                        if(ic==jc)
                            m+=-beta*(rho_phi[istart+2]*dE[i*c_dim+ic]+rho_phi[istart+1]*dE[j*c_dim+jc]+rho_phi[istart]);
                        else
                            m-=-beta*(rho_phi[istart+2]*dE[i*c_dim+ic]+rho_phi[istart+1]*dE[j*c_dim+jc]+rho_phi[istart]);
                    }
                
                istart-=3*c_dim*c_dim;
                
                m*=oo;
                type0 delta=m;
                if(delta>delta_max)
                    delta=delta_max;
                if(delta<delta_min)
                    delta=delta_min;
                
                
                if(delta!=0.0)
                {
                    

                    B_pair[iistart]=c[i*c_dim]*c[j*c_dim]+delta;
                    B_pair[iistart+1]=c[i*c_dim]*(1.0-c[j*c_dim])-delta;
                    B_pair[iistart+2]=c[j*c_dim]*(1.0-c[i*c_dim])-delta;
                    B_pair[iistart+3]=(1.0-c[i*c_dim])*(1.0-c[j*c_dim])+delta;
                        
                    
                    if(j<natms_lcl)
                    {
                        if(oo!=0.0)
                            ent_corr_lcl+=0.5*oo*log(1.0+delta*delta/(oo*oo));
                            
                        for(int ic=0;ic<c_dim;ic++)
                            for(int jc=0;jc<c_dim;jc++,istart+=3,iistart++)
                            {
                                __vec_lcl[0]+=B_pair[iistart]*rho_phi[istart];
                                if(c[i*c_dim+ic]>0.0)
                                    rho[i*c_dim+ic]+=B_pair[iistart]*rho_phi[istart+2]/c[i*c_dim+ic];
                                if(c[j*c_dim+jc]>0.0)
                                    rho[j*c_dim+jc]+=B_pair[iistart]*rho_phi[istart+1]/c[j*c_dim+jc];
                            }
                    }
                    else
                    {
                        if(oo!=0.0)
                            ent_corr_lcl+=0.25*oo*log(1.0+delta*delta/(oo*oo));
                        
                        for(int ic=0;ic<c_dim;ic++)
                            for(int jc=0;jc<c_dim;jc++,istart+=3,iistart++)
                            {
                                __vec_lcl[0]+=0.5*B_pair[iistart]*rho_phi[istart];
                                if(c[i*c_dim+ic]>0.0)
                                    rho[i*c_dim+ic]+=B_pair[iistart]*rho_phi[istart+2]/c[i*c_dim+ic];
                            }
                    }
                    
                }
                else
                {
                    B_pair[iistart]=c[i*c_dim]*c[j*c_dim];
                    B_pair[iistart+1]=c[i*c_dim]*c[j*c_dim+1];
                    B_pair[iistart+2]=c[i*c_dim+1]*c[j*c_dim];
                    B_pair[iistart+3]=c[i*c_dim+1]*c[j*c_dim+1];
                    
                    
                    if(j<natms_lcl)
                    {
                        for(int ic=0;ic<c_dim;ic++)
                            for(int jc=0;jc<c_dim;jc++,istart+=3,iistart++)
                            {
                                __vec_lcl[0]+=c[i*c_dim+ic]*c[j*c_dim+jc]*rho_phi[istart];
                                rho[i*c_dim+ic]+=c[j*c_dim+jc]*rho_phi[istart+2];
                                rho[j*c_dim+jc]+=c[i*c_dim+ic]*rho_phi[istart+1];
                                
                            }
                    }
                    else{
                        for(int ic=0;ic<c_dim;ic++)
                            for(int jc=0;jc<c_dim;jc++,istart+=3,iistart++)
                            {
                                __vec_lcl[0]+=0.5*c[i*c_dim+ic]*c[j*c_dim+jc]*rho_phi[istart];
                                rho[i*c_dim+ic]+=c[j*c_dim+jc]*rho_phi[istart+2];
                            }
                    }
                }
                
                
            }
            
            for(int ic=0;ic<c_dim;ic++)
            {
                
                p=rho[i*c_dim+ic]*drho_inv;
                m=static_cast<size_t>(p);
                m=MIN(m,nrho-2);
                p-=m;
                p=MIN(p,1.0);
                coef=F_arr[ic][m];
                
                tmp0=(((coef[4]*p+coef[3])*p+coef[2])*p+coef[1])*p+coef[0];
                tmp1=(((4.0*coef[4]*p+3.0*coef[3])*p+2.0*coef[2])*p+coef[1])*drho_inv;
                
                if(rho[i*c_dim+ic]>rho_max) tmp0+=tmp1*(rho[i*c_dim+ic]-rho_max);
                
                dE[i*c_dim+ic]=tmp1;
                __vec_lcl[0]+=c[i*c_dim+ic]*tmp0;
                
            }
        }
        
        
        
        
        __vec_lcl[1+__nvoigt__]+=ent_corr_lcl;
        __vec_lcl[0]+=ent_corr_lcl*kbT;
        
        MPI_Allreduce(__vec_lcl,&curr_en,1,Vec<type0>::MPI_T,MPI_SUM,world);
        en_diff=fabs(curr_en-prev_en);
        //printf("%e %e\n",curr_en,prev_en);
        prev_en=curr_en;
        dynamic->update(dE_ptr);
    }
    
    
    atoms->BB=tt;
}
/*--------------------------------------------
 
 --------------------------------------------*/
void ForceFieldEAMDMDCLUSTER::___sc_loop()
{
    
    
    type0* dE=dE_ptr->begin();
    type0* rho=rho_ptr->begin();
    
    type0* coef;
    
    size_t m,iistart,istart;
    type0 ent_corr_lcl=0.0,prev_en,curr_en;
    MPI_Allreduce(__vec_lcl,&prev_en,1,Vec<type0>::MPI_T,MPI_SUM,world);
    type0 p,en_diff=1.0,tmp0,tmp1;
    int natms_lcl=atoms->natms_lcl;
    int** neighbor_list=neighbor->neighbor_list;
    int* neighbor_list_size=neighbor->neighbor_list_size;
    type0 const* c=atoms->c->begin();
    type0* x=atoms->x->begin();
    
    type0 tt=0.0;
    while(en_diff>TOL)
    {
        for(int i=0;i<natms_lcl*c_dim;i++) rho[i]=0.0;
        
        __vec_lcl[0]=__vec_lcl[1+__nvoigt__]=0.0;
        istart=0;
        iistart=0;
        ent_corr_lcl=0.0;
        for(int i=0;i<natms_lcl;i++)
        {
            const int neigh_sz=neighbor_list_size[i];
            for(int j,__j=0;__j<neigh_sz;__j++)
            {
                j=neighbor_list[i][__j];
                
                type0 m=0.0,u;
                for(int ic=0;ic<c_dim;ic++)
                    for(int jc=0;jc<c_dim;jc++,istart+=3)
                    {
                        if(ic==jc)
                            m-=-beta*(rho_phi[istart+2]*dE[i*c_dim+ic]+rho_phi[istart+1]*dE[j*c_dim+jc]+rho_phi[istart]);
                        else
                            m+=-beta*(rho_phi[istart+2]*dE[i*c_dim+ic]+rho_phi[istart+1]*dE[j*c_dim+jc]+rho_phi[istart]);
                    }
                istart-=3*c_dim*c_dim;
                
                
                
                type0 rsq=Algebra::RSQ<__dim__>(x+i*__dim__,x+j*__dim__)/(atoms->H[0][0]*atoms->H[0][0]);
                if(m!=0.0 && rsq<.9)
                {
                    
                    type0 __a,__b,__c,__delta;
                    
                    if((c[i*c_dim]==0.0 && c[j*c_dim]==0.0) || (c[i*c_dim]==1.0 && c[j*c_dim]==1.0))
                    {
                        u=c[i*c_dim];
                    }
                    else
                    {
                        __a=0.0;
                        
                        if(m<0.0)
                        {
                            __a=exp(m)-1.0;
                            __b=-1.0-(c[j*c_dim]+c[i*c_dim])*__a;
                            __c=c[j*c_dim]*c[i*c_dim]*(__a+1.0);
                        }
                        else
                        {
                            __a=1.0-exp(-m);
                            __b=-(c[j*c_dim]+c[i*c_dim]-1.0)*__a-1.0;
                            __c=c[j*c_dim]*c[i*c_dim];
                        }
                        //printf("%e\n",__a);
                        /*
                         if(__a==-1.0 || __a==1.0)
                         printf("%e %e\n",__a,m);
                         */
                        
                        type0 __x0=0.0,__x1=0.0;
                        if(__a==0.0)
                        {
                            
                            
                            u=-__c/__b;
                        }
                        else
                        {
                            __delta=sqrt(__b*__b-4.0*__a*__c);
                            
                            __x1=0.5*(-__b+__delta)/__a;
                            __x0=0.5*(-__b-__delta)/__a;
                            if(__x0<=c[i*c_dim] && __x0<=c[j*c_dim] && __x0>=0.0)
                                u=__x0;
                            else
                                u=__x1;
                        }
                        
                        /*
                         
                         if(c[0]==1.0-0.031250000)
                         printf(">>>> %e %e %e | %e %e\n",u,c[0],m,__x0,__x1);
                         */
                        
                    }
                    
                    if(i==0 && j==1)
                    {
                        if(c[i*c_dim]==1.0)
                            tt=1.0;
                        else
                            tt=(1.0-c[i*c_dim]-c[j*c_dim]+u)/((1.0-c[i*c_dim])*(1.0-c[i*c_dim]));
                        
                    }
                    
                    B_pair[iistart]=u;
                    B_pair[iistart+1]=c[i*c_dim]-u;
                    B_pair[iistart+2]=c[j*c_dim]-u;
                    B_pair[iistart+3]=1.0-c[i*c_dim]-c[j*c_dim]+u;
                    
                    
                    if(j<natms_lcl)
                    {
                        for(int ic=0;ic<c_dim;ic++)
                            for(int jc=0;jc<c_dim;jc++,istart+=3,iistart++)
                            {
                                __vec_lcl[0]+=B_pair[iistart]*rho_phi[istart];
                                if(c[i*c_dim+ic]>0.0)
                                    rho[i*c_dim+ic]+=B_pair[iistart]*rho_phi[istart+2]/c[i*c_dim+ic];
                                if(c[j*c_dim+jc]>0.0)
                                    rho[j*c_dim+jc]+=B_pair[iistart]*rho_phi[istart+1]/c[j*c_dim+jc];
                                
                                
                                //ent_corr_lcl+=4.0*(c[i*c_dim+ic]*c[j*c_dim+jc]-B_pair[iistart]);
                                
                                
                                if(B_pair[iistart]>0.0 && c[i*c_dim+ic]*c[j*c_dim+jc]>0.0)
                                {
                                    ent_corr_lcl-=B_pair[iistart]*log(B_pair[iistart]/(c[i*c_dim+ic]*c[j*c_dim+jc]));
                                }
                                
                            }
                    }
                    else
                    {
                        for(int ic=0;ic<c_dim;ic++)
                            for(int jc=0;jc<c_dim;jc++,istart+=3,iistart++)
                            {
                                __vec_lcl[0]+=0.5*B_pair[iistart]*rho_phi[istart];
                                if(c[i*c_dim+ic]>0.0)
                                    rho[i*c_dim+ic]+=B_pair[iistart]*rho_phi[istart+2]/c[i*c_dim+ic];
                                
                                //ent_corr_lcl+=2.0*(c[i*c_dim+ic]*c[j*c_dim+jc]-B_pair[iistart]);
                                
                                
                                if(B_pair[iistart]>0.0 && c[i*c_dim+ic]*c[j*c_dim+jc]>0.0)
                                {
                                    ent_corr_lcl-=0.5*B_pair[iistart]*log(B_pair[iistart]/(c[i*c_dim+ic]*c[j*c_dim+jc]));
                                }
                            }
                    }
                    
                }
                else
                {
                    B_pair[iistart]=c[i*c_dim]*c[j*c_dim];
                    B_pair[iistart+1]=c[i*c_dim]*c[j*c_dim+1];
                    B_pair[iistart+2]=c[i*c_dim+1]*c[j*c_dim];
                    B_pair[iistart+3]=c[i*c_dim+1]*c[j*c_dim+1];
                    
                    
                    if(j<natms_lcl)
                    {
                        for(int ic=0;ic<c_dim;ic++)
                            for(int jc=0;jc<c_dim;jc++,istart+=3,iistart++)
                            {
                                __vec_lcl[0]+=c[i*c_dim+ic]*c[j*c_dim+jc]*rho_phi[istart];
                                rho[i*c_dim+ic]+=c[j*c_dim+jc]*rho_phi[istart+2];
                                rho[j*c_dim+jc]+=c[i*c_dim+ic]*rho_phi[istart+1];
                                
                            }
                    }
                    else{
                        for(int ic=0;ic<c_dim;ic++)
                            for(int jc=0;jc<c_dim;jc++,istart+=3,iistart++)
                            {
                                __vec_lcl[0]+=0.5*c[i*c_dim+ic]*c[j*c_dim+jc]*rho_phi[istart];
                                rho[i*c_dim+ic]+=c[j*c_dim+jc]*rho_phi[istart+2];
                            }
                    }
                }
                
                
            }
            
            for(int ic=0;ic<c_dim;ic++)
            {
                
                p=rho[i*c_dim+ic]*drho_inv;
                m=static_cast<size_t>(p);
                m=MIN(m,nrho-2);
                p-=m;
                p=MIN(p,1.0);
                coef=F_arr[ic][m];
                
                tmp0=(((coef[4]*p+coef[3])*p+coef[2])*p+coef[1])*p+coef[0];
                tmp1=(((4.0*coef[4]*p+3.0*coef[3])*p+2.0*coef[2])*p+coef[1])*drho_inv;
                
                if(rho[i*c_dim+ic]>rho_max) tmp0+=tmp1*(rho[i*c_dim+ic]-rho_max);
                
                dE[i*c_dim+ic]=tmp1;
                __vec_lcl[0]+=c[i*c_dim+ic]*tmp0;
                
            }
        }
        
        
        
        
        __vec_lcl[1+__nvoigt__]+=ent_corr_lcl;
        __vec_lcl[0]+=ent_corr_lcl*kbT;
        
        MPI_Allreduce(__vec_lcl,&curr_en,1,Vec<type0>::MPI_T,MPI_SUM,world);
        en_diff=fabs(curr_en-prev_en);
        //printf("%e %e\n",curr_en,prev_en);
        prev_en=curr_en;
        dynamic->update(dE_ptr);
    }
    
    
    atoms->BB=tt;
    
    
    
}
/*--------------------------------------------
 force calculation
 --------------------------------------------*/
type0 ForceFieldEAMDMDCLUSTER::prep(VecTens<type0,2>& f)
{ return 0.0;}
/*--------------------------------------------
 
 --------------------------------------------*/
void ForceFieldEAMDMDCLUSTER::J(VecTens<type0,2>& Dx,VecTens<type0,2>& ADx)
{}
/*--------------------------------------------
 init
 --------------------------------------------*/
void ForceFieldEAMDMDCLUSTER::init()
{
    pre_init();
    set_temp();
    
    
    mu_ptr=new DMDVec<type0>(atoms,0.0,"mu");
    d_ptr=new Vec<type0>(atoms,1);
    b_ptr=new Vec<type0>(atoms,1);
    theta_ptr=new Vec<type0>(atoms,1);
    E_ptr=new Vec<type0>(atoms,c_dim);
    dE_ptr=new Vec<type0>(atoms,c_dim);
    rho_ptr=new Vec<type0>(atoms,c_dim);
}
/*--------------------------------------------
 fin
 --------------------------------------------*/
void ForceFieldEAMDMDCLUSTER::fin()
{
    
    Memory::dealloc(rho_phi);
    Memory::dealloc(drho_phi_dr);
    Memory::dealloc(drho_phi_dalpha);
    Memory::dealloc(B_pair);
    Memory::dealloc(is_pair);
    max_pairs=0;
    
    delete rho_ptr;
    delete dE_ptr;
    delete E_ptr;
    delete theta_ptr;
    delete b_ptr;
    delete d_ptr;
    delete mu_ptr;
    
    dE_ptr=rho_ptr=d_ptr=vec0=vec1=vec2=vec3=NULL;
    mu_ptr=NULL;
    post_fin();
}
/*--------------------------------------------
 set the temperature in the simulation
 --------------------------------------------*/
void ForceFieldEAMDMDCLUSTER::set_temp()
{
    type0 T=atoms->temp;
    type0 kb=atoms->kB;
    type0 hP=atoms->hP;
    type0 mass;
    type0 deb_l;
    
    for(size_t i=0;i<nelems;i++)
    {
        mass=atoms->elements.masses[i];
        c_1[i]=sqrt(0.5*kb*T/mass)/M_PI;
        deb_l=hP*hP/(2.0*M_PI*M_PI*mass*kb*T);
        c_0[i]=1.5*kb*T*(log(deb_l)-1.0);
    }
    
    kbT=kb*T;
    beta=1.0/kbT;
}
/*--------------------------------------------
 return M_{ij}^{\alpha}
 --------------------------------------------*/
inline type0 ForceFieldEAMDMDCLUSTER::calc_ent(type0 x)
{
    if(x==0.0) return 0.0;
    return x*log(x);
}

/*--------------------------------------------
 python constructor
 --------------------------------------------*/
void ForceFieldEAMDMDCLUSTER::ml_new(PyMethodDef& method_0,PyMethodDef& method_1,PyMethodDef& method_2)
{
    method_0.ml_flags=METH_VARARGS | METH_KEYWORDS;
    method_0.ml_name="ff_eam_funcfl_cluster";
    method_0.ml_meth=(PyCFunction)(PyCFunctionWithKeywords)
    [](PyObject* self,PyObject* args,PyObject* kwds)->PyObject*
    {
        AtomsDMD::Object* __self=reinterpret_cast<AtomsDMD::Object*>(self);
        FuncAPI<std::string*,type0*,type0*,std::string*> f("ff_eam_funcfl_cluster",{"funcfl_files","r_crd","C","elems"});
        f.noptionals=1;
        f.logics<1>()[0]=VLogics("gt",0.0);
        f.logics<2>()[0]=VLogics("gt",0.0);
        
        const std::string* names=__self->atoms->elements.names;
        const size_t nelems=__self->atoms->elements.nelems;
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
        __self->ff=new ForceFieldEAMDMDCLUSTER(__self->atoms,dr,drho,nr,nrho,std::move(r_phi),std::move(rho),std::move(F),std::move(r_c),f.mov<1>(),f.mov<2>());
        Py_RETURN_NONE;
    };
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
    method_1.ml_name="ff_eam_setfl_cluster";
    method_1.ml_meth=(PyCFunction)(PyCFunctionWithKeywords)
    [](PyObject* self,PyObject* args,PyObject* kwds)->PyObject*
    {
        AtomsDMD::Object* __self=reinterpret_cast<AtomsDMD::Object*>(self);
        FuncAPI<std::string,type0*,type0*,std::string*> f("ff_eam_setfl_cluster",{"setfl_file","r_crd","C","elems"});
        f.noptionals=1;
        f.logics<1>()[0]=VLogics("gt",0.0);
        f.logics<2>()[0]=VLogics("gt",0.0);
        
        const std::string* names=__self->atoms->elements.names;
        const size_t nelems=__self->atoms->elements.nelems;
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
        __self->ff=new ForceFieldEAMDMDCLUSTER(__self->atoms,dr,drho,nr,nrho,std::move(r_phi),std::move(rho),std::move(F),std::move(r_c),f.mov<1>(),f.mov<2>());
        Py_RETURN_NONE;
    };
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
    method_2.ml_name="ff_eam_fs_cluster";
    method_2.ml_meth=(PyCFunction)(PyCFunctionWithKeywords)
    [](PyObject* self,PyObject* args,PyObject* kwds)->PyObject*
    {
        AtomsDMD::Object* __self=reinterpret_cast<AtomsDMD::Object*>(self);
        FuncAPI<std::string,type0*,type0*,std::string*> f("ff_eam_fs_cluster",{"fs_file","r_crd","C","elems"});
        f.noptionals=1;
        f.logics<1>()[0]=VLogics("gt",0.0);
        f.logics<2>()[0]=VLogics("gt",0.0);
        
        const std::string* names=__self->atoms->elements.names;
        const size_t nelems=__self->atoms->elements.nelems;
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
        __self->ff=new ForceFieldEAMDMDCLUSTER(__self->atoms,dr,drho,nr,nrho,std::move(r_phi),std::move(rho),std::move(F),std::move(r_c),f.mov<1>(),f.mov<2>());
        Py_RETURN_NONE;
    };
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


