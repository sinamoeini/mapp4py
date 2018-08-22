#include "ff_eam_fit.h"
#include "neighbor_md.h"
#include "atoms_md.h"
#include "elements.h"
#include "memory.h"
#include "dynamic_md.h"
using namespace MAPP_NS;
/*--------------------------------------------
 This is for my personal and fitting of iron
 and hydrogen
 --------------------------------------------*/
/*--------------------------------------------
 constructor
 --------------------------------------------*/
ForceFieldEAMFit::ForceFieldEAMFit(AtomsMD* __atoms,
type0(***&& __rho_AR)[2],
size_t**&& __rho_sz,
type0(***&& __phi_AR)[2],
size_t**&& __phi_sz,
type0(*&& __F_A)[3]
):
ForceFieldMD(__atoms),
rho_sz(__rho_sz),
phi_sz(__phi_sz),
alloc(true)
{
    
    __rho_sz=NULL;
    __phi_sz=NULL;
    
    size_t nelems=2;
    size_t tot_sz=0;
    for(size_t ielem=0;ielem<nelems;ielem++)
        for(size_t jelem=0;jelem<nelems;jelem++)
        {

            tot_sz+=rho_sz[ielem][jelem];
            tot_sz+=phi_sz[ielem][jelem];
        }
    
    type0* __data=NULL;
    Memory::alloc(__data,2*tot_sz+6);
    
    Memory::alloc(rho_A,nelems,nelems);
    Memory::alloc(phi_A,nelems,nelems);
    Memory::alloc(F_A,nelems);
    Memory::alloc(rho_R,nelems,nelems);
    Memory::alloc(phi_R,nelems,nelems);
    
    size_t __i=0;
    for(size_t ielem=0;ielem<nelems;ielem++)
        for(size_t jelem=0;jelem<nelems;jelem++)
        {
            rho_A[ielem][jelem]=__data+__i;
            for(size_t i=0;i<rho_sz[ielem][jelem];i++)
                __data[__i++]=__rho_AR[ielem][jelem][i][0];
        }

    for(size_t ielem=0;ielem<nelems;ielem++)
        for(size_t jelem=0;jelem<nelems;jelem++)
        {
            phi_A[ielem][jelem]=__data+__i;
            for(size_t i=0;i<phi_sz[ielem][jelem];i++)
                __data[__i++]=__phi_AR[ielem][jelem][i][0];
        }
    
    for(size_t ielem=0;ielem<nelems;ielem++)
    {
        F_A[ielem]=__data+__i;
        for(size_t i=0;i<3;i++)
            __data[__i++]=__F_A[ielem][i];
    }
    
    for(size_t ielem=0;ielem<nelems;ielem++)
        for(size_t jelem=0;jelem<nelems;jelem++)
        {
            rho_R[ielem][jelem]=__data+__i;
            for(size_t i=0;i<rho_sz[ielem][jelem];i++)
                __data[__i++]=__rho_AR[ielem][jelem][i][1];
        }

    for(size_t ielem=0;ielem<nelems;ielem++)
        for(size_t jelem=0;jelem<nelems;jelem++)
        {
            phi_R[ielem][jelem]=__data+__i;
            for(size_t i=0;i<phi_sz[ielem][jelem];i++)
                __data[__i++]=__phi_AR[ielem][jelem][i][1];
        }
    
    
    Memory::dealloc(__rho_AR);
    __rho_AR=NULL;
    Memory::dealloc(__phi_AR);
    __phi_AR=NULL;
    Memory::dealloc(__F_A);
    __F_A=NULL;
    

    
    Memory::alloc(phi_cut,nelems,nelems);
    Memory::alloc(rho_cut,nelems,nelems);
    
    for(size_t ielem=0;ielem<nelems;ielem++)
        for(size_t jelem=0;jelem<nelems;jelem++)
            cut[ielem][jelem]=0.0;
    
    for(size_t ielem=0;ielem<nelems;ielem++)
        for(size_t jelem=0;jelem<nelems;jelem++)
        {
            if(rho_sz[ielem][jelem])
            {
                sort_AR_ij(rho_A[ielem][jelem],rho_R[ielem][jelem],rho_sz[ielem][jelem]);
                rho_cut[ielem][jelem]=rho_R[ielem][jelem][0];
            }
            
            
            if(phi_sz[ielem][jelem])
            {
                sort_AR_ij(phi_A[ielem][jelem],phi_R[ielem][jelem],phi_sz[ielem][jelem]);
                phi_cut[jelem][ielem]=phi_cut[ielem][jelem]=phi_R[ielem][jelem][0];
            }
        }
    phi_cut[1][1]=MAX(2.8,rho_cut[1][1]);
    
    for(size_t ielem=0;ielem<nelems;ielem++)
        for(size_t jelem=ielem;jelem<nelems;jelem++)
        {
            cut[ielem][jelem]=cut[jelem][ielem]=MAX(phi_cut[ielem][jelem],MAX(rho_cut[ielem][jelem],rho_cut[jelem][ielem]));
            cut_sq[ielem][jelem]=cut_sq[jelem][ielem]=cut[ielem][jelem]*cut[ielem][jelem];
        }

/*
    for(size_t ielem=0;ielem<nelems;ielem++)
        for(size_t jelem=0;jelem<nelems;jelem++)
        {
            for(size_t i=0;i<rho_sz[ielem][jelem];i++)
                printf("%s %s %lf %lf\n",atoms->elements.names[ielem].c_str(),atoms->elements.names[jelem].c_str(),rho_A[ielem][jelem][i],rho_R[ielem][jelem][i]);
            printf("\n");
        }
    
    for(size_t ielem=0;ielem<nelems;ielem++)
        for(size_t jelem=0;jelem<nelems;jelem++)
        {
            for(size_t i=0;i<phi_sz[ielem][jelem];i++)
                printf("%s %s %lf %lf\n",atoms->elements.names[ielem].c_str(),atoms->elements.names[jelem].c_str(),phi_A[ielem][jelem][i][0],phi_R[ielem][jelem][i]);
            printf("\n");
        }
  */
    

}
/*--------------------------------------------
 constructor
 --------------------------------------------*/
ForceFieldEAMFit::ForceFieldEAMFit(
AtomsMD* __atoms,
type0***& __rho_A,
type0***& __rho_R,
size_t**& __rho_sz,
type0***& __phi_A,
type0***& __phi_R,
size_t**& __phi_sz,
type0**& __F_A):
ForceFieldMD(__atoms),
rho_A(__rho_A),
rho_R(__rho_R),
rho_sz(__rho_sz),
phi_A(__phi_A),
phi_R(__phi_R),
phi_sz(__phi_sz),
F_A(__F_A),
alloc(false)
{

    Memory::alloc(phi_cut,nelems,nelems);
    Memory::alloc(rho_cut,nelems,nelems);
    
    for(size_t ielem=0;ielem<nelems;ielem++)
        for(size_t jelem=0;jelem<nelems;jelem++)
            cut[ielem][jelem]=0.0;
    
    for(size_t ielem=0;ielem<nelems;ielem++)
        for(size_t jelem=0;jelem<nelems;jelem++)
        {
            if(rho_sz[ielem][jelem])
            {
                sort_AR_ij(rho_A[ielem][jelem],rho_R[ielem][jelem],rho_sz[ielem][jelem]);
                rho_cut[ielem][jelem]=rho_R[ielem][jelem][0];
            }
            
            
            if(phi_sz[ielem][jelem])
            {
                sort_AR_ij(phi_A[ielem][jelem],phi_R[ielem][jelem],phi_sz[ielem][jelem]);
                phi_cut[jelem][ielem]=phi_cut[ielem][jelem]=phi_R[ielem][jelem][0];
            }
        }
    phi_cut[1][1]=MAX(2.8,rho_cut[1][1]);
    
    for(size_t ielem=0;ielem<nelems;ielem++)
        for(size_t jelem=ielem;jelem<nelems;jelem++)
        {
            cut[ielem][jelem]=cut[jelem][ielem]=MAX(phi_cut[ielem][jelem],MAX(rho_cut[ielem][jelem],rho_cut[jelem][ielem]));
            cut_sq[ielem][jelem]=cut_sq[jelem][ielem]=cut[ielem][jelem]*cut[ielem][jelem];
        }
}
/*--------------------------------------------
 destructor
 --------------------------------------------*/
ForceFieldEAMFit::~ForceFieldEAMFit()
{
    if(alloc)
    {
        **phi_R=NULL;
        **rho_R=NULL;
        *F_A=NULL;
        **phi_A=NULL;
        
        Memory::dealloc(phi_R);
        Memory::dealloc(rho_R);
        Memory::dealloc(F_A);
        Memory::dealloc(phi_A);
        Memory::dealloc(rho_A);
        
    }
    
    Memory::dealloc(rho_cut);
    Memory::dealloc(phi_cut);
}
/*--------------------------------------------
 destructor
 --------------------------------------------*/
void ForceFieldEAMFit::sort_AR_ij(type0*& A_ij,type0*& R_ij,size_t sz)
{
    
    size_t* key=NULL;
    Memory::alloc(key,sz);
    for(size_t i=0;i<sz;i++) key[i]=i;
    
    XMath::quicksort(key,key+sz,
    [&R_ij](size_t* rank_i,size_t* rank_j)
    {return (R_ij[*rank_i]>R_ij[*rank_j]);},
    [&R_ij,&A_ij](size_t* rank_i,size_t* rank_j)
    {
        std::swap(R_ij[*rank_i],R_ij[*rank_j]);
        std::swap(A_ij[*rank_i],A_ij[*rank_j]);
    });
    
    Memory::dealloc(key);
    
}
/*--------------------------------------------
 force and energy calculation
 --------------------------------------------*/
void ForceFieldEAMFit::__force_calc()
{
    type0* xvec=atoms->x->begin();
    type0* fvec=f->begin();
    type0* rho=rho_ptr->begin();
    elem_type* evec=atoms->elem->begin();
    
    int iatm,jatm;
    elem_type ielem,jelem;
    type0 r,rsq;
    type0 fpair;
    type0 dx_ij[__dim__];
    int** neighbor_list=neighbor->neighbor_list;
    int* neighbor_list_size=neighbor->neighbor_list_size;
    const int natms_lcl=atoms->natms_lcl;
    for(int i=0;i<natms_lcl;i++)
        rho[i]=0.0;
    
    for(iatm=0;iatm<natms_lcl;iatm++)
    {
        //printf("%e %e %e\n",xvec[3*iatm],xvec[3*iatm+1],xvec[3*iatm+2]);
        ielem=evec[iatm];
        for(int j=0;j<neighbor_list_size[iatm];j++)
        {
            jatm=neighbor_list[iatm][j];
            jelem=evec[jatm];
            rsq=Algebra::RSQ<__dim__>(xvec+iatm*__dim__,xvec+jatm*__dim__);
            if(rsq>=cut_sq[ielem][jelem]) continue;
            r=sqrt(rsq);
            

            rho[iatm]+=calc_rho(jelem,ielem,r);
            if(jatm<natms_lcl)
            {
                rho[jatm]+=calc_rho(ielem,jelem,r);
                __vec_lcl[0]+=calc_phi(ielem,jelem,r);;
            }
            else
                __vec_lcl[0]+=0.5*calc_phi(ielem,jelem,r);;
        }
        // add the embedded energy here
        __vec_lcl[0]+=calc_F(ielem,rho[iatm]);
    }
    
    /*
    for(int i=0;i<natms_lcl;i++)
    {
        printf("%.16lf %.16lf %.16lf %.16lf\n",xvec[3*i],xvec[3*i+1],xvec[3*i+2],rho[i]);
    }*/

    update(rho_ptr);


    //printf("------------------------------------------\n");
    for(iatm=0;iatm<natms_lcl;iatm++)
    {
        ielem=evec[iatm];

        for(int j=0;j<neighbor_list_size[iatm];j++)
        {
            jatm=neighbor_list[iatm][j];
            jelem=evec[jatm];

            rsq=Algebra::DX_RSQ<__dim__>(xvec+iatm*__dim__,xvec+jatm*__dim__,dx_ij);
            if(rsq>=cut_sq[ielem][jelem]) continue;
            
            r=sqrt(rsq);
            fpair=-(calc_dF(ielem,rho[iatm])*calc_drho(jelem,ielem,r)+
                    calc_dF(jelem,rho[jatm])*calc_drho(ielem,jelem,r)+
                    calc_dphi(ielem,jelem,r))/r;
            if(fpair==0.0) continue;
            
            Algebra::V_add_x_mul_V<__dim__>(fpair,dx_ij,fvec+iatm*__dim__);
            if(jatm<natms_lcl)
                Algebra::V_add_x_mul_V<__dim__>(-fpair,dx_ij,fvec+jatm*__dim__);
            else
                fpair*=0.5;
            
            Algebra::DyadicV<__dim__>(-fpair,dx_ij,&__vec_lcl[1]);
        }
    }
    type0 f_sum_lcl[__dim__];
    Algebra::zero<__dim__>(f_sum_lcl);
    for(int i=0;i<natms_lcl;i++)
        Algebra::V_add<__dim__>(fvec+__dim__*i,f_sum_lcl);

    type0 f_corr[__dim__];
    MPI_Allreduce(f_sum_lcl,f_corr,__dim__,Vec<type0>::MPI_T,MPI_SUM,world);
    type0 a=-1.0/static_cast<type0>(atoms->natms);
    Algebra::Do<__dim__>::func([&a,&f_corr](int i){f_corr[i]*=a;});
    for(int i=0;i<natms_lcl;i++)
        Algebra::V_add<__dim__>(f_corr,fvec+__dim__*i);
    
}
/*--------------------------------------------
 only energy calculation this is useful for
 minimization/linesearch methods that do not
 use derivatives of energy
 --------------------------------------------*/
void ForceFieldEAMFit::__energy_calc()
{
    type0* xvec=atoms->x->begin();
    type0* rho=rho_ptr->begin();
    elem_type* evec=atoms->elem->begin();
    
    int iatm,jatm;
    elem_type ielem,jelem;
    type0 r,rsq;
    
    int** neighbor_list=neighbor->neighbor_list;
    int* neighbor_list_size=neighbor->neighbor_list_size;
    const int natms_lcl=atoms->natms_lcl;
    for(int i=0;i<natms_lcl;i++)
        rho[i]=0.0;
    
    for(iatm=0;iatm<natms_lcl;iatm++)
    {
        ielem=evec[iatm];
        for(int j=0;j<neighbor_list_size[iatm];j++)
        {
            jatm=neighbor_list[iatm][j];
            jelem=evec[jatm];
            rsq=Algebra::RSQ<__dim__>(xvec+iatm*__dim__,xvec+jatm*__dim__);
            if(rsq>=cut_sq[ielem][jelem]) continue;
            r=sqrt(rsq);
            
            
            rho[iatm]+=calc_rho(jelem,ielem,r);
            if(jatm<natms_lcl)
            {
                rho[jatm]+=calc_rho(ielem,jelem,r);
                __vec_lcl[0]+=calc_phi(ielem,jelem,r);
            }
            else
                __vec_lcl[0]+=0.5*calc_phi(ielem,jelem,r);
        }
        // add the embedded energy here
        __vec_lcl[0]+=calc_F(ielem,rho[iatm]);
    }
    /*
    for(int i=0;i<natms_lcl;i++)
    {
        printf("EEE %.16lf %.16lf %.16lf %.16lf\n",xvec[3*i],xvec[3*i+1],xvec[3*i+2],rho[i]);
    }*/
}
/*--------------------------------------------
 only energy calculation this is useful for
 minimization/linesearch methods that do not
 use derivatives of energy
 --------------------------------------------*/
void ForceFieldEAMFit::prep4deriv()
{
    type0* xvec=atoms->x->begin();
    type0* fvec=f->begin();
    type0* rho=rho_ptr->begin();
    elem_type* evec=atoms->elem->begin();
    
    int iatm,jatm;
    elem_type ielem,jelem;
    type0 r,rsq;
    type0 fpair;
    type0 dx_ij[__dim__];
    int** neighbor_list=neighbor->neighbor_list;
    int* neighbor_list_size=neighbor->neighbor_list_size;
    const int natms_lcl=atoms->natms_lcl;
    for(int i=0;i<natms_lcl;i++)
        rho[i]=0.0;
    
    for(iatm=0;iatm<natms_lcl;iatm++)
    {
        ielem=evec[iatm];
        for(int j=0;j<neighbor_list_size[iatm];j++)
        {
            jatm=neighbor_list[iatm][j];
            jelem=evec[jatm];
            rsq=Algebra::RSQ<__dim__>(xvec+iatm*__dim__,xvec+jatm*__dim__);
            if(rsq>=cut_sq[ielem][jelem]) continue;
            r=sqrt(rsq);
            

            rho[iatm]+=calc_rho(jelem,ielem,r);
            
            if(jatm<natms_lcl)
                rho[jatm]+=calc_rho(ielem,jelem,r);
            
        }
    }
    
    update(rho_ptr);

    
    type0* S=S_ptr->begin();
    for(int i=0;i<natms_lcl*__nvoigt__;i++)
        S[i]=0.0;
    type0 drho_ji,drho_ij;
    for(iatm=0;iatm<natms_lcl;iatm++)
    {
        ielem=evec[iatm];
        
        for(int j=0;j<neighbor_list_size[iatm];j++)
        {
            jatm=neighbor_list[iatm][j];
            jelem=evec[jatm];
            
            rsq=Algebra::DX_RSQ<__dim__>(xvec+iatm*__dim__,xvec+jatm*__dim__,dx_ij);
            if(rsq>=cut_sq[ielem][jelem]) continue;
            
            r=sqrt(rsq);
            drho_ji=calc_drho(jelem,ielem,r);
            drho_ij=calc_drho(ielem,jelem,r);
            fpair=-(calc_dF(ielem,rho[iatm])*drho_ji+calc_dF(jelem,rho[jatm])*drho_ij+calc_dphi(ielem,jelem,r))/r;
            if(fpair==0.0) continue;
            
            Algebra::V_add_x_mul_V<__dim__>(fpair,dx_ij,fvec+iatm*__dim__);
            Algebra::DyadicV<__dim__>(drho_ji/r,dx_ij,S+__nvoigt__*iatm);
            if(jatm<natms_lcl)
            {
                Algebra::V_add_x_mul_V<__dim__>(-fpair,dx_ij,fvec+jatm*__dim__);
                Algebra::DyadicV<__dim__>(drho_ij/r,dx_ij,S+__nvoigt__*jatm);
            }
            else
                fpair*=0.5;
            

        }
    }

    

    
}
/*--------------------------------------------
 initiate before a run
 --------------------------------------------*/
void ForceFieldEAMFit::init()
{
    pre_init();
    rho_ptr=new Vec<type0>(atoms,1,"rho");
    S_ptr=new Vec<type0>(atoms,1+__nvoigt__);
}
/*--------------------------------------------
 after a run
 --------------------------------------------*/
void ForceFieldEAMFit::fin()
{
    delete S_ptr;
    delete rho_ptr;
    post_fin();
}
/*--------------------------------------------
 init xchng
 --------------------------------------------*/
void ForceFieldEAMFit::init_xchng()
{
}
/*--------------------------------------------
 fin xchng
 --------------------------------------------*/
void ForceFieldEAMFit::fin_xchng()
{
}
/*--------------------------------------------
 pre xchng energy
 --------------------------------------------*/
void ForceFieldEAMFit::pre_xchng_energy(GCMC*)
{
}
/*--------------------------------------------
 xchng energy
 --------------------------------------------*/
type0 ForceFieldEAMFit::xchng_energy(GCMC*)
{
    return 0.0;
}
/*--------------------------------------------
 post xchng energy
 --------------------------------------------*/
void ForceFieldEAMFit::post_xchng_energy(GCMC*)
{
}
/*--------------------------------------------
 post xchng energy
 --------------------------------------------*/
void ForceFieldEAMFit::drho_A(elem_type i,elem_type j,type0 A,type0 R,type0(&tmp_vec)[1+__nvoigt__])
{
    Vec<type0>* Drho_ptr=new Vec<type0>(atoms,1);
    
    type0* xvec=atoms->x->begin();
    elem_type* evec=atoms->elem->begin();
    type0* Drho=Drho_ptr->begin();
    type0* rho=rho_ptr->begin();
    type0* S=S_ptr->begin();
    
    int iatm,jatm;
    elem_type ielem,jelem;
    type0 r;
    type0 Drho_ij,Ddrho_ij;
    type0 dx_ij[__dim__];
    
    int** neighbor_list=neighbor->neighbor_list;
    int* neighbor_list_size=neighbor->neighbor_list_size;
    const int natms_lcl=atoms->natms_lcl;
    for(int i=0;i<natms_lcl;i++)
        Drho[i]=0.0;
    
    type0 tmp_vec_lcl[1+__nvoigt__];
    Algebra::zero<1+__nvoigt__>(tmp_vec_lcl);

    
    for(iatm=0;iatm<natms_lcl;iatm++)
    {
        ielem=evec[iatm];
        for(int __j=0;__j<neighbor_list_size[iatm];__j++)
        {
            jatm=neighbor_list[iatm][__j];
            jelem=evec[jatm];
            
            if(!((ielem==i && jelem==j) || (ielem==j && jelem==i))) continue;
            r=sqrt(Algebra::DX_RSQ<__dim__>(xvec+iatm*__dim__,xvec+jatm*__dim__,dx_ij));
            if(r>=R) continue;
                
            Drho_ij=(R-r)*(R-r)*(R-r);
            Ddrho_ij=-3.0*(R-r)*(R-r);
            
            if(ielem==j)
            {
                Drho[iatm]+=Drho_ij;
                Algebra::DyadicV<__dim__>(Ddrho_ij*calc_dF(j,rho[iatm])/r,dx_ij,&tmp_vec_lcl[1]);
            }

            if(jatm<natms_lcl && jelem==j)
            {
                Drho[jatm]+=Drho_ij;
                Algebra::DyadicV<__dim__>(Ddrho_ij*calc_dF(j,rho[jatm])/r,dx_ij,&tmp_vec_lcl[1]);
            }
            
        }
        if(j==ielem)
        {
            Algebra::V_add_x_mul_V<__nvoigt__>(Drho[iatm]*calc_ddF(j,rho[iatm]),S+iatm*__nvoigt__,&tmp_vec_lcl[1]);
            tmp_vec_lcl[0]+=calc_dF(j,rho[iatm])*Drho[iatm];
        }
    }

    
    /*
    if(i==1 && j==1)
    {
        type0 rsq,rho_ij,drho_ij,dFH,ddFH,Dphi_ij,Ddphi_ij;
        
        //it should be completed 
        for(iatm=0;iatm<natms_lcl;iatm++)
        {
            ielem=evec[iatm];
            if(ielem!=1) continue;
            for(int __j=0;__j<neighbor_list_size[iatm];__j++)
            {
                jatm=neighbor_list[iatm][__j];
                jelem=evec[jatm];
                
                if(jelem!=1) continue;
                
                rsq=Algebra::DX_RSQ<__dim__>(xvec+iatm*__dim__,xvec+jatm*__dim__,dx_ij);
                r=sqrt(rsq);
                if(r>=R) continue;
                
                rho_ij=calc_rho(1,1,r);
                drho_ij=calc_drho(1,1,r);
                dFH=calc_dF(i,rho_ij);
                ddFH=calc_ddF(i,rho_ij);
                Drho_ij=(R-r)*(R-r)*(R-r);
                Ddrho_ij=-3.0*(R-r)*(R-r);
                
                Dphi_ij=-2.0*dFH*Drho_ij;
                Ddphi_ij=-2.0*ddFH*Drho_ij*drho_ij-2.0*dFH*Ddrho_ij;
                
                
                if(jatm<natms_lcl)
                {
                    Algebra::DyadicV<__dim__>(Ddphi_ij/r,dx_ij,&tmp_vec_lcl[1]);
                    tmp_vec_lcl[0]+=Dphi_ij;
                }
                else
                {
                    Algebra::DyadicV<__dim__>(0.5*Ddphi_ij/r,dx_ij,&tmp_vec_lcl[1]);
                    tmp_vec_lcl[0]+=0.5*Dphi_ij;
                }                
            }
           
        }
    }*/
    
    MPI_Allreduce(tmp_vec_lcl,tmp_vec,__nvoigt__+1,Vec<type0>::MPI_T,MPI_SUM,world);
    const type0 vol=atoms->vol;
    Algebra::Do<__nvoigt__>::func([&tmp_vec,&vol](int i){tmp_vec[i+1]/=vol;});
}
/*--------------------------------------------
 post xchng energy
 --------------------------------------------*/
void ForceFieldEAMFit::dphi_A(elem_type i,elem_type j,type0 A,type0 R,type0(&tmp_vec)[1+__nvoigt__])
{
    type0* xvec=atoms->x->begin();
    elem_type* evec=atoms->elem->begin();
    
    
    int iatm,jatm;
    elem_type ielem,jelem;
    type0 r;
    type0 Dphi_ij,Ddphi_ij;
    type0 dx_ij[__dim__];
    
    int** neighbor_list=neighbor->neighbor_list;
    int* neighbor_list_size=neighbor->neighbor_list_size;
    const int natms_lcl=atoms->natms_lcl;
    type0 tmp_vec_lcl[1+__nvoigt__];
    Algebra::zero<1+__nvoigt__>(tmp_vec_lcl);
    
    for(iatm=0;iatm<natms_lcl;iatm++)
    {
        ielem=evec[iatm];
        for(int __j=0;__j<neighbor_list_size[iatm];__j++)
        {
            jatm=neighbor_list[iatm][__j];
            jelem=evec[jatm];
            
            if(!((ielem==i && jelem==j) || (ielem==j && jelem==i))) continue;
            r=sqrt(Algebra::DX_RSQ<__dim__>(xvec+iatm*__dim__,xvec+jatm*__dim__,dx_ij));
            if(r>=R) continue;
            
            Dphi_ij=(R-r)*(R-r)*(R-r);
            Ddphi_ij=-3.0*(R-r)*(R-r);
        
            
            if(jatm<natms_lcl)
            {
                tmp_vec_lcl[0]+=Dphi_ij;
                Algebra::DyadicV<__dim__>(Ddphi_ij/r,dx_ij,&tmp_vec_lcl[1]);
            }
            else
            {
                tmp_vec_lcl[0]+=0.5*Dphi_ij;
                Algebra::DyadicV<__dim__>(0.5*Ddphi_ij/r,dx_ij,&tmp_vec_lcl[1]);
            }

        }
    }
    
    
    MPI_Allreduce(tmp_vec_lcl,tmp_vec,__nvoigt__+1,Vec<type0>::MPI_T,MPI_SUM,world);
    const type0 vol=atoms->vol;
    Algebra::Do<__nvoigt__>::func([&tmp_vec,&vol](int i){tmp_vec[i+1]/=vol;});

}
/*--------------------------------------------
 post xchng energy
 --------------------------------------------*/
void ForceFieldEAMFit::dFH_A(type0(&tmp_A_vec)[1+__nvoigt__],type0(&tmp_rho0_vec)[1+__nvoigt__],type0(&tmp_alpha_vec)[1+__nvoigt__])
{
    
    type0* S=S_ptr->begin();
    
    //type0* xvec=atoms->x->begin();
    elem_type* evec=atoms->elem->begin();
    type0* rho=rho_ptr->begin();
    
    //int iatm,jatm;
    //elem_type ielem,jelem;
    int iatm;
    elem_type ielem;
    type0 DFH_DA,DFH_Drho0,DFH_Dalpha;
    type0 DdFH_DA,DdFH_Drho0,DdFH_Dalpha;
    //type0 rho_ij,drho_ij,r,rsq;
    //type0 dx_ij[__dim__];
    
    //int** neighbor_list=neighbor->neighbor_list;
    //int* neighbor_list_size=neighbor->neighbor_list_size;
    const int natms_lcl=atoms->natms_lcl;
    
    type0 tmp_A_vec_lcl[1+__nvoigt__];
    Algebra::zero<1+__nvoigt__>(tmp_A_vec_lcl);
    type0 tmp_rho0_vec_lcl[1+__nvoigt__];
    Algebra::zero<1+__nvoigt__>(tmp_rho0_vec_lcl);
    type0 tmp_alpha_vec_lcl[1+__nvoigt__];
    Algebra::zero<1+__nvoigt__>(tmp_alpha_vec_lcl);
    
    
    for(iatm=0;iatm<natms_lcl;iatm++)
    {
        ielem=evec[iatm];
        if(ielem!=1) continue;
        
        calc_DFH(rho[iatm],DFH_DA,DFH_Drho0,DFH_Dalpha);
        calc_DdFH(rho[iatm],DdFH_DA,DdFH_Drho0,DdFH_Dalpha);
        tmp_A_vec_lcl[0]+=DFH_DA;
        Algebra::V_add_x_mul_V<__nvoigt__>(DdFH_DA,S+iatm*__nvoigt__,&tmp_A_vec_lcl[1]);
        tmp_rho0_vec_lcl[0]+=DFH_Drho0;
        Algebra::V_add_x_mul_V<__nvoigt__>(DdFH_Drho0,S+iatm*__nvoigt__,&tmp_rho0_vec_lcl[1]);
        tmp_alpha_vec_lcl[0]+=DFH_Dalpha;
        Algebra::V_add_x_mul_V<__nvoigt__>(DdFH_Dalpha,S+iatm*__nvoigt__,&tmp_alpha_vec_lcl[1]);
        
        
        /*
        for(int __j=0;__j<neighbor_list_size[iatm];__j++)
        {
            jatm=neighbor_list[iatm][__j];
            jelem=evec[jatm];
            if(jelem!=1) continue;
            rsq=Algebra::DX_RSQ<__dim__>(xvec+iatm*__dim__,xvec+jatm*__dim__,dx_ij);
            r=sqrt(rsq);
            rho_ij=drho_ij=0.0;
            rho_ij=calc_rho(1,1,r);
            drho_ij=calc_drho(1,1,r);
            calc_DFH(rho_ij,DFH_DA,DFH_Drho0,DFH_Dalpha);
            calc_DdFH(rho_ij,DdFH_DA,DdFH_Drho0,DdFH_Dalpha);
            
            if(jatm<natms_lcl)
            {
                DFH_DA*=2.0;
                DFH_Drho0*=2.0;
                DFH_Dalpha*=2.0;
                DdFH_DA*=2.0;
                DdFH_Drho0*=2.0;
                DdFH_Dalpha*=2.0;
            }

            tmp_A_vec_lcl[0]-=DFH_DA;
            Algebra::DyadicV<__dim__>(-DdFH_DA*drho_ij/r,dx_ij,&tmp_A_vec_lcl[1]);
            tmp_rho0_vec_lcl[0]-=DFH_Drho0;
            Algebra::DyadicV<__dim__>(-DdFH_Drho0*drho_ij/r,dx_ij,&tmp_rho0_vec_lcl[1]);
            tmp_alpha_vec_lcl[0]-=DFH_Dalpha;
            Algebra::DyadicV<__dim__>(-DdFH_Dalpha*drho_ij/r,dx_ij,&tmp_alpha_vec_lcl[1]);
        }
         */
    }
    
    MPI_Allreduce(tmp_A_vec_lcl,tmp_A_vec,__nvoigt__+1,Vec<type0>::MPI_T,MPI_SUM,world);
    MPI_Allreduce(tmp_rho0_vec_lcl,tmp_rho0_vec,__nvoigt__+1,Vec<type0>::MPI_T,MPI_SUM,world);
    MPI_Allreduce(tmp_alpha_vec_lcl,tmp_alpha_vec,__nvoigt__+1,Vec<type0>::MPI_T,MPI_SUM,world);
    

    const type0 vol=atoms->vol;
    Algebra::Do<__nvoigt__>::func([&tmp_A_vec,&tmp_rho0_vec,&tmp_alpha_vec,&vol](int i)
    {tmp_A_vec[i+1]/=vol; tmp_rho0_vec[i+1]/=vol; tmp_alpha_vec[i+1]/=vol;});
    
}
/*--------------------------------------------
 
 --------------------------------------------*/
void ForceFieldEAMFit::calc_deriv(
bool***& Arho_dof,type0 (***& DArho)[1+__nvoigt__],
bool***& Aphi_dof,type0 (***& DAphi)[1+__nvoigt__],
bool**& AF_dof,type0 (**& DAF)[1+__nvoigt__])
{
    
    
    prep4deriv();
    elem_type nelems=2;
    for(elem_type ielem=0;ielem<nelems;ielem++)
        for(elem_type jelem=0;jelem<nelems;jelem++)
            for(size_t i=0;i<rho_sz[ielem][jelem];i++)
            {
                if(Arho_dof[ielem][jelem][i])
                    drho_A(ielem,jelem,rho_A[ielem][jelem][i],rho_R[ielem][jelem][i],DArho[ielem][jelem][i]);
                else
                    Algebra::zero<1+__nvoigt__>(DArho[ielem][jelem][i]);
            }
    
    for(elem_type ielem=0;ielem<nelems;ielem++)
        for(elem_type jelem=0;jelem<ielem+1;jelem++)
            for(size_t i=0;i<phi_sz[ielem][jelem];i++)
            {
                if(Aphi_dof[ielem][jelem][i])
                    dphi_A(ielem,jelem,phi_A[ielem][jelem][i],phi_R[ielem][jelem][i],DAphi[ielem][jelem][i]);
                else
                    Algebra::zero<1+__nvoigt__>(DAphi[ielem][jelem][i]);
                if(ielem!=jelem)
                    Algebra::V_eq<1+__nvoigt__>(DAphi[ielem][jelem][i],DAphi[jelem][ielem][i]);
            }

    for(size_t i=0;i<3;i++)
    {
        Algebra::zero<1+__nvoigt__>(DAF[0][i]);
        Algebra::zero<1+__nvoigt__>(DAF[1][i]);
    }
    if(AF_dof[1][0] || AF_dof[1][1] || AF_dof[1][2])
    {
        dFH_A(DAF[1][0],DAF[1][1],DAF[1][2]);
        
        for(size_t i=0;i<3;i++)
            if(!AF_dof[1][i])
                Algebra::zero<1+__nvoigt__>(DAF[1][i]);
    
    }
}
/*--------------------------------------------
 
 --------------------------------------------*/
size_t ForceFieldEAMFit::get_rFeH(type0*& Rs,int*& Ns)
{
    
    
    type0* xvec=atoms->x->begin();
    type0* rho=rho_ptr->begin();
    elem_type* evec=atoms->elem->begin();
    
    int iatm,jatm;
    elem_type ielem,jelem;
    type0 r,rsq;
    
    int** neighbor_list=neighbor->neighbor_list;
    int* neighbor_list_size=neighbor->neighbor_list_size;
    const int natms_lcl=atoms->natms_lcl;
    for(int i=0;i<natms_lcl;i++)
    rho[i]=0.0;
    
    type0* rs_lcl=NULL;
    Memory::alloc(rs_lcl,natms_lcl);
    int n_lcl=0;
    
    for(iatm=0;iatm<natms_lcl;iatm++)
    {
        ielem=evec[iatm];
        for(int j=0;j<neighbor_list_size[iatm];j++)
        {
            jatm=neighbor_list[iatm][j];
            jelem=evec[jatm];
            if(!((ielem==1&&jelem==0) || ((ielem==0&&jelem==1)&&jatm<natms_lcl))) continue;
            
            
            
            
            
            rsq=Algebra::RSQ<__dim__>(xvec+iatm*__dim__,xvec+jatm*__dim__);
            if(rsq>=cut_sq[ielem][jelem]) continue;
            r=sqrt(rsq);
            
            rs_lcl[n_lcl++]=r;

        }

    }
    
    int n;
    MPI_Allreduce(&n_lcl,&n,1,Vec<int>::MPI_T,MPI_SUM,world);
    type0* rs=NULL;
    Memory::alloc(rs,n);
    int __n_lcl;
    type0* __rs=rs;
    for(int i=0;i<atoms->comm_size;i++)
    {
        __n_lcl=n_lcl;
        MPI_Bcast(&__n_lcl,1,Vec<int>::MPI_T,i,world);
        if(atoms->comm_rank==i)
            memcpy(__rs,rs_lcl,__n_lcl*sizeof(type0));
        
        MPI_Bcast(__rs,__n_lcl,Vec<type0>::MPI_T,i,world);
        __rs+=__n_lcl;
    }
    Memory::dealloc(rs_lcl);
    
    XMath::quicksort(rs,rs+n,
    [](type0* r_i,type0* r_j)
    {return (*r_i<*r_j);},
    [](type0* r_i,type0* r_j)
    {
        std::swap(*r_i,*r_j);
    });
    
    type0 tol=1.0e-5;
    int i0=0,i1=0;
    type0 r0,n0;
    size_t sz=0;
    
    type0* __Rs=NULL;
    Memory::alloc(__Rs,n);
    int* __Ns=NULL;
    Memory::alloc(__Ns,n);
    
    while(i0<n)
    {
        r0=0.0;
        n0=0.0;
        while (i1<n && fabs(rs[i1]-rs[i0])<tol)
        {
            r0+=rs[i1];
            n0++;
            i1++;
        }
        r0/=n0;
        
        __Rs[sz]=r0;
        __Ns[sz]=static_cast<int>(n0);
        sz++;
        i0=i1;
    }
    Memory::dealloc(rs);
    
    Rs=NULL;
    Memory::alloc(Rs,sz);
    Ns=NULL;
    Memory::alloc(Ns,sz);
    memcpy(Rs,__Rs,sz*sizeof(type0));
    memcpy(Ns,__Ns,sz*sizeof(int));
    
    Memory::dealloc(__Rs);
    Memory::dealloc(__Ns);
    
    
    return sz;
    
    
    
}







/*--------------------------------------------
 
 --------------------------------------------*/
type0 ForceFieldEAMFit::calc_F(elem_type i,type0 rho)
{
    if(i==0)
        return F_A[0][0]*sqrt(rho)+F_A[0][1]*rho*rho+F_A[0][2]*rho*rho*rho*rho;
    
    type0 rhob=rho/F_A[1][1];
    type0 eta=pow(rhob,F_A[1][2]);
    type0 falpha=pow((1.0-3.0/(2.0*F_A[1][2]-1.0))/4.0,1.0-2.0/F_A[1][2]);
    return F_A[1][0]*(eta*(eta-2.0)+falpha*rhob*rhob);
    
}
/*--------------------------------------------
 
 --------------------------------------------*/
type0 ForceFieldEAMFit::calc_dF(elem_type i,type0 rho)
{
    if(i==0)
        return 0.5*F_A[0][0]/sqrt(rho)+2.0*F_A[0][1]*rho+4.0*F_A[0][2]*rho*rho*rho;
    if(rho==0.0) return 0.0;
    type0 rhob=rho/F_A[1][1];
    type0 eta=pow(rhob,F_A[1][2]);
    type0 falpha=pow((1.0-3.0/(2.0*F_A[1][2]-1.0))/4.0,1.0-2.0/F_A[1][2]);
    return 2.0*F_A[1][0]*(F_A[1][2]*eta*(eta-1.0)+falpha*rhob*rhob)/rho;
}
/*--------------------------------------------
 
 --------------------------------------------*/
type0 ForceFieldEAMFit::calc_ddF(elem_type i,type0 rho)
{
    if(i==0)
        return -0.25*F_A[0][0]/sqrt(rho*rho*rho)+2.0*F_A[0][1]+12.0*F_A[0][2]*rho*rho;
    if(rho==0.0) return 0.0;
    type0 rhob=rho/F_A[1][1];
    type0 eta=pow(rhob,F_A[1][2]);
    type0 falpha=pow((1.0-3.0/(2.0*F_A[1][2]-1.0))/4.0,1.0-2.0/F_A[1][2]);
    return 2.0*F_A[1][0]*(F_A[1][2]*eta*((2.0*F_A[1][2]-1.0)*eta+1.0-F_A[1][2])+falpha*rhob*rhob)/(rho*rho);
    
}
/*--------------------------------------------
 
 --------------------------------------------*/
void ForceFieldEAMFit::calc_DFH(type0& rho,type0& DF_DA,type0& DF_Drho0,type0& DF_Dalpha)
{
    if(rho==0.0)
    {
        DF_DA=0.0;
        DF_Drho0=0.0;
        DF_Dalpha=0.0;
        return;
    }
    type0 rhob=rho/F_A[1][1];
    type0 eta=pow(rhob,F_A[1][2]);
    type0 falpha=pow((1.0-3.0/(2.0*F_A[1][2]-1.0))/4.0,1.0-2.0/F_A[1][2]);
    type0 dlogfalpha=(2.0*log((1.0-3.0/(2.0*F_A[1][2]-1.0))/4.0)/F_A[1][2]+3.0/(2.0*F_A[1][2]-1.0))/F_A[1][2];
    DF_DA=eta*(eta-2.0)+falpha*rhob*rhob;
    DF_Drho0=-2.0*F_A[1][0]*(F_A[1][2]*eta*(eta-1.0)+falpha*rhob*rhob)/F_A[1][1];
    DF_Dalpha=F_A[1][0]*(2.0*eta*(eta-1.0)*log(rhob)+dlogfalpha*falpha*rhob*rhob);
}
/*--------------------------------------------
 
 --------------------------------------------*/
void ForceFieldEAMFit::calc_DdFH(type0& rho,type0& DdF_DA,type0& DdF_Drho0,type0& DdF_Dalpha)
{
    if(rho==0.0)
    {
        DdF_DA=0.0;
        DdF_Drho0=0.0;
        DdF_Dalpha=0.0;
        return;
    }
    type0 rhob=rho/F_A[1][1];
    type0 eta=pow(rhob,F_A[1][2]);
    type0 falpha=pow((1.0-3.0/(2.0*F_A[1][2]-1.0))/4.0,1.0-2.0/F_A[1][2]);
    type0 dlogfalpha=(2.0*log((1.0-3.0/(2.0*F_A[1][2]-1.0))/4.0)/F_A[1][2]+3.0/(2.0*F_A[1][2]-1.0))/F_A[1][2];
    DdF_DA=2.0*(F_A[1][2]*eta*(eta-1.0)+falpha*rhob*rhob)/rho;
    DdF_Drho0=-2.0*F_A[1][0]*(F_A[1][2]*F_A[1][2]*eta*(2.0*eta-1.0)
                              +2.0*falpha*rhob*rhob)/(F_A[1][1]*rho);
    DdF_Dalpha=2.0*F_A[1][0]*(eta*((1.0+F_A[1][2]*log(rhob))*(2.0*eta-1.0)-eta)+dlogfalpha*falpha*rhob*rhob)/rho;
}
/*--------------------------------------------
 
 --------------------------------------------*/
type0 ForceFieldEAMFit::EHH(type0 r)
{
    if(r>=2.8) return 0.0;
    
    return 1.2105506230313432 + (33.34234844583974 - 97.08416920979838*r)*exp(-2.5871227340727088*r)- 0.3714821016657931*r;
}
/*--------------------------------------------
 
 --------------------------------------------*/
type0 ForceFieldEAMFit::dEHH(type0 r)
{
    if(r>=2.8) return 0.0;
    
    return -0.3714821016657931 + (-183.34491688140423 + 251.16866128123107*r)*exp(-2.5871227340727088*r);
}
/*--------------------------------------------
 
 --------------------------------------------*/
type0 ForceFieldEAMFit::ddEHH(type0 r)
{
    if(r>=2.8) return 0.0;
    
    return (725.5044639217831 - 649.8041536872805*r)*exp(-2.5871227340727088*r);
}
/*--------------------------------------------
 
 --------------------------------------------*/
type0 ForceFieldEAMFit::calc_rho(elem_type i,elem_type j,type0 r)
{
    if(r>=rho_cut[i][j]) return 0.0;
    type0 ans=0.0;
    for(size_t __i=0;__i<rho_sz[i][j] && r<rho_R[i][j][__i];__i++)
        ans+=rho_A[i][j][__i]*(rho_R[i][j][__i]-r)*(rho_R[i][j][__i]-r)*(rho_R[i][j][__i]-r);
    return ans;
}
/*--------------------------------------------
 
 --------------------------------------------*/
type0 ForceFieldEAMFit::calc_drho(elem_type i,elem_type j,type0 r)
{
    if(r>=rho_cut[i][j]) return 0.0;
    type0 ans=0.0;
    for(size_t __i=0;__i<rho_sz[i][j] && r<rho_R[i][j][__i];__i++)
        ans-=3.0*rho_A[i][j][__i]*(rho_R[i][j][__i]-r)*(rho_R[i][j][__i]-r);
    return ans;
}
/*--------------------------------------------
 
 --------------------------------------------*/
type0 ForceFieldEAMFit::calc_phi(elem_type i,elem_type j,type0 r)
{
    if(r>=phi_cut[i][j]) return 0.0;
    /*
    if(i==1 && j==1)
    {
        type0 rho=calc_rho(1,1,r);
        if(rho<0.0) rho=0.0;
        return EHH(r)-2.0*calc_F(1,rho);
    }*/
    type0 ans=0.0;
    for(size_t __i=0;__i<phi_sz[i][j] && r<phi_R[i][j][__i];__i++)
        ans+=phi_A[i][j][__i]*(phi_R[i][j][__i]-r)*(phi_R[i][j][__i]-r)*(phi_R[i][j][__i]-r);
    return ans;
    
}
/*--------------------------------------------
 
 --------------------------------------------*/
type0 ForceFieldEAMFit::calc_dphi(elem_type i,elem_type j,type0 r)
{
    if(r>=phi_cut[i][j]) return 0.0;
    /*
    if(i==1 && j==1)
    {
        type0 rho=calc_rho(1,1,r);
        if(rho<0.0) rho=0.0;
        return dEHH(r)-2.0*calc_dF(1,rho)*calc_drho(1,1,r);
    }
    */
    type0 ans=0.0;
    for(size_t __i=0;__i<phi_sz[i][j] && r<phi_R[i][j][__i];__i++)
        ans-=3.0*phi_A[i][j][__i]*(phi_R[i][j][__i]-r)*(phi_R[i][j][__i]-r);
    return ans;
    
}





























/*--------------------------------------------
 python constructor
 --------------------------------------------*/
void ForceFieldEAMFit::ml_new(PyMethodDef& tp_methods)
{
    tp_methods.ml_flags=METH_VARARGS | METH_KEYWORDS;
    tp_methods.ml_name="ff_eam_fit";
    
    tp_methods.ml_meth=(PyCFunction)(PyCFunctionWithKeywords)(
    [](PyObject* self,PyObject* args,PyObject* kwds)->PyObject*
    {
        AtomsMD::Object* __self=reinterpret_cast<AtomsMD::Object*>(self);
        
        FuncAPI<type0(***)[2],symm<type0(***)[2]>,type0(*)[3]>f("ff_eam_fit",{"rho_AR","phi_AR","F_A"});
        
        
        const size_t nelems=__self->atoms->elements.nelems;
        if(f(args,kwds)) return NULL;

        size_t** rho_sz=NULL;
        Memory::alloc(rho_sz,nelems,nelems);
        size_t** phi_sz=NULL;
        Memory::alloc(phi_sz,nelems,nelems);
        for(size_t ielem=0;ielem<nelems;ielem++)
            for(size_t jelem=0;jelem<nelems;jelem++)
            {
                rho_sz[ielem][jelem]=f.v<0>()[ielem][jelem].size;
                phi_sz[ielem][jelem]=f.v<1>()[ielem][jelem].size;
            }
        
        
        delete __self->ff;
        __self->ff=new ForceFieldEAMFit(__self->atoms,f.mov<0>(),std::move(rho_sz),f.mov<1>(),std::move(phi_sz),f.mov<2>());

        Py_RETURN_NONE;
    });
    
    tp_methods.ml_doc=(char*)"";
}


