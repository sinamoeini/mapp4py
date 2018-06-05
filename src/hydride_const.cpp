#include "hydride_const.h"
#include "memory.h"
#include "xmath.h"
using namespace MAPP_NS;
/*--------------------------------------------
 
 --------------------------------------------*/
HydrideConst::HydrideConst(PotFit* __potfit,
type0(* phi_AR_FeH)[2],size_t phi_FeH_sz,
type0(* rho_AR_H)[2],size_t rho_H_sz,
size_t*&& __phi_const_list,size_t __phi_const_list_sz,
size_t*&& __phi_nconst_list,size_t __phi_nconst_list_sz,
type0*&& __phi_const0,type0**&& __A,type0*&& __rho0s,type0**&& __B,type0**&& __L,type0**&& __D0,type0***&& __C):
potfit(__potfit)
{
    size_t* R_FeH_key;
    Memory::alloc(R_FeH_key,phi_FeH_sz);
    for(size_t i=0;i<phi_FeH_sz;i++)
        R_FeH_key[i]=i;
    
    size_t* R_H_key;
    Memory::alloc(R_H_key,rho_H_sz);
    for(size_t i=0;i<rho_H_sz;i++)
        R_H_key[i]=i;
        
    XMath::quicksort(
    R_FeH_key,R_FeH_key+phi_FeH_sz,
    [&phi_AR_FeH](size_t* rank_i,size_t* rank_j){return (phi_AR_FeH[*rank_i][1]>phi_AR_FeH[*rank_j][1]);},
    [](size_t* rank_i,size_t* rank_j){std::swap(*rank_i,*rank_j);}
    );
    
    XMath::quicksort(
    R_H_key,R_H_key+rho_H_sz,
    [&rho_AR_H](size_t* rank_i,size_t* rank_j){return (rho_AR_H[*rank_i][1]>rho_AR_H[*rank_j][1]);},
    [](size_t* rank_i,size_t* rank_j){std::swap(*rank_i,*rank_j);}
    );
    
    
    phi_const_sz=__phi_const_list_sz;
    phi_nconst_sz=__phi_nconst_list_sz;
    Memory::alloc(phi_const_list,phi_const_sz);
    Memory::alloc(phi_nconst_list,phi_nconst_sz);
    
    for(size_t i=0;i<phi_const_sz;i++) phi_const_list[i]=R_FeH_key[__phi_const_list[i]];
    for(size_t i=0;i<phi_nconst_sz;i++) phi_nconst_list[i]=R_FeH_key[__phi_nconst_list[i]];
    
    Memory::dealloc(__phi_const_list);
    Memory::dealloc(__phi_nconst_list);
    __phi_const_list=__phi_nconst_list=NULL;
    
    rho_sz=rho_H_sz;
    FH_sz=3;
    nFe=3;
    natms=5;
    nH=natms-nFe;
    
    phi_const0=__phi_const0;
    __phi_const0=NULL;
    
    A=__A;
    __A=NULL;
    
    Memory::alloc(B,natms,rho_sz);
    for(int j=0;j<natms;j++)
        for(size_t k=0;k<rho_sz;k++)
            B[j][R_H_key[k]]=__B[j][k];
    
    Memory::dealloc(__B);
    __B=NULL;
    
    L=__L;
    __L=NULL;
    
    D0=__D0;
    __D0=NULL;
    
    Memory::alloc(C,phi_const_sz,natms,rho_sz);
    for(size_t i=0;i<phi_const_sz;i++)
        for(int j=0;j<natms;j++)
            for(size_t k=0;k<rho_sz;k++)
                C[i][j][R_H_key[k]]=__C[i][j][k];
    Memory::dealloc(__C);
    __C=NULL;
    
    Memory::alloc(rho0s,5*natms);
    rhos=rho0s+natms;
    Fs=rhos+natms;
    dFs=Fs+natms;
    ddFs=dFs+natms;
    memcpy(rho0s,__rho0s,natms*sizeof(type0));
    Memory::dealloc(__rho0s);
    __rho0s=NULL;
    
    Memory::dealloc(R_FeH_key);
    Memory::dealloc(R_H_key);
    Memory::alloc(D,phi_const_sz,natms);
    Memory::alloc(Dphi_Drho,phi_const_sz,rho_sz);
    Memory::alloc(Dphi_DFH,phi_const_sz,FH_sz);
}
/*--------------------------------------------
 
 --------------------------------------------*/
HydrideConst::~HydrideConst()
{
    Memory::dealloc(Dphi_DFH);
    Memory::dealloc(Dphi_Drho);
    Memory::dealloc(D);
    Memory::dealloc(L);
    Memory::dealloc(D0);
    Memory::dealloc(C);
    Memory::dealloc(B);
    Memory::dealloc(A);
    Memory::dealloc(phi_const0);
    Memory::dealloc(rho0s);
    Memory::dealloc(phi_const_list);
    Memory::dealloc(phi_nconst_list);
    
}
/*--------------------------------------------
 
 --------------------------------------------*/
void HydrideConst::prep()
{
    for(int i=0;i<natms;i++)
    {
        rhos[i]=rho0s[i];
        for(size_t j=0;j<rho_sz;j++)
            rhos[i]+=B[i][j]*potfit->rho_A[1][0][j];
        for(int j=0;j<phi_const_sz;j++)
        {
            D[j][i]=D0[j][i];
            for(size_t k=0;k<rho_sz;k++)
                D[j][i]+=C[j][i][k]*potfit->rho_A[1][0][k];
        }
        
    }
    

    
    
    
    
    for(int i=0;i<nFe;i++)
    {
        Fs[i]=potfit->ff->calc_F(0,rhos[i]);
        dFs[i]=potfit->ff->calc_dF(0,rhos[i]);
    }
    

    for(int i=nFe;i<natms;i++)
    {
        Fs[i]=potfit->ff->calc_F(1,rhos[i]);
        dFs[i]=potfit->ff->calc_dF(1,rhos[i]);
    }
    
    
   
    

}
/*--------------------------------------------
 
 --------------------------------------------*/
void HydrideConst::adj_deriv()
{
    for(int i=0;i<nFe;i++) ddFs[i]=potfit->ff->calc_ddF(0,rhos[i]);
    
    type0 DFH_DA[3],DdFH_DA[3];
    for(int j=0;j<phi_const_sz;j++)
        for(int k=0;k<FH_sz;k++)Dphi_DFH[j][k]=0.0;
    for(int i=nFe;i<natms;i++)
    {
        ddFs[i]=potfit->ff->calc_ddF(1,rhos[i]);
        potfit->ff->calc_DFH(rhos[i],DFH_DA[0],DFH_DA[1],DFH_DA[2]);
        potfit->ff->calc_DdFH(rhos[i],DdFH_DA[0],DdFH_DA[1],DdFH_DA[2]);
        for(int j=0;j<phi_const_sz;j++)
        for(int k=0;k<FH_sz;k++)
        Dphi_DFH[j][k]+=L[j][i]*DFH_DA[k]+D[j][i]*DdFH_DA[k];
        
    }

    for(int i=0;i<phi_const_sz;i++)
        for(int j=0;j<rho_sz;j++)
        {
            Dphi_Drho[i][j]=0.0;
            for(int k=0;k<natms;k++)
                Dphi_Drho[i][j]+=(B[k][j]*L[i][k]+C[i][k][j])*dFs[k]+B[k][j]*D[i][k]*ddFs[k];
        }
    
    
    
    
    
    for(int i=0;i<rho_sz;i++)
        if(potfit->rho_A_dof[1][0][i])
        {
            for(int j=0;j<phi_const_sz;j++)
                potfit->rho_A_f[1][0][i]+=potfit->phi_A_f[1][0][phi_const_list[j]]*Dphi_Drho[j][i];
        }
    
    for(int i=0;i<FH_sz;i++)
        if(potfit->F_A_dof[1][i])
        {
            for(int j=0;j<phi_const_sz;j++)
                potfit->F_A_f[1][i]+=potfit->phi_A_f[1][0][phi_const_list[j]]*Dphi_DFH[j][i];
        }
    
    for(int i=0;i<phi_nconst_sz;i++)
        if(potfit->phi_A_dof[1][0][phi_nconst_list[i]])
        {
            for(int j=0;j<phi_const_sz;j++)
                potfit->phi_A_f[1][0][phi_nconst_list[i]]+=potfit->phi_A_f[1][0][phi_const_list[j]]*A[j][i];
        }
    
    
    for(int j=0;j<phi_const_sz;j++)
        potfit->phi_A_f[1][0][phi_const_list[j]]=0.0;
    
    
    
    for(int i=0;i<rho_sz;i++)
        potfit->rho_A_f[1][1][i]=potfit->rho_A_f[1][0][i];
    for(int i=0;i<phi_nconst_sz+phi_const_sz;i++)
        potfit->phi_A_f[0][1][i]=potfit->phi_A_f[1][0][i];
}
/*--------------------------------------------
 
 --------------------------------------------*/
void HydrideConst::update()
{

    
    for(int i=0;i<phi_const_sz;i++)
    {
        potfit->phi_A[1][0][phi_const_list[i]]=phi_const0[i];
        
        for(int j=0;j<phi_nconst_sz;j++)
            potfit->phi_A[1][0][phi_const_list[i]]+=A[i][j]*potfit->phi_A[1][0][phi_nconst_list[j]];
         
        for(int j=0;j<natms;j++)
            potfit->phi_A[1][0][phi_const_list[i]]+=L[i][j]*Fs[j]+D[i][j]*dFs[j];
        
        potfit->phi_A[0][1][phi_const_list[i]]=potfit->phi_A[1][0][phi_const_list[i]];
    }

    
}
/*--------------------------------------------
 
 --------------------------------------------*/
void HydrideConst::readj_dof()
{
    for(int i=0;i<phi_const_sz;i++)
        potfit->phi_A_dof[0][1][phi_const_list[i]]=potfit->phi_A_dof[1][0][phi_const_list[i]]=true;
    
}



