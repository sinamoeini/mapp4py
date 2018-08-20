#include "ff_eam_fit_o.h"
#include "neighbor_md.h"
#include "atoms_md.h"
#include "elements.h"
#include "memory.h"
#include "dynamic_md.h"
#define rho_Fe_gauge 26.0
type0 ForceFieldEAMFitO::rc_rho_H=2.8;
type0 ForceFieldEAMFitO::rc_rho_Fe=4.2;
type0 ForceFieldEAMFitO::rc_phi_FeFe=5.3;
type0 ForceFieldEAMFitO::rc_phi_FeH=4.2;
type0 ForceFieldEAMFitO::rc_phi_HH=2.8;

using namespace MAPP_NS;
type0 ForceFieldEAMFitO::r0_FeFe=1.0;
type0 ForceFieldEAMFitO::r1_FeFe=2.05;
type0 ForceFieldEAMFitO::AR_rho_Fe[3][2]=
{{11.686859407970/rho_Fe_gauge,2.4},
    {-0.01471074009883/rho_Fe_gauge,3.2},
    {0.47193527075943/rho_Fe_gauge,4.2}};

type0 ForceFieldEAMFitO::AR_phi_FeFe[13][2]=
{{-27.444805994228,2.2},
    {15.738054058489,2.3},
    {2.2077118733936,2.4},
    {-2.4989799053251,2.5},
    {4.2099676494795,2.6},
    {-0.77361294129713,2.7},
    {0.80656414937789,2.8},
    {-2.3194358924605,3.0},
    {2.6577406128280,3.3},
    {-1.0260416933564,3.7},
    {0.35018615891957,4.2},
    {-0.058531821042271,4.7},
    {-0.0030458824556234,5.3}};
type0 ForceFieldEAMFitO::B_phi_FeFe[4]={7.4122709384068,-0.64180690713367,-2.6043547961722,0.6262539393123};
type0 ForceFieldEAMFitO::A_F_Fe[3]={-1.0*sqrt(rho_Fe_gauge),-6.7314115586063e-4*rho_Fe_gauge*rho_Fe_gauge,7.6514905604792e-8*rho_Fe_gauge*rho_Fe_gauge*rho_Fe_gauge*rho_Fe_gauge};
#ifdef FeH_SPLINE
type0 ForceFieldEAMFitO::R_phi_FeH[nphi_FeH]={1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,2.0,2.2,2.4,2.6,2.8,3.0,3.3,3.6,3.9,4.2};
//type0 ForceFieldEAMFitO::R_phi_FeH[nphi_FeH]={1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2., 2.1, 2.2, 2.3,2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 3., 3.1, 3.2, 3.3, 3.4, 3.5, 3.6, 3.7,3.8, 3.9, 4., 4.1, 4.2};
#endif
/*--------------------------------------------
 This is for my personal and fitting of iron
 and hydrogen
 --------------------------------------------*/
/*--------------------------------------------
 constructor
 --------------------------------------------*/
ForceFieldEAMFitO::ForceFieldEAMFitO(AtomsMD* __atoms,
type0(&__A_rho_H)[nrho_H],type0(&__A_phi_FeH)[nphi_FeH],type0(&__A_phi_HH)[nphi_HH],type0(&__A_F_H)[nF_H]):
ForceFieldMD(__atoms)
{
    A_rho_H=v+rho_H_offset;
    A_phi_FeH=v+phi_FeH_offset;
    A_phi_HH=v+phi_HH_offset;
    A_F_H=v+F_H_offset;
    
    dA_rho_H=dv+rho_H_offset;
    dA_phi_FeH=dv+phi_FeH_offset;
    dA_phi_HH=dv+phi_HH_offset;
    dA_F_H=dv+F_H_offset;
    
    
    dA_rho_H_lcl=dv_lcl+rho_H_offset;
    dA_phi_FeH_lcl=dv_lcl+phi_FeH_offset;
    dA_phi_HH_lcl=dv_lcl+phi_HH_offset;
    dA_F_H_lcl=dv_lcl+F_H_offset;
    
    
    
    Algebra::V_eq<nrho_H>(__A_rho_H,A_rho_H);
    Algebra::V_eq<nphi_FeH>(__A_phi_FeH,A_phi_FeH);
    Algebra::V_eq<nphi_HH>(__A_phi_HH,A_phi_HH);
    Algebra::V_eq<nF_H>(__A_F_H,A_F_H);
    
    if(nelems==2)
    {
        cut[0][0]=rc_phi_FeFe;
        cut[0][1]=cut[1][0]=rc_phi_FeH;
        cut[1][1]=rc_phi_HH;
    }
    else
        cut[0][0]=rc_phi_FeFe;
    
    for(size_t ielem=0;ielem<nelems;ielem++)
        for(size_t jelem=ielem;jelem<nelems;jelem++)
            cut_sq[ielem][jelem]=cut_sq[jelem][ielem]=cut[ielem][jelem]*cut[ielem][jelem];
    
    
    Algebra::Do<nHvars>::func([this](int i){dof[i]=true;});
    
}
/*--------------------------------------------
 destructor
 --------------------------------------------*/
ForceFieldEAMFitO::~ForceFieldEAMFitO()
{

}
/*--------------------------------------------
 force and energy calculation
 --------------------------------------------*/
void ForceFieldEAMFitO::force_calc()
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
            

            rho[iatm]+=calc_rho(jelem,r);
            if(jatm<natms_lcl)
            {
                rho[jatm]+=calc_rho(ielem,r);
                __vec_lcl[0]+=calc_phi(ielem,jelem,r);;
            }
            else
                __vec_lcl[0]+=0.5*calc_phi(ielem,jelem,r);;
        }
        __vec_lcl[0]+=calc_F(ielem,rho[iatm]);
    }
    
#ifdef NEW_UPDATE
    update(rho_ptr);
#else
    dynamic->update(rho_ptr);
#endif

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
            fpair=-(calc_dF(ielem,rho[iatm])*calc_drho(jelem,r)+
                    calc_dF(jelem,rho[jatm])*calc_drho(ielem,r)+
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
void ForceFieldEAMFitO::energy_calc()
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
            
            
            rho[iatm]+=calc_rho(jelem,r);
            if(jatm<natms_lcl)
            {
                rho[jatm]+=calc_rho(ielem,r);
                __vec_lcl[0]+=calc_phi(ielem,jelem,r);
            }
            else
                __vec_lcl[0]+=0.5*calc_phi(ielem,jelem,r);
        }
        // add the embedded energy here
        __vec_lcl[0]+=calc_F(ielem,rho[iatm]);
    }
}
/*--------------------------------------------
 
 --------------------------------------------*/
type0 ForceFieldEAMFitO::mean_rho_H()
{
    type0* rho=rho_ptr->begin();
    elem_type* evec=atoms->elem->begin();
    const int natms_lcl=atoms->natms_lcl;
    type0 data_lcl[2]={0.0,0.0};
    type0 data[2];
    for(int i=0;i<natms_lcl;i++)
        if(evec[i]==1)
        {
            data_lcl[0]+=rho[i];
            data_lcl[1]++;
        }
    MPI_Allreduce(data_lcl,data,2,Vec<type0>::MPI_T,MPI_SUM,world);
    if(data[1]==0.0) return 0.0;
    return data[0]/data[1];
}
/*--------------------------------------------

 --------------------------------------------*/
void ForceFieldEAMFitO::gradient()
{
    Algebra::zero<nHvars>(dv_lcl);
    type0* xvec=atoms->x->begin();
    type0* rho=rho_ptr->begin();
    elem_type* evec=atoms->elem->begin();
    
    int iatm,jatm;
    elem_type ielem,jelem;
    type0 r,rsq,dFi;
    
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
            
            if(((ielem==0 && jelem==1) || (ielem==1 && jelem==0)) && r<rc_phi_FeH)
                calc_Dphi_FeH(r,jatm<natms_lcl);
            
            if((ielem==1 && jelem==1)  && r<rc_phi_FeH)
                calc_Dphi_HH(r,jatm<natms_lcl);
            
            rho[iatm]+=calc_rho(jelem,r);
            if(jatm<natms_lcl)
                rho[jatm]+=calc_rho(ielem,r);
        }

        if(ielem==1) calc_DF_H(rho[iatm]);
        else calc_DF_Fe(rho[iatm]);
    }
    
    
    
    type0 dFj;
    
    for(iatm=0;iatm<natms_lcl;iatm++)
    {
        ielem=evec[iatm];
        dFi=calc_dF(ielem,rho[iatm]);
        for(int j=0;j<neighbor_list_size[iatm];j++)
        {
            jatm=neighbor_list[iatm][j];
            jelem=evec[jatm];
            rsq=Algebra::RSQ<__dim__>(xvec+iatm*__dim__,xvec+jatm*__dim__);
            if(rsq>=cut_sq[ielem][jelem]) continue;
            r=sqrt(rsq);
           
            
            dFj=jatm<natms_lcl ? calc_dF(jelem,rho[jatm]):0.0;
            
            if(ielem==0 && jelem==0)
            {
                calc_Drho_Fe(dFi+dFj,r);
            }
            else if(ielem==1 && jelem==1)
            {
                calc_Drho_H(dFi+dFj,r);
            }
            else if(ielem==1 && jelem==0)
            {
                calc_Drho_Fe(dFi,r);
                calc_Drho_H(dFj,r);
            }
            else
            {
                calc_Drho_H(dFi,r);
                calc_Drho_Fe(dFj,r);
            }
            
            
        }
    }
    
    MPI_Allreduce(dv_lcl,dv,nHvars,Vec<type0>::MPI_T,MPI_SUM,world);
}
/*--------------------------------------------
 
 --------------------------------------------*/
void ForceFieldEAMFitO::f_sq_gradient()
{
    Algebra::zero<nHvars>(dv_lcl);
    type0* xvec=atoms->x->begin();
    type0* rho=rho_ptr->begin();
    elem_type* evec=atoms->elem->begin();
    
    int iatm,jatm;
    elem_type ielem,jelem;
    type0 r,rsq,dFi,fpair;
    
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
            
            if(((ielem==0 && jelem==1) || (ielem==1 && jelem==0)) && r<rc_phi_FeH)
            calc_Dphi_FeH(r,jatm<natms_lcl);
            
            if((ielem==1 && jelem==1)  && r<rc_phi_FeH)
            calc_Dphi_HH(r,jatm<natms_lcl);
            
            rho[iatm]+=calc_rho(jelem,r);
            if(jatm<natms_lcl)
            rho[jatm]+=calc_rho(ielem,r);
        }
        
        if(ielem==1) calc_DF_H(rho[iatm]);
        else calc_DF_Fe(rho[iatm]);
    }
    
    
    
    type0 dFj;
    
    for(iatm=0;iatm<natms_lcl;iatm++)
    {
        ielem=evec[iatm];
        dFi=calc_dF(ielem,rho[iatm]);
        for(int j=0;j<neighbor_list_size[iatm];j++)
        {
            jatm=neighbor_list[iatm][j];
            jelem=evec[jatm];
            rsq=Algebra::RSQ<__dim__>(xvec+iatm*__dim__,xvec+jatm*__dim__);
            if(rsq>=cut_sq[ielem][jelem]) continue;
            r=sqrt(rsq);
            
            
            fpair=-(calc_dF(ielem,rho[iatm])*calc_drho(jelem,r)+
                    calc_dF(jelem,rho[jatm])*calc_drho(ielem,r)+
                    calc_dphi(ielem,jelem,r));
            
            dFj=jatm<natms_lcl ? calc_dF(jelem,rho[jatm]):0.0;
            
            if(ielem==0 && jelem==0)
            {
                calc_Drho_Fe(dFi+dFj,r);
            }
            else if(ielem==1 && jelem==1)
            {
                calc_Drho_H(dFi+dFj,r);
            }
            else if(ielem==1 && jelem==0)
            {
                calc_Drho_Fe(dFi,r);
                calc_Drho_H(dFj,r);
            }
            else
            {
                calc_Drho_H(dFi,r);
                calc_Drho_Fe(dFj,r);
            }
            
            
        }
    }
    
    MPI_Allreduce(dv_lcl,dv,nHvars,Vec<type0>::MPI_T,MPI_SUM,world);
}
/*--------------------------------------------
 
 --------------------------------------------*/
size_t ForceFieldEAMFitO::get_rFeH(type0*& Rs,type0*& Fs,int*& Ns)
{
    
    
    type0* xvec=atoms->x->begin();
    type0* rho=rho_ptr->begin();
    elem_type* evec=atoms->elem->begin();
    
    int iatm,jatm;
    elem_type ielem,jelem;
    type0 r,rsq,fpair;
    
    int** neighbor_list=neighbor->neighbor_list;
    int* neighbor_list_size=neighbor->neighbor_list_size;
    const int natms_lcl=atoms->natms_lcl;
    type0* rs_lcl=NULL;
    type0* fs_lcl=NULL;
    int sz_lcl=natms_lcl;
    Memory::alloc(rs_lcl,sz_lcl);
    Memory::alloc(fs_lcl,sz_lcl);
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
            fpair=-(calc_dF(ielem,rho[iatm])*calc_drho(jelem,r)+
                    calc_dF(jelem,rho[jatm])*calc_drho(ielem,r)+
                    calc_dphi(ielem,jelem,r));
            if(n_lcl+1>=sz_lcl)
            {
                Memory::grow(rs_lcl,sz_lcl,sz_lcl+100);
                Memory::grow(fs_lcl,sz_lcl,sz_lcl+100);
                sz_lcl+=100;
            }
            fs_lcl[n_lcl]=fpair;
            rs_lcl[n_lcl]=r;
            ++n_lcl;
        }
        
    }
    
    int n;
    MPI_Allreduce(&n_lcl,&n,1,Vec<int>::MPI_T,MPI_SUM,world);
    type0* rs=NULL;
    Memory::alloc(rs,n);
    type0* fs=NULL;
    Memory::alloc(fs,n);
    int __n_lcl;
    type0* __rs=rs;
    type0* __fs=fs;
    for(int i=0;i<atoms->comm_size;i++)
    {
        __n_lcl=n_lcl;
        MPI_Bcast(&__n_lcl,1,Vec<int>::MPI_T,i,world);
        if(atoms->comm_rank==i)
        {
            memcpy(__rs,rs_lcl,__n_lcl*sizeof(type0));
            memcpy(__fs,fs_lcl,__n_lcl*sizeof(type0));
        }
        
        MPI_Bcast(__rs,__n_lcl,Vec<type0>::MPI_T,i,world);
        MPI_Bcast(__fs,__n_lcl,Vec<type0>::MPI_T,i,world);
        __rs+=__n_lcl;
        __fs+=__n_lcl;
    }
    Memory::dealloc(rs_lcl);
    Memory::dealloc(fs_lcl);
    
    size_t* key=NULL;
    Memory::alloc(key,n);
    for(size_t i=0;i<n;i++) key[i]=i;
    
    XMath::quicksort(key,key+n,
    [&rs](size_t* rank_i,size_t* rank_j){return (rs[*rank_i]<rs[*rank_j]);},
    [&rs,&fs](size_t* rank_i,size_t* rank_j)
    {
        std::swap(rs[*rank_i],rs[*rank_j]);
        std::swap(fs[*rank_i],fs[*rank_j]);
    });
    
    Memory::dealloc(key);
    
    
    type0 tol=1.0e-5;
    int i0=0,i1=0;
    type0 r0,f0,n0;
    size_t sz=0;
    
    type0* __Rs=NULL;
    Memory::alloc(__Rs,n);
    type0* __Fs=NULL;
    Memory::alloc(__Fs,n);
    int* __Ns=NULL;
    Memory::alloc(__Ns,n);
    
    while(i0<n)
    {
        r0=0.0;
        f0=0.0;
        n0=0.0;
        while (i1<n && fabs(rs[i1]-rs[i0])<tol)
        {
            r0+=rs[i1];
            f0+=fs[i1];
            n0++;
            i1++;
        }
        r0/=n0;
        f0/=n0;
        
        __Rs[sz]=r0;
        __Fs[sz]=f0;
        __Ns[sz]=static_cast<int>(n0);
        sz++;
        i0=i1;
    }
    Memory::dealloc(rs);
    Memory::dealloc(fs);
    
    Rs=NULL;
    Memory::alloc(Rs,sz);
    Fs=NULL;
    Memory::alloc(Fs,sz);
    Ns=NULL;
    Memory::alloc(Ns,sz);
    memcpy(Rs,__Rs,sz*sizeof(type0));
    memcpy(Fs,__Fs,sz*sizeof(type0));
    memcpy(Ns,__Ns,sz*sizeof(int));
    
    Memory::dealloc(__Rs);
    Memory::dealloc(__Fs);
    Memory::dealloc(__Ns);
    
    
    return sz;
    
}
/*--------------------------------------------
 initiate before a run
 --------------------------------------------*/
void ForceFieldEAMFitO::init()
{
    pre_init();
    rho_ptr=new Vec<type0>(atoms,1,"rho");
    //S_ptr=new Vec<type0>(atoms,1+__nvoigt__);
}
/*--------------------------------------------
 after a run
 --------------------------------------------*/
void ForceFieldEAMFitO::fin()
{
    //delete S_ptr;
    delete rho_ptr;
    post_fin();
}
/*--------------------------------------------
 init xchng
 --------------------------------------------*/
void ForceFieldEAMFitO::init_xchng()
{
}
/*--------------------------------------------
 fin xchng
 --------------------------------------------*/
void ForceFieldEAMFitO::fin_xchng()
{
}
/*--------------------------------------------
 pre xchng energy
 --------------------------------------------*/
void ForceFieldEAMFitO::pre_xchng_energy(GCMC*)
{
}
/*--------------------------------------------
 xchng energy
 --------------------------------------------*/
type0 ForceFieldEAMFitO::xchng_energy(GCMC*)
{
    return 0.0;
}
/*--------------------------------------------
 post xchng energy
 --------------------------------------------*/
void ForceFieldEAMFitO::post_xchng_energy(GCMC*)
{
}
/*--------------------------------------------
 E_HH
 --------------------------------------------*/
type0 ForceFieldEAMFitO::E_HH(type0 r)
{
    return (38.00339647547726-100.67818306817757*r)*exp(-2.758422844154626*r);
}
/*--------------------------------------------
 E_HH
 --------------------------------------------*/
type0 ForceFieldEAMFitO::dE_HH(type0 r)
{
    return (-205.50762006159945+277.7130000832425*r)*exp(-2.758422844154626*r);
}
/*--------------------------------------------
 E_HH
 --------------------------------------------*/
type0 ForceFieldEAMFitO::E_HH_(type0 r)
{
    //return (38.00339647547726-100.67818306817757*r)*exp(-2.758422844154626*r);
    
    return exp(1.0/(r-2.8 )-2.4809870256282154*r)*(48.796680008387405-131.20995880102427*r);
}
/*--------------------------------------------
 E_HH
 --------------------------------------------*/
type0 ForceFieldEAMFitO::dE_HH_(type0 r)
{
    return (exp(1.0/(r-2.8 )-2.4809870256282154*r)*(-2026.6239681656182 + r*(4096.100546537651 + r*(-2075.2430391394673 + 325.5302054185539*r))))/((r-2.8)*(r-2.8));
}
/*--------------------------------------------
 calc F
 --------------------------------------------*/
type0 ForceFieldEAMFitO::calc_F(elem_type i,type0 rho)
{
    if(i==0)
    {
        rho/=A_F_H[4];
        return A_F_Fe[0]*sqrt(rho)+A_F_Fe[1]*rho*rho+A_F_Fe[2]*rho*rho*rho*rho;
    }

    if(rho==0.0) return 0.0;
    type0 rhob=rho/A_F_H[3];
    return (A_F_H[0]*rhob+A_F_H[1])*rhob+A_F_H[2]*pow(rhob,2.0/3.0)*(rhob/(1.0+rhob));
    
    /*
    type0 rhob=rho/A_F_H[1];
    type0 eta=pow(rhob,A_F_H[2]);
    type0 falpha=pow((1.0-3.0/(2.0*A_F_H[2]-1.0))/4.0,1.0-2.0/A_F_H[2]);
    return A_F_H[0]*(eta*(eta-2.0)+falpha*rhob*rhob);
     */
}
/*--------------------------------------------
 calc dF
 --------------------------------------------*/
type0 ForceFieldEAMFitO::calc_dF(elem_type i,type0 rho)
{
    if(i==0)
    {
        rho/=A_F_H[4];
        return (0.5*A_F_Fe[0]/sqrt(rho)+2.0*A_F_Fe[1]*rho+4.0*A_F_Fe[2]*rho*rho*rho)/A_F_H[4];
    }
    
    if(rho==0.0) return 0.0;
    type0 rhob=rho/A_F_H[3];
    return (2.0*A_F_H[0]*rhob+A_F_H[1]+A_F_H[2]*pow(rhob,2.0/3.0)*(5.0+2.0*rhob)/(3.0*(1.0+rhob)*(1.0+rhob)))/A_F_H[3];

    
    
    /*
    type0 rhob=rho/A_F_H[1];
    type0 eta=pow(rhob,A_F_H[2]);
    type0 falpha=pow((1.0-3.0/(2.0*A_F_H[2]-1.0))/4.0,1.0-2.0/A_F_H[2]);
    return 2.0*A_F_H[0]*(A_F_H[2]*eta*(eta-1.0)+falpha*rhob*rhob)/rho;
    */
    
}
/*--------------------------------------------
 calc ddF
 --------------------------------------------*/
type0 ForceFieldEAMFitO::calc_ddF(elem_type i,type0 rho)
{
    if(i==0)
    {
        rho/=A_F_H[4];
        return (-0.25*A_F_Fe[0]/sqrt(rho*rho*rho)+2.0*A_F_Fe[1]+12.0*A_F_Fe[2]*rho*rho)/(A_F_H[4]*A_F_H[4]);
    }
    
    if(rho==0.0) return 0.0;
    type0 rhob=rho/A_F_H[3];
    return (2.0*A_F_H[0]-A_F_H[2]*2.0*((rhob+5.0)*rhob-5.0)/(9.0*pow(rhob,1.0/3.0)*(1.0+rhob)*(1.0+rhob)*(1.0+rhob)))/(A_F_H[3]*A_F_H[3]);
     
    
    /*
    type0 rhob=rho/A_F_H[1];
    type0 eta=pow(rhob,A_F_H[2]);
    type0 falpha=pow((1.0-3.0/(2.0*A_F_H[2]-1.0))/4.0,1.0-2.0/A_F_H[2]);
    return 2.0*A_F_H[0]*(A_F_H[2]*eta*((2.0*A_F_H[2]-1.0)*eta+1.0-A_F_H[2])+falpha*rhob*rhob)/(rho*rho);
    */
}
/*--------------------------------------------
 clac rho
 --------------------------------------------*/
type0 ForceFieldEAMFitO::calc_rho(elem_type ielem,type0 r)
{
    
    if(ielem==0)
    {
        if(r>=rc_rho_Fe) return 0.0;
        return A_F_H[4]*spline(AR_rho_Fe,r);
    }
    if(r>=rc_rho_H) return 0.0;
    
    type0 y=exp(-A_rho_H[1]*r);
    return A_F_H[3]*A_rho_H[0]*exp(1.0/(r-rc_rho_H))*Algebra::pow<6>(r)*y*(1.0+512.0*y);
}

/*--------------------------------------------
 clac drho
 --------------------------------------------*/
type0 ForceFieldEAMFitO::calc_drho(elem_type ielem,type0 r)
{
    
    if(ielem==0)
    {
        if(r>=rc_rho_Fe) return 0.0;
        return A_F_H[4]*dspline(AR_rho_Fe,r);
    }
    
    if(r>=rc_rho_H) return 0.0;
    type0 y=exp(-A_rho_H[1]*r);
    type0 rho=A_F_H[3]*A_rho_H[0]*exp(1.0/(r-rc_rho_H))*Algebra::pow<6>(r)*y*(1.0+512.0*y);
    return rho*(6.0/r-1.0/((r-rc_rho_H)*(r-rc_rho_H))-A_rho_H[1]*(1.0+y/(0.001953125+y)));
}
/*--------------------------------------------
 clac phi
 --------------------------------------------*/
type0 ForceFieldEAMFitO::calc_phi(elem_type ielem,elem_type jelem,type0 r)
{
    
    if(ielem==0 && jelem==0)
    {
        if(r>=rc_phi_FeFe) return 0.0;
        if(r>=r1_FeFe) return spline(AR_phi_FeFe,r);
        if(r>=r0_FeFe) return exp(B_phi_FeFe[0]+r*(B_phi_FeFe[1]+r*(B_phi_FeFe[2]+r*B_phi_FeFe[3])));
        return xi(r);
    }
    if(ielem==1 && jelem==1)
    {
        
        type0 c0=tanh(A_phi_HH[2]*(r-A_phi_HH[3]));
        type0 rhoH=calc_rho(1,r);

        return
        0.5*(1.0-c0)*(E_HH(r)-2.0*calc_F(1,rhoH))+
        0.5*(1.0+c0)*(A_phi_HH[0]*exp(1.0/(r-rc_phi_HH))+(A_phi_HH[1]/A_F_H[3])*rhoH);
    }
    
#ifdef FeH_SPLINE
    if(r>=rc_phi_FeH) return 0.0;
    return spline(A_phi_FeH,R_phi_FeH,r);
#else
    if(r>=rc_phi_FeH) return 0.0;
    type0 y=exp(-A_phi_FeH[1]*(r-A_phi_FeH[2]));
    return exp(1.0/(r-rc_phi_FeH))*(A_phi_FeH[0]*y*(y-2.0))+A_phi_FeH[3]*spline(AR_rho_Fe,r);
#endif
}
/*--------------------------------------------
 clac phi
 --------------------------------------------*/
type0 ForceFieldEAMFitO::calc_dphi(elem_type ielem,elem_type jelem,type0 r)
{
    if(ielem==0 && jelem==0)
    {
        if(r>=rc_phi_FeFe) return 0.0;
        if(r>=r1_FeFe) return dspline(AR_phi_FeFe,r);
        if(r>=r0_FeFe) return (B_phi_FeFe[1]+r*(2.0*B_phi_FeFe[2]+r*3.0*B_phi_FeFe[3]))
            *exp(B_phi_FeFe[0]+r*(B_phi_FeFe[1]+r*(B_phi_FeFe[2]+r*B_phi_FeFe[3])));
        
        return dxi(r);
    }
    if(ielem==1 && jelem==1)
    {
        
        type0 c0=tanh(A_phi_HH[2]*(r-A_phi_HH[3]));
        type0 dc0=A_phi_HH[2]*(1.0-c0*c0);
        type0 rhoH=calc_rho(1,r);
        type0 drhoH=calc_drho(1,r);
        type0 fc=exp(1.0/(r-rc_phi_HH));
        type0 dfc=-fc/((r-rc_phi_HH)*(r-rc_phi_HH));
        
        
        return
        -0.5*dc0*(E_HH(r)-2.0*calc_F(1,rhoH))+
        0.5*(1.0-c0)*(dE_HH(r)-2.0*drhoH*calc_dF(1,rhoH))+
        0.5*dc0*(A_phi_HH[0]*fc+(A_phi_HH[1]/A_F_H[3])*rhoH)+
        0.5*(1.0+c0)*(A_phi_HH[0]*dfc+(A_phi_HH[1]/A_F_H[3])*drhoH);

    }
#ifdef FeH_SPLINE
    if(r>=rc_phi_FeH) return 0.0;
    return dspline(A_phi_FeH,R_phi_FeH,r);
#else
    if(r>=rc_phi_FeH) return 0.0;
    type0 y=exp(-A_phi_FeH[1]*(r-A_phi_FeH[2]));
    type0 c1=1.0/((r-rc_phi_FeH)*(r-rc_phi_FeH));
    type0 c0=c1+A_phi_FeH[1];
    return exp(1.0/(r-rc_phi_FeH))*(A_phi_FeH[0]*y*(-(A_phi_FeH[1]+c0)*y+2.0*c0))+A_phi_FeH[3]*dspline(AR_rho_Fe,r);
#endif
}
/*--------------------------------------------

--------------------------------------------*/
void ForceFieldEAMFitO::calc_Dphi_FeH(type0 r,bool cond)
{
    type0 coef= cond ? 1.0:0.5;
#ifdef FeH_SPLINE
    Dspline(coef,dA_phi_FeH_lcl,R_phi_FeH,r);
#else
    type0 y=exp(-A_phi_FeH[1]*(r-A_phi_FeH[2]));
    type0 fc=exp(1.0/(r-rc_phi_FeH));
    dA_phi_FeH_lcl[0]+=coef*fc*y*(y-2.0);
    dA_phi_FeH_lcl[1]+=-coef*A_phi_FeH[0]*fc*2.0*y*(y-1.0)*(r-A_phi_FeH[2]);
    dA_phi_FeH_lcl[2]+=coef*A_phi_FeH[0]*fc*2.0*y*(y-1.0)*A_phi_FeH[1];
    dA_phi_FeH_lcl[3]+=coef*spline(AR_rho_Fe,r);
#endif
}
/*--------------------------------------------
 
 --------------------------------------------*/
void ForceFieldEAMFitO::calc_Dphi_HH(type0 r,bool cond)
{
    type0 coef= cond ? 1.0:0.5;
    type0 c0=tanh(A_phi_HH[2]*(r-A_phi_HH[3]));
    type0 rhoH=calc_rho(1,r);
    type0 fc=exp(1.0/(r-rc_phi_HH));
    type0 c1=A_phi_HH[0]*fc+(A_phi_HH[1]/A_F_H[3])*rhoH-E_HH(r)+2.0*calc_F(1,rhoH);
    dA_phi_HH_lcl[0]+=coef*0.5*(1.0+c0)*fc;
    dA_phi_HH_lcl[1]+=coef*0.5*(1.0+c0)*rhoH/A_F_H[3];
    dA_F_H_lcl[3]+=-A_phi_HH[1]*coef*0.5*(1.0+c0)*rhoH/(A_F_H[3]*A_F_H[3]);
    dA_phi_HH_lcl[2]+=coef*0.5*(r-A_phi_HH[3])*(1.0-c0*c0)*c1;
    dA_phi_HH_lcl[3]+=-coef*0.5*A_phi_HH[2]*(1.0-c0*c0)*c1;
    
    c1=0.5*(1.0+c0)*A_phi_HH[1]/A_F_H[3];
    c0=c0-1.0;
    // c0*calc_F(1,rhoH)+c1*rhoH
    //(c0*calc_dF(1,rhoH)+c1)
    
    type0 rhob=rhoH/A_F_H[3];
    type0 rhob_2_3=pow(rhob,2.0/3.0);
    
    type0 dFH=((2.0*A_F_H[0]*rhob+A_F_H[1]+A_F_H[2]*rhob_2_3*(5.0+2.0*rhob)/(3.0*(1.0+rhob)*(1.0+rhob)))/A_F_H[3]);
    calc_Drho_H(coef*(c0*dFH+c1),r);
    dA_F_H_lcl[0]+=c0*coef*rhob*rhob;
    dA_F_H_lcl[1]+=c0*coef*rhob;
    dA_F_H_lcl[2]+=c0*coef*rhob_2_3*(rhob/(1.0+rhob));
    dA_F_H_lcl[3]+=-c0*coef*rhob*dFH;
     
    /*
    type0 rhob=rhoH/A_F_H[1];
    type0 eta=pow(rhob,A_F_H[2]);
    type0 falpha=pow((1.0-3.0/(2.0*A_F_H[2]-1.0))/4.0,1.0-2.0/A_F_H[2]);
    type0 dFH=2.0*A_F_H[0]*(A_F_H[2]*eta*(eta-1.0)+falpha*rhob*rhob)/rhoH;
    calc_Drho_H((c0*dFH+c1),r);
    type0 dlogfalpha=(2.0*log((1.0-3.0/(2.0*A_F_H[2]-1.0))/4.0)/A_F_H[2]+3.0/(2.0*A_F_H[2]-1.0))/A_F_H[2];
    dA_F_H_lcl[0]+=coef*eta*(eta-2.0)+falpha*rhob*rhob;
    dA_F_H_lcl[1]+=-coef*2.0*A_F_H[0]*(A_F_H[2]*eta*(eta-1.0)+falpha*rhob*rhob)/A_F_H[1];
    dA_F_H_lcl[2]+=coef*A_F_H[0]*(2.0*eta*(eta-1.0)*log(rhob)+dlogfalpha*falpha*rhob*rhob);
     */
}
/*--------------------------------------------
 
 --------------------------------------------*/
void ForceFieldEAMFitO::calc_DF_H(type0 rho)
{
    type0 rhob=rho/A_F_H[3];
    type0 rhob_2_3=pow(rhob,2.0/3.0);
    dA_F_H_lcl[0]+=rhob*rhob;
    dA_F_H_lcl[1]+=rhob;
    dA_F_H_lcl[2]+=rhob_2_3*(rhob/(1.0+rhob));
    dA_F_H_lcl[3]+=-rhob*((2.0*A_F_H[0]*rhob+A_F_H[1]+A_F_H[2]*rhob_2_3*(5.0+2.0*rhob)/(3.0*(1.0+rhob)*(1.0+rhob)))/A_F_H[3]);
    
    /*
    if(rho==0.0) return;
    type0 rhob=rho/A_F_H[1];
    type0 eta=pow(rhob,A_F_H[2]);
    type0 falpha=pow((1.0-3.0/(2.0*A_F_H[2]-1.0))/4.0,1.0-2.0/A_F_H[2]);
    type0 dlogfalpha=(2.0*log((1.0-3.0/(2.0*A_F_H[2]-1.0))/4.0)/A_F_H[2]+3.0/(2.0*A_F_H[2]-1.0))/A_F_H[2];
    dA_F_H_lcl[0]+=eta*(eta-2.0)+falpha*rhob*rhob;
    dA_F_H_lcl[1]+=-2.0*A_F_H[0]*(A_F_H[2]*eta*(eta-1.0)+falpha*rhob*rhob)/A_F_H[1];
    dA_F_H_lcl[2]+=A_F_H[0]*(2.0*eta*(eta-1.0)*log(rhob)+dlogfalpha*falpha*rhob*rhob);
     */
}
/*--------------------------------------------
 
 --------------------------------------------*/
void ForceFieldEAMFitO::calc_DF_Fe(type0 rho)
{
    rho/=A_F_H[4];
    dA_F_H_lcl[4]+=-rho*(0.5*A_F_Fe[0]/sqrt(rho)+2.0*A_F_Fe[1]*rho+4.0*A_F_Fe[2]*rho*rho*rho)/A_F_H[4];
}
/*--------------------------------------------
 
 --------------------------------------------*/
void ForceFieldEAMFitO::calc_Drho_H(type0 coef,type0 r)
{
    if(r>rc_rho_H) return;
    type0 y=exp(-A_rho_H[1]*r);
    type0 c0=exp(1.0/(r-rc_rho_H))*Algebra::pow<6>(r);
    dA_rho_H_lcl[0]+=coef*A_F_H[3]*c0*y*(1.0+512.0*y);
    dA_F_H_lcl[3]+=coef*A_rho_H[0]*c0*y*(1.0+512.0*y);
    dA_rho_H_lcl[1]+=-coef*A_rho_H[0]*A_F_H[3]*r*c0*y*(1.0+1024.0*y);
}
/*--------------------------------------------
 
 --------------------------------------------*/
void ForceFieldEAMFitO::calc_Drho_Fe(type0 coef,type0 r)
{
    if(r>rc_rho_Fe) return;
    dA_F_H_lcl[4]+=coef*spline(AR_rho_Fe,r);
}
/*--------------------------------------------
 python constructor
 --------------------------------------------*/
void ForceFieldEAMFitO::ml_new(PyMethodDef& tp_methods)
{
    tp_methods.ml_flags=METH_VARARGS | METH_KEYWORDS;
    tp_methods.ml_name="ff_eam_fit_o";
    tp_methods.ml_meth=(PyCFunction)(PyCFunctionWithKeywords)(
    [](PyObject* self,PyObject* args,PyObject* kwds)->PyObject*
    {
        AtomsMD::Object* __self=reinterpret_cast<AtomsMD::Object*>(self);
        
        FuncAPI<type0[nrho_H],type0[nphi_FeH],type0[nphi_HH],type0[nF_H]> f("ff_eam_fit_o",{"A_rho_H","A_phi_FeH","A_phi_HH","A_F_H"});
        if(f(args,kwds)) return NULL;

        
        delete __self->ff;
        __self->ff=new ForceFieldEAMFitO(__self->atoms,f.val<0>(),f.val<1>(),f.val<2>(),f.val<3>());
        Py_RETURN_NONE;
    });

    tp_methods.ml_doc="";
}
















