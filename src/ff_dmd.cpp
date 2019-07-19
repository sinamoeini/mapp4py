#include "ff_dmd.h"
#include "atoms_dmd.h"
#include "neighbor_dmd.h"
#include "xmath.h"
#include "memory.h"
using namespace MAPP_NS;
/*--------------------------------------------
 
 --------------------------------------------*/
ForceFieldDMD::ForceFieldDMD(AtomsDMD* __atoms):
ForceField(__atoms),
atoms(__atoms),
rsq_crd(NULL),
r_crd(NULL),
ave_mu(NULL)
{
    Memory::alloc(lambda,nelems);
    Memory::alloc(mu_0,nelems);
    Memory::alloc(ave_mu,nelems);
    Memory::alloc(cut_sk,nelems,nelems);
    Memory::alloc(rsq_crd,nelems);
    Memory::alloc(r_crd,nelems);
    Memory::alloc(ext_mu,nelems);
    neighbor=new NeighborDMD(__atoms,cut_sk,rsq_crd);
    f=new Vec<type0>(atoms,__dim__,"f");
    f_alpha=new DMDVec<type0>(atoms,0.0,"f_alpha");
    c_d=new DMDVec<type0>(atoms,0.0,"c_d");
    f_c=new DMDVec<type0>(atoms,0.0,"mu");
    for(size_t i=0;i<nelems;i++) ext_mu[i]=0.0;
}
/*--------------------------------------------
 
 --------------------------------------------*/
ForceFieldDMD::~ForceFieldDMD()
{
    delete f_c;
    delete c_d;
    delete f_alpha;
    delete f;
    delete neighbor;
    Memory::dealloc(ext_mu);
    Memory::dealloc(r_crd);
    Memory::dealloc(rsq_crd);
    Memory::dealloc(cut_sk);
    Memory::dealloc(ave_mu);
    Memory::dealloc(mu_0);
    Memory::dealloc(lambda);
}
/*--------------------------------------------
 
 --------------------------------------------*/
void ForceFieldDMD::calc_ndof()
{
    nx_dof=atoms->natms*__dim__;
    if(!atoms->x_dof->is_empty())
    {
        const int n=atoms->natms_lcl*__dim__;
        int nx_dof_lcl=n;
        bool* x_dof=atoms->x_dof->begin();
        for(int i=0;i<n;i++) if(!x_dof[i]) nx_dof_lcl--;
        MPI_Allreduce(&nx_dof_lcl,&nx_dof,1,Vec<int>::MPI_T,MPI_SUM,world);
    }
    
    nalpha_dof=atoms->natms*atoms->c_dim;
    if(!atoms->alpha_dof->is_empty())
    {
        const int n=atoms->natms_lcl*atoms->c_dim;
        int nalpha_dof_lcl=n;
        bool* alpha_dof=atoms->alpha_dof->begin();
        for(int i=0;i<n;i++) if(!alpha_dof[i]) nalpha_dof_lcl--;
        MPI_Allreduce(&nalpha_dof_lcl,&nalpha_dof,1,Vec<int>::MPI_T,MPI_SUM,world);
    }
    
    
    nc_dof=atoms->natms*atoms->c_dim;
    if(!atoms->c_dof->is_empty())
    {
        const int n=atoms->natms_lcl*atoms->c_dim;
        int nc_dof_lcl=n;
        bool* c_dof=atoms->c_dof->begin();
        for(int i=0;i<n;i++) if(!c_dof[i]) nc_dof_lcl--;
        MPI_Allreduce(&nc_dof_lcl,&nc_dof,1,Vec<int>::MPI_T,MPI_SUM,world);
    }
}
/*--------------------------------------------
 
 --------------------------------------------*/
void ForceFieldDMD::reset()
{
    type0* __f=f->begin();
    const int n=atoms->natms_lcl*__dim__;
    for(int i=0;i<n;i++) __f[i]=0.0;
    __f=f_alpha->begin();
    const int m=atoms->natms_lcl*c_dim;
    for(int i=0;i<m;i++) __f[i]=0.0;
    for(int i=0;i<__nvoigt__+2;i++) __vec_lcl[i]=0.0;
}
/*--------------------------------------------
 
 --------------------------------------------*/
void ForceFieldDMD::impose_dof(type0* __f,type0* __f_alpha)
{
    if(!dof_empty)
    {
        bool* dof=atoms->x_dof->begin();
        const int n=atoms->natms_lcl*__dim__;
        for(int i=0;i<n;i++) __f[i]=dof[i] ? __f[i]:0.0;
    }
    
    if(!dof_alpha_empty)
    {
        bool* dof=atoms->alpha_dof->begin();
        const int n=atoms->natms_lcl*c_dim;
        for(int i=0;i<n;i++) __f_alpha[i]=dof[i] ? __f_alpha[i]:0.0;
    }
}
/*--------------------------------------------
 
 --------------------------------------------*/
void ForceFieldDMD::impose_c_dof(type0* __mu)
{
    if(!dof_c_empty)
    {
        bool* dof=atoms->c_dof->begin();
        const int n=atoms->natms_lcl*c_dim;
        for(int i=0;i<n;i++) __mu[i]=dof[i] ? __mu[i]:0.0;
    }
}
/*--------------------------------------------
 
 --------------------------------------------*/
void ForceFieldDMD::norm_sq(type0* __f,type0& __err_sq_x,type0* __f_alpha,type0& __err_sq_alpha)
{
    type0 err_lcl=0.0;
    int n=atoms->natms_lcl*__dim__;
    for(int i=0;i<n;i++) err_lcl+=__f[i]*__f[i];
    MPI_Allreduce(&err_lcl,&__err_sq_x,1,Vec<type0>::MPI_T,MPI_SUM,world);
    
    err_lcl=0.0;
    n=atoms->natms_lcl*c_dim;
    for(int i=0;i<n;i++) err_lcl+=__f_alpha[i]*__f_alpha[i];
    
    MPI_Allreduce(&err_lcl,&__err_sq_alpha,1,Vec<type0>::MPI_T,MPI_SUM,world);
}
/*--------------------------------------------
 
 --------------------------------------------*/
void ForceFieldDMD::reset_c_d()
{
    type0* __c_d=c_d->begin();
    const int n=atoms->natms_lcl*c_dim;
    for(int i=0;i<n;i++) __c_d[i]=0.0;
}
/*--------------------------------------------
 
 --------------------------------------------*/
void ForceFieldDMD::pre_init()
{
    dof_empty=atoms->x_dof->is_empty();
    dof_alpha_empty=atoms->alpha_dof->is_empty();
    dof_c_empty=atoms->c_dof->is_empty();
    
    type0 tmp;
    max_cut=0.0;
    kBT=atoms->kB*atoms->temp;
    beta=1.0/kBT;
    for(size_t i=0;i<nelems;i++)
    {
        lambda[i]=atoms->hP/sqrt(2.0*M_PI*kBT*atoms->elements.masses[i]);
        mu_0[i]=kBT*1.5*(2.0*log(lambda[i])-1.0-log(M_PI));
        for(size_t j=0;j<i+1;j++)
        {
            tmp=cut[i][j];
            max_cut=MAX(max_cut,tmp);
            cut_sq[i][j]=cut_sq[j][i]=tmp*tmp;
            tmp+=atoms->comm.skin;
            cut_sk[i][j]=cut_sk[j][i]=tmp;
            cut_sk_sq[i][j]=cut_sk_sq[j][i]=tmp*tmp;
        }
        rsq_crd[i]=r_crd[i]*r_crd[i];
    }
    c_dim=atoms->c->dim;
    
    
#ifdef CONSERVATIVE
    Memory::alloc(tot_mu_lcl,nelems);
    Memory::alloc(tot_mu,nelems);
    Memory::alloc(tot_nelms_lcl,nelems);
    Memory::alloc(tot_nelms,nelems);
    
    
    for(int i=0;i<c_dim;i++) tot_nelms_lcl[i]=0.0;
    
    elem_type* __elem=atoms->elem->begin();
    type0* c=atoms->c->begin();
    int n=atoms->natms_lcl*c_dim;
    for(int i=0;i<n;i++)
    {
        if(c[i]<0.0) continue;
        ++tot_nelms_lcl[__elem[i]];
    }
    MPI_Allreduce(tot_nelms_lcl,tot_nelms,static_cast<int>(nelems),Vec<type0>::MPI_T,MPI_SUM,world);
#endif
}
/*--------------------------------------------
 
 --------------------------------------------*/
void ForceFieldDMD::post_fin()
{
#ifdef CONSERVATIVE
    Memory::dealloc(tot_mu_lcl);
    Memory::dealloc(tot_mu);
    Memory::dealloc(tot_nelms_lcl);
    Memory::dealloc(tot_nelms);
#endif
}
/*--------------------------------------------
 
 --------------------------------------------*/
type0 ForceFieldDMD::value()
{
    //timer->start(FORCE_TIME_mode);
    __vec_lcl[0]=__vec_lcl[1+__nvoigt__]=0.0;
    __energy_calc();
    type0 en;
    MPI_Allreduce(&__vec_lcl[0],&en,1,Vec<type0>::MPI_T,MPI_SUM,world);
    //timer->stop(FORCE_TIME_mode);
    return en;
}
/*--------------------------------------------
 
 --------------------------------------------*/
type0* ForceFieldDMD::derivative()
{
    reset();
    __force_calc();
    MPI_Allreduce(__vec_lcl,__vec,__nvoigt__+2,Vec<type0>::MPI_T,MPI_SUM,world);
    Algebra::Do<__nvoigt__>::func([this](int i){__vec[i+1]*=-1.0;});
    Algebra::DyadicV_2_MSY(__vec+1,F_H);
    atoms->fe=__vec[0];
    type0 vol_neg=-atoms->vol;
    Algebra::Do<__nvoigt__>::func([this,&vol_neg](int i){__vec[i+1]/=vol_neg;});
    Algebra::DyadicV_2_MSY(__vec+1,atoms->S_fe);
    atoms->s=__vec[1+__nvoigt__];

    impose_dof(f->begin(),f_alpha->begin());
    return __vec;
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<>
type0 ForceFieldDMD::new_value<false>()
{
    __vec_lcl[0]=__vec_lcl[1+__nvoigt__]=0.0;
    __energy_calc();
    type0 en;
    MPI_Allreduce(&__vec_lcl[0],&en,1,Vec<type0>::MPI_T,MPI_SUM,world);
    return en;
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<>
type0 ForceFieldDMD::new_value<true>()
{
    type0 en=new_value<false>();
    type0* c=atoms->c->begin();
    elem_type* __elem=atoms->elem->begin();
    int n=atoms->natms_lcl*c_dim;
    type0 tot_cmu_lcl=0.0,tot_cmu;
    for(int i=0;i<n;i++)
    {
        if(c[i]<0.0) continue;
        tot_cmu_lcl+=c[i]*ext_mu[__elem[i]];
    }
    
    MPI_Allreduce(&tot_cmu_lcl,&tot_cmu,1,Vec<type0>::MPI_T,MPI_SUM,world);
    
    return en-tot_cmu;
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<>
type0* ForceFieldDMD::new_derivative<false>()
{
    reset();
    __force_calc();
    MPI_Allreduce(__vec_lcl,__vec,__nvoigt__+2,Vec<type0>::MPI_T,MPI_SUM,world);
    Algebra::Do<__nvoigt__>::func([this](int i){__vec[i+1]*=-1.0;});
    Algebra::DyadicV_2_MSY(__vec+1,F_H);
    atoms->fe=__vec[0];
    type0 vol_neg=-atoms->vol;
    Algebra::Do<__nvoigt__>::func([this,&vol_neg](int i){__vec[i+1]/=vol_neg;});
    Algebra::DyadicV_2_MSY(__vec+1,atoms->S_fe);
    atoms->s=__vec[1+__nvoigt__];
    impose_dof(f->begin(),f_alpha->begin());
    return __vec;
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<>
type0* ForceFieldDMD::new_derivative<true>()
{
    new_derivative<false>();
    int natms_lcl=atoms->natms_lcl;
    type0* c=atoms->c->begin();
    type0* __mu=f_c->begin();
    type0 cv;
    type0 logcv;
    type0 kBT=atoms->temp*atoms->kB;
    elem_type* __elem=atoms->elem->begin();
    type0* __alpha=atoms->alpha->begin();
    type0 tot_cmu_lcl=0.0,tot_cmu;
    for(int i=0;i<natms_lcl;i++,c+=c_dim,__mu+=c_dim,__elem+=c_dim,__alpha+=c_dim)
    {
        cv=1.0;
        for(int j=0;j<c_dim;j++)
        {
            if(c[j]>=0.0)
            cv-=c[j];
        }
        logcv=Algebra::mod_log(cv);
        for(int j=0;j<c_dim;j++)
        {
            if(c[j]>=0.0)
            {
                
                __mu[j]=-__mu[j]-kBT*(Algebra::mod_log(c[j])-logcv-3.0*log(__alpha[j]))-mu_0[__elem[j]]+ext_mu[__elem[j]];
                tot_cmu_lcl+=c[j]*ext_mu[__elem[j]];
            }
        }
    }
    MPI_Allreduce(&tot_cmu_lcl,&tot_cmu,1,Vec<type0>::MPI_T,MPI_SUM,world);
    __vec[0]-=tot_cmu;
    impose_c_dof(f_c->begin());

    
    
#ifdef CONSERVATIVE
    __mu=f_c->begin();
    __elem=atoms->elem->begin();
    c=atoms->c->begin();
    int n=natms_lcl*c_dim;
    for(int i=0;i<c_dim;i++) tot_mu_lcl[i]=0.0;
    for(int i=0;i<n;i++)
    {
        if(c[i]<0.0) continue;
        tot_mu_lcl[__elem[i]]+=__mu[i];
    }
    
    MPI_Allreduce(tot_mu_lcl,tot_mu,static_cast<int>(nelems),Vec<type0>::MPI_T,MPI_SUM,world);
    for(int i=0;i<c_dim;i++) tot_mu[i]=tot_nelms[i]==0.0 ? 0.0: tot_mu[i]/tot_nelms[i];
    for(int i=0;i<n;i++)
    for(int i=0;i<n;i++)
    {
        if(c[i]<0.0) continue;
        __mu[i]-=tot_mu[__elem[i]];
    }
    
#endif
    
    return __vec;
}
/*--------------------------------------------
 
 --------------------------------------------*/
type0* ForceFieldDMD::derivative_gp()
{
    reset();
    __force_calc_gp();
    MPI_Allreduce(__vec_lcl,__vec,__nvoigt__+3,Vec<type0>::MPI_T,MPI_SUM,world);
    Algebra::Do<__nvoigt__>::func([this](int i){__vec[i+1]*=-1.0;});
    Algebra::DyadicV_2_MSY(__vec+1,F_H);
    atoms->gp=__vec[0];
    atoms->fe=__vec[1+__nvoigt__];
    atoms->pe=__vec[2+__nvoigt__];
    type0 vol_neg=-atoms->vol;
    Algebra::Do<__nvoigt__>::func([this,&vol_neg](int i){__vec[i+1]/=vol_neg;});
    Algebra::DyadicV_2_MSY(__vec+1,atoms->S_fe);
    impose_dof(f->begin(),f_alpha->begin());
    return __vec;
}
/*--------------------------------------------
 
 --------------------------------------------*/
type0 ForceFieldDMD::value_gp()
{
    __vec_lcl[0]=0.0;
    __energy_calc_gp();
    type0 en;
    MPI_Allreduce(&__vec_lcl[0],&en,1,Vec<type0>::MPI_T,MPI_SUM,world);
    return en;
}
#include "dynamic_dmd.h"
/*--------------------------------------------
 
 --------------------------------------------*/
void ForceFieldDMD::calc_thermo()
{
    DynamicDMD __dynamic(atoms,this,false,{},{atoms->x_dof,atoms->alpha_dof,atoms->c_dof},{});
    __dynamic.init();
    calc_ndof();
    derivative();
    
    
    
    type0* mu_sum_lcl=NULL;
    Memory::alloc(mu_sum_lcl,nelems);
    for(size_t i=0;i<nelems;i++) mu_sum_lcl[i]=ave_mu[i]=0.0;
    
    type0* sum_lcl=NULL;
    type0* sum=NULL;
    Memory::alloc(sum_lcl,nelems);
    Memory::alloc(sum,nelems);
    for(size_t i=0;i<nelems;i++) sum_lcl[i]=sum[i]=0.0;
    
    
    type0* mu_vec=f_c->begin();
    type0 const* __c=atoms->c->begin();
    elem_type const* __elem_vec=atoms->elem->begin();
    type0 const* __alpha=atoms->alpha->begin();
    const int n=atoms->natms_lcl;
    type0 kBT=atoms->temp*atoms->kB;
    type0 __cv;
    type0 __logcv;
    for(int i=0;i<n;i++,__c+=c_dim,__elem_vec+=c_dim,mu_vec+=c_dim,__alpha+=c_dim)
    {
        __cv=1.0;
        for(int j=0;j<c_dim;j++)
        {
            if(__c[j]<0.0) continue;
            __cv-=__c[j];
            
        }
        __logcv=Algebra::mod_log(__cv);
        for(int j=0;j<c_dim;j++)
        {
            mu_sum_lcl[__elem_vec[j]]+=mu_vec[j]+kBT*(Algebra::mod_log(__c[j])-__logcv-3.0*log(__alpha[j]))+mu_0[__elem_vec[j]];
            ++sum_lcl[__elem_vec[j]];
        }
    }
    
    MPI_Allreduce(mu_sum_lcl,ave_mu,static_cast<int>(nelems),Vec<type0>::MPI_T,MPI_SUM,world);
    MPI_Allreduce(sum_lcl,sum,static_cast<int>(nelems),Vec<type0>::MPI_T,MPI_SUM,world);
    for(size_t i=0;i<nelems;i++) if(sum[i]!=0.0) ave_mu[i]/=sum[i];
    
    Memory::dealloc(sum_lcl);
    Memory::dealloc(sum);
    Memory::dealloc(mu_sum_lcl);
    
    __dynamic.fin();


}
/*--------------------------------------------
 
 --------------------------------------------*/
void ForceFieldDMD::force_calc_static()
{
    reset();
    __force_calc_static();
    MPI_Allreduce(__vec_lcl,__vec,__nvoigt__+2,Vec<type0>::MPI_T,MPI_SUM,world);
    Algebra::Do<__nvoigt__>::func([this](int i){__vec[i+1]*=-1.0;});
    Algebra::DyadicV_2_MSY(__vec+1,F_H);
    atoms->fe=__vec[0];
    type0 vol_neg=-atoms->vol;
    Algebra::Do<__nvoigt__>::func([this,&vol_neg](int i){__vec[i+1]/=vol_neg;});
    Algebra::DyadicV_2_MSY(__vec+1,atoms->S_fe);
    atoms->s=__vec[1+__nvoigt__];
    
    impose_dof(f->begin(),f_alpha->begin());
}
/*--------------------------------------------
 
 --------------------------------------------*/
void ForceFieldDMD::c_d_calc()
{
    reset_c_d();
    __c_d_calc();
    
    impose_dof(f->begin(),f_alpha->begin());
    norm_sq(f->begin(),err_sq_x,f_alpha->begin(),err_sq_alpha);
}
/*--------------------------------------------
 
 --------------------------------------------*/
type0 ForceFieldDMD::c_dd_norm()
{
    reset_c_d();
    __c_d_calc();
    __J(c_d,f_alpha);
    int n=atoms->natms_lcl*c_dim;
    type0* __f_alpha=f_alpha->begin();
    type0 norm_lcl=0.0;
    for(int i=0;i<n;i++)
        norm_lcl+=__f_alpha[i]*__f_alpha[i];
    type0 norm;
    MPI_Allreduce(&norm_lcl,&norm,1,Vec<type0>::MPI_T,MPI_SUM,world);
    return sqrt(norm);
}
/*--------------------------------------------
 
 --------------------------------------------*/
void ForceFieldDMD::J(Vec<type0>* x,Vec<type0>* Jx)
{
    __J(x,Jx);
}
/*--------------------------------------------
 
 --------------------------------------------*/
type0* ForceFieldDMD::J(Vec<type0>* x_ptr,Vec<type0>* alpha_ptr,Vec<type0>* Jx_ptr,Vec<type0>* Jalpha_ptr)
{
    Algebra::zero<__nvoigt__+2>(__vec_lcl);
    __J(x_ptr,alpha_ptr,Jx_ptr,Jalpha_ptr);
    impose_dof(Jx_ptr->begin(),Jalpha_ptr->begin());
    MPI_Allreduce(__vec_lcl,__vec,__nvoigt__,Vec<type0>::MPI_T,MPI_SUM,world);
    return __vec;
}
/*--------------------------------------------
 
 --------------------------------------------*/
void ForceFieldDMD::prepJ_n_res(Vec<type0>* x_ptr,Vec<type0>* alpha_ptr)
{
    Algebra::zero<__nvoigt__+2>(__vec_lcl);
    __prepJ_n_res(x_ptr,alpha_ptr);
    MPI_Allreduce(__vec_lcl,__vec,__nvoigt__+2,Vec<type0>::MPI_T,MPI_SUM,world);
    Algebra::Do<__nvoigt__>::func([this](int i){__vec[i+1]*=-1.0;});
    Algebra::DyadicV_2_MSY(__vec+1,F_H);
    atoms->fe=__vec[0];
    type0 vol_neg=-atoms->vol;
    Algebra::Do<__nvoigt__>::func([this,&vol_neg](int i){__vec[i+1]/=vol_neg;});
    Algebra::DyadicV_2_MSY(__vec+1,atoms->S_fe);
    atoms->s=__vec[1+__nvoigt__];
    
    impose_dof(x_ptr->begin(),alpha_ptr->begin());
    norm_sq(x_ptr->begin(),err_sq_x,alpha_ptr->begin(),err_sq_alpha);
    
}





