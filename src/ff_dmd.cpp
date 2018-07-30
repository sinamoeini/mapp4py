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
r_crd(NULL)
{
    Memory::alloc(cut_sk,nelems,nelems);
    Memory::alloc(rsq_crd,nelems);
    Memory::alloc(r_crd,nelems);
    neighbor=new NeighborDMD(__atoms,cut_sk,rsq_crd);
    f=new Vec<type0>(atoms,__dim__,"f");
    f_alpha=new DMDVec<type0>(atoms,0.0,"f_alpha");
    c_d=new DMDVec<type0>(atoms,0.0,"c_d");
}
/*--------------------------------------------
 
 --------------------------------------------*/
ForceFieldDMD::~ForceFieldDMD()
{
    delete c_d;
    delete f_alpha;
    delete f;
    delete neighbor;
    Memory::dealloc(r_crd);
    Memory::dealloc(rsq_crd);
    Memory::dealloc(cut_sk);
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
    dof_empty=atoms->dof->is_empty();
    dof_alpha_empty=atoms->dof_alpha->is_empty();
    dof_c_empty=atoms->dof_c->is_empty();
    
    type0 tmp;
    max_cut=0.0;
    for(size_t i=0;i<nelems;i++)
    {
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
}
/*--------------------------------------------
 
 --------------------------------------------*/
void ForceFieldDMD::post_fin()
{
}
/*--------------------------------------------
 
 --------------------------------------------*/
type0 ForceFieldDMD::value_timer()
{
    //timer->start(FORCE_TIME_mode);
    __vec_lcl[0]=__vec_lcl[1+__nvoigt__]=0.0;
    energy_calc();
    type0 en;
    MPI_Allreduce(&__vec_lcl[0],&en,1,Vec<type0>::MPI_T,MPI_SUM,world);
    //timer->stop(FORCE_TIME_mode);
    return en;
}
/*--------------------------------------------
 
 --------------------------------------------*/
type0* ForceFieldDMD::derivative_timer()
{
    reset();
    force_calc();
    MPI_Allreduce(__vec_lcl,__vec,__nvoigt__+2,Vec<type0>::MPI_T,MPI_SUM,world);
    Algebra::Do<__nvoigt__>::func([this](int i){__vec[i+1]/=atoms->vol;});
    
    if(!dof_empty)
    {
        type0* fvec=f->begin();
        bool* dof=atoms->dof->begin();
        const int n=atoms->natms_lcl*__dim__;
        for(int i=0;i<n;i++) fvec[i]=dof[i] ? fvec[i]:0.0;
    }
    
    if(!dof_alpha_empty)
    {
        type0* fvec=f_alpha->begin();
        bool* dof=atoms->dof_alpha->begin();
        const int n=atoms->natms_lcl*c_dim;
        for(int i=0;i<n;i++) fvec[i]=dof[i] ? fvec[i]:0.0;
    }
    
    atoms->fe=__vec[0];
    Algebra::DyadicV_2_MSY(__vec+1,atoms->S_fe);
    atoms->s=__vec[1+__nvoigt__];
    return __vec;
}
/*--------------------------------------------
 this does not sound right hs to be check later
 --------------------------------------------*/
//void ForceFieldDMD::derivative_timer(type0(*&S)[__dim__])
//{
//    reset();
//    force_calc();
//    type0* fvec=f->begin();
//    type0* xvec=atoms->x->begin();
//    if(dof_empty)
//    {
//        const int natms_lcl=atoms->natms_lcl;
//        for(int i=0;i<natms_lcl;i++,fvec+=__dim__,xvec+=__dim__)
//            Algebra::DyadicV<__dim__>(xvec,fvec,__vec_lcl+1);
//    }
//    else
//    {
//        bool* dof=atoms->dof->begin();
//        const int natms_lcl=atoms->natms_lcl;
//        for(int i=0;i<natms_lcl;i++,fvec+=__dim__,xvec+=__dim__,dof+=__dim__)
//        {
//            Algebra::Do<__dim__>::func([&dof,&fvec](int i){fvec[i]=dof[i] ? fvec[i]:0.0;});
//            Algebra::DyadicV<__dim__>(xvec,fvec,__vec_lcl+1);
//        }
//    }
//    
//    
//    MPI_Allreduce(__vec_lcl,__vec,__nvoigt__+2,Vec<type0>::MPI_T,MPI_SUM,world);
//    Algebra::Do<__nvoigt__>::func([this](int i){__vec[i+1]*=-1.0;});
//    Algebra::NONAME_DyadicV_mul_MLT(__vec+1,atoms->B,S);
//    const type0 vol=atoms->vol;
//    Algebra::Do<__nvoigt__>::func([this,&vol](int i){__vec[i+1]/=-vol;});
//    
//    if(!dof_alpha_empty)
//    {
//        type0* fvec=f_alpha->begin();
//        bool* dof=atoms->dof_alpha->begin();
//        const int n=atoms->natms_lcl*c_dim;
//        for(int i=0;i<n;i++) fvec[i]=dof[i] ? fvec[i]:0.0;
//    }
//    
//    atoms->fe=__vec[0];
//    Algebra::DyadicV_2_MSY(__vec+1,atoms->S_fe);
//    atoms->s=__vec[1+__nvoigt__];
//}
/*--------------------------------------------
 this does not sound right hs to be check later
 --------------------------------------------*/
//void ForceFieldDMD::derivative_timer(bool affine,type0(*&S)[__dim__])
//{
//    reset();
//    force_calc();
//    if(!affine)
//    {
//        type0* fvec=f->begin();
//        type0* xvec=atoms->x->begin();
//        if(dof_empty)
//        {
//            const int natms_lcl=atoms->natms_lcl;
//            for(int i=0;i<natms_lcl;i++,fvec+=__dim__,xvec+=__dim__)
//                Algebra::DyadicV<__dim__>(xvec,fvec,__vec_lcl+1);
//        }
//        else
//        {
//            bool* dof=atoms->dof->begin();
//            const int natms_lcl=atoms->natms_lcl;
//            for(int i=0;i<natms_lcl;i++,fvec+=__dim__,xvec+=__dim__,dof+=__dim__)
//            {
//                Algebra::Do<__dim__>::func([&dof,&fvec](int i){fvec[i]=dof[i] ? fvec[i]:0.0;});
//                Algebra::DyadicV<__dim__>(xvec,fvec,__vec_lcl+1);
//            }
//        }
//    }
//    
//    MPI_Allreduce(__vec_lcl,__vec,__nvoigt__+2,Vec<type0>::MPI_T,MPI_SUM,world);
//    Algebra::Do<__nvoigt__>::func([this](int i){__vec[i+1]*=-1.0;});
//    Algebra::NONAME_DyadicV_mul_MLT(__vec+1,atoms->B,S);
//    const type0 vol=atoms->vol;
//    Algebra::Do<__nvoigt__>::func([this,&vol](int i){__vec[i+1]/=-vol;});
//    
//    if(!dof_alpha_empty)
//    {
//        type0* fvec=f_alpha->begin();
//        bool* dof=atoms->dof_alpha->begin();
//        const int n=atoms->natms_lcl*c_dim;
//        for(int i=0;i<n;i++) fvec[i]=dof[i] ? fvec[i]:0.0;
//    }
//    
//    atoms->fe=__vec[0];
//    Algebra::DyadicV_2_MSY(__vec+1,atoms->S_fe);
//    atoms->s=__vec[1+__nvoigt__];
//}
/*--------------------------------------------
 
 --------------------------------------------*/
void ForceFieldDMD::force_calc_static_timer()
{
    reset();
    force_calc_static();
    MPI_Allreduce(__vec_lcl,__vec,__nvoigt__+2,Vec<type0>::MPI_T,MPI_SUM,world);
    Algebra::Do<__nvoigt__>::func([this](int i){__vec[i+1]/=atoms->vol;});
    
    atoms->fe=__vec[0];
    Algebra::DyadicV_2_MSY(__vec+1,atoms->S_fe);
    atoms->s=__vec[1+__nvoigt__];
}
/*--------------------------------------------
 
 --------------------------------------------*/
void ForceFieldDMD::c_d_calc_timer()
{
    reset_c_d();
    c_d_calc();
}
/*--------------------------------------------
 
 --------------------------------------------*/
type0 ForceFieldDMD::c_dd_norm_timer()
{
    reset_c_d();
    c_d_calc();
    J(c_d,f_alpha);
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
void ForceFieldDMD::J_timer(Vec<type0>* x,Vec<type0>* Jx)
{
    J(x,Jx);
}
/*--------------------------------------------
 
 --------------------------------------------*/
void ForceFieldDMD::J_timer(VecTens<type0,2>& x,VecTens<type0,2>& Jx)
{
    Algebra::zero<__nvoigt__+2>(__vec_lcl);
    J(x,Jx);
    if(x.chng_box)
    {
        MPI_Allreduce(__vec_lcl,__vec,__nvoigt__,Vec<type0>::MPI_T,MPI_SUM,world);
        Algebra::DyadicV_2_MLT(__vec,Jx.A);
    }
}
/*--------------------------------------------
 
 --------------------------------------------*/
type0 ForceFieldDMD::prep_timer(VecTens<type0,2>& f)
{
    Algebra::zero<__nvoigt__+2>(__vec_lcl);
    type0 err_sq=prep(f);
    MPI_Allreduce(__vec_lcl,__vec,__nvoigt__+1,Vec<type0>::MPI_T,MPI_SUM,world);
    Algebra::Do<__nvoigt__>::func([this](int i){__vec[i+1]/=atoms->vol;});
    
    atoms->fe=__vec[0];
    Algebra::DyadicV_2_MSY(__vec+1,atoms->S_fe);
    atoms->s=__vec[1+__nvoigt__];
    return sqrt(err_sq);
}
/*--------------------------------------------
 
 --------------------------------------------*/
type0 ForceFieldDMD::prep_timer(VecTens<type0,2>& f,type0 (&S)[__dim__][__dim__])
{
    Algebra::zero<__nvoigt__+2>(__vec_lcl);
    type0 err_sq=prep(f);
    MPI_Allreduce(__vec_lcl,__vec,__nvoigt__+1,Vec<type0>::MPI_T,MPI_SUM,world);
    type0 vol=atoms->vol;
    Algebra::DoLT<__dim__>::func([&S,&err_sq,&f,this,&vol](int i,int j)
    {
        
        int k=j*(__dim__-1)-j*(j-1)/2+i+1;
        
        if(!std::isnan(S[i][j]))
        {
            f.A[i][j]=S[i][j]*vol-__vec[k];
            err_sq+=f.A[i][j]*f.A[i][j];
        }
        else
            f.A[i][j]=0.0;
        
        __vec[k]/=vol;
        
    });
    
    atoms->fe=__vec[0];
    Algebra::DyadicV_2_MSY(__vec+1,atoms->S_fe);
    atoms->s=__vec[1+__nvoigt__];
    return sqrt(err_sq);
}




