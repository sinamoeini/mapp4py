#include "ff_md.h"
#include "neighbor_md.h"
#include "atoms_md.h"
#include "xmath.h"
#include "memory.h"
using namespace MAPP_NS;
/*--------------------------------------------
 
 --------------------------------------------*/
ForceFieldMD::ForceFieldMD(AtomsMD* __atoms):
ForceField(__atoms),
f_alloc(false),
atoms(__atoms)
{
    neighbor=new NeighborMD(__atoms,cut_sk_sq);
    vec* __f=atoms->find_vec("f");
    if(__f) f=reinterpret_cast<Vec<type0>*>(__f);
    else
    {
        f_alloc=true;
        f=new Vec<type0>(atoms,__dim__,"f");
    }
}
/*--------------------------------------------
 
 --------------------------------------------*/
ForceFieldMD::~ForceFieldMD()
{
    if(f_alloc) delete f;
    delete neighbor;
}
/*--------------------------------------------
 
 --------------------------------------------*/
ForceFieldZero::ForceFieldZero(class AtomsMD* atoms,type0**&& __cut):
ForceFieldMD(atoms)
{
    for(size_t i=0;i<nelems;i++)
    for(size_t j=0;j<i+1;j++)
    {
        cut[i][j]=cut[j][i]=__cut[i][j];
        cut_sq[i][j]=cut_sq[j][i]=__cut[i][j]*__cut[i][j];
    }
    Memory::dealloc(__cut);
}
/*--------------------------------------------
 
 --------------------------------------------*/
ForceFieldZero::~ForceFieldZero()
{
}
/*--------------------------------------------
 
 --------------------------------------------*/
void ForceFieldMD::calc_ndof()
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
}
/*--------------------------------------------
 
 --------------------------------------------*/
void ForceFieldMD::reset()
{
    type0* __f=f->begin();
    const int n=f->vec_sz*__dim__;
    for(int i=0;i<n;i++) __f[i]=0.0;
    for(int i=0;i<__nvoigt__+1;i++) __vec_lcl[i]=0.0;
}
/*--------------------------------------------
 
 --------------------------------------------*/
void ForceFieldMD::pre_init()
{
    dof_empty=atoms->x_dof->is_empty();
    type0 tmp;
    max_cut=0.0;
    for(size_t i=0;i<nelems;i++)
        for(size_t j=0;j<i+1;j++)
        {
            tmp=cut[i][j];
            max_cut=MAX(max_cut,tmp);
            cut_sq[i][j]=cut_sq[j][i]=tmp*tmp;
            tmp+=atoms->comm.skin;
            cut_sk_sq[i][j]=cut_sk_sq[j][i]=tmp*tmp;
        }
}
/*--------------------------------------------
 
 --------------------------------------------*/
void ForceFieldMD::post_fin()
{
}
/*--------------------------------------------
 
 --------------------------------------------*/
void ForceFieldMD::pre_xchng_energy_timer(GCMC* gcmc)
{
    pre_xchng_energy(gcmc);
}
/*--------------------------------------------
 
 --------------------------------------------*/
type0 ForceFieldMD::xchng_energy_timer(GCMC* gcmc)
{
    return xchng_energy(gcmc);
}
/*--------------------------------------------
 
 --------------------------------------------*/
void ForceFieldMD::post_xchng_energy_timer(GCMC* gcmc)
{
    post_xchng_energy(gcmc);
}
/*--------------------------------------------
 
 --------------------------------------------*/
void ForceFieldMD::force_calc()
{
    reset();
    __force_calc();
    MPI_Allreduce(__vec_lcl,__vec,__nvoigt__+1,Vec<type0>::MPI_T,MPI_SUM,world);
    Algebra::Do<__nvoigt__>::func([this](int i){__vec[i+1]*=-1.0;});
    Algebra::DyadicV_2_MSY(__vec+1,F_H);
    atoms->pe=__vec[0];
    type0 vol_neg=-atoms->vol;
    Algebra::Do<__nvoigt__>::func([this,&vol_neg](int i){__vec[i+1]/=vol_neg;});
    Algebra::DyadicV_2_MSY(__vec+1,atoms->S_pe);
    
    if(!dof_empty)
    {
        type0* fvec=f->begin();
        bool* dof=atoms->x_dof->begin();
        const int n=atoms->natms_lcl*__dim__;
        for(int i=0;i<n;i++) fvec[i]=dof[i] ? fvec[i]:0.0;
    }
}
/*--------------------------------------------
 
 --------------------------------------------*/
type0 ForceFieldMD::value()
{
    __vec_lcl[0]=0.0;
    __energy_calc();
    type0 en;
    MPI_Allreduce(&__vec_lcl[0],&en,1,Vec<type0>::MPI_T,MPI_SUM,world);
    return en;
}
/*--------------------------------------------
 
 --------------------------------------------*/
type0* ForceFieldMD::derivative()
{
    reset();
    __force_calc();
    MPI_Allreduce(__vec_lcl,__vec,__nvoigt__+1,Vec<type0>::MPI_T,MPI_SUM,world);
    Algebra::Do<__nvoigt__>::func([this](int i){__vec[i+1]*=-1.0;});
    Algebra::DyadicV_2_MSY(__vec+1,F_H);
    atoms->pe=__vec[0];
    type0 vol_neg=-atoms->vol;
    Algebra::Do<__nvoigt__>::func([this,&vol_neg](int i){__vec[i+1]/=vol_neg;});
    Algebra::DyadicV_2_MSY(__vec+1,atoms->S_pe);
    
    if(!dof_empty)
    {
        type0* fvec=f->begin();
        bool* dof=atoms->x_dof->begin();
        const int n=atoms->natms_lcl*__dim__;
        for(int i=0;i<n;i++) fvec[i]=dof[i] ? fvec[i]:0.0;
    }
    return __vec;
}
/*--------------------------------------------
 this does not sound right hs to be check later
 --------------------------------------------*/
//void ForceFieldMD::derivative_timer(type0(*&S)[__dim__])
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
//    MPI_Allreduce(__vec_lcl,__vec,__nvoigt__+1,Vec<type0>::MPI_T,MPI_SUM,world);
//    Algebra::Do<__nvoigt__>::func([this](int i){__vec[i+1]*=-1.0;});
//    Algebra::NONAME_DyadicV_mul_MLT(__vec+1,atoms->B,S);
//    const type0 vol=atoms->vol;
//    Algebra::Do<__nvoigt__>::func([this,&vol](int i){__vec[i+1]/=-vol;});
//    atoms->pe=__vec[0];
//    Algebra::DyadicV_2_MSY(__vec+1,atoms->S_pe);
//}
/*--------------------------------------------
 this does not sound right hs to be check later
 --------------------------------------------*/
//void ForceFieldMD::derivative_timer(bool affine,type0(*&S)[__dim__])
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
//    MPI_Allreduce(__vec_lcl,__vec,__nvoigt__+1,Vec<type0>::MPI_T,MPI_SUM,world);
//    Algebra::Do<__nvoigt__>::func([this](int i){__vec[i+1]*=-1.0;});
//    Algebra::NONAME_DyadicV_mul_MLT(__vec+1,atoms->B,S);
//    const type0 vol=atoms->vol;
//    Algebra::Do<__nvoigt__>::func([this,&vol](int i){__vec[i+1]/=-vol;});
//    atoms->pe=__vec[0];
//    Algebra::DyadicV_2_MSY(__vec+1,atoms->S_pe);
//}









