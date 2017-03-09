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
    f_alpha=new Vec<type0>(atoms,atoms->c_dim,"f_alpha");
}
/*--------------------------------------------
 
 --------------------------------------------*/
ForceFieldDMD::~ForceFieldDMD()
{
    delete f_alpha;
    delete f;
    delete neighbor;
    Memory::dealloc(r_crd);
    Memory::dealloc(rsq_crd);
    Memory::dealloc(cut_sk);
}
/*--------------------------------------------
 
 --------------------------------------------*/
void ForceFieldDMD::pre_init()
{
    dof_empty=atoms->dof->is_empty();
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
    nrgy_strss_lcl[0]=0.0;
    energy_calc();
    type0 en;
    MPI_Allreduce(&nrgy_strss_lcl[0],&en,1,Vec<type0>::MPI_T,MPI_SUM,world);
    //timer->stop(FORCE_TIME_mode);
    return en;
}
/*--------------------------------------------
 
 --------------------------------------------*/
void ForceFieldDMD::derivative_timer()
{
    force_calc();
    MPI_Allreduce(nrgy_strss_lcl,nrgy_strss,__nvoigt__+1,Vec<type0>::MPI_T,MPI_SUM,world);
    Algebra::Do<__nvoigt__>::func([this](int i){nrgy_strss[i+1]/=atoms->vol;});
    if(dof_empty) return;
    type0* fvec=f->begin();
    bool* dof=atoms->dof->begin();
    const int n=atoms->natms_lcl*__dim__;
    for(int i=0;i<n;i++) fvec[i]*=dof[i];
}
/*--------------------------------------------
 
 --------------------------------------------*/
void ForceFieldDMD::derivative_timer(type0(*&S)[__dim__])
{
    force_calc();
    type0* fvec=f->begin();
    type0* xvec=atoms->x->begin();
    if(dof_empty)
    {
        const int natms_lcl=atoms->natms_lcl;
        for(int i=0;i<natms_lcl;i++,fvec+=__dim__,xvec+=__dim__)
            Algebra::DyadicV<__dim__>(xvec,fvec,nrgy_strss_lcl+1);
    }
    else
    {
        bool* dof=atoms->dof->begin();
        const int natms_lcl=atoms->natms_lcl;
        for(int i=0;i<natms_lcl;i++,fvec+=__dim__,xvec+=__dim__)
        {
            Algebra::Do<__dim__>::func([&dof,&fvec](int i){fvec[i]*=dof[i];});
            Algebra::DyadicV<__dim__>(xvec,fvec,nrgy_strss_lcl+1);
        }
    }
    
    
    MPI_Allreduce(nrgy_strss_lcl,nrgy_strss,__nvoigt__+1,Vec<type0>::MPI_T,MPI_SUM,world);
    Algebra::Do<__nvoigt__>::func([this](int i){nrgy_strss[i+1]*=-1.0;});
    Algebra::NONAME_DyadicV_mul_MLT(nrgy_strss+1,atoms->B,S);
    Algebra::Do<__nvoigt__>::func([this](int i){nrgy_strss[i+1]*=-1.0/atoms->vol;});
    
}
/*--------------------------------------------
 
 --------------------------------------------*/
void ForceFieldDMD::force_calc_static_timer()
{
    force_calc_static();
    MPI_Allreduce(nrgy_strss_lcl,nrgy_strss,__nvoigt__+1,Vec<type0>::MPI_T,MPI_SUM,world);
    Algebra::Do<__nvoigt__>::func([this](int i){nrgy_strss[i+1]/=atoms->vol;});
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
    for(int i=0;i<__nvoigt__+1;i++) nrgy_strss_lcl[i]=0.0;
}


