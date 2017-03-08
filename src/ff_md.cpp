#include "ff_md.h"
#include "timer.h"
#include "neighbor_md.h"
#include "atoms_md.h"
#include "xmath.h"
using namespace MAPP_NS;
/*--------------------------------------------
 
 --------------------------------------------*/
ForceFieldMD::ForceFieldMD(AtomsMD* __atoms):
ForceField(__atoms),
atoms(__atoms),
elem(__atoms->elem)
{
    neighbor=new NeighborMD(__atoms,cut_sk_sq);
}
/*--------------------------------------------
 
 --------------------------------------------*/
ForceFieldMD::~ForceFieldMD()
{
    delete neighbor;
}
/*--------------------------------------------
 
 --------------------------------------------*/
void ForceFieldMD::setup()
{
    dof_empty=atoms->dof->is_empty();
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
void ForceFieldMD::force_calc_timer()
{
    force_calc();
    MPI_Allreduce(nrgy_strss_lcl,nrgy_strss,__nvoigt__+1,Vec<type0>::MPI_T,MPI_SUM,world);
    const type0 vol=atoms->vol;
    Algebra::Do<__nvoigt__>::func([this,&vol](int i){nrgy_strss[i+1]/=vol;});
}
/*--------------------------------------------
 
 --------------------------------------------*/
type0 ForceFieldMD::energy_calc_timer()
{
    energy_calc();
    MPI_Allreduce(nrgy_strss_lcl,nrgy_strss,1,Vec<type0>::MPI_T,MPI_SUM,world);
    return nrgy_strss[0];
}
/*--------------------------------------------
 
 --------------------------------------------*/
type0 ForceFieldMD::value_timer()
{
    nrgy_strss_lcl[0]=0.0;
    energy_calc();
    type0 en;
    MPI_Allreduce(&nrgy_strss_lcl[0],&en,1,Vec<type0>::MPI_T,MPI_SUM,world);
    return en;
}
/*--------------------------------------------
 
 --------------------------------------------*/
void ForceFieldMD::derivative_timer()
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
void ForceFieldMD::derivative_timer(type0(*&S)[__dim__])
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
    const type0 vol=atoms->vol;
    Algebra::Do<__nvoigt__>::func([this,&vol](int i){nrgy_strss[i+1]*=-1.0/vol;});

}










