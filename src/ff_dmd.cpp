#include "ff_dmd.h"
#include "atoms_dmd.h"
#include "neighbor_dmd.h"
#include "xmath.h"
#include "memory.h"
using namespace MAPP_NS;
/*--------------------------------------------
 
 --------------------------------------------*/
ForceFieldDMD::ForceFieldDMD(AtomsDMD*& __atoms):
ForceField(__atoms),
atoms(__atoms),
rsq_crd(NULL),
r_crd(NULL),
f_alpha(new Vec<type0>(__atoms,__atoms->c->dim))
{
    Memory::alloc(cut_sk,nelems,nelems);
    if(nelems)
    {
        rsq_crd=new type0[nelems];
        r_crd=new type0[nelems];
    }
    neighbor=neighbor_dmd=new NeighborDMD(__atoms,cut_sk,rsq_crd);
    
}
/*--------------------------------------------
 
 --------------------------------------------*/
ForceFieldDMD::~ForceFieldDMD()
{
    Memory::dealloc(cut_sk);
    delete f_alpha;
    delete neighbor_dmd;
    delete [] rsq_crd;
    delete [] r_crd;
}
/*--------------------------------------------
 
 --------------------------------------------*/
void ForceFieldDMD::setup()
{
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
    x_dim=atoms->x->dim;
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
    if(!atoms->dof) return;
    type0* fvec=f->begin();
    bool* dof=atoms->dof->begin();
    const int n=atoms->natms*__dim__;
    for(int i=0;i<n;i++) fvec[i]*=dof[i];
}
/*--------------------------------------------
 
 --------------------------------------------*/
void ForceFieldDMD::derivative_timer(type0(*&S)[__dim__])
{
    force_calc();
    type0* fvec=f->begin();
    type0* xvec=atoms->x->begin();
    if(!atoms->dof)
    {
        const int natms=atoms->natms;
        for(int i=0;i<natms;i++,fvec+=x_dim,xvec+=x_dim)
            Algebra::DyadicV<__dim__>(xvec,fvec,nrgy_strss_lcl+1);
    }
    else
    {
        bool* dof=atoms->dof->begin();
        const int natms=atoms->natms;
        for(int i=0;i<natms;i++,fvec+=x_dim,xvec+=x_dim)
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
void ForceFieldDMD::reset()
{
    type0* __f=f->begin();
    const int n=atoms->natms*x_dim;
    for(int i=0;i<n;i++) __f[i]=0.0;
    __f=f_alpha->begin();
    const int m=atoms->natms*c_dim;
    for(int i=0;i<m;i++) __f[i]=0.0;
    for(int i=0;i<__nvoigt__+1;i++) nrgy_strss_lcl[i]=0.0;
}


