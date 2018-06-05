#include "min_cg_fit.h"
#include "potfit.h"
#include <stdlib.h>
using namespace MAPP_NS;
/*--------------------------------------------
 constructor
 --------------------------------------------*/
MinCGFit::MinCGFit(type0 __e_tol,
bool(&__H_dof)[__dim__][__dim__],bool __affine,type0 __max_dx,LineSearch* __ls,vec* __ext_vec_0,vec* __ext_vec_1):
Min(__e_tol,__H_dof,__affine,__max_dx,__ls),
atoms(NULL),
ff(NULL),
xprt(NULL),
ext_vec_0(__ext_vec_0),
ext_vec_1(__ext_vec_1)
{
}
/*--------------------------------------------
 destructor
 --------------------------------------------*/
MinCGFit::~MinCGFit()
{
    ext_vec_0=ext_vec_1=NULL;
    atoms=NULL;
    ff=NULL;
}
/*--------------------------------------------
 
 --------------------------------------------*/
void MinCGFit::force_calc()
{
    ff->derivative_timer();
    /*
    if(chng_box)
    {
        Algebra::MLT_mul_MLT(atoms->S_pe,atoms->B,f.A);
        type0 neg_v=-atoms->vol;
        Algebra::DoLT<__dim__>::func([this,&neg_v](int i,int j)
        {f.A[i][j]=H_dof[i][j] ? f.A[i][j]*neg_v:0.0;});
    }*/
    
    
    if(chng_box)
    {
        type0 (&S_pe)[__dim__][__dim__]=atoms->S_pe;
        type0 neg_v=-atoms->vol;
        Algebra::DoLT<__dim__>::func([this,&S_pe,&neg_v](int i,int j)
        {S_tmp[i][j]=H_dof[i][j] ? S_pe[i][j]*neg_v:0.0;});
        Algebra::MLT_mul_MLT(S_tmp,atoms->B,f.A);
    }
    
}
/*--------------------------------------------
 
 --------------------------------------------*/
void MinCGFit::prep()
{
    if(chng_box)
    {
        const int natms_lcl=atoms->natms_lcl;
        type0 MLT[__dim__][__dim__];
        Algebra::MLT_mul_MLT(atoms->B,h.A,MLT);
        type0* xvec=x0.vecs[0]->begin();
        type0* hvec=h.vecs[0]->begin();
        type0* x_dvec=x_d.vecs[0]->begin();
        Algebra::DoLT<__dim__>::func([this](int i,int j)
        {
            x_d.A[i][j]=h.A[i][j];
        });
        
        if(affine)
        {
            for(int iatm=0;iatm<natms_lcl;iatm++,xvec+=__dim__,x_dvec+=__dim__)
                Algebra::V_mul_MLT(xvec,MLT,x_dvec);
        }
        else
        {
            for(int iatm=0;iatm<natms_lcl;iatm++,xvec+=__dim__,x_dvec+=__dim__,hvec+=__dim__)
            {
                Algebra::V_mul_MLT(xvec,MLT,x_dvec);
                Algebra::V_add<__dim__>(hvec,x_dvec);
            }
        }
    }
    else
        x_d=h;
    
}
/*--------------------------------------------
 init before a run
 --------------------------------------------*/
void MinCGFit::init()
{
    x.~VecTens();
    new (&x) VecTens<type0,1>(atoms,chng_box,atoms->H,atoms->x);
    f.~VecTens();
    new (&f) VecTens<type0,1>(atoms,chng_box,ff->f);
    h.~VecTens();
    new (&h) VecTens<type0,1>(atoms,chng_box,__dim__);
    x0.~VecTens();
    new (&x0) VecTens<type0,1>(atoms,chng_box,__dim__);
    x_d.~VecTens();
    new (&x_d) VecTens<type0,1>(atoms,chng_box,__dim__);
    f0.~VecTens();
    new (&f0) VecTens<type0,1>(atoms,chng_box,__dim__);
    
    dynamic=new DynamicMD(atoms,ff,chng_box,{},{atoms->dof,h.vecs[0],x0.vecs[0],x_d.vecs[0],f0.vecs[0],ext_vec_0,ext_vec_1},{});
    dynamic->init();
    
    if(xprt)
    {
        try
        {
            xprt->atoms=atoms;
            xprt->init();
        }
        catch(std::string& err_msg)
        {
            fin();
            throw err_msg;
        }
    }
}
/*--------------------------------------------
 finishing minimization
 --------------------------------------------*/
void MinCGFit::fin()
{
    if(xprt)
    {
        xprt->fin();
        xprt->atoms=NULL;
    }
    
    dynamic->fin();
    delete dynamic;
    dynamic=NULL;
    
    f0.~VecTens();
    x_d.~VecTens();
    x0.~VecTens();
    h.~VecTens();
    f.~VecTens();
    x.~VecTens();
}
/*--------------------------------------------
 min
 --------------------------------------------*/
void MinCGFit::run(int nsteps)
{
    if(dynamic_cast<LineSearchGoldenSection*>(ls))
        return run(dynamic_cast<LineSearchGoldenSection*>(ls),nsteps);
    
    if(dynamic_cast<LineSearchBrent*>(ls))
        return run(dynamic_cast<LineSearchBrent*>(ls),nsteps);
    
    if(dynamic_cast<LineSearchBackTrack*>(ls))
        return run(dynamic_cast<LineSearchBackTrack*>(ls),nsteps);
}
/*--------------------------------------------
 
 --------------------------------------------*/
type0 MinCGFit::F(type0 alpha)
{
    x=x0+alpha*x_d;
    if(chng_box)
        atoms->update_H();    
    
    dynamic->update(atoms->x);
    return ff->value_timer();
}
/*--------------------------------------------
 inner product of f and h
 --------------------------------------------*/
type0 MinCGFit::dF(type0 alpha,type0& drev)
{
    x=x0+alpha*x_d;
    if(chng_box)
        atoms->update_H();
    
    dynamic->update(atoms->x);
    force_calc();
    
    drev=-(f*h);
    return atoms->pe;
}
/*--------------------------------------------
 find maximum h
 lets find the next sensible number 
 
 x=x_0+h*alpha
 (x-x0)/alpha=sqrt(eps)/alpha

 --------------------------------------------*/
void MinCGFit::ls_prep(type0& dfa,type0& h_norm,type0& max_a)
{
    
    h_norm=h*h;
    
    dfa=-f_h;
    
    if(h_norm==0.0)
    {
        max_a=0.0;
        dfa=0.0;
        return;
    }
    
    if(dfa>=0.0)
    {
        max_a=0.0;
        dfa=1.0;
        return;
    }
    
    h_norm=sqrt(h_norm);
    
    type0 max_x_d_lcl=0.0;
    type0 max_x_d;
    type0* x_dvec=x_d.vecs[0]->begin();
    const int n=atoms->natms_lcl*__dim__;
    for(int i=0;i<n;i++)
        max_x_d_lcl=MAX(max_x_d_lcl,fabs(x_dvec[i]));
    
    MPI_Allreduce(&max_x_d_lcl,&max_x_d,1,Vec<type0>::MPI_T,MPI_MAX,atoms->world);
    max_a=fabs(max_dx/max_x_d);
    
}
/*--------------------------------------------
 reset to initial position
 --------------------------------------------*/
void MinCGFit::F_reset()
{
    x=x0;
    if(chng_box) atoms->update_H();
    dynamic->update(atoms->x);
}

