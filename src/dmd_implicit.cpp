#include "dmd_implicit.h"
#include "atoms_dmd.h"
#include "ff_dmd.h"
#include "dynamic_dmd.h"
#include "gmres.h"
using namespace MAPP_NS;
/*--------------------------------------------
 
 --------------------------------------------*/
DMDImplicit::DMDImplicit():
DMD(),
max_niters_nonlin(10),
max_niters_lin(10),
y_0(NULL),
c_0(NULL),
del_c(NULL),
F(NULL),
a(NULL),
gmres(NULL)
{
}
/*--------------------------------------------
 
 --------------------------------------------*/
DMDImplicit::~DMDImplicit()
{
    delete [] y_0;
    y_0=c_0=del_c=F=a=NULL;
}
/*--------------------------------------------
 
 --------------------------------------------*/
void DMDImplicit::init_static()
{
    DMD::init_static();
    y_0=new type0[5*ncs];
    c_0=y_0+ncs;
    del_c=y_0+2*ncs;
    F=y_0+3*ncs;
    a=y_0+4*ncs;
    gmres=new GMRES<type0,ForceFieldDMD>(atoms,max_niters_lin,c_dim,*ff);
}
/*--------------------------------------------
 
 --------------------------------------------*/
void DMDImplicit::fin_static()
{
    delete gmres;
    gmres=NULL;
    delete [] y_0;
    y_0=c_0=del_c=F=a=NULL;
    DMD::fin_static();
}
/*--------------------------------------------
 
 --------------------------------------------*/
inline type0 DMDImplicit::update_c()
{
    type0 tmp;
    type0 r,r_lcl=1.0;
    for(int i=0;i<ncs;i++)
    {
        if(c[i]>=0.0)
        {
            tmp=c[i]+r_lcl*del_c[i];
            if(tmp>1.0)
            {
                r_lcl=(1.0-c[i])/del_c[i];
                while(c[i]+r_lcl*del_c[i]>1.0)
                    r_lcl=nextafter(r_lcl,0.0);
            }
            if(tmp<0.0)
            {
                r_lcl=-c[i]/del_c[i];
                while(c[i]+r_lcl*del_c[i]<0.0)
                    r_lcl=nextafter(r_lcl,0.0);
            }
        }
        else
            del_c[i]=0.0;
    }
    
    
    MPI_Allreduce(&r_lcl,&r,1,MPI_TYPE0,MPI_MIN,atoms->world);
    
    volatile type0 c0;
    for(int i=0;i<ncs;i++)
    {
        c0=c[i]+r*del_c[i];
        --++c0;
        c[i]=c0;
    }
    return r;
}
/*--------------------------------------------
 solve the implicit equation
 
 nc_dofs
 err_fac
 --------------------------------------------*/
bool DMDImplicit::nonlin()
{
    
    type0 res_tol=0.005*a_tol*sqrt(nc_dofs)/err_fac;
    type0 denom=1.0*a_tol*sqrt(nc_dofs)/err_fac;
    type0 cost,cost_p;
    type0 norm=1.0,del=0.0,delp=0.0;
    int iter=0;
    int solver_iter;


    
    memcpy(c_0,c,ncs*sizeof(type0));
    memcpy(c,y_0,ncs*sizeof(type0));
    dynamic->update(atoms->c);
    cost=cost_p=ff->update_J(beta,a,F)/res_tol;
    
    while(cost>=1.0 && iter<max_niters_nonlin)
    {
        for(int i=0;i<ncs;i++) del_c[i]=0.0;
        
        gmres->solve(res_tol,F,del_c,solver_iter,norm);
        
        del=fabs(update_c()*norm/denom);

        dynamic->update(atoms->c);
        cost_p=ff->update_J(beta,a,F)/res_tol;
        cost=MIN(cost_p,del*err_fac*10.0);
        delp=del;
        iter++;
    }
    
    

    if(cost<1.0)
    {
        if(iter) nnonlin_acc++;
        memcpy(c,c_0,ncs*sizeof(type0));
        return true;
    }
    
    nnonlin_rej++;
    return false;
}
