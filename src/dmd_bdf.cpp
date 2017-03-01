#include "dmd_bdf.h"
#include "atoms_dmd.h"
#include "ff_dmd.h"
#include "dynamic_dmd.h"
#include "neighbor_dmd.h"
#include "memory.h"
#include <limits>
using namespace MAPP_NS;
/*--------------------------------------------
 
 --------------------------------------------*/
DMDBDF::DMDBDF():
DMDImplicit(),
t{[0 ... max_q]=0.0},
l{[0 ... max_q]=0.0},
A{[0 ... max_q]={[0 ... max_q]= 0.0}},
A_bar{[0 ... max_q]={[0 ... max_q]= 0.0}},
z(NULL),
dy(NULL)
{
    for(int i=0;i<max_q+1;i++)
        A[0][i]=A[i][i]=1.0;
    
    for(int i=2;i<max_q+1;i++)
        for(int j=1;j<i;j++)
            A[j][i]=A[j][i-1]+A[j-1][i-1];

}
/*--------------------------------------------
 
 --------------------------------------------*/
DMDBDF::~DMDBDF()
{
    delete dy;
    delete z;
}
/*--------------------------------------------
 
 --------------------------------------------*/
void DMDBDF::init_static()
{
    DMDImplicit::init_static();
    dy=new type0[(max_q+2)*ncs];
    z=dy+ncs;
    
    dynamic=new DynamicDMD(atoms,ff,false,{atoms->elem,atoms->c},{},{});
    dynamic->init();
}
/*--------------------------------------------
 
 --------------------------------------------*/
void DMDBDF::fin_static()
{
    delete dynamic;
    dynamic=NULL;
    delete [] dy;
    dy=z=NULL;
    DMDImplicit::fin_static();
}
/*--------------------------------------------
 
 --------------------------------------------*/
void DMDBDF::run_static(type0 t_tot)
{
    ff->derivative_timer();
    ff->neighbor->init_static();
    ff->neighbor->create_2nd_list();
    
    t_cur=0.0;
    t_fin=t_tot;
    nconst_q=0;
    nconst_dt=0;
    
    type0 dt_prev;
    int q_prev;
    int istep=0;
    for(;istep<max_nsteps;istep++)
    {
        dt_prev=dt;
        q_prev=q;
        
        while(!integrate()) integrate_fail();
        
        
        if(q_prev!=q) nconst_q=1;
        else nconst_q++;
        
        if(dt_prev!=dt) nconst_dt=1;
        else nconst_dt++;
        
        
        
        
        dt_prev=dt;
        q_prev=q;
        prep_for_next();
        
        if(dt_prev!=dt)
            nconst_dt=0;
        if(q_prev!=q)
            nconst_q=0;
        t_cur+=dt_prev;
    }
    
    ff->neighbor->fin_static();
}
/*--------------------------------------------
 
 --------------------------------------------*/
bool DMDBDF::integrate()
{
    
    while (true)
    {
        /*
         do interpolation
         
         if jumps out of boundary do
         a fail safe interpolation
         */
        
        if(!interpolate())
        {
            interpolate_fail();
            continue;
        }
        
            
        
        /*
         solve the nonlinear equation
         if unable to solve
         do somthing and start over (continue;)
         possible things to do in case of failure:
         reduce order
         reduce time step
         */
        
        if(!nonlin())
        {
            nonlin_fail();
            continue;
        }
        
        break;
    }
    
    
    
    /*
     here we have the solution
     now check for the integration error
     if it is acceptable we are done
     and continue, if not do somethoing
     
     possible things to do in case of failure:
     reduce order
     reduce time step
     */
    
    type0 norm_lcl=0.0,norm;
    for(int i=0;i<ncs;i++)
    {
        dy[i]=c[i]-y_0[i];
        norm_lcl+=dy[i]*dy[i];
    }
    MPI_Allreduce(&norm_lcl,&norm,1,MPI_TYPE0,MPI_SUM,atoms->world);
    err=err_fac*sqrt(norm/nc_dofs)/a_tol;
    
    if(err>1.0)
    {
        ninteg_rej++;
        return false;
    }
    
    ninteg_acc++;
    return true;
}
/*--------------------------------------------
 we want to solve:
 c - y_0 = beta * (c_d - y_d_0)
 or 
 c - beta * c_d + a = 0
 where 
 a = beta * y_d_0 - y_0
 
 but we have to be careful not to exceed the
 boundary i.e.
 0 <= y_0 <= 1
 
 --------------------------------------------*/
bool DMDBDF::interpolate()
{
    type0 __dt;
    for(int i=0;i<q+1;i++)
    {
        __dt=1.0;
        for(int j=i;j<q+1;j++)
        {
            A_bar[i][j]=A[i][j]*__dt;
            __dt*=dt;
        }
    }
    
    /*
     to unroll the loops
     otherwise it'll be shitty
     */
    bool in_domain=true;
    switch(q)
    {
        case 1: in_domain=DMDBDFMath::interpolate<1>(ncs,beta,A[0],&A[1][1],z,y_0,a);
            break;
        case 2: in_domain=DMDBDFMath::interpolate<2>(ncs,beta,A[0],&A[1][1],z,y_0,a);
            break;
        case 3: in_domain=DMDBDFMath::interpolate<3>(ncs,beta,A[0],&A[1][1],z,y_0,a);
            break;
        case 4: in_domain=DMDBDFMath::interpolate<4>(ncs,beta,A[0],&A[1][1],z,y_0,a);
            break;
        default: in_domain=DMDBDFMath::interpolate<max_q>(ncs,beta,A[0],&A[1][1],z,y_0,a);
    }
    
    int domain_err,domain_err_lcl=in_domain ? 0:1;
    MPI_Allreduce(&domain_err_lcl,&domain_err,1,MPI_INT,MPI_MAX,atoms->world);
    
    if(domain_err)
    {
        nintpol_rej++;
        return false;
    }
    
    err_fac=DMDBDFMath::err_fac_calc(q,dt,t,lo_err_fac,hi_err_fac);
    
    nintpol_acc++;
    return true;
}
/*--------------------------------------------
 here A and l are already updated
 --------------------------------------------*/
void DMDBDF::prep_for_next()
{
    //first decide the new order and new time step
    
    
    type0 tmp,norm,iq=static_cast<type0>(q);
    eta[1]=pow(0.5/err,1.0/(iq+1.0));
    eta[0]=eta[2]=0.0;
    
    if(nconst_q>q && q<max_q)
    {
        type0 norm_lcl=0.0;
        for(int i=0;i<ncs;i++)
        {
            if(c[i]>=0.0)
            {
                tmp=c[i]-y_0[i]+hi_err_fac[1]*dy[i];
                norm_lcl+=tmp*tmp;
            }
        }
        
        MPI_Allreduce(&norm_lcl,&tmp,1,MPI_TYPE0,MPI_SUM,atoms->world);
        norm=hi_err_fac[0]*sqrt(tmp/nc_dofs)/a_tol;
        eta[2]=pow(0.5/norm,1.0/(iq+2.0));
    }
    
    if(q>1)
    {
        type0* __z=z;
        type0 norm_lcl=0.0;
        for(int i=0;i<ncs;i++)
        {
            if(c[i]>=0.0)
            {
                tmp=c[i]-y_0[i]+lo_err_fac[1]*__z[q];
                norm_lcl+=tmp*tmp;
            }
            
            __z+=max_q+1;
        }
        MPI_Allreduce(&norm_lcl,&tmp,1,MPI_TYPE0,MPI_SUM,atoms->world);
        norm=lo_err_fac[0]*sqrt(tmp/nc_dofs)/a_tol;
        eta[0]=pow(0.5/norm,1.0/(iq));
    }
    
    
    type0 dt_r;
    if(eta[0]>eta[1] && eta[0]>eta[2])
    {
        dq=-1;
        dt_r=eta[0];
    }
    else if(eta[2]>eta[0] && eta[2]>eta[1])
    {
        dq=1;
        dt_r=eta[2];
    }
    else
    {
        dq=0;
        dt_r=eta[1];
    }
    
    
    if(10.0<dt_r)
        dt_r=10.0;
    else if(1.0<=dt_r && dt_r<1.5)
        dt_r=1.0;
    else if(dt_r<=0.5)
        dt_r=0.5;
    
    type0 new_dt;
    if(dt_r*dt>t_fin-t_cur)
        new_dt=t_fin-t_cur;
    else
    {
        if(t_fin-t_cur<=2.0*dt_min)
            new_dt=2.0*dt_min;
        else
        {
            if(dt_r*dt<dt_min)
                new_dt=dt_min;
            else if(dt_r*dt>=t_fin-t_cur-dt_min)
                new_dt=t_fin-t_cur-dt_min;
            else
                new_dt=dt_r*dt;
        }
    }
    
    
    update_z();
    
    for(int i=max_q;i>0;i--)
        t[i]=t[i-1]-dt;
    
    dt=new_dt;
    q+=dq;
}
/*--------------------------------------------
 
 --------------------------------------------*/
void DMDBDF::update_z()
{
    DMDBDFMath::prep_A_bar_l(q,dt,t,dq,l,A_bar);
    
    if(dq==1)
    {
        switch(q)
        {
            case 1: DMDBDFMath::update_z_inc<1>(ncs,A_bar[0],z,dy,l);
                break;
            case 2: DMDBDFMath::update_z_inc<2>(ncs,A_bar[0],z,dy,l);
                break;
            case 3: DMDBDFMath::update_z_inc<3>(ncs,A_bar[0],z,dy,l);
                break;
            case 4: DMDBDFMath::update_z_inc<4>(ncs,A_bar[0],z,dy,l);
                break;
            default: DMDBDFMath::update_z_inc<max_q>(ncs,A_bar[0],z,dy,l);
        }
    }
    else
    {
        switch(q+dq)
        {
            case 1: DMDBDFMath::update_z<1>(ncs,A_bar[0],z,dy,l);
                break;
            case 2: DMDBDFMath::update_z<2>(ncs,A_bar[0],z,dy,l);
                break;
            case 3: DMDBDFMath::update_z<3>(ncs,A_bar[0],z,dy,l);
                break;
            case 4: DMDBDFMath::update_z<4>(ncs,A_bar[0],z,dy,l);
                break;
            default: DMDBDFMath::update_z<max_q>(ncs,A_bar[0],z,dy,l);

        }
    }
}
/*--------------------------------------------
 run
 --------------------------------------------*/
void DMDBDF::nonlin_fail()
{
    if(t_fin-t_cur<=2.0*dt_min)
    {
        if(q>1)
            q--;
        /*
        else
            Error::abort("reached minimum order & del_t (%e) at line %d",dt,__LINE__);
         */
    }
    else
    {
        if(dt==dt_min)
        {
            if(q>1)
                q--;
            /*
            else
                Error::abort("reached minimum order & del_t (%e) at line %d",dt,__LINE__);
             */
        }
        else
        {
            type0 r=0.25;
            if(r*dt<dt_min)
                dt=dt_min;
            else if(r*dt>t_fin-t_cur-dt_min)
                dt=t_fin-t_cur-dt_min;
            else
                dt*=r;
        }
    }
}
/*--------------------------------------------
 
 --------------------------------------------*/
void DMDBDF::interpolate_fail()
{

    auto get_max_dt=[](type0 x,type0 x_d)->type0
    {
        
        if(x_d>0.0)
        {
            type0 ans=(1.0-x)/x_d;
            while(x+x_d*ans>1.0)
                ans=std::nextafter(ans,-1.0);
            return ans;
        }
        else if(x_d<0.0)
        {
            type0 ans=-x/x_d;
            while(x+x_d*ans<0.0)
                ans=std::nextafter(ans,1.0);
            return ans;
        }
        else
            return std::numeric_limits<type0>::infinity();
    };
    
    
    type0 max_dt_lcl=std::numeric_limits<type0>::infinity();
    type0* __z=z;
    for(int i=0;i<ncs;i++,__z+=max_q+1)
    {
        if(c[i]>=0.0)
        {
            max_dt_lcl=MIN(max_dt_lcl,get_max_dt(c[i],c_d[i]));
            __z[0]=c[i];
            __z[1]=c_d[i];
        }
        else
            __z[0]=__z[1]=0.0;
    }
    
    type0 max_dt;
    MPI_Allreduce(&max_dt_lcl,&max_dt,1,MPI_TYPE0,MPI_MIN,atoms->world);
    
    // ddc^2
    type0 norm=ff->ddc_norm()/sqrt(nc_dofs);
    q=1;
    type0 __dt=MIN(sqrt(2.0*a_tol/norm),max_dt);
    if(__dt<min_dt)
    {
        
    }
}
/*--------------------------------------------
 run
 --------------------------------------------*/
void DMDBDF::integrate_fail()
{
    if(t_fin-t_cur<=2.0*dt_min)
    {
        if(q>1)
            q--;
        /*
        else
            Error::abort("reached minimum order & del_t (%e) at line %d",dt,__LINE__);
         */
    }
    else
    {
        if(dt==dt_min)
        {
            if(q>1)
                q--;
            /*
            else
                Error::abort("reached minimum order & del_t (%e) at line %d",dt,__LINE__);
             */
        }
        else
        {
            type0 r=pow(0.5/err,1.0/static_cast<type0>(q+1));
            
            if(r*dt<dt_min)
                dt=dt_min;
            else if(r*dt>t_fin-t_cur-dt_min)
                dt=t_fin-t_cur-dt_min;
            else
                dt*=r;
        }
    }
}
/*------------------------------------------------------------------------------------------------------------------------------------
 
 ------------------------------------------------------------------------------------------------------------------------------------*/
PyObject* DMDBDF::__new__(PyTypeObject* type,PyObject* args,PyObject* kwds)
{
    Object* __self=reinterpret_cast<Object*>(type->tp_alloc(type,0));
    PyObject* self=reinterpret_cast<PyObject*>(__self);
    return self;
}
/*--------------------------------------------
 
 --------------------------------------------*/
int DMDBDF::__init__(PyObject* self,PyObject* args,PyObject* kwds)
{
    FuncAPI<> f("__init__");

    
    if(f(args,kwds)==-1) return -1;
    
    
    return 0;
}
/*--------------------------------------------
 
 --------------------------------------------*/
PyObject* DMDBDF::__alloc__(PyTypeObject* type,Py_ssize_t)
{
    Object* __self=new Object;
    __self->ob_type=type;
    __self->ob_refcnt=1;
    __self->dmd=NULL;
    return reinterpret_cast<PyObject*>(__self);
}
/*--------------------------------------------
 
 --------------------------------------------*/
void DMDBDF::__dealloc__(PyObject* self)
{
    Object* __self=reinterpret_cast<Object*>(self);
    delete __self->dmd;
    __self->dmd=NULL;
    delete __self;
}





