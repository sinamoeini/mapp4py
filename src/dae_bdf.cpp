#include "dae_bdf.h"
#include "atoms_dmd.h"
#include "ff_dmd.h"
#include "dynamic_dmd.h"
#include "neighbor_dmd.h"
#include "memory.h"
#include "thermo_dynamics.h"
#include "MAPP.h"
#include <limits>
using namespace MAPP_NS;
/*--------------------------------------------
 
 --------------------------------------------*/
DAEBDF::DAEBDF():
DAEImplicit(),
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
DAEBDF::~DAEBDF()
{
    delete dy;
    delete z;
}
/*--------------------------------------------
 
 --------------------------------------------*/
void DAEBDF::init_static()
{
    DAEImplicit::init_static();
    Memory::alloc(dy,(max_q+2)*ncs);
    z=dy+ncs;
}
/*--------------------------------------------
 
 --------------------------------------------*/
void DAEBDF::run_static(type0 t_tot)
{
    init_static();
    ff->reset();
    ff->derivative_timer();
    ff->neighbor->init_static();
    ff->neighbor->create_2nd_list();
    ff->init_static();
    /*
    type0* alpha=atoms->alpha->begin();
    for(int i=0;i<ncs;i++)
        printf("%d\t%e\n",i,alpha[i]);
    */
    
    t_cur=0.0;
    t_fin=t_tot;
    nconst_q=nconst_dt=nnonlin_acc=nnonlin_rej=ninteg_acc=ninteg_rej=nintpol_acc=nintpol_rej=0;
    int __max_q=1;
    type0 __max_dt=0.0,__min_dt=std::numeric_limits<type0>::infinity();
    reset();

    
    type0 S[__dim__][__dim__];
    
    ThermoDynamics thermo(6,
    "Time",t_cur,
    "FE",ff->nrgy_strss[0],
    "S[0][0]",S[0][0],
    "S[1][1]",S[1][1],
    "S[2][2]",S[2][2],
    "S[1][2]",S[2][1],
    "S[2][0]",S[2][0],
    "S[0][1]",S[1][0]);
    thermo.init();
    Algebra::DyadicV_2_MLT(ff->nrgy_strss+1,S);
    thermo.print(0);
    
    type0 dt_prev;
    int q_prev;
    int istep=0;
    for(;istep<max_nsteps && t_cur<t_fin;istep++)
    {
        dt_prev=dt;
        q_prev=q;
        
        
        
        while(!integrate())
            integrate_fail();
        //printf("step %d %e %e %e\n",istep,t_cur+dt,c[0],c_d[0]);
        
        __max_q=MAX(q,__max_q);
        __max_dt=MAX(dt,__max_dt);
        __min_dt=MIN(dt,__min_dt);
        
        
        if(q_prev!=q) nconst_q=1;
        else nconst_q++;
        
        if(dt_prev!=dt) nconst_dt=1;
        else nconst_dt++;
        
        
        
        
        dt_prev=dt;
        q_prev=q;
        t_cur+=dt_prev;
        prep_for_next();
        
        if(dt_prev!=dt)
            nconst_dt=0;
        if(q_prev!=q)
            nconst_q=0;
        
        if((istep+1)%ntally==0)
        {
            ff->reset();
            ff->force_calc_static_timer();
            Algebra::DyadicV_2_MLT(ff->nrgy_strss+1,S);
            thermo.print(istep+1);
        }
        
    }
    
    if(istep%ntally)
    {
        ff->reset();
        ff->force_calc_static_timer();
        Algebra::DyadicV_2_MLT(ff->nrgy_strss+1,S);
        thermo.print(istep);
    }
    
    thermo.fin();
    
    /*
    //type0* alpha=atoms->alpha->begin();
    for(int i=0;i<ncs;i++)
        printf("%d\t%e\t%e\n",i,c[i],c_d[i]);
    */
    
    fprintf(MAPP::mapp_out,"nonlin: accepted = %d rejected = %d\n",nnonlin_acc,nnonlin_rej);
    fprintf(MAPP::mapp_out,"intrtp: accepetd = %d rejected = %d\n",nintpol_acc,nintpol_rej);
    fprintf(MAPP::mapp_out,"integr: accepetd = %d rejected = %d\n",ninteg_acc,ninteg_rej);
    fprintf(MAPP::mapp_out,"maximum order: %d\n",__max_q);
    fprintf(MAPP::mapp_out,"maximum timestep: %e\n",__max_dt);
    fprintf(MAPP::mapp_out,"minimum timestep: %e\n",__min_dt);
    

    fin_static();
}
/*--------------------------------------------
 
 --------------------------------------------*/
void DAEBDF::fin_static()
{
    ff->fin_static();
    ff->neighbor->fin_static();
    Memory::dealloc(dy);
    z=NULL;
    DAEImplicit::fin_static();
}
/*--------------------------------------------
 
 --------------------------------------------*/
bool DAEBDF::integrate()
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
        norm_lcl+=(c[i]-y_0[i])*(c[i]-y_0[i]);
    MPI_Allreduce(&norm_lcl,&norm,1,Vec<type0>::MPI_T,MPI_SUM,atoms->world);
    err=err_fac*sqrt(norm)/a_tol_sqrt_nc_dofs;
    
    //printf("Error: %e\n",err);
    if(err<1.0)
    {
        ninteg_acc++;
        return true;
    }

    memcpy(c,c_0,ncs*sizeof(type0));
    dynamic->update(atoms->c);
    ninteg_rej++;
    return false;
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
bool DAEBDF::interpolate()
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
    err_fac=VC_y::err_fac_calc(q,dt,t,lo_err_fac,hi_err_fac,beta);
    
    /*
     to unroll the loops
     otherwise it'll be shitty
     */
    
    bool in_domain=true;
    switch(q)
    {
        case 1: in_domain=DAEBDFMath::interpolate<1>(ncs,beta,A_bar[0],&A_bar[1][1],z,y_0,a);
            break;
        case 2: in_domain=DAEBDFMath::interpolate<2>(ncs,beta,A_bar[0],&A_bar[1][1],z,y_0,a);
            break;
        case 3: in_domain=DAEBDFMath::interpolate<3>(ncs,beta,A_bar[0],&A_bar[1][1],z,y_0,a);
            break;
        case 4: in_domain=DAEBDFMath::interpolate<4>(ncs,beta,A_bar[0],&A_bar[1][1],z,y_0,a);
            break;
        default: in_domain=DAEBDFMath::interpolate<max_q>(ncs,beta,A_bar[0],&A_bar[1][1],z,y_0,a);
    }
    
    int domain_err,domain_err_lcl=in_domain ? 0:1;
    MPI_Allreduce(&domain_err_lcl,&domain_err,1,MPI_INT,MPI_MAX,atoms->world);
    
    if(domain_err)
    {
        nintpol_rej++;
        return false;
    }
    
    
    
    nintpol_acc++;
    return true;
}
/*--------------------------------------------
 here A and l are already updated
 --------------------------------------------*/
void DAEBDF::prep_for_next()
{
    //first decide the new order and new time step
    
    
    type0 tmp,norm,iq=static_cast<type0>(q);
    eta[1]=pow(0.5/err,1.0/(iq+1.0));
    eta[0]=eta[2]=0.0;
    
    if(q<max_q)
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
        
        MPI_Allreduce(&norm_lcl,&tmp,1,Vec<type0>::MPI_T,MPI_SUM,atoms->world);
        norm=hi_err_fac[0]*sqrt(tmp)/a_tol_sqrt_nc_dofs;
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
        MPI_Allreduce(&norm_lcl,&tmp,1,Vec<type0>::MPI_T,MPI_SUM,atoms->world);
        norm=lo_err_fac[0]*sqrt(tmp)/a_tol_sqrt_nc_dofs;
        eta[0]=pow(0.5/norm,1.0/(iq));
    }
    
    
    for(int i=0;i<ncs;i++)
    {
        dy[i]=c[i]-y_0[i];
    }
    
    
    type0 dt_r=1.0;
    dq=0;
    
    
    if(nconst_q>2+q && nconst_dt>2+q)
    {
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
    }
    
    if(2.0<dt_r)
        dt_r=2.0;
    else if(dt_r<=0.5)
        dt_r=0.5;
    
    type0 new_dt;
    if(dt_r*dt>t_fin-t_cur)
    {
        
        new_dt=t_fin-t_cur;
    }
    else
    {
        if(t_fin-t_cur<=2.0*min_dt)
            new_dt=2.0*min_dt;
        else
        {
            if(dt_r*dt<min_dt)
                new_dt=min_dt;
            else if(dt_r*dt>=t_fin-t_cur-min_dt)
                new_dt=t_fin-t_cur-min_dt;
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
void DAEBDF::update_z()
{
    VC_y::prep_A_bar_l(q,dt,t,dq,l,A_bar);
    
    if(dq==1)
    {
        switch(q)
        {
            case 1: DAEBDFMath::update_z_inc<1>(ncs,A_bar[0],z,dy,l);
                break;
            case 2: DAEBDFMath::update_z_inc<2>(ncs,A_bar[0],z,dy,l);
                break;
            case 3: DAEBDFMath::update_z_inc<3>(ncs,A_bar[0],z,dy,l);
                break;
            case 4: DAEBDFMath::update_z_inc<4>(ncs,A_bar[0],z,dy,l);
                break;
            default: DAEBDFMath::update_z_inc<max_q>(ncs,A_bar[0],z,dy,l);
        }
    }
    else
    {
        switch(q+dq)
        {
            case 1: DAEBDFMath::update_z<1>(ncs,A_bar[0],z,dy,l);
                break;
            case 2: DAEBDFMath::update_z<2>(ncs,A_bar[0],z,dy,l);
                break;
            case 3: DAEBDFMath::update_z<3>(ncs,A_bar[0],z,dy,l);
                break;
            case 4: DAEBDFMath::update_z<4>(ncs,A_bar[0],z,dy,l);
                break;
            default: DAEBDFMath::update_z<max_q>(ncs,A_bar[0],z,dy,l);

        }
    }
}
/*--------------------------------------------
 run
 --------------------------------------------*/
void DAEBDF::nonlin_fail()
{
    if(t_fin-t_cur<=2.0*min_dt)
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
        if(dt==min_dt)
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
            if(r*dt<min_dt)
                dt=min_dt;
            else if(r*dt>t_fin-t_cur-min_dt)
                dt=t_fin-t_cur-min_dt;
            else
                dt*=r;
        }
    }
}
/*--------------------------------------------
 
 --------------------------------------------*/
void DAEBDF::reset()
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
                ans=std::nextafter(ans,-1.0);
            return ans;
        }
        else
            return std::numeric_limits<type0>::infinity();
    };
    
    // ddc^2
    type0 norm=ff->ddc_norm()/a_tol_sqrt_nc_dofs;
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
    MPI_Allreduce(&max_dt_lcl,&max_dt,1,Vec<type0>::MPI_T,MPI_MIN,atoms->world);

    
    q=1;
    dt=MAX(MIN(MIN(sqrt(2.0/norm),(t_fin-t_cur)*0.001),max_dt),min_dt);
    Algebra::zero<max_q+1>(t);
    memset(dy,0,ncs*sizeof(type0));
}
/*--------------------------------------------
 
 --------------------------------------------*/
void DAEBDF::interpolate_fail()
{
    reset();
}
/*--------------------------------------------
 run
 --------------------------------------------*/
void DAEBDF::integrate_fail()
{
    if(t_fin-t_cur<=2.0*min_dt)
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
        if(dt==min_dt)
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
            
            
            
            
            if(r*dt<min_dt)
                dt=min_dt;
            else if(r*dt>t_fin-t_cur-min_dt)
                dt=t_fin-t_cur-min_dt;
            else
                dt*=r;
        }
    }
}
/*------------------------------------------------------------------------------------------------------------------------------------
 
 ------------------------------------------------------------------------------------------------------------------------------------*/
PyObject* DAEBDF::__new__(PyTypeObject* type,PyObject* args,PyObject* kwds)
{
    Object* __self=reinterpret_cast<Object*>(type->tp_alloc(type,0));
    PyObject* self=reinterpret_cast<PyObject*>(__self);
    return self;
}
/*--------------------------------------------
 
 --------------------------------------------*/
int DAEBDF::__init__(PyObject* self,PyObject* args,PyObject* kwds)
{
    FuncAPI<> f("__init__");

    
    if(f(args,kwds)==-1) return -1;
    Object* __self=reinterpret_cast<Object*>(self);
    __self->dae=new DAEBDF();
    return 0;
}
/*--------------------------------------------
 
 --------------------------------------------*/
PyObject* DAEBDF::__alloc__(PyTypeObject* type,Py_ssize_t)
{
    Object* __self=new Object;
    __self->ob_type=type;
    __self->ob_refcnt=1;
    __self->dae=NULL;
    return reinterpret_cast<PyObject*>(__self);
}
/*--------------------------------------------
 
 --------------------------------------------*/
void DAEBDF::__dealloc__(PyObject* self)
{
    Object* __self=reinterpret_cast<Object*>(self);
    delete __self->dae;
    __self->dae=NULL;
    delete __self;
}
/*--------------------------------------------*/
PyTypeObject DAEBDF::TypeObject={PyObject_HEAD_INIT(NULL)};
/*--------------------------------------------*/
void DAEBDF::setup_tp()
{
    TypeObject.tp_name="mapp.dmd.bdf";
    TypeObject.tp_doc="chemical integration";
    
    TypeObject.tp_flags=Py_TPFLAGS_DEFAULT;
    TypeObject.tp_basicsize=sizeof(Object);
    
    TypeObject.tp_new=__new__;
    TypeObject.tp_init=__init__;
    TypeObject.tp_alloc=__alloc__;
    TypeObject.tp_dealloc=__dealloc__;
    setup_tp_methods();
    TypeObject.tp_methods=methods;
    setup_tp_getset();
    TypeObject.tp_getset=getset;
}
/*--------------------------------------------*/
PyGetSetDef DAEBDF::getset[]={[0 ... 6]={NULL,NULL,NULL,NULL,NULL}};
/*--------------------------------------------*/
void DAEBDF::setup_tp_getset()
{
    getset_a_tol(getset[0]);
    getset_max_nsteps(getset[1]);
    getset_min_dt(getset[2]);
    getset_max_ngmres_iters(getset[3]);
    getset_max_nnewton_iters(getset[4]);
    getset_ntally(getset[5]);
}
/*--------------------------------------------*/
PyMethodDef DAEBDF::methods[]={[0 ... 1]={NULL}};
/*--------------------------------------------*/
void DAEBDF::setup_tp_methods()
{
    ml_run(methods[0]);
}    
/*--------------------------------------------
 
 --------------------------------------------*/
void DAEBDF::ml_run(PyMethodDef& tp_methods)
{
    tp_methods.ml_flags=METH_VARARGS | METH_KEYWORDS;
    tp_methods.ml_name="run";
    tp_methods.ml_doc="run chemical integration";
    
    tp_methods.ml_meth=(PyCFunction)(PyCFunctionWithKeywords)
    [](PyObject* self,PyObject* args,PyObject* kwds)->PyObject*
    {
        Object* __self=reinterpret_cast<Object*>(self);
        FuncAPI<OP<AtomsDMD>,type0> f("run",{"atoms","t"});
        f.logics<1>()[0]=VLogics("ge",0.0);
        if(f(args,kwds)) return NULL;
        
        AtomsDMD* __atoms=reinterpret_cast<AtomsDMD::Object*>(f.val<0>().ob)->atoms;
        ForceFieldDMD* __ff=reinterpret_cast<AtomsDMD::Object*>(f.val<0>().ob)->ff;

        try
        {
            __self->dae->pre_run_chk(__atoms,__ff);
        }
        catch(std::string err_msg)
        {
            PyErr_SetString(PyExc_TypeError,err_msg.c_str());
            return NULL;
        }
        
        __self->dae->atoms=
        __atoms;
        __self->dae->ff=__ff;
        
        __self->dae->run_static(f.val<1>());
        
        __self->dae->ff=NULL;
        __self->dae->atoms=NULL;
        
        Py_RETURN_NONE;
    };
}
