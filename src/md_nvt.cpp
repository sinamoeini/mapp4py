#include "md_nvt.h"
#include "xmath.h"
#include "atoms_md.h"
#include "ff_md.h"
#include "elements.h"
#include "dynamic_md.h"
#include "MAPP.h"
#include "thermo_dynamics.h"
#include <math.h>
using namespace MAPP_NS;
/*--------------------------------------------
 
 --------------------------------------------*/
MDNVT::MDNVT(type0 __T,type0 __dt):
dynamic(NULL),
atoms(NULL),
ff(NULL),
xprt(NULL),
dt(__dt),
T(__T),
dt2(__dt/2.0),
ntally(10000),
thermo_part(__dt/2.0,100.0*__dt,3,1)
{
}
/*--------------------------------------------
 
 --------------------------------------------*/
MDNVT::~MDNVT()
{
}
/*--------------------------------------------
 
 --------------------------------------------*/
void MDNVT::change_dt(type0 __dt)
{
    if(__dt==dt) return;
    dt=__dt;
    dt2=__dt/2;
    
    int L=thermo_part.L;
    int niters=thermo_part.niters;
    type0 t_relax=thermo_part.t_relax;
    
    thermo_part.~ThermostatNHC();
    new (&thermo_part) ThermostatNHC(__dt/2.0,t_relax,L,niters);
}
/*--------------------------------------------
 
 --------------------------------------------*/
void MDNVT::update_x_d_final(type0 fac_x_d)
{
    const int n=atoms->natms_lcl*__dim__;
    type0* x_d=atoms->x_d->begin();
    for(int i=0;i<n;i++) x_d[i]*=fac_x_d;
}
/*--------------------------------------------
 
 --------------------------------------------*/
void MDNVT::update_x_d_final_w_dof(type0 fac_x_d)
{
    const int n=atoms->natms_lcl*__dim__;
    type0* x_d=atoms->x_d->begin();
    bool* dof=atoms->x_dof->begin();
    for(int i=0;i<n;i++) if(dof[i]) x_d[i]*=fac_x_d;
}
/*--------------------------------------------
 
 --------------------------------------------*/
void MDNVT::update_x_d__x__x_d(type0 fac_x_d)
{
    type0* x=atoms->x->begin();
    type0* f=ff->f->begin();
    type0* x_d=atoms->x_d->begin();
    elem_type* elem=atoms->elem->begin();
    type0* m=atoms->elements.masses;
    type0 m_i;
    type0 dx_lcl[__dim__]={DESIG(__dim__,0.0)};
    const int natms0=atoms->natms_lcl;
    for(int i=0;i<natms0;++i)
    {
        m_i=m[*elem];
        Algebra::Do<__dim__>::func([&dx_lcl,&x_d,&x,&f,&m_i,&fac_x_d,this](const int j)
        {
            x_d[j]=x_d[j]*fac_x_d+f[j]*dt2/m_i;
            dx_lcl[j]+=x_d[j]*dt;
            x[j]+=x_d[j]*dt;
        });
        
        f+=__dim__;
        x_d+=__dim__;
        x+=__dim__;
        ++elem;
    }
    type0 dx[__dim__]={DESIG(__dim__,0.0)};
    MPI_Allreduce(dx_lcl,dx,__dim__,Vec<type0>::MPI_T,MPI_SUM,atoms->world);
    type0 natms=static_cast<type0>(atoms->natms);
    Algebra::Do<__dim__>::func([&dx,natms](const int i){dx[i]/=natms;});
    x=atoms->x->begin();
    for(int i=0;i<natms0;++i,x+=__dim__)
        Algebra::Do<__dim__>::func([&dx,&x](const int j){x[j]-=dx[j];});
    
    
    
    dynamic->update(atoms->x);
    ff->force_calc_timer();
    
    f=ff->f->begin();
    x_d=atoms->x_d->begin();
    elem=atoms->elem->begin();
    Algebra::zero<__nvoigt__>(__vec_lcl);
    const int natms1=atoms->natms_lcl;

    for(int i=0;i<natms1;++i)
    {
        m_i=m[*elem];
        Algebra::Do<__dim__>::func([&x_d,&f,&m_i,this](const int j){x_d[j]+=f[j]*dt2/m_i;});
        Algebra::DyadicV<__dim__>(m_i,x_d,__vec_lcl);
        
        f+=__dim__;
        x_d+=__dim__;
        ++elem;
    }
    MPI_Allreduce(__vec_lcl,mvv,__nvoigt__,Vec<type0>::MPI_T,MPI_SUM,atoms->world);
    T_part=Algebra::Tr_DyadicV<__dim__>(mvv)/(ndof_part*kB);
}
/*--------------------------------------------
 
 --------------------------------------------*/
void MDNVT::update_x_d__x__x_d_w_dof(type0 fac_x_d)
{
    type0* x=atoms->x->begin();
    type0* f=ff->f->begin();
    type0* x_d=atoms->x_d->begin();
    elem_type* elem=atoms->elem->begin();
    type0* m=atoms->elements.masses;
    type0 m_i;
    type0 dx_lcl[__dim__]={DESIG(__dim__,0.0)};
    const int natms0=atoms->natms_lcl;
    bool* dof=atoms->x_dof->begin();
    for(int i=0;i<natms0;++i)
    {
        m_i=m[*elem];
        Algebra::Do<__dim__>::func([&dof,&dx_lcl,&x_d,&x,&f,&m_i,&fac_x_d,this](const int j)
        {
            if(dof[j])
            {
                x_d[j]=x_d[j]*fac_x_d+f[j]*dt2/m_i;
                dx_lcl[j]+=x_d[j]*dt;
            }
            
            x[j]+=x_d[j]*dt;
        });
        
        f+=__dim__;
        x_d+=__dim__;
        x+=__dim__;
        dof+=__dim__;
        ++elem;
    }
    type0 dx[__dim__]={DESIG(__dim__,0.0)};
    MPI_Allreduce(dx_lcl,dx,__dim__,Vec<type0>::MPI_T,MPI_SUM,atoms->world);
    
    Algebra::Do<__dim__>::func([&dx,this](const int i){dx[i]/=Ndof_part[i];});
    x=atoms->x->begin();
    dof=atoms->x_dof->begin();
    for(int i=0;i<natms0;++i,x+=__dim__,dof+=__dim__)
        Algebra::Do<__dim__>::func([&dx,&x,&dof](const int j){if(dof[j]) x[j]-=dx[j];});
    
    
    
    dynamic->update(atoms->x);
    ff->force_calc_timer();
    
    f=ff->f->begin();
    x_d=atoms->x_d->begin();
    elem=atoms->elem->begin();
    dof=atoms->x_dof->begin();
    Algebra::zero<__nvoigt__>(__vec_lcl);
    type0 __x_d[__dim__];
    const int natms1=atoms->natms_lcl;

    for(int i=0;i<natms1;++i)
    {
        m_i=m[*elem];
        Algebra::Do<__dim__>::func([&x_d,&f,&m_i,&dof,&__x_d,this](const int j)
        {
            x_d[j]+=f[j]*dt2/m_i;
            if(dof[j])
                __x_d[j]=x_d[j];
            else
                __x_d[j]=0.0;
        });
        Algebra::DyadicV<__dim__>(m_i,__x_d,__vec_lcl);
        
        f+=__dim__;
        x_d+=__dim__;
        dof+=__dim__;
        ++elem;
    }
    MPI_Allreduce(__vec_lcl,mvv,__nvoigt__,Vec<type0>::MPI_T,MPI_SUM,atoms->world);
    T_part=Algebra::Tr_DyadicV<__dim__>(mvv)/(ndof_part*kB);
}
/*--------------------------------------------
 pre run check it throw excepctions
 --------------------------------------------*/
void MDNVT::pre_run_chk(AtomsMD* atoms,ForceFieldMD* ff)
{
    //check if configuration is loaded
    if(!atoms)
        throw std::string("cannot start md without initial conditions");
    
    //check if force field is loaded
    if(!ff)
        throw std::string("cannot start md without governing equations (force field)");
    
    if(std::isnan(atoms->kB))
        throw std::string("boltzmann constant should be set prior to MD");
}
/*--------------------------------------------
 
 --------------------------------------------*/
void MDNVT::pre_init()
{
    dof_empty=atoms->x_dof->is_empty();
    kB=atoms->kB;
    /*
     calculating the number of degress of freedom
     */
    
    int natms=atoms->natms;
    ndof_part=atoms->natms*__dim__;
    Algebra::Do<__dim__>::func([this,&natms](int i){Ndof_part[i]=static_cast<type0>(natms);});
    
    if(!dof_empty)
    {
        bool* dof=atoms->x_dof->begin();
        int Ndof_red_lcl[__dim__]={DESIG(__dim__,0)};
        int natms_lcl=atoms->natms_lcl;
        
        for(int i=0;i<natms_lcl;i++,dof+=__dim__)
            Algebra::Do<__dim__>::func([&dof,&Ndof_red_lcl](int j){if(!dof[j])Ndof_red_lcl[j]++;});
        
        int Ndof_red[__dim__]={DESIG(__dim__,0)};
        
        MPI_Allreduce(Ndof_red_lcl,Ndof_red,__dim__,MPI_INT,MPI_SUM,atoms->world);

        Algebra::Do<__dim__>::func([&Ndof_red,this](int i)
        {
            ndof_part-=Ndof_red[i];
            Ndof_part[i]-=static_cast<type0>(Ndof_red[i]);
        });
    }
    
    /*
     temperature and kinetic energy tensor
     */
    if(atoms->x_d->is_empty())
    {
        atoms->x_d->fill();
        Algebra::Do<__nvoigt__>::func([this](int i){mvv[i]=0.0;});
        T_part=0.0;
    }
    else
    {
        type0* x_d=atoms->x_d->begin();
        type0* m=atoms->elements.masses;
        type0 m_i;
        elem_type* elem=atoms->elem->begin();
        Algebra::zero<__nvoigt__>(__vec_lcl);
        const int natms_lcl=atoms->natms_lcl;
        if(dof_empty)
        {
            for(int i=0;i<natms_lcl;i++)
            {
                m_i=m[*elem];
                Algebra::DyadicV<__dim__>(m_i,x_d,__vec_lcl);
                x_d+=__dim__;
                ++elem;
            }
        }
        else
        {
            bool* dof=atoms->x_dof->begin();
            type0 __x_d[__dim__];
            for(int i=0;i<natms_lcl;i++)
            {
                Algebra::V_eq<__dim__>(x_d,__x_d);
                Algebra::Do<__dim__>::func([&dof,&x_d,&__x_d](int i){__x_d[i]=dof[i] ? x_d[i]:0.0;});
                m_i=m[*elem];
                Algebra::DyadicV<__dim__>(m_i,__x_d,__vec_lcl);
                
                dof+=__dim__;
                x_d+=__dim__;
                ++elem;
            }
        }
        
        MPI_Allreduce(__vec_lcl,mvv,__nvoigt__,Vec<type0>::MPI_T,MPI_SUM,atoms->world);
        T_part=Algebra::Tr_DyadicV<__dim__>(mvv)/(ndof_part*kB);
    }
}
/*--------------------------------------------
 
 --------------------------------------------*/
void MDNVT::init()
{
    pre_init();
    
    dynamic=new DynamicMD(atoms,ff,false,{},{atoms->x_d,atoms->x_dof},{});
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
 
 --------------------------------------------*/
void MDNVT::fin()
{
    if(xprt)
    {
        xprt->fin();
        xprt->atoms=NULL;
    }
    
    dynamic->fin();
    delete dynamic;
    dynamic=NULL;
}
/*--------------------------------------------
 
 --------------------------------------------*/
void MDNVT::run(int nsteps)
{
    int step=atoms->step;
    
    ff->force_calc_timer();
    
    int nevery_xprt=xprt==NULL ? 0:xprt->nevery;
    if(nevery_xprt) xprt->write(step);
    
    ThermoDynamics thermo(6,"T",T_part,"PE",atoms->pe,
    "S[0][0]",S_part[0][0],
    "S[1][1]",S_part[1][1],
    "S[2][2]",S_part[2][2],
    "S[1][2]",S_part[2][1],
    "S[2][0]",S_part[2][0],
    "S[0][1]",S_part[1][0]);
    
    if(ntally) thermo.init();
    Algebra::DoLT<__dim__>::func([this](const int i,const int j)
    {
        S_part[i][j]=atoms->S_pe[i][j]-mvv[i+j*__dim__-j*(j+1)/2]/atoms->vol;
    });
    
    
    if(ntally) thermo.print(step);
    
    type0 fac,fac_x_d=1.0;
    if(dof_empty)
    {
        for(int istep=0;istep<nsteps;istep++)
        {
            // particle thermostat
            fac_x_d*=fac=thermo_part(T_part/T,ndof_part);
            fac*=fac;
            Algebra::Do<__nvoigt__>::func([this,&fac](int i){mvv[i]*=fac;});
            T_part*=fac;
            
            update_x_d__x__x_d(fac_x_d);
            
            
            // particle thermostat
            fac_x_d=fac=thermo_part(T_part/T,ndof_part);
            fac*=fac;
            Algebra::Do<__nvoigt__>::func([this,&fac](int i){mvv[i]*=fac;});
            T_part*=fac;
            
            
            Algebra::DoLT<__dim__>::func([this](const int i,const int j)
            {
                S_part[i][j]=atoms->S_pe[i][j]-mvv[i+j*__dim__-j*(j+1)/2]/atoms->vol;
            });
            
            if(ntally && (istep+1)%ntally==0) thermo.print(step+istep+1);
            if(nevery_xprt && (istep+1)%nevery_xprt==0) xprt->write(step+istep+1);
        }
    }
    else
    {
        for(int istep=0;istep<nsteps;istep++)
        {
            // particle thermostat
            fac_x_d*=fac=thermo_part(T_part/T,ndof_part);
            fac*=fac;
            Algebra::Do<__nvoigt__>::func([this,&fac](int i){mvv[i]*=fac;});
            T_part*=fac;
            
            update_x_d__x__x_d_w_dof(fac_x_d);
            
            
            // particle thermostat
            fac_x_d=fac=thermo_part(T_part/T,ndof_part);
            fac*=fac;
            Algebra::Do<__nvoigt__>::func([this,&fac](int i){mvv[i]*=fac;});
            T_part*=fac;
            
            
            Algebra::DoLT<__dim__>::func([this](const int i,const int j)
            {
                S_part[i][j]=atoms->S_pe[i][j]-mvv[i+j*__dim__-j*(j+1)/2]/atoms->vol;
            });
            
            if(ntally && (istep+1)%ntally==0) thermo.print(step+istep+1);
            if(nevery_xprt && (istep+1)%nevery_xprt==0) xprt->write(step+istep+1);
        }
    }
    
    if(ntally && nsteps%ntally) thermo.print(step+nsteps);
    if(nevery_xprt && nsteps%nevery_xprt) xprt->write(step+nsteps);
    
    if(dof_empty)
        update_x_d_final(fac_x_d);
    else
        update_x_d_final_w_dof(fac_x_d);
    
    if(ntally) thermo.fin();
    
    atoms->step+=nsteps;
}
/*------------------------------------------------------------------------------------------------------------------------------------
 
 ------------------------------------------------------------------------------------------------------------------------------------*/
PyObject* MDNVT::__new__(PyTypeObject* type,PyObject* args,PyObject* kwds)
{
    Object* __self=reinterpret_cast<Object*>(type->tp_alloc(type,0));
    PyObject* self=reinterpret_cast<PyObject*>(__self);
    return self;
}
/*--------------------------------------------
 
 --------------------------------------------*/
int MDNVT::__init__(PyObject* self,PyObject* args,PyObject* kwds)
{
    FuncAPI<type0,type0> f("__init__",{"T","dt"});
    
    f.logics<0>()[0]=VLogics("gt",0.0);
    f.logics<1>()[0]=VLogics("gt",0.0);
    
    if(f(args,kwds)==-1) return -1;
    Object* __self=reinterpret_cast<Object*>(self);
    __self->md=new MDNVT(f.val<0>(),f.val<1>());
    __self->xprt=NULL;
    
    return 0;
}
/*--------------------------------------------
 
 --------------------------------------------*/
PyObject* MDNVT::__alloc__(PyTypeObject* type,Py_ssize_t)
{
    Object* __self=new Object;
    Py_TYPE(__self)=type;
    Py_REFCNT(__self)=1;
    __self->md=NULL;
    __self->xprt=NULL;
    return reinterpret_cast<PyObject*>(__self);
}
/*--------------------------------------------
 
 --------------------------------------------*/
void MDNVT::__dealloc__(PyObject* self)
{
    Object* __self=reinterpret_cast<Object*>(self);
    delete __self->md;
    __self->md=NULL;
    if(__self->xprt) Py_DECREF(__self->xprt);
    __self->xprt=NULL;
    delete __self;
}
/*--------------------------------------------*/
PyTypeObject MDNVT::TypeObject ={PyObject_HEAD_INIT(NULL)};
/*--------------------------------------------*/
int MDNVT::setup_tp()
{
    TypeObject.tp_name="mapp.md.nvt";
    TypeObject.tp_doc=R"---(
    __init__(T,dt)
    
    :math:`NVT` ensemble
    
    Molecular dynamics of canonical ensemble
        
    Parameters
    ----------
    T : double
       Temperature of the ensemble
    dt : double
       Time step for simulation
    
    Notes
    -----
       * Thermostat Details
          Nose-Hoover chain
    
    )---";
    
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
    
    int ichk=PyType_Ready(&TypeObject);
    if(ichk<0) return ichk;
    Py_INCREF(&TypeObject);
    GET_WRAPPER_DOC(TypeObject,__init__)=NULL;
    return ichk;
}
/*--------------------------------------------*/
PyGetSetDef MDNVT::getset[]=EmptyPyGetSetDef(8);
/*--------------------------------------------*/
void MDNVT::setup_tp_getset()
{
    getset_T(getset[0]);
    getset_dt(getset[1]);
    getset_niters(getset[2]);
    getset_L(getset[3]);
    getset_t_relax(getset[4]);
    getset_export(getset[5]);
    getset_ntally(getset[6]);
}
/*--------------------------------------------*/
PyMethodDef MDNVT::methods[]=EmptyPyMethodDef(2);
/*--------------------------------------------*/
void MDNVT::setup_tp_methods()
{
    ml_run(methods[0]);
}
/*--------------------------------------------
 
 --------------------------------------------*/
void MDNVT::getset_T(PyGetSetDef& getset)
{
    getset.name=(char*)"T";
    getset.doc=(char*)R"---(
    (double) temperature
    
    Temperature of the ensemble
    )---";
    getset.get=[](PyObject* self,void*)->PyObject*
    {
        return var<type0>::build(reinterpret_cast<Object*>(self)->md->T);
    };
    getset.set=[](PyObject* self,PyObject* op,void*)->int
    {
        VarAPI<type0> T("T");
        T.logics[0]=VLogics("gt",0.0);
        int ichk=T.set(op);
        if(ichk==-1) return -1;
        reinterpret_cast<Object*>(self)->md->T=T.val;
        return 0;
    };
}
/*--------------------------------------------
 
 --------------------------------------------*/
void MDNVT::getset_dt(PyGetSetDef& getset)
{
    getset.name=(char*)"dt";
    getset.doc=(char*)R"---(
    (double) dt
    
    Time step for simulation
    )---";
    getset.get=[](PyObject* self,void*)->PyObject*
    {
        return var<type0>::build(reinterpret_cast<Object*>(self)->md->dt);
    };
    getset.set=[](PyObject* self,PyObject* op,void*)->int
    {
        VarAPI<type0> dt("dt");
        dt.logics[0]=VLogics("gt",0.0);
        int ichk=dt.set(op);
        if(ichk==-1) return -1;
        
        reinterpret_cast<Object*>(self)->md->change_dt(dt.val);
        return 0;
    };
}
/*--------------------------------------------
 
 --------------------------------------------*/
void MDNVT::getset_L(PyGetSetDef& getset)
{
    getset.name=(char*)"L";
    getset.doc=(char*)R"---(
    (int) length of particle NHC
    
    Number of links in NHC of particle thermostat, see :ref:`nhc-ref`
    )---";
    getset.get=[](PyObject* self,void*)->PyObject*
    {
        int L=reinterpret_cast<Object*>(self)->md->thermo_part.L;
        return var<int>::build(L);
    };
    getset.set=[](PyObject* self,PyObject* op,void*)->int
    {
        VarAPI<int> L("L");
        L.logics[0]=VLogics("gt",0);
        int ichk=L.set(op);
        if(ichk==-1) return -1;
        if(reinterpret_cast<Object*>(self)->md->thermo_part.L==L.val)
            return 0;
        
        ThermostatNHC& thermo_part=reinterpret_cast<Object*>(self)->md->thermo_part;
        int niters=thermo_part.niters;
        type0 t_relax=thermo_part.t_relax;
        type0 __dt=reinterpret_cast<Object*>(self)->md->dt2;
        thermo_part.~ThermostatNHC();
        new (&thermo_part) ThermostatNHC(__dt,t_relax,L.val,niters);
        return 0;
    };
}
/*--------------------------------------------
 
 --------------------------------------------*/
void MDNVT::getset_niters(PyGetSetDef& getset)
{
    getset.name=(char*)"niters";
    getset.doc=(char*)R"---(
    (int) number of iterations in NHC for particles
    
    Number of iterations in particle thermostat per each step of evolution, see :ref:`nhc-ref`
    )---";
    getset.get=[](PyObject* self,void*)->PyObject*
    {
        int niters=reinterpret_cast<Object*>(self)->md->thermo_part.niters;
        return var<int>::build(niters);
    };
    getset.set=[](PyObject* self,PyObject* op,void*)->int
    {
        VarAPI<int> niters("niters");
        niters.logics[0]=VLogics("gt",0);
        int ichk=niters.set(op);
        if(ichk==-1) return -1;
        if(reinterpret_cast<Object*>(self)->md->thermo_part.niters==niters.val)
            return 0;
        
        ThermostatNHC& thermo_part=reinterpret_cast<Object*>(self)->md->thermo_part;
        int L=thermo_part.L;
        type0 t_relax=thermo_part.t_relax;
        type0 __dt=reinterpret_cast<Object*>(self)->md->dt2;
        thermo_part.~ThermostatNHC();
        new (&thermo_part) ThermostatNHC(__dt,t_relax,L,niters.val);
        return 0;
    };
}
/*--------------------------------------------
 
 --------------------------------------------*/
void MDNVT::getset_t_relax(PyGetSetDef& getset)
{
    getset.name=(char*)"t_relax";
    getset.doc=(char*)"thermostat parameter, relaxation time";
    getset.doc=(char*)R"---(
    (double) NHC relaxation parameter for particles
    
    NHC relaxation parameter (unit of time) for partiacle thermostat. Roughly equivalent to relaxation time of particles with thermal bath, see :ref:`nhc-ref`
    )---";
    getset.get=[](PyObject* self,void*)->PyObject*
    {
        type0 t_relax=reinterpret_cast<Object*>(self)->md->thermo_part.t_relax;
        return var<type0>::build(t_relax);
    };
    getset.set=[](PyObject* self,PyObject* op,void*)->int
    {
        VarAPI<type0> t_relax("t_relax");
        t_relax.logics[0]=VLogics("gt",0.0);
        int ichk=t_relax.set(op);
        if(ichk==-1) return -1;
        if(reinterpret_cast<Object*>(self)->md->thermo_part.t_relax==t_relax.val)
            return 0;
        
        ThermostatNHC& thermo_part=reinterpret_cast<Object*>(self)->md->thermo_part;
        int L=thermo_part.L;
        int niters=thermo_part.niters;
        type0 __dt=reinterpret_cast<Object*>(self)->md->dt2;
        thermo_part.~ThermostatNHC();
        new (&thermo_part) ThermostatNHC(__dt,t_relax.val,L,niters);
        return 0;
    };
}
/*--------------------------------------------
 
 --------------------------------------------*/
void MDNVT::getset_export(PyGetSetDef& getset)
{
    getset.name=(char*)"export";
    getset.doc=(char*)R"---(
    (mapp.md.export) export object
    
    Export object to record the snapshots of the system while running
    )---";
    getset.get=[](PyObject* self,void*)->PyObject*
    {
        ExportMD::Object* xprt=reinterpret_cast<Object*>(self)->xprt;
        if(!xprt) Py_RETURN_NONE;
        Py_INCREF(xprt);
        return reinterpret_cast<PyObject*>(xprt);
    };
    getset.set=[](PyObject* self,PyObject* op,void*)->int
    {
        VarAPI<OP<ExportMD>> xprt("export");
        int ichk=xprt.set(op);
        if(ichk==-1) return -1;
        if(reinterpret_cast<Object*>(self)->xprt) Py_DECREF(reinterpret_cast<Object*>(self)->xprt);
        Py_INCREF(xprt.val.ob);
        reinterpret_cast<Object*>(self)->xprt=reinterpret_cast<ExportMD::Object*>(xprt.val.ob);
        return 0;
    };
}
/*--------------------------------------------
 
 --------------------------------------------*/
void MDNVT::getset_ntally(PyGetSetDef& getset)
{
    getset.name=(char*)"ntally";
    getset.doc=(char*)R"---(
    (int) thermodynamic tallying period
    
    Number of steps to be taken from one thermodynamics output to the next.
    )---";
    getset.get=[](PyObject* self,void*)->PyObject*
    {
        return var<int>::build(reinterpret_cast<Object*>(self)->md->ntally);
    };
    getset.set=[](PyObject* self,PyObject* op,void*)->int
    {
        VarAPI<int> ntally("ntally");
        ntally.logics[0]=VLogics("ge",0);
        int ichk=ntally.set(op);
        if(ichk==-1) return -1;
        reinterpret_cast<Object*>(self)->md->ntally=ntally.val;
        return 0;
    };
}
/*--------------------------------------------
 
 --------------------------------------------*/
void MDNVT::ml_run(PyMethodDef& tp_methods)
{
    tp_methods.ml_flags=METH_VARARGS | METH_KEYWORDS;
    tp_methods.ml_name="run";
    tp_methods.ml_meth=(PyCFunction)(PyCFunctionWithKeywords)(
    [](PyObject* self,PyObject* args,PyObject* kwds)->PyObject*
    {
        Object* __self=reinterpret_cast<Object*>(self);
        FuncAPI<OP<AtomsMD>,int> f("run",{"atoms","nsteps"});
        f.logics<1>()[0]=VLogics("ge",0);
        if(f(args,kwds)) return NULL;
        
        AtomsMD* __atoms=reinterpret_cast<AtomsMD::Object*>(f.val<0>().ob)->atoms;
        ForceFieldMD* __ff=reinterpret_cast<AtomsMD::Object*>(f.val<0>().ob)->ff;
        ExportMD* __xprt=__self->xprt==NULL ? NULL:__self->xprt->xprt;
        
        try
        {
            __self->md->pre_run_chk(__atoms,__ff);
        }
        catch(std::string& err_msg)
        {
            PyErr_SetString(PyExc_TypeError,err_msg.c_str());
            return NULL;
        }
        
        __self->md->atoms=__atoms;
        __self->md->ff=__ff;
        __self->md->xprt=__xprt;
        
        try
        {
            __self->md->init();
        }
        catch(std::string& err_msg)
        {
            __self->md->xprt=NULL;
            __self->md->ff=NULL;
            __self->md->atoms=NULL;
            PyErr_SetString(PyExc_TypeError,err_msg.c_str());
            return NULL;
        }
        
        __self->md->run(f.val<1>());
        
        __self->md->fin();
        
        __self->md->xprt=NULL;
        __self->md->ff=NULL;
        __self->md->atoms=NULL;
        
        
        Py_RETURN_NONE;
    });
    
    tp_methods.ml_doc=(char*)R"---(
    run(atoms,nsteps)
   
    Execute molecular dynamics
    
    This method starts the molecular dynamics for a given atoms object and number of steps.
    
    Parameters
    ----------
    atoms : mapp.md.atoms
        System of interest
    nsteps : int
        Number of steps to perform molecular dynamics
        
    Returns
    -------
    None

    )---";
}

















