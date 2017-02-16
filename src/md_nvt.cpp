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
MDNVT::MDNVT(AtomsMD*& __atoms,ForceFieldMD*& __ff,
type0 __T,type0 __dt):
dynamic(NULL),
atoms(__atoms),
ff(__ff),
world(__atoms->comm.world),
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
    
    int nchains=thermo_part.nchains;
    int niters=thermo_part.niters;
    type0 t_relax=thermo_part.t_relax;
    
    thermo_part.~ThermostatNHC();
    new (&thermo_part) ThermostatNHC(__dt/2.0,t_relax,nchains,niters);
}
/*--------------------------------------------
 
 --------------------------------------------*/
void MDNVT::update_x()
{
    type0* x=atoms->x->begin();
    type0* x_d=atoms->x_d->begin();
    const int natms=atoms->natms;
    for(int i=0;i<natms;++i)
    {
        Algebra::Do<__dim__>::func(
        [&x,&x_d,this](const int j){x[j]+=x_d[j]*dt;});
        x+=__dim__;
        x_d+=__dim__;
    }
}
/*--------------------------------------------
 
 --------------------------------------------*/
void MDNVT::update_x_d()
{
    
    type0* f=ff->f->begin();
    type0* x_d=atoms->x_d->begin();
    elem_type* elem=atoms->elem->begin();
    type0* m=atoms->elements->masses;
    type0 m_i;
    Algebra::zero(__vec_lcl);
    const int natms=atoms->natms;
    if(atoms->dof)
    {
        bool* dof=atoms->dof->begin();
        for(int i=0;i<natms;++i)
        {
            m_i=m[*elem];
            Algebra::Do<__dim__>::func([&dof,&x_d,&f,&m_i,this](const int j){if(dof[j])x_d[j]+=f[j]*dt2/m_i;});
            Algebra::DyadicV<__dim__>(m_i,x_d,__vec_lcl);
            
            dof+=__dim__;
            f+=__dim__;
            x_d+=__dim__;
            ++elem;
        }
    }
    else
        for(int i=0;i<natms;++i)
        {
            m_i=m[*elem];
            Algebra::Do<__dim__>::func([&x_d,&f,&m_i,this](const int j){x_d[j]+=f[j]*dt2/m_i;});
            Algebra::DyadicV<__dim__>(m_i,x_d,__vec_lcl);
            f+=__dim__;
            x_d+=__dim__;
            ++elem;
        }
    
    MPI_Allreduce(__vec_lcl,mvv,__nvoigt__,Vec<type0>::MPI_T,MPI_SUM,world);
    T_part=Algebra::Tr_DyadicV<__dim__>(mvv)/(ndof_part*kB);
}
/*--------------------------------------------
 
 --------------------------------------------*/
void MDNVT::update_x_d_final(type0 fac_x_d)
{
    const int n=atoms->natms*__dim__;
    type0* x_d=atoms->x_d->begin();
    for(int i=0;i<n;i++) x_d[i]*=fac_x_d;
}
/*--------------------------------------------
 
 --------------------------------------------*/
void MDNVT::update_x_d(type0 fac_x_d)
{
    
    type0* f=ff->f->begin();
    type0* x_d=atoms->x_d->begin();
    elem_type* elem=atoms->elem->begin();
    type0* m=atoms->elements->masses;
    type0 m_i;
    Algebra::zero(__vec_lcl);
    const int natms=atoms->natms;
    if(atoms->dof)
    {
        bool* dof=atoms->dof->begin();
        for(int i=0;i<natms;++i)
        {
            m_i=m[*elem];
            Algebra::Do<__dim__>::func([&dof,&x_d,&f,&m_i,&fac_x_d,this](const int j){if(dof[j])x_d[j]=x_d[j]*fac_x_d+f[j]*dt2/m_i;});
            Algebra::DyadicV<__dim__>(m_i,x_d,__vec_lcl);
            
            dof+=__dim__;
            f+=__dim__;
            x_d+=__dim__;
            ++elem;
        }
    }
    else
        for(int i=0;i<natms;++i)
        {
            m_i=m[*elem];
            Algebra::Do<__dim__>::func([&x_d,&f,&m_i,&fac_x_d,this](const int j){x_d[j]=x_d[j]*fac_x_d+f[j]*dt2/m_i;});
            Algebra::DyadicV<__dim__>(m_i,x_d,__vec_lcl);
            f+=__dim__;
            x_d+=__dim__;
            ++elem;
        }
    
    MPI_Allreduce(__vec_lcl,mvv,__nvoigt__,Vec<type0>::MPI_T,MPI_SUM,world);
    T_part=Algebra::Tr_DyadicV<__dim__>(mvv)/(ndof_part*kB);
}
/*--------------------------------------------
 
 --------------------------------------------*/
void MDNVT::update_x_d__x__x_d(type0 fac_x_d)
{
    type0* x=atoms->x->begin();
    type0* f=ff->f->begin();
    type0* x_d=atoms->x_d->begin();
    elem_type* elem=atoms->elem->begin();
    type0* m=atoms->elements->masses;
    type0 m_i;
    const int natms0=atoms->natms;
    for(int i=0;i<natms0;++i)
    {
        m_i=m[*elem];
        Algebra::Do<__dim__>::func([&x_d,&x,&f,&m_i,&fac_x_d,this](const int j)
        {
            x_d[j]=x_d[j]*fac_x_d+f[j]*dt2/m_i;
            x[j]+=x_d[j]*dt;
        });
        
        f+=__dim__;
        x_d+=__dim__;
        x+=__dim__;
        ++elem;
    }
    
    
    dynamic->update(atoms->x);
    ff->reset();
    ff->force_calc_timer();
    
    f=ff->f->begin();
    x_d=atoms->x_d->begin();
    elem=atoms->elem->begin();
    Algebra::zero(__vec_lcl);
    const int natms1=atoms->natms;
    for(int i=0;i<natms1;++i)
    {
        m_i=m[*elem];
        Algebra::Do<__dim__>::func([&x_d,&f,&m_i,this](const int j){x_d[j]+=f[j]*dt2/m_i;});
        Algebra::DyadicV<__dim__>(m_i,x_d,__vec_lcl);
        
        f+=__dim__;
        x_d+=__dim__;
        ++elem;
    }
    MPI_Allreduce(__vec_lcl,mvv,__nvoigt__,Vec<type0>::MPI_T,MPI_SUM,world);
    T_part=Algebra::Tr_DyadicV<__dim__>(mvv)/(ndof_part*kB);
}
/*--------------------------------------------
 
 --------------------------------------------*/
void MDNVT::init()
{
    kB=atoms->kB;
    /*
     calculating the number of degress of freedom
     */
    ndof_part=atoms->tot_natms*__dim__;
    if(atoms->dof)
    {
        bool* dof=atoms->dof->begin();
        int ndof_red_lcl=0;
        const int n=atoms->natms*__dim__;
        for(int i=0;i<n;i++)
            if(!dof[i]) ndof_red_lcl++;
        int ndof_red;
        MPI_Allreduce(&ndof_red_lcl,&ndof_red,1,MPI_INT,MPI_SUM,world);
        ndof_part-=ndof_red;
    }
    
    /*
     temperature and kinetic energy tensor
     */
    if(atoms->x_d)
    {
        type0* x_d=atoms->x_d->begin();
        type0* m=atoms->elements->masses;
        type0 m_i;
        elem_type* elem=atoms->elem->begin();
        Algebra::zero(__vec_lcl);
        const int natms=atoms->natms;
        for(int i=0;i<natms;i++)
        {
            m_i=m[*elem];
            Algebra::DyadicV<__dim__>(m_i,x_d,__vec_lcl);
            x_d+=__dim__;
            ++elem;
        }
        
        MPI_Allreduce(__vec_lcl,mvv,__nvoigt__,Vec<type0>::MPI_T,MPI_SUM,world);
        T_part=Algebra::Tr_DyadicV<__dim__>(mvv)/(ndof_part*kB);
    }
    else
    {
        atoms->x_d=new Vec<type0>(atoms,__dim__,"x_d");
        type0* x_d=atoms->x_d->begin();
        const int n=atoms->natms*__dim__;
        for(int i=0;i<n;i++) x_d[i]=0.0;
        Algebra::Do<__nvoigt__>::func([this](int i){mvv[i]=0.0;});
        T_part=0.0;
    }

}
/*--------------------------------------------
 
 --------------------------------------------*/
void MDNVT::fin()
{
}
/*--------------------------------------------
 
 --------------------------------------------*/
void MDNVT::run(int nsteps)
{    
    init();
    
    
    if(atoms->dof)
        dynamic=new DynamicMD(atoms,ff,false,{atoms->elem},{atoms->x_d,atoms->dof});
    else
        dynamic=new DynamicMD(atoms,ff,false,{atoms->elem},{atoms->x_d});
    
    dynamic->init();
    
    
    ff->reset();
    ff->force_calc_timer();
    ThermoDynamics thermo(6,"T",T_part,"PE",ff->nrgy_strss[0],
    "S[0][0]",S_part[0][0],
    "S[1][1]",S_part[1][1],
    "S[2][2]",S_part[2][2],
    "S[1][2]",S_part[2][1],
    "S[2][0]",S_part[2][0],
    "S[0][1]",S_part[1][0]);
    
    
    thermo.init();
    Algebra::DoLT<__dim__>::func([this](const int i,const int j)
    {
        S_part[i][j]=ff->nrgy_strss[1+i+j*__dim__-j*(j+1)/2]-mvv[i+j*__dim__-j*(j+1)/2]/atoms->vol;
    });
    
    
    thermo.print(0);
    
    type0 fac,fac_x_d=1.0;
    for(int i=0;i<nsteps;i++)
    {
        // particle thermostat
        fac_x_d*=fac=thermo_part(T_part/T,ndof_part);
        fac*=fac;
        Algebra::Do<__nvoigt__>::func([this,&fac](int i){mvv[i]*=fac;});
        T_part*=fac;
        
        
        update_x_d__x__x_d(fac_x_d);
        
        
        /*
        update_x_d(fac_x_d);
        update_x();
        update_x_d();
        */
        
        // particle thermostat
        fac_x_d=fac=thermo_part(T_part/T,ndof_part);
        fac*=fac;
        Algebra::Do<__nvoigt__>::func([this,&fac](int i){mvv[i]*=fac;});
        T_part*=fac;
        
        
        Algebra::DoLT<__dim__>::func([this](const int i,const int j)
        {
            S_part[i][j]=ff->nrgy_strss[1+i+j*__dim__-j*(j+1)/2]-mvv[i+j*__dim__-j*(j+1)/2]/atoms->vol;
        });
        
        if((i+1)%ntally==0) thermo.print(i+1);
    }
    
    if(nsteps%ntally) thermo.print(nsteps);
    update_x_d_final(fac_x_d);
    
    thermo.fin();
    
    dynamic->fin();
    delete dynamic;
    dynamic=NULL;
    fin();

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
    FuncAPI<OP<AtomsMD>,type0,type0> f("__init__",{"atoms_md","T","dt"});
    
    f.logics<1>()[0]=VLogics("gt",0.0);
    f.logics<2>()[0]=VLogics("gt",0.0);
    
    if(f(args,kwds)==-1) return -1;
    Object* __self=reinterpret_cast<Object*>(self);
    AtomsMD::Object* atoms_md=reinterpret_cast<AtomsMD::Object*>(f.val<0>().ob);
    __self->md=new MDNVT(atoms_md->atoms,atoms_md->ff,f.val<1>(),f.val<2>());
    __self->atoms_md=atoms_md;
    Py_INCREF(atoms_md);
    return 0;
}
/*--------------------------------------------
 
 --------------------------------------------*/
PyObject* MDNVT::__alloc__(PyTypeObject* type,Py_ssize_t)
{
    Object* __self=new Object;
    __self->ob_type=type;
    __self->ob_refcnt=1;
    __self->md=NULL;
    __self->atoms_md=NULL;
    return reinterpret_cast<PyObject*>(__self);
}
/*--------------------------------------------
 
 --------------------------------------------*/
void MDNVT::__dealloc__(PyObject* self)
{
    Object* __self=reinterpret_cast<Object*>(self);
    delete __self->md;
    __self->md=NULL;
    if(__self->atoms_md) Py_DECREF(__self->atoms_md);
    __self->atoms_md=NULL;
    delete __self;
}
/*--------------------------------------------*/
PyMethodDef MDNVT::methods[]={[0 ... 1]={NULL}};
/*--------------------------------------------*/
void MDNVT::setup_tp_methods()
{
    ml_run(methods[0]);
}
/*--------------------------------------------*/
PyTypeObject MDNVT::TypeObject ={PyObject_HEAD_INIT(NULL)};
/*--------------------------------------------*/
void MDNVT::setup_tp()
{
    TypeObject.tp_name="nvt";
    TypeObject.tp_doc="MD of canonical ensemble";
    
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
PyGetSetDef MDNVT::getset[]={[0 ... 6]={NULL,NULL,NULL,NULL,NULL}};
/*--------------------------------------------*/
void MDNVT::setup_tp_getset()
{
    getset_niters(getset[0]);
    getset_nchains(getset[1]);
    getset_T(getset[2]);
    getset_dt(getset[3]);
    getset_t_relax(getset[4]);
    getset_ntally(getset[5]);
}
/*--------------------------------------------
 
 --------------------------------------------*/
void MDNVT::getset_T(PyGetSetDef& getset)
{
    getset.name=(char*)"T";
    getset.doc=(char*)"external temperature of simulation";
    getset.get=[](PyObject* self,void*)->PyObject*
    {
        return var<type0>::build(reinterpret_cast<Object*>(self)->md->T,NULL);
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
void MDNVT::getset_nchains(PyGetSetDef& getset)
{
    getset.name=(char*)"nchains";
    getset.doc=(char*)"number of chains in Nose-Hoover chain thermostat";
    getset.get=[](PyObject* self,void*)->PyObject*
    {
        int nchains=reinterpret_cast<Object*>(self)->md->thermo_part.nchains;
        return var<int>::build(nchains,NULL);
    };
    getset.set=[](PyObject* self,PyObject* op,void*)->int
    {
        VarAPI<int> nchains("nchains");
        nchains.logics[0]=VLogics("gt",0);
        int ichk=nchains.set(op);
        if(ichk==-1) return -1;
        if(reinterpret_cast<Object*>(self)->md->thermo_part.nchains==nchains.val)
            return 0;
        
        ThermostatNHC& thermo_part=reinterpret_cast<Object*>(self)->md->thermo_part;
        int niters=thermo_part.niters;
        type0 t_relax=thermo_part.t_relax;
        type0 __dt=reinterpret_cast<Object*>(self)->md->dt2;
        thermo_part.~ThermostatNHC();
        new (&thermo_part) ThermostatNHC(__dt,t_relax,nchains.val,niters);
        return 0;
    };
}
/*--------------------------------------------
 
 --------------------------------------------*/
void MDNVT::getset_niters(PyGetSetDef& getset)
{
    getset.name=(char*)"niters";
    getset.doc=(char*)"number of iterations in Nose-Hoover chain thermostat";
    getset.get=[](PyObject* self,void*)->PyObject*
    {
        int niters=reinterpret_cast<Object*>(self)->md->thermo_part.niters;
        return var<int>::build(niters,NULL);
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
        int nchains=thermo_part.nchains;
        type0 t_relax=thermo_part.t_relax;
        type0 __dt=reinterpret_cast<Object*>(self)->md->dt2;
        thermo_part.~ThermostatNHC();
        new (&thermo_part) ThermostatNHC(__dt,t_relax,nchains,niters.val);
        return 0;
    };
}
/*--------------------------------------------
 
 --------------------------------------------*/
void MDNVT::getset_t_relax(PyGetSetDef& getset)
{
    getset.name=(char*)"t_relax";
    getset.doc=(char*)"thermostat parameter, relaxation time";
    getset.get=[](PyObject* self,void*)->PyObject*
    {
        type0 t_relax=reinterpret_cast<Object*>(self)->md->thermo_part.t_relax;
        return var<type0>::build(t_relax,NULL);
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
        int nchains=thermo_part.nchains;
        int niters=thermo_part.niters;
        type0 __dt=reinterpret_cast<Object*>(self)->md->dt2;
        thermo_part.~ThermostatNHC();
        new (&thermo_part) ThermostatNHC(__dt,t_relax.val,nchains,niters);
        return 0;
    };
}
/*--------------------------------------------
 
 --------------------------------------------*/
void MDNVT::getset_dt(PyGetSetDef& getset)
{
    getset.name=(char*)"dt";
    getset.doc=(char*)"timestep of MD simulation";
    getset.get=[](PyObject* self,void*)->PyObject*
    {
        return var<type0>::build(reinterpret_cast<Object*>(self)->md->dt,NULL);
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
void MDNVT::getset_ntally(PyGetSetDef& getset)
{
    getset.name=(char*)"ntally";
    getset.doc=(char*)"tally thermodynamic quantities every ntally steps";
    getset.get=[](PyObject* self,void*)->PyObject*
    {
        return var<int>::build(reinterpret_cast<Object*>(self)->md->ntally,NULL);
    };
    getset.set=[](PyObject* self,PyObject* op,void*)->int
    {
        VarAPI<int> ntally("ntally");
        ntally.logics[0]=VLogics("gt",0);
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
    tp_methods.ml_doc="run simulation for n steps";
    
    tp_methods.ml_meth=(PyCFunction)(PyCFunctionWithKeywords)
    [](PyObject* self,PyObject* args,PyObject* kwds)->PyObject*
    {
        Object* __self=reinterpret_cast<Object*>(self);
        FuncAPI<int> f("run",{"nsteps"});
        f.logics<0>()[0]=VLogics("ge",0);
        if(f(args,kwds)) return NULL;
        
        try
        {
            __self->md->dof_consistency();
        }
        catch (std::string& err_msg)
        {
            PyErr_SetString(PyExc_TypeError,err_msg.c_str());
            return NULL;
        }
        
        __self->md->run(f.val<0>());
        
        Py_RETURN_NONE;
    };
}

















