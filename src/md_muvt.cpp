#include "md_muvt.h"
#include "thermo_dynamics.h"
#include "dynamic_md.h"
#include "elements.h"
#include "ff_md.h"
#include "pgcmc.h"
#include "xmath.h"
/*--------------------------------------------
 
 --------------------------------------------*/
MDMuVT::MDMuVT(type0 __mu,type0 __T,type0 __dt,std::string __gas_elem_name,int __seed):
MDNVT(__T,__dt),
seed(__seed),
mu(__mu),
nevery(1000),
nattempts(1000),
gas_elem_name(__gas_elem_name)
{
}
/*--------------------------------------------
 
 --------------------------------------------*/
MDMuVT::~MDMuVT()
{
}
/*--------------------------------------------
 
 --------------------------------------------*/
void MDMuVT::update_x_d__x(type0 fac_x_d)
{
    type0* x=atoms->x_d->begin();
    type0* f=ff->f->begin();
    type0* x_d=atoms->x_d->begin();
    elem_type* elem=atoms->elem->begin();
    type0* m=atoms->elements.masses;
    type0 m_i;
    
    const int natms_lcl=atoms->natms_lcl;
    for(int i=0;i<natms_lcl;++i)
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
}
/*--------------------------------------------
 pre run check it throw excepctions
 --------------------------------------------*/
void MDMuVT::pre_run_chk(AtomsMD* atoms,ForceFieldMD* ff)
{
    //check if configuration is loaded
    if(!atoms)
        throw std::string("cannot start md without initial conditions");
    
    //check if force field is loaded
    if(!ff)
        throw std::string("cannot start md without governing equations (force field)");
    
    elem_type ielem;
    try
    {
        ielem=atoms->elements.find(gas_elem_name.c_str());
    }
    catch(int)
    {
        throw "atom "+gas_elem_name+" is not assigned to configuration";
    }
    
}
/*--------------------------------------------
 
 --------------------------------------------*/
void MDMuVT::pre_init()
{
    MDNVT::pre_init();
}
/*--------------------------------------------
 
 --------------------------------------------*/
void MDMuVT::init()
{
    pre_init();
    
    dynamic=new DynamicMD(atoms,ff,false,{},{atoms->x_d,atoms->dof},{});
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
void MDMuVT::fin()
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
void MDMuVT::run(int nsteps)
{
    int step=atoms->step;
    
    PGCMC gcmc(atoms,ff,dynamic,1,gas_elem,mu,T,seed);
    gcmc.init();
               
    
    ff->force_calc_timer();
    
    int nevery_xprt=xprt==NULL ? 0:xprt->nevery;
    if(nevery_xprt) xprt->write(step);
    
    ThermoDynamics thermo(6,"T",T_part,"PE",ff->nrgy_strss[0],
    "S[0][0]",S_part[0][0],
    "S[1][1]",S_part[1][1],
    "S[2][2]",S_part[2][2],
    "S[1][2]",S_part[2][1],
    "S[2][0]",S_part[2][0],
    "S[0][1]",S_part[1][0]);
    
    if(ntally) thermo.init();
    Algebra::DoLT<__dim__>::func([this](const int i,const int j)
    {
        S_part[i][j]=ff->nrgy_strss[1+i+j*__dim__-j*(j+1)/2]-mvv[i+j*__dim__-j*(j+1)/2]/atoms->vol;
    });
    
    
    if(ntally) thermo.print(step);
    
    type0 fac,fac_x_d=1.0;
    
    for(int istep=0;istep<nsteps;istep++)
    {
        // particle thermostat
        fac_x_d*=fac=thermo_part(T_part/T,ndof_part);
        fac*=fac;
        Algebra::Do<__nvoigt__>::func([this,&fac](int i){mvv[i]*=fac;});
        T_part*=fac;
        
        update_x_d__x(fac_x_d);
        
        if((istep+1)%nevery)
            dynamic->update(atoms->x);
        else
        {
            gcmc.xchng(false,nattempts);
            ndof_part+=static_cast<type0>(gcmc.dof_diff);
        }
        
        ff->force_calc_timer();
        
        update_x_d();
        
        // particle thermostat
        fac_x_d=fac=thermo_part(T_part/T,ndof_part);
        fac*=fac;
        Algebra::Do<__nvoigt__>::func([this,&fac](int i){mvv[i]*=fac;});
        T_part*=fac;
        
        
        Algebra::DoLT<__dim__>::func([this](const int i,const int j)
        {
            S_part[i][j]=ff->nrgy_strss[1+i+j*__dim__-j*(j+1)/2]-mvv[i+j*__dim__-j*(j+1)/2]/atoms->vol;
        });
        
        if(ntally && (istep+1)%ntally==0) thermo.print(step+istep+1);
        if(nevery_xprt && (istep+1)%nevery_xprt==0) xprt->write(step+istep+1);
    }
    
    if(ntally && nsteps%ntally) thermo.print(step+nsteps);
    if(nevery_xprt && nsteps%nevery_xprt) xprt->write(step+nsteps);
    
    update_x_d_final(fac_x_d);
    if(ntally) thermo.fin();
    
    gcmc.fin();
    
    atoms->step+=nsteps;
}
/*------------------------------------------------------------------------------------------------------------------------------------
 
 ------------------------------------------------------------------------------------------------------------------------------------*/
PyObject* MDMuVT::__new__(PyTypeObject* type,PyObject* args,PyObject* kwds)
{
    Object* __self=reinterpret_cast<Object*>(type->tp_alloc(type,0));
    PyObject* self=reinterpret_cast<PyObject*>(__self);
    return self;
}
/*--------------------------------------------
 
 --------------------------------------------*/
int MDMuVT::__init__(PyObject* self,PyObject* args,PyObject* kwds)
{
    FuncAPI<type0,type0,type0,std::string,int> f("__init__",{"mu","T","dt","gas_element","seed"});
    
    f.logics<1>()[0]=VLogics("gt",0.0);
    f.logics<2>()[0]=VLogics("gt",0.0);
    f.logics<4>()[0]=VLogics("gt",0);
    if(f(args,kwds)==-1) return -1;
    
    Object* __self=reinterpret_cast<Object*>(self);
    __self->md=new MDMuVT(f.val<0>(),f.val<1>(),f.val<2>(),f.val<3>(),f.val<4>());
    __self->xprt=NULL;
    
    return 0;
}
/*--------------------------------------------
 
 --------------------------------------------*/
PyObject* MDMuVT::__alloc__(PyTypeObject* type,Py_ssize_t)
{
    Object* __self=new Object;
    __self->ob_type=type;
    __self->ob_refcnt=1;
    __self->md=NULL;
        __self->xprt=NULL;
    return reinterpret_cast<PyObject*>(__self);
}
/*--------------------------------------------
 
 --------------------------------------------*/
void MDMuVT::__dealloc__(PyObject* self)
{
    Object* __self=reinterpret_cast<Object*>(self);
    delete __self->md;
    __self->md=NULL;
    if(__self->xprt) Py_DECREF(__self->xprt);
    __self->xprt=NULL;
    delete __self;
}
/*--------------------------------------------*/
PyTypeObject MDMuVT::TypeObject ={PyObject_HEAD_INIT(NULL)};
/*--------------------------------------------*/
void MDMuVT::setup_tp()
{
    TypeObject.tp_name="mapp.md.muvt";
    TypeObject.tp_doc="MD of canonical ensemble";
    
    TypeObject.tp_flags=Py_TPFLAGS_DEFAULT;
    TypeObject.tp_basicsize=sizeof(Object);
    
    TypeObject.tp_new=__new__;
    TypeObject.tp_init=__init__;
    TypeObject.tp_alloc=__alloc__;
    TypeObject.tp_dealloc=__dealloc__;
    setup_tp_getset();
    TypeObject.tp_getset=getset;
    
    TypeObject.tp_base=&MDNVT::TypeObject;
}
/*--------------------------------------------*/
PyGetSetDef MDMuVT::getset[]={[0 ... 4]={NULL,NULL,NULL,NULL,NULL}};
/*--------------------------------------------*/
void MDMuVT::setup_tp_getset()
{
    getset_nevery(getset[0]);
    getset_nattempts(getset[1]);
    getset_seed(getset[2]);
    getset_gas_element(getset[3]);
}
/*--------------------------------------------
 
 --------------------------------------------*/
void MDMuVT::getset_nevery(PyGetSetDef& getset)
{
    getset.name=(char*)"nevery";
    getset.doc=(char*)"perform deletion/insertion attempts every nevery steps";
    getset.get=[](PyObject* self,void*)->PyObject*
    {
        return var<int>::build(reinterpret_cast<Object*>(self)->md->nevery,NULL);
    };
    getset.set=[](PyObject* self,PyObject* op,void*)->int
    {
        VarAPI<int> nevery("nevery");
        nevery.logics[0]=VLogics("gt",0);
        int ichk=nevery.set(op);
        if(ichk==-1) return -1;
        
        reinterpret_cast<Object*>(self)->md->nevery=nevery.val;
        return 0;
    };
}
/*--------------------------------------------
 
 --------------------------------------------*/
void MDMuVT::getset_nattempts(PyGetSetDef& getset)
{
    getset.name=(char*)"nattempts";
    getset.doc=(char*)"number of deletion/insertion attempts";
    getset.get=[](PyObject* self,void*)->PyObject*
    {
        return var<int>::build(reinterpret_cast<Object*>(self)->md->nattempts,NULL);
    };
    getset.set=[](PyObject* self,PyObject* op,void*)->int
    {
        VarAPI<int> nattempts("nattempts");
        nattempts.logics[0]=VLogics("gt",0);
        int ichk=nattempts.set(op);
        if(ichk==-1) return -1;
        
        reinterpret_cast<Object*>(self)->md->nattempts=nattempts.val;
        return 0;
    };
}
/*--------------------------------------------
 
 --------------------------------------------*/
void MDMuVT::getset_seed(PyGetSetDef& getset)
{
    getset.name=(char*)"seed";
    getset.doc=(char*)"gcmc random seed";
    getset.get=[](PyObject* self,void*)->PyObject*
    {
        return var<int>::build(reinterpret_cast<Object*>(self)->md->seed,NULL);
    };
    getset.set=[](PyObject* self,PyObject* op,void*)->int
    {
        VarAPI<int> seed("seed");
        seed.logics[0]=VLogics("gt",0);
        int ichk=seed.set(op);
        if(ichk==-1) return -1;
        
        reinterpret_cast<Object*>(self)->md->seed=seed.val;
        return 0;
    };
}
/*--------------------------------------------
 
 --------------------------------------------*/
void MDMuVT::getset_gas_element(PyGetSetDef& getset)
{
    getset.name=(char*)"gas_element";
    getset.doc=(char*)"gcmc element";
    getset.get=[](PyObject* self,void*)->PyObject*
    {
        elem_type gas_elem=reinterpret_cast<Object*>(self)->md->gas_elem;
        std::string gas_element(reinterpret_cast<Object*>(self)->md->atoms->elements.names[gas_elem]);
        
        return var<std::string>::build(gas_element,NULL);
    };
    getset.set=[](PyObject* self,PyObject* op,void*)->int
    {
        VarAPI<std::string> gas_element("gas_element");
        int ichk=gas_element.set(op);
        if(ichk==-1) return -1;
        reinterpret_cast<Object*>(self)->md->gas_elem_name=gas_element.val;
        return 0;
    };
}







