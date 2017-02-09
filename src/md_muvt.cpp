#include "md_muvt.h"
#include "thermo_dynamics.h"
#include "dynamic_md.h"
#include "elements.h"
#include "ff_md.h"
#include "pgcmc.h"
#include "xmath.h"
/*--------------------------------------------
 
 --------------------------------------------*/
MDMuVT::MDMuVT(AtomsMD*& __atoms,ForceFieldMD*& __ff,
type0 __mu,type0 __T,type0 __dt,elem_type __gas_elem,int __seed):
MDNVT(__atoms,__ff,__T,__dt),
seed(__seed),
mu(__mu),
nevery(1000),
nattempts(1000),
gas_elem(__gas_elem)
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
    type0* m=atoms->elements->masses;
    type0 m_i;
    
    const int natms=atoms->natms;
    for(int i=0;i<natms;++i)
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
 
 --------------------------------------------*/
void MDMuVT::run(int nsteps)
{
    
    init();
    
    
    
    if(atoms->dof)
        dynamic=new DynamicMD(atoms,ff,false,{atoms->elem},{atoms->x_d,atoms->dof});
    else
        dynamic=new DynamicMD(atoms,ff,false,{atoms->elem},{atoms->x_d});
    
    dynamic->init();
    
    PGCMC gcmc(atoms,ff,dynamic,1,gas_elem,mu,T,seed);
    gcmc.init();
               
    
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
        
        update_x_d__x(fac_x_d);
        
        if((i+1)%nevery)
            dynamic->update(atoms->x);
        else
        {
            gcmc.xchng(false,nattempts);
            ndof_part+=static_cast<type0>(gcmc.dof_diff);
        }
        
        ff->reset();
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
        
        if((i+1)%ntally==0) thermo.print(i+1);
    }
    
    if(nsteps%ntally) thermo.print(nsteps);
    update_x_d_final(fac_x_d);
    thermo.fin();
    
    gcmc.fin();
    
    dynamic->fin();
    delete dynamic;
    dynamic=NULL;
    fin();
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
    FuncAPI<OP<AtomsMD>,type0,type0,type0,std::string,int> f("__init__",{"atoms_md","mu","T","dt","gas_element","seed"});
    
    f.logics<2>()[0]=VLogics("gt",0.0);
    f.logics<3>()[0]=VLogics("gt",0.0);
    f.logics<5>()[0]=VLogics("gt",0);
    if(f(args,kwds)==-1) return -1;
    
    Object* __self=reinterpret_cast<Object*>(self);
    AtomsMD::Object* atoms_md=reinterpret_cast<AtomsMD::Object*>(f.val<0>().ob);
    
    elem_type ielem;
    try
    {
        ielem=atoms_md->atoms->elements->find(f.val<4>().c_str());
    }
    catch(const char* err_msg)
    {
        PyErr_SetString(PyExc_TypeError,err_msg);
        return -1;
    }
    
    
    __self->md=new MDMuVT(atoms_md->atoms,atoms_md->ff,f.val<1>(),f.val<2>(),f.val<3>(),ielem,f.val<5>());
    __self->atoms_md=atoms_md;
    Py_INCREF(atoms_md);
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
    __self->atoms_md=NULL;
    return reinterpret_cast<PyObject*>(__self);
}
/*--------------------------------------------
 
 --------------------------------------------*/
void MDMuVT::__dealloc__(PyObject* self)
{
    Object* __self=reinterpret_cast<Object*>(self);
    delete __self->md;
    __self->md=NULL;
    if(__self->atoms_md) Py_DECREF(__self->atoms_md);
    __self->atoms_md=NULL;
    delete __self;
}
/*--------------------------------------------*/
PyTypeObject MDMuVT::TypeObject ={PyObject_HEAD_INIT(NULL)};
/*--------------------------------------------*/
void MDMuVT::setup_tp()
{
    TypeObject.tp_name="md_muvt";
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
        std::string gas_element(reinterpret_cast<Object*>(self)->md->atoms->elements->names[gas_elem]);
        
        return var<std::string>::build(gas_element,NULL);
    };
    getset.set=[](PyObject* self,PyObject* op,void*)->int
    {
        VarAPI<std::string> gas_element("gas_element");
        int ichk=gas_element.set(op);
        if(ichk==-1) return -1;
        elem_type ielem;
        try
        {
            ielem=reinterpret_cast<Object*>(self)->md->atoms->elements->find(gas_element.val.c_str());
        }
        catch(const char* err_msg)
        {
            PyErr_SetString(PyExc_TypeError,err_msg);
            return -1;
        }
        
        
        reinterpret_cast<Object*>(self)->md->gas_elem=ielem;
        return 0;
    };
}







