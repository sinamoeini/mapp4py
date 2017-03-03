#include "dmd.h"
#include "atoms_dmd.h"
#include "dynamic_dmd.h"
using namespace MAPP_NS;
/*--------------------------------------------
 
 --------------------------------------------*/
DMD::DMD():
max_nsteps(1000),
a_tol(sqrt(std::numeric_limits<type0>::epsilon())),
min_dt(std::numeric_limits<type0>::epsilon()),
c(NULL),
c_d(NULL)
{
}
/*--------------------------------------------
 
 --------------------------------------------*/
DMD::~DMD()
{
}
/*--------------------------------------------
 
 --------------------------------------------*/
void DMD::init_static()
{
    c_dim=atoms->c_dim;
    ncs=atoms->natms*c_dim;
    // for now
    nc_dofs=static_cast<type0>(ncs);
    if(!atoms->c_d) atoms->c_d=new Vec<type0>(atoms,c_dim);
    dynamic=new DynamicDMD(atoms,ff,false,{atoms->elem,atoms->c},{},{});
    dynamic->init();
    c=atoms->c->begin();
    c_d=atoms->c_d->begin();
}
/*--------------------------------------------
 
 --------------------------------------------*/
void DMD::fin_static()
{
    dynamic->fin();
    delete dynamic;
    dynamic=NULL;
    c=c_d=NULL;
}
/*------------------------------------------------------------------------------------------------------------------------------------
 
 ------------------------------------------------------------------------------------------------------------------------------------*/
PyObject* DMD::__new__(PyTypeObject* type,PyObject* args,PyObject* kwds)
{
    Object* __self=reinterpret_cast<Object*>(type->tp_alloc(type,0));
    PyObject* self=reinterpret_cast<PyObject*>(__self);
    return self;
}
/*--------------------------------------------
 
 --------------------------------------------*/
int DMD::__init__(PyObject* self,PyObject* args,PyObject* kwds)
{
    FuncAPI<> f("__init__");

    
    if(f(args,kwds)==-1) return -1;
    Object* __self=reinterpret_cast<Object*>(self);
    __self->dmd=new DMD();
    return 0;
}
/*--------------------------------------------
 
 --------------------------------------------*/
PyObject* DMD::__alloc__(PyTypeObject* type,Py_ssize_t)
{
    Object* __self=new Object;
    __self->ob_type=type;
    __self->ob_refcnt=1;
    __self->dmd=NULL;
    return reinterpret_cast<PyObject*>(__self);
}
/*--------------------------------------------
 
 --------------------------------------------*/
void DMD::__dealloc__(PyObject* self)
{
    Object* __self=reinterpret_cast<Object*>(self);
    delete __self->dmd;
    __self->dmd=NULL;
    delete __self;
}
/*--------------------------------------------*/
PyMethodDef DMD::methods[]={[0 ... 0]={NULL}};
/*--------------------------------------------*/
void DMD::setup_tp_methods()
{
}
/*--------------------------------------------*/
PyTypeObject DMD::TypeObject={PyObject_HEAD_INIT(NULL)};
/*--------------------------------------------*/
void DMD::setup_tp()
{
    TypeObject.tp_name="mapp.dmd.dmd";
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
PyGetSetDef DMD::getset[]={[0 ... 3]={NULL,NULL,NULL,NULL,NULL}};
/*--------------------------------------------*/
void DMD::setup_tp_getset()
{
    getset_a_tol(getset[0]);
    getset_max_nsteps(getset[1]);
    getset_min_dt(getset[2]);
}
/*--------------------------------------------
 
 --------------------------------------------*/
void DMD::getset_a_tol(PyGetSetDef& getset)
{
    getset.name=(char*)"a_tol";
    getset.doc=(char*)"absolute error tolerence in LTE";
    getset.get=[](PyObject* self,void*)->PyObject*
    {
        return var<type0>::build(reinterpret_cast<Object*>(self)->dmd->a_tol,NULL);
    };
    getset.set=[](PyObject* self,PyObject* op,void*)->int
    {
        VarAPI<type0> a_tol("a_tol");
        a_tol.logics[0]=VLogics("gt",0.0)*VLogics("lt",1.0);
        int ichk=a_tol.set(op);
        if(ichk==-1) return -1;
        reinterpret_cast<Object*>(self)->dmd->a_tol=a_tol.val;
        return 0;
    };
}
/*--------------------------------------------
 
 --------------------------------------------*/
void DMD::getset_max_nsteps(PyGetSetDef& getset)
{
    getset.name=(char*)"max_nsteps";
    getset.doc=(char*)"maximum number of steps";
    getset.get=[](PyObject* self,void*)->PyObject*
    {
        return var<int>::build(reinterpret_cast<Object*>(self)->dmd->max_nsteps,NULL);
    };
    getset.set=[](PyObject* self,PyObject* op,void*)->int
    {
        VarAPI<int> max_nsteps("max_nsteps");
        max_nsteps.logics[0]=VLogics("ge",0);
        int ichk=max_nsteps.set(op);
        if(ichk==-1) return -1;
        reinterpret_cast<Object*>(self)->dmd->max_nsteps=max_nsteps.val;
        return 0;
    };
}
/*--------------------------------------------
 
 --------------------------------------------*/
void DMD::getset_min_dt(PyGetSetDef& getset)
{
    getset.name=(char*)"min_dt";
    getset.doc=(char*)"minimum time step";
    getset.get=[](PyObject* self,void*)->PyObject*
    {
        return var<type0>::build(reinterpret_cast<Object*>(self)->dmd->min_dt,NULL);
    };
    getset.set=[](PyObject* self,PyObject* op,void*)->int
    {
        VarAPI<type0> min_dt("min_dt");
        min_dt.logics[0]=VLogics("gt",0.0);
        int ichk=min_dt.set(op);
        if(ichk==-1) return -1;
        reinterpret_cast<Object*>(self)->dmd->min_dt=min_dt.val;
        return 0;
    };
}














