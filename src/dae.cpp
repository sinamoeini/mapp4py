#include "dae.h"
#include "atoms_dmd.h"
#include "dynamic_dmd.h"
#include "ff_dmd.h"
#include "memory.h"
using namespace MAPP_NS;
/*--------------------------------------------
 
 --------------------------------------------*/
DAE::DAE():
max_nsteps(1000),
a_tol(sqrt(std::numeric_limits<type0>::epsilon())),
min_dt(std::numeric_limits<type0>::epsilon()),
c_init(NULL),
c(NULL),
c_d(NULL),
xprt(NULL),
ncs(0),
ntally(1000),
max_dc_rel(std::numeric_limits<type0>::infinity())
{
}
/*--------------------------------------------
 
 --------------------------------------------*/
DAE::~DAE()
{
}
/*--------------------------------------------
 
 --------------------------------------------*/
void DAE::pre_run_chk(AtomsDMD* __atoms, ForceFieldDMD* __ff)
{
    
}
/*--------------------------------------------
 
 --------------------------------------------*/
int DAE::calc_ndofs(AtomsDMD* __atoms)
{
    int ndof_lcl=__atoms->c_dim*__atoms->natms_lcl;
    int n=ndof_lcl;
    type0* c=__atoms->c->begin();
    if(!__atoms->dof_c->is_empty())
    {
        bool* c_dof=__atoms->dof_c->begin();
        for(int i=0;i<n;i++)
            if(!c_dof[i] || c[i]<0.0) ndof_lcl--;
    }
    else
        for(int i=0;i<n;i++)
            if(c[i]<0.0) ndof_lcl--;
    int ndof;
    MPI_Allreduce(&ndof_lcl,&ndof,1,MPI_INT,MPI_SUM,__atoms->world);
    return ndof;
}
/*--------------------------------------------
 
 --------------------------------------------*/
void DAE::init_static()
{
    c_dim=atoms->c_dim;
    ncs=atoms->natms_lcl*c_dim;
    a_tol_sqrt_nc_dofs=a_tol*sqrt(static_cast<type0>(calc_ndofs(atoms)));
    ff->c_d->fill();
    dynamic=new DynamicDMD(atoms,ff,false,{},{},{});
    dynamic->init();
    c=atoms->c->begin();
    c_d=ff->c_d->begin();
    Memory::alloc(c_init,ncs);
    memcpy(c_init,c,ncs*sizeof(type0));
}
/*--------------------------------------------
 
 --------------------------------------------*/
void DAE::fin_static()
{
    Memory::dealloc(c_init);
    c=c_d=NULL;
    dynamic->fin();
    delete dynamic;
    dynamic=NULL;
}
/*------------------------------------------------------------------------------------------------------------------------------------
 
 ------------------------------------------------------------------------------------------------------------------------------------*/
PyObject* DAE::__new__(PyTypeObject* type,PyObject* args,PyObject* kwds)
{
    Object* __self=reinterpret_cast<Object*>(type->tp_alloc(type,0));
    PyObject* self=reinterpret_cast<PyObject*>(__self);
    return self;
}
/*--------------------------------------------
 
 --------------------------------------------*/
int DAE::__init__(PyObject* self,PyObject* args,PyObject* kwds)
{
    FuncAPI<> f("__init__");
    
    if(f(args,kwds)==-1) return -1;
    Object* __self=reinterpret_cast<Object*>(self);
    __self->dae=new DAE();
    __self->xprt=NULL;
    
    return 0;
}
/*--------------------------------------------
 
 --------------------------------------------*/
PyObject* DAE::__alloc__(PyTypeObject* type,Py_ssize_t)
{
    Object* __self=new Object;
    __self->ob_type=type;
    __self->ob_refcnt=1;
    __self->dae=NULL;
    __self->xprt=NULL;
    return reinterpret_cast<PyObject*>(__self);
}
/*--------------------------------------------
 
 --------------------------------------------*/
void DAE::__dealloc__(PyObject* self)
{
    Object* __self=reinterpret_cast<Object*>(self);
    delete __self->dae;
    __self->dae=NULL;
    if(__self->xprt) Py_DECREF(__self->xprt);
    __self->xprt=NULL;
    delete __self;
}
/*--------------------------------------------*/
PyTypeObject DAE::TypeObject={PyObject_HEAD_INIT(NULL)};
/*--------------------------------------------*/
int DAE::setup_tp()
{
    TypeObject.tp_name="mapp.dmd.dae";
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
    
    int ichk=PyType_Ready(&TypeObject);
    if(ichk<0) return ichk;
    Py_INCREF(&TypeObject);
    return ichk;
    
}
/*--------------------------------------------*/
PyGetSetDef DAE::getset[]={[0 ... 6]={NULL,NULL,NULL,NULL,NULL}};
/*--------------------------------------------*/
void DAE::setup_tp_getset()
{
    getset_a_tol(getset[0]);
    getset_max_nsteps(getset[1]);
    getset_min_dt(getset[2]);
    getset_max_dc_rel(getset[3]);
    getset_ntally(getset[4]);
    getset_export(getset[5]);
}
/*--------------------------------------------*/
PyMethodDef DAE::methods[]={[0 ... 0]={NULL}};
/*--------------------------------------------*/
void DAE::setup_tp_methods()
{
}
/*--------------------------------------------
 
 --------------------------------------------*/
void DAE::getset_a_tol(PyGetSetDef& getset)
{
    getset.name=(char*)"a_tol";
    getset.doc=(char*)R"---(
    (double) LTE tolerance
    
    Absolute error tolerence in local trucation error
    )---";
    getset.get=[](PyObject* self,void*)->PyObject*
    {
        return var<type0>::build(reinterpret_cast<Object*>(self)->dae->a_tol,NULL);
    };
    getset.set=[](PyObject* self,PyObject* op,void*)->int
    {
        VarAPI<type0> a_tol("a_tol");
        a_tol.logics[0]=VLogics("gt",0.0)*VLogics("lt",1.0);
        int ichk=a_tol.set(op);
        if(ichk==-1) return -1;
        reinterpret_cast<Object*>(self)->dae->a_tol=a_tol.val;
        return 0;
    };
}
/*--------------------------------------------
 
 --------------------------------------------*/
void DAE::getset_max_nsteps(PyGetSetDef& getset)
{
    getset.name=(char*)"max_nsteps";
    getset.doc=(char*)R"---(
    (int) maximum number of steps
    
    Maximum number of steps to achieve energy minimization
    )---";
    getset.get=[](PyObject* self,void*)->PyObject*
    {
        return var<int>::build(reinterpret_cast<Object*>(self)->dae->max_nsteps,NULL);
    };
    getset.set=[](PyObject* self,PyObject* op,void*)->int
    {
        VarAPI<int> max_nsteps("max_nsteps");
        max_nsteps.logics[0]=VLogics("ge",0);
        int ichk=max_nsteps.set(op);
        if(ichk==-1) return -1;
        reinterpret_cast<Object*>(self)->dae->max_nsteps=max_nsteps.val;
        return 0;
    };
}
/*--------------------------------------------
 
 --------------------------------------------*/
void DAE::getset_min_dt(PyGetSetDef& getset)
{
    getset.name=(char*)"min_dt";
    getset.doc=(char*)R"---(
    (double) minimum time step
    
    Minimum time step
    )---";
    getset.get=[](PyObject* self,void*)->PyObject*
    {
        return var<type0>::build(reinterpret_cast<Object*>(self)->dae->min_dt,NULL);
    };
    getset.set=[](PyObject* self,PyObject* op,void*)->int
    {
        VarAPI<type0> min_dt("min_dt");
        min_dt.logics[0]=VLogics("gt",0.0);
        int ichk=min_dt.set(op);
        if(ichk==-1) return -1;
        reinterpret_cast<Object*>(self)->dae->min_dt=min_dt.val;
        return 0;
    };
}
/*--------------------------------------------
 
 --------------------------------------------*/
void DAE::getset_max_dc_rel(PyGetSetDef& getset)
{
    getset.name=(char*)"max_dc_rel";
    getset.doc=(char*)R"---(
    (double) maximum relative c variation
    
    Maximum relative c variation
    )---";
    getset.get=[](PyObject* self,void*)->PyObject*
    {
        return var<type0>::build(reinterpret_cast<Object*>(self)->dae->max_dc_rel,NULL);
    };
    getset.set=[](PyObject* self,PyObject* op,void*)->int
    {
        VarAPI<type0> max_dc_rel("max_dc_rel");
        max_dc_rel.logics[0]=VLogics("gt",0.0);
        int ichk=max_dc_rel.set(op);
        if(ichk==-1) return -1;
        reinterpret_cast<Object*>(self)->dae->max_dc_rel=max_dc_rel.val;
        return 0;
    };
}
/*--------------------------------------------
 
 --------------------------------------------*/
void DAE::getset_ntally(PyGetSetDef& getset)
{
    getset.name=(char*)"ntally";
    getset.doc=(char*)R"---(
    (int) thermodynamic tallying period
    
    Number of steps to be taken from one thermodynamics output to the next.
    )---";
    getset.get=[](PyObject* self,void*)->PyObject*
    {
        return var<int>::build(reinterpret_cast<Object*>(self)->dae->ntally,NULL);
    };
    getset.set=[](PyObject* self,PyObject* op,void*)->int
    {
        VarAPI<int> ntally("ntally");
        ntally.logics[0]=VLogics("ge",0);
        int ichk=ntally.set(op);
        if(ichk==-1) return -1;
        reinterpret_cast<Object*>(self)->dae->ntally=ntally.val;
        return 0;
    };
}
/*--------------------------------------------
 
 --------------------------------------------*/
void DAE::getset_export(PyGetSetDef& getset)
{
    getset.name=(char*)"export";
    getset.doc=(char*)R"---(
    (mapp.dmd.export) export object
    
    Export object to record the snapshots of the system while minimizing
    )---";
    getset.get=[](PyObject* self,void*)->PyObject*
    {
        ExportDMD::Object* xprt=reinterpret_cast<Object*>(self)->xprt;
        if(!xprt) Py_RETURN_NONE;
        Py_INCREF(xprt);
        return reinterpret_cast<PyObject*>(xprt);
    };
    getset.set=[](PyObject* self,PyObject* op,void*)->int
    {
        VarAPI<OP<ExportDMD>> xprt("export");
        int ichk=xprt.set(op);
        if(ichk==-1) return -1;
        if(reinterpret_cast<Object*>(self)->xprt) Py_DECREF(reinterpret_cast<Object*>(self)->xprt);
        Py_INCREF(xprt.val.ob);
        reinterpret_cast<Object*>(self)->xprt=reinterpret_cast<ExportDMD::Object*>(xprt.val.ob);
        return 0;
    };
}
/*--------------------------------------------
 
 --------------------------------------------*/
void DAE::ml_run(PyMethodDef& tp_methods)
{
    tp_methods.ml_flags=METH_VARARGS | METH_KEYWORDS;
    tp_methods.ml_name="run";    
    tp_methods.ml_meth=(PyCFunction)(PyCFunctionWithKeywords)
    [](PyObject* self,PyObject* args,PyObject* kwds)->PyObject*
    {
        Object* __self=reinterpret_cast<Object*>(self);
        FuncAPI<OP<AtomsDMD>,type0> f("run",{"atoms","t"});
        f.logics<1>()[0]=VLogics("ge",0.0);
        if(f(args,kwds)) return NULL;
        
        AtomsDMD* __atoms=reinterpret_cast<AtomsDMD::Object*>(f.val<0>().ob)->atoms;
        ForceFieldDMD* __ff=reinterpret_cast<AtomsDMD::Object*>(f.val<0>().ob)->ff;
        ExportDMD* __xprt=__self->xprt==NULL ? NULL:__self->xprt->xprt;
        try
        {
            __self->dae->pre_run_chk(__atoms,__ff);
        }
        catch(std::string& err_msg)
        {
            PyErr_SetString(PyExc_TypeError,err_msg.c_str());
            return NULL;
        }
        
        __self->dae->atoms=__atoms;
        __self->dae->ff=__ff;
        __self->dae->xprt=__xprt;
        
        try
        {
            __self->dae->init_static();
        }
        catch(std::string& err_msg)
        {
            __self->dae->xprt=NULL;
            __self->dae->ff=NULL;
            __self->dae->atoms=NULL;
            PyErr_SetString(PyExc_TypeError,err_msg.c_str());
            return NULL;
        }
        
        __self->dae->run_static(f.val<1>());
        
        __self->dae->fin_static();
        
        __self->dae->xprt=NULL;
        __self->dae->ff=NULL;
        __self->dae->atoms=NULL;
        
        Py_RETURN_NONE;
    };
    
    tp_methods.ml_doc=(char*)R"---(
    run(atoms,t)
   
    Execute DEA
    
    This method starts differential-algebraic equations for a given atoms object and number of time.
    
    Parameters
    ----------
    atoms : mapp.dmd.atoms
        System of interest
    t : double
        Desired time
        
    Returns
    -------
    None

    )---";
}


