#include "min_cg.h"
using namespace MAPP_NS;


/*--------------------------------------------
 constructor
 --------------------------------------------*/
MinCG::MinCG(type0 __e_tol,
bool(&__H_dof)[__dim__][__dim__],bool __affine,type0 __max_dx,LineSearch* __ls):
Min(__e_tol,__H_dof,__affine,__max_dx,__ls),
X_DOF(true),
atoms(NULL),
ff(NULL),
xprt(NULL)
{
}
/*--------------------------------------------
 destructor
 --------------------------------------------*/
MinCG::~MinCG()
{
    atoms=NULL;
    ff=NULL;
}
/*--------------------------------------------
 pre run check it throw excepctions
 --------------------------------------------*/
void MinCG::pre_run_chk(AtomsMD* __atoms,ForceFieldMD* __ff)
{
    //check if configuration is loaded
    if(!__atoms)
        throw std::string("cannot start minimization without initial conditions");
    
    //check if force field is loaded
    if(!__ff)
        throw std::string("cannot start minimization without governing equations (force field)");
}
/*--------------------------------------------
lss=['LineSearchBrent','LineSearchGoldenSection','LineSearchBackTrack']
cases=['true','false']
tab='    '
for i in cases:
    for j in cases:
        for ls in lss:
            if i=='true' or j=='true':
                print(tab+'if(dynamic_cast<%s*>(ls)!=NULL && chng_box==%s &&  X_DOF==%s)' %(ls,i,j))
                print(tab+tab+'return run<%s,%s>(dynamic_cast<%s*>(ls),nsteps);' %(i,j,ls))
 --------------------------------------------*/
/*--------------------------------------------
 
 --------------------------------------------*/
void MinCG::init()
{
    chng_box=false;
    Algebra::DoLT<__dim__>::func([this](int i,int j)
    {
        if(H_dof[i][j]) chng_box=true;
    });
    
    try
    {
        MinHelper::CondB<>::init(*this,chng_box,X_DOF);
    }
    catch(std::string& err_msg)
    {
        throw err_msg;
    }
}
/*--------------------------------------------
 
 --------------------------------------------*/
void MinCG::run(int nsteps)
{
    MinHelper::CondLS<LineSearchBrent,LineSearchGoldenSection,LineSearchBackTrack>::run(*this,nsteps,ls,chng_box,X_DOF);
}
/*--------------------------------------------
 
 --------------------------------------------*/
void MinCG::fin()
{
    MinHelper::CondB<>::fin(*this,chng_box,X_DOF);
}
/*------------------------------------------------------------------------------------------------------------------------------------
 
 ------------------------------------------------------------------------------------------------------------------------------------*/
PyObject* MinCG::__new__(PyTypeObject* type,PyObject* args,PyObject* kwds)
{
    Object* __self=reinterpret_cast<Object*>(type->tp_alloc(type,0));
    PyObject* self=reinterpret_cast<PyObject*>(__self);
    return self;
}
/*--------------------------------------------
 
 --------------------------------------------*/
int MinCG::__init__(PyObject* self,PyObject* args,PyObject* kwds)
{
    FuncAPI<type0,symm<bool[__dim__][__dim__]>,bool,type0,OP<LineSearch>> f("__init__",{"e_tol","H_dof","affine","max_dx","ls"});
    f.noptionals=5;
    f.logics<0>()[0]=VLogics("ge",0.0);
    f.logics<3>()[0]=VLogics("gt",0.0);
    
    //set the defualts
    f.val<0>()=sqrt(std::numeric_limits<type0>::epsilon());
    for(int i=0;i<__dim__;i++) for(int j=0;j<__dim__;j++)f.val<1>()[i][j]=false;
    f.val<2>()=false;
    f.val<3>()=1.0;
    PyObject* empty_tuple=PyTuple_New(0);
    PyObject* empty_dict=PyDict_New();
    PyObject* __ls=LineSearchBackTrack::__new__(&LineSearchBackTrack::TypeObject,empty_tuple,empty_dict);
    LineSearchBackTrack::__init__(__ls,empty_tuple,empty_dict);
    Py_DECREF(empty_dict);
    Py_DECREF(empty_tuple);
    f.val<4>().ob=__ls;
    
    
    if(f(args,kwds)==-1) return -1;
    
    
    
    Object* __self=reinterpret_cast<Object*>(self);
    Py_INCREF(f.val<4>().ob);
    __self->ls=reinterpret_cast<LineSearch::Object*>(f.val<4>().ob);
    __self->min=new MinCG(f.val<0>(),f.val<1>(),f.val<2>(),f.val<3>(),&(__self->ls->ls));
    __self->xprt=NULL;
    
    return 0;
}
/*--------------------------------------------
 
 --------------------------------------------*/
PyObject* MinCG::__alloc__(PyTypeObject* type,Py_ssize_t)
{
    Object* __self=new Object;
    Py_TYPE(__self)=type;
    Py_REFCNT(__self)=1;
    __self->min=NULL;
    __self->ls=NULL;
    __self->xprt=NULL;
    return reinterpret_cast<PyObject*>(__self);
}
/*--------------------------------------------
 
 --------------------------------------------*/
void MinCG::__dealloc__(PyObject* self)
{
    Object* __self=reinterpret_cast<Object*>(self);
    delete __self->min;
    __self->min=NULL;
    if(__self->ls) Py_DECREF(__self->ls);
    __self->ls=NULL;
    if(__self->xprt) Py_DECREF(__self->xprt);
    __self->xprt=NULL;
    delete __self;
}
/*--------------------------------------------*/
PyTypeObject MinCG::TypeObject={PyObject_HEAD_INIT(NULL)};
/*--------------------------------------------*/
int MinCG::setup_tp()
{
    TypeObject.tp_name="mapp4py.md.min_cg_new";
    TypeObject.tp_doc=R"---(
    __init__(e_tol=1.0e-8,H_dof=[[False],[False,False],[False,False,False]],affine=False,max_dx=1.0,ls=mapp4py.ls_bt())
    
    CG minimization algorithm
    
    Parameters
    ----------
    e_tol : double
       Energy tolerance criterion for stopping minimization
    H_dof : symm<bool[dim][dim]>
       Unitcell degrees of freedom during minimization, here dim is the dimension of simulation
    affine : bool
       If set to True atomic displacements would be affine
    max_dx : double
       Maximum displacement of any atom in one step of minimization
    ls : mapp4py.ls
       Line search method

    Notes
    -----
    Cojugate Gradient (CG) algorithm for minimization.

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
PyGetSetDef MinCG::getset[]=EmptyPyGetSetDef(8);
/*--------------------------------------------*/
void MinCG::setup_tp_getset()
{
    getset_e_tol(getset[0]);
    getset_H_dof(getset[1]);
    getset_max_dx(getset[2]);
    getset_ls(getset[3]);
    getset_ntally(getset[4]);
    getset_export(getset[5]);
    getset_x_dof(getset[6]);
    
}
/*--------------------------------------------*/
PyMethodDef MinCG::methods[]=EmptyPyMethodDef(2);
/*--------------------------------------------*/
void MinCG::setup_tp_methods()
{
    ml_run(methods[0]);
}
/*--------------------------------------------
 
 --------------------------------------------*/
void MinCG::getset_export(PyGetSetDef& getset)
{
    getset.name=(char*)"export";
    getset.doc=(char*)R"---(
    (mapp4py.md.export) export object
    
    Export object to record the snapshots of the system while minimizing
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
void MinCG::getset_x_dof(PyGetSetDef& getset)
{
    getset.name=(char*)"x_dof";
    getset.doc=(char*)R"---(
    (bool) if set true position of atoms will considered as degrees of freedom
    
    If set true position of atoms will considered as degrees of freedom
    )---";
    getset.get=[](PyObject* self,void*)->PyObject*
    {
        return var<bool>::build(reinterpret_cast<Object*>(self)->min->X_DOF);
    };
    getset.set=[](PyObject* self,PyObject* op,void*)->int
    {
        VarAPI<bool> x_dof("x_dof");
        if(x_dof.set(op)==-1) return -1;
        reinterpret_cast<Object*>(self)->min->X_DOF=x_dof.val;
        return 0;
    };
}
/*--------------------------------------------
 
 --------------------------------------------*/
void MinCG::ml_run(PyMethodDef& tp_methods)
{
    tp_methods.ml_flags=METH_VARARGS | METH_KEYWORDS;
    tp_methods.ml_name="run";
    tp_methods.ml_meth=(PyCFunction)(PyCFunctionWithKeywords)(
    [](PyObject* self,PyObject* args,PyObject* kwds)->PyObject*
    {
        Object* __self=reinterpret_cast<Object*>(self);
        FuncAPI<OP<AtomsMD>,int> f("run",{"atoms","max_nsteps"});
        f.logics<1>()[0]=VLogics("ge",0);
        
        
        
        if(f(args,kwds)) return NULL;
        
        
        AtomsMD* __atoms=reinterpret_cast<AtomsMD::Object*>(f.val<0>().ob)->atoms;
        ForceFieldMD* __ff=reinterpret_cast<AtomsMD::Object*>(f.val<0>().ob)->ff;
        ExportMD* __xprt=__self->xprt==NULL ? NULL:__self->xprt->xprt;
        

        
        try
        {
            __self->min->pre_run_chk(__atoms,__ff);
        }
        catch(std::string& err_msg)
        {
            PyErr_SetString(PyExc_TypeError,err_msg.c_str());
            return NULL;
        }
        
        __self->min->atoms=__atoms;
        __self->min->ff=__ff;
        __self->min->xprt=__xprt;
        try
        {
            __self->min->init();
        }
        catch(std::string& err_msg)
        {
            __self->min->xprt=NULL;
            __self->min->ff=NULL;
            __self->min->atoms=NULL;
            PyErr_SetString(PyExc_TypeError,err_msg.c_str());
            return NULL;
        }
        
        __self->min->run(f.val<1>());
        
        __self->min->fin();
        
        __self->min->xprt=NULL;
        __self->min->ff=NULL;
        __self->min->atoms=NULL;
        
        Py_RETURN_NONE;
    });
    
    tp_methods.ml_doc=(char*)R"---(
    run(atoms,max_nsteps)
   
    Execute minimization
    
    This method starts the energy minimization for a given atoms object and maximum number of steps.
    
    Parameters
    ----------
    atoms : mapp4py.md.atoms
        System of interest
    max_nsteps : int
        Maximum number of steps to achieve energy minimization
        
    Returns
    -------
    None

    )---";
}
