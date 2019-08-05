#include "min_cg_dmd.h"
#include <stdlib.h>
using namespace MAPP_NS;
/*--------------------------------------------
 constructor
 --------------------------------------------*/
MinCGDMD::MinCGDMD(type0 __e_tol,
bool(&__H_dof)[__dim__][__dim__],bool __affine,type0 __max_dx,type0 __max_dalpha,LineSearch* __ls):
Min(__e_tol,__H_dof,__affine,__max_dx,__ls),
ALPHA_DOF(true),
C_DOF(false),
max_dalpha(__max_dalpha),
atoms(NULL),
ff(NULL),
xprt(NULL)
{
}
/*--------------------------------------------
 destructor
 --------------------------------------------*/
MinCGDMD::~MinCGDMD()
{
    atoms=NULL;
    ff=NULL;
}
/*--------------------------------------------
 pre run check it throw excepctions
 --------------------------------------------*/
void MinCGDMD::pre_run_chk(AtomsDMD* __atoms,ForceFieldDMD* __ff)
{
    //check if configuration is loaded
    if(!__atoms)
        throw std::string("cannot start minimization without initial conditions");
    
    //check if force field is loaded
    if(!__ff)
        throw std::string("cannot start minimization without governing equations (force field)");
    
    if(std::isnan(__atoms->kB))
        throw std::string("boltzmann constant should be set prior to minimizatiom");
    
    if(std::isnan(__atoms->hP))
        throw std::string("planck constant should be set prior to minimizatiom");

    
    if(std::isnan(__atoms->temp))
        throw std::string("temperature should be set prior to minimizatiom");
    
    if(C_DOF)
    {
        int err;
        const int c_dim=__atoms->c_dim;
        type0* c=__atoms->c->begin();
        int natms_lcl=__atoms->natms_lcl;
        bool chk=true;
        type0 cv;
        if(__atoms->c_dof->is_empty())
        {
            for(int i=0;i<natms_lcl && chk;i++,c+=c_dim)
            {
                cv=1.0;
                for(int j=0;j<c_dim;j++)
                {
                    if(c[j]<0.0) continue;
                    cv-=c[j];
                    if(c[j]==0.0) chk=false;
                }
                if(cv<=0.0) chk=false;
            }
        }
        else
        {
            bool* c_dof=__atoms->c_dof->begin();
            bool cv_dof;
            for(int i=0;i<natms_lcl && chk;i++,c+=c_dim,c_dof+=c_dim)
            {
                cv=1.0;
                cv_dof=false;
                for(int j=0;j<c_dim;j++)
                {
                    if(c[j]<0.0) continue;
                    cv-=c[j];
                    if(c_dof[j]==false) continue;
                    cv_dof=true;
                    if(c[j]==0.0) chk=false;

                }
                if(cv_dof && cv<=0.0) chk=false;
            }
            
        }
        
        
        
        int err_lcl=chk ? 0:1;
        
        MPI_Allreduce(&err_lcl,&err,1,Vec<int>::MPI_T,MPI_MAX,__atoms->world);
        if(err)
            throw std::string("for energy minimization no c can be either 0.0 or 1.0");
    }
}
/*--------------------------------------------
 init before a run
 --------------------------------------------*/
void MinCGDMD::init()
{
    chng_box=false;
    Algebra::DoLT<__dim__>::func([this](int i,int j)
    {
        if(H_dof[i][j]) chng_box=true;
    });
    
    try
    {
        MinHelper::CondB<>::init(*this,chng_box,!affine,ALPHA_DOF,C_DOF);
    }
    catch(std::string& err_msg)
    {
        throw err_msg;
    }
}
/*--------------------------------------------
 min
 --------------------------------------------*/
void MinCGDMD::run(int nsteps)
{
    MinHelper::CondLS<LineSearchBrent,LineSearchGoldenSection,LineSearchBackTrack>::run(*this,nsteps,ls,chng_box,!affine,ALPHA_DOF,C_DOF
    ,ntally!=0,xprt!=NULL);
}
/*--------------------------------------------
 finishing minimization
 --------------------------------------------*/
void MinCGDMD::fin()
{
    MinHelper::CondB<>::fin(*this,chng_box,!affine,ALPHA_DOF,C_DOF);
}
/*--------------------------------------------
 
 --------------------------------------------*/
void MinCGDMD::create_thermo()
{
    thermo=new(&thermo_mem) ThermoDynamics(6,
    "FE",atoms->fe,
    "S[0][0]",atoms->S_fe[0][0],
    "S[1][1]",atoms->S_fe[1][1],
    "S[2][2]",atoms->S_fe[2][2],
    "S[1][2]",atoms->S_fe[2][1],
    "S[2][0]",atoms->S_fe[2][0],
    "S[0][1]",atoms->S_fe[1][0]);
}
/*--------------------------------------------
 
 --------------------------------------------*/
void MinCGDMD::destroy_thermo()
{
    thermo->~ThermoDynamics();
}
/*------------------------------------------------------------------------------------------------------------------------------------
 
 ------------------------------------------------------------------------------------------------------------------------------------*/
PyObject* MinCGDMD::__new__(PyTypeObject* type,PyObject* args,PyObject* kwds)
{
    Object* __self=reinterpret_cast<Object*>(type->tp_alloc(type,0));
    PyObject* self=reinterpret_cast<PyObject*>(__self);
    return self;
}
/*--------------------------------------------
 
 --------------------------------------------*/
int MinCGDMD::__init__(PyObject* self,PyObject* args,PyObject* kwds)
{
    FuncAPI<type0,symm<bool[__dim__][__dim__]>,bool,type0,type0,OP<LineSearch>> f("__init__",{"e_tol","H_dof","affine","max_dx","max_dalpha","ls"});
    f.noptionals=6;
    f.logics<0>()[0]=VLogics("ge",0.0);
    f.logics<3>()[0]=VLogics("gt",0.0);
    f.logics<4>()[0]=VLogics("gt",0.0);
    
    //set the defualts
    f.val<0>()=sqrt(std::numeric_limits<type0>::epsilon());
    for(int i=0;i<__dim__;i++) for(int j=0;j<__dim__;j++)f.val<1>()[i][j]=false;
    f.val<2>()=false;
    f.val<3>()=1.0;
    f.val<4>()=0.1;
    PyObject* empty_tuple=PyTuple_New(0);
    PyObject* empty_dict=PyDict_New();
    PyObject* __ls=LineSearchBackTrack::__new__(&LineSearchBackTrack::TypeObject,empty_tuple,empty_dict);
    LineSearchBackTrack::__init__(__ls,empty_tuple,empty_dict);
    Py_DECREF(empty_dict);
    Py_DECREF(empty_tuple);
    f.val<5>().ob=__ls;
    
    
    if(f(args,kwds)==-1) return -1;
    
    
    
    Object* __self=reinterpret_cast<Object*>(self);
    Py_INCREF(f.val<5>().ob);
    __self->ls=reinterpret_cast<LineSearch::Object*>(f.val<5>().ob);
    __self->min=new MinCGDMD(f.val<0>(),f.val<1>(),f.val<2>(),f.val<3>(),f.val<4>(),&(__self->ls->ls));
    __self->xprt=NULL;
    
    return 0;
}
/*--------------------------------------------
 
 --------------------------------------------*/
PyObject* MinCGDMD::__alloc__(PyTypeObject* type,Py_ssize_t)
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
void MinCGDMD::__dealloc__(PyObject* self)
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
PyTypeObject MinCGDMD::TypeObject={PyObject_HEAD_INIT(NULL)};
/*--------------------------------------------*/
int MinCGDMD::setup_tp()
{
    TypeObject.tp_name="mapp4py.dmd.min_cg_new";
    TypeObject.tp_doc=R"---(
    __init__(e_tol=1.0e-8,H_dof=[[False],[False,False],[False,False,False]],affine=False,max_dx=1.0,max_dalpha=0.1,ls=mapp4py.dmd.ls_bt())
    
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
    max_dalpha : double
       Maximum change in alpha component of any atom in one step of minimization
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
PyGetSetDef MinCGDMD::getset[]=EmptyPyGetSetDef(11);
/*--------------------------------------------*/
void MinCGDMD::setup_tp_getset()
{
    getset_e_tol(getset[0]);
    getset_H_dof(getset[1]);
    getset_max_dx(getset[2]);
    getset_max_dalpha(getset[3]);
    getset_ls(getset[4]);
    getset_ntally(getset[5]);
    getset_export(getset[6]);
    getset_affine(getset[7]);
    getset_alpha_dof(getset[8]);
    getset_c_dof(getset[9]);
}
/*--------------------------------------------*/
PyMethodDef MinCGDMD::methods[]=EmptyPyMethodDef(2);
/*--------------------------------------------*/
void MinCGDMD::setup_tp_methods()
{
    ml_run(methods[0]);
}
/*--------------------------------------------
 
 --------------------------------------------*/
void MinCGDMD::getset_max_dalpha(PyGetSetDef& getset)
{
    getset.name=(char*)"max_dalpha";
    getset.doc=(char*)R"---(
    (double) mximum alpha change
    
    Maximum change in alpha component of any atom in one step of minimization
    )---";
    getset.get=[](PyObject* self,void*)->PyObject*
    {
        return var<type0>::build(reinterpret_cast<Object*>(self)->min->max_dx);
    };
    getset.set=[](PyObject* self,PyObject* op,void*)->int
    {
        VarAPI<type0> max_dalpha("max_dalpha");
        max_dalpha.logics[0]=VLogics("gt",0.0);
        int ichk=max_dalpha.set(op);
        if(ichk==-1) return -1;
        reinterpret_cast<Object*>(self)->min->max_dalpha=max_dalpha.val;
        return 0;
    };
}
/*--------------------------------------------
 
 --------------------------------------------*/
void MinCGDMD::getset_export(PyGetSetDef& getset)
{
    getset.name=(char*)"export";
    getset.doc=(char*)R"---(
    (mapp4py.dmd.export) export object
    
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
void MinCGDMD::getset_alpha_dof(PyGetSetDef& getset)
{
    getset.name=(char*)"alpha_dof";
    getset.doc=(char*)R"---(
    (bool) if set true alpha of atoms will considered as degrees of freedom
    
    If set true alpha of atoms will considered as degrees of freedom
    )---";
    getset.get=[](PyObject* self,void*)->PyObject*
    {
        return var<bool>::build(reinterpret_cast<Object*>(self)->min->ALPHA_DOF);
    };
    getset.set=[](PyObject* self,PyObject* op,void*)->int
    {
        VarAPI<bool> alpha_dof("alpha_dof");
        if(alpha_dof.set(op)==-1) return -1;
        reinterpret_cast<Object*>(self)->min->ALPHA_DOF=alpha_dof.val;
        return 0;
    };
}
/*--------------------------------------------
 
 --------------------------------------------*/
void MinCGDMD::getset_c_dof(PyGetSetDef& getset)
{
    getset.name=(char*)"c_dof";
    getset.doc=(char*)R"---(
    (bool) if set true c of atoms will considered as degrees of freedom
        
        If set true c of atoms will considered as degrees of freedom
        )---";
        getset.get=[](PyObject* self,void*)->PyObject*
    {
        return var<bool>::build(reinterpret_cast<Object*>(self)->min->C_DOF);
    };
    getset.set=[](PyObject* self,PyObject* op,void*)->int
    {
        VarAPI<bool> c_dof("c_dof");
        if(c_dof.set(op)==-1) return -1;
        reinterpret_cast<Object*>(self)->min->C_DOF=c_dof.val;
        return 0;
    };
}
/*--------------------------------------------
 
 --------------------------------------------*/
void MinCGDMD::ml_run(PyMethodDef& tp_methods)
{
    tp_methods.ml_flags=METH_VARARGS | METH_KEYWORDS;
    tp_methods.ml_name="run";
    tp_methods.ml_meth=(PyCFunction)(PyCFunctionWithKeywords)(
    [](PyObject* self,PyObject* args,PyObject* kwds)->PyObject*
    {
        Object* __self=reinterpret_cast<Object*>(self);
        FuncAPI<OP<AtomsDMD>,int> f("run",{"atoms","max_nsteps"});
        f.logics<1>()[0]=VLogics("ge",0);
        
        
        
        if(f(args,kwds)) return NULL;
        
        
        AtomsDMD* __atoms=reinterpret_cast<AtomsDMD::Object*>(f.val<0>().ob)->atoms;
        ForceFieldDMD* __ff=reinterpret_cast<AtomsDMD::Object*>(f.val<0>().ob)->ff;
        ExportDMD* __xprt=__self->xprt==NULL ? NULL:__self->xprt->xprt;
        

        
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



