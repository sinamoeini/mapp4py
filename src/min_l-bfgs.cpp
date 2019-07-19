/*--------------------------------------------
 Created by Sina on 07/02/14.
 Copyright (c) 2013 MIT. All rights reserved.
 
 L-BFGS minimization is written based on
 Numerical Optimization written by Nocedal & 
 Wright, second edition, pages 177-179, 
 Algorithm 7.4 & 7.5 Equation (7.20)
 
 with respect to notations in Nocedal:
 new_y_i=y_i
 new_rho_i=rho_i
 new_alpha_i=alpha_i
 new_beta=beta
 --------------------------------------------*/
#include "min_l-bfgs.h"
#include <stdlib.h>
#include "ff.h"
#include "thermo_dynamics.h"
#include "ls.h"
#include "dynamic_md.h"
#include "memory.h"
#include "atoms_md.h"
#include "ff_styles.h"
#include "MAPP.h"
using namespace MAPP_NS;
/*--------------------------------------------
 constructor
 --------------------------------------------*/
MinLBFGS::MinLBFGS(int __m,type0 __e_tol,
bool(&__H_dof)[__dim__][__dim__],bool __affine,type0 __max_dx,LineSearch* __ls):
MinCG(__e_tol,__H_dof,__affine,__max_dx,__ls),
m(__m),
rho(NULL),
alpha(NULL),
s(NULL),
y(NULL)
{
}
/*--------------------------------------------
 destructor
 --------------------------------------------*/
MinLBFGS::~MinLBFGS()
{
}
/*--------------------------------------------
 init before a run
 --------------------------------------------*/
void MinLBFGS::init()
{
    x.~VecTens();
    new (&x) VecTens<type0,1>(atoms,chng_box,atoms->H,atoms->x);
    f.~VecTens();
    new (&f) VecTens<type0,1>(atoms,chng_box,ff->F_H,ff->f);
    h.~VecTens();
    new (&h) VecTens<type0,1>(atoms,chng_box,__dim__);
    x0.~VecTens();
    new (&x0) VecTens<type0,1>(atoms,chng_box,__dim__);
    x_d.~VecTens();
    new (&x_d) VecTens<type0,1>(atoms,chng_box,__dim__);
    f0.~VecTens();
    new (&f0) VecTens<type0,1>(atoms,chng_box,__dim__);
    
    if(m)
    {
        s=new VecTens<type0,1>[m];
        y=new VecTens<type0,1>[m];
        for(int i=0;i<m;i++)
        {
            s[i].~VecTens();
            new (s+i) VecTens<type0,1>(atoms,chng_box,__dim__);
            y[i].~VecTens();
            new (y+i) VecTens<type0,1>(atoms,chng_box,__dim__);
        }
    }
    
    dynamic=new DynamicMD(atoms,ff,chng_box,{},{atoms->x_dof,h.vecs[0],x0.vecs[0],f0.vecs[0]},{atoms->x_d});
    
    for(int i=0;i<m;i++)
    {
        dynamic->add_xchng(s[i].vecs[0]);
        dynamic->add_xchng(y[i].vecs[0]);
    }
    
    dynamic->init();
    
    Memory::alloc(rho,m);
    Memory::alloc(alpha,m);
    
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
 fin after a run
 --------------------------------------------*/
void MinLBFGS::fin()
{
    if(xprt)
    {
        xprt->fin();
        xprt->atoms=NULL;
    }
    
    Memory::dealloc(alpha);
    Memory::dealloc(rho);
    rho=alpha=NULL;
    
    dynamic->fin();
    delete dynamic;
    dynamic=NULL;
    
    delete [] s;
    delete [] y;
    s=y=NULL;
    
    f0.~VecTens();
    x_d.~VecTens();
    x0.~VecTens();
    h.~VecTens();
    f.~VecTens();
    x.~VecTens();
}
/*--------------------------------------------
 run
 --------------------------------------------*/
void MinLBFGS::run(int nsteps)
{
    if(dynamic_cast<LineSearchGoldenSection*>(ls))
        return run(dynamic_cast<LineSearchGoldenSection*>(ls),nsteps);
    
    if(dynamic_cast<LineSearchBrent*>(ls))
        return run(dynamic_cast<LineSearchBrent*>(ls),nsteps);
    
    if(dynamic_cast<LineSearchBackTrack*>(ls))
        return run(dynamic_cast<LineSearchBackTrack*>(ls),nsteps);
}
/*------------------------------------------------------------------------------------------------------------------------------------
 
 ------------------------------------------------------------------------------------------------------------------------------------*/
PyObject* MinLBFGS::__new__(PyTypeObject* type,PyObject* args,PyObject* kwds)
{
    Object* __self=reinterpret_cast<Object*>(type->tp_alloc(type,0));
    PyObject* self=reinterpret_cast<PyObject*>(__self);
    return self;
}
/*--------------------------------------------
 
 --------------------------------------------*/
int MinLBFGS::__init__(PyObject* self,PyObject* args,PyObject* kwds)
{
    
    FuncAPI<int,type0,symm<bool[__dim__][__dim__]>,bool,type0,OP<LineSearch>> f("__init__",{"m","e_tol","H_dof","affine","max_dx","ls"});
    f.noptionals=6;
    f.logics<0>()[0]=VLogics("ge",0);
    f.logics<1>()[0]=VLogics("ge",0.0);
    f.logics<4>()[0]=VLogics("gt",0.0);
    
    //set the defualts
    f.val<0>()=2;
    f.val<1>()=sqrt(std::numeric_limits<type0>::epsilon());
    for(int i=0;i<__dim__;i++) for(int j=0;j<__dim__;j++)f.val<2>()[i][j]=false;
    f.val<3>()=false;
    f.val<4>()=1.0;
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
    __self->min=new MinLBFGS(f.val<0>(),f.val<1>(),f.val<2>(),f.val<3>(),f.val<4>(),&(__self->ls->ls));
    __self->xprt=NULL;
    
    return 0;
}
/*--------------------------------------------
 
 --------------------------------------------*/
PyObject* MinLBFGS::__alloc__(PyTypeObject* type,Py_ssize_t)
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
void MinLBFGS::__dealloc__(PyObject* self)
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
PyTypeObject MinLBFGS::TypeObject ={PyObject_HEAD_INIT(NULL)};
/*--------------------------------------------*/
int MinLBFGS::setup_tp()
{
    TypeObject.tp_name="mapp.md.min_lbfgs";
    TypeObject.tp_doc=R"---(
    __init__(m=2,e_tol=1.0e-8,H_dof=[[False],[False,False],[False,False,False]],affine=False,max_dx=1.0,ls=mapp.md.ls_bt())
    
    L-BFGS minimization algorithm
    
    Parameters
    ----------
    m : int
       Maximum number of vectors to be stored in memory
    e_tol : double
       Energy tolerance criterion for stopping minimization
    H_dof : symm<bool[dim][dim]>
       Unitcell degrees of freedom during minimization, here dim is the dimension of simulation
    affine : bool
       If set to True atomic displacements would be affine
    max_dx : double
       Maximum displacement of any atom in one step of minimization
    ls : mapp.ls
       Line search method
    
    Notes
    -----
    Limited memory Broyden–Fletcher–Goldfarb–Shanno (L-BFGS) algorithm for minimization, see :cite:`nocedal_numerical_2006`.
    
    References
    ----------
    .. bibliography:: ../refs.bib
       :filter: docname in docnames
       :style: unsrt
    
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
    
    TypeObject.tp_base=&MinCG::TypeObject;
    
    int ichk=PyType_Ready(&TypeObject);
    if(ichk<0) return ichk;
    Py_INCREF(&TypeObject);
    GET_WRAPPER_DOC(TypeObject,__init__)=NULL;
    return ichk;
}
/*--------------------------------------------*/
PyMethodDef MinLBFGS::methods[]=EmptyPyMethodDef(1);
/*--------------------------------------------*/
void MinLBFGS::setup_tp_methods()
{
}
/*--------------------------------------------*/
PyGetSetDef MinLBFGS::getset[]=EmptyPyGetSetDef(2);
/*--------------------------------------------*/
void MinLBFGS::setup_tp_getset()
{
    getset_m(getset[0]);
}
/*--------------------------------------------
 
 --------------------------------------------*/
void MinLBFGS::getset_m(PyGetSetDef& getset)
{
    getset.name=(char*)"m";
    getset.doc=(char*)R"---(
    (int) Maximum No. of vectors
    
    Maximum number of vectors to be stored in memory
    )---";
    getset.get=[](PyObject* self,void*)->PyObject*
    {
        return var<int>::build(reinterpret_cast<Object*>(self)->min->m);
    };
    getset.set=[](PyObject* self,PyObject* op,void*)->int
    {
        VarAPI<int> m("m");
        m.logics[0]=VLogics("ge",0);
        int ichk=m.set(op);
        if(ichk==-1) return -1;
        reinterpret_cast<Object*>(self)->min->m=m.val;
        return 0;
    };
}
