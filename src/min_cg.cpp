#include "min_cg.h"
#include <stdlib.h>
using namespace MAPP_NS;
/*--------------------------------------------
 constructor
 --------------------------------------------*/
MinCG::MinCG(type0 __e_tol,
bool(&__H_dof)[__dim__][__dim__],bool __affine,type0 __max_dx,LineSearch* __ls):
Min(__e_tol,__H_dof,__affine,__max_dx,__ls),
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
 
 --------------------------------------------*/
type0 MinCG::calc_ndofs()
{
    return 0.0;
}
/*--------------------------------------------
 
 --------------------------------------------*/
void MinCG::force_calc()
{
    ff->derivative_timer();
    if(chng_box)
        Algebra::DoLT<__dim__>::func([this](int i,int j){f.A[i][j]=H_dof[i][j] ? f.A[i][j]:0.0;});
}
/*--------------------------------------------
 
 --------------------------------------------*/
void MinCG::prep()
{
    x_d=h;
    if(!chng_box) return;
    
    const int natms_lcl=atoms->natms_lcl;
    type0* xvec=x0.vecs[0]->begin();
    type0* hvec=h.vecs[0]->begin();
    type0* x_dvec=x_d.vecs[0]->begin();
    Algebra::MLT_mul_MLT(atoms->H,h.A,x_d.A);
    
    
    if(affine)
    {
        for(int iatm=0;iatm<natms_lcl;iatm++,xvec+=__dim__,x_dvec+=__dim__)
            Algebra::V_mul_MLT(xvec,h.A,x_dvec);
    }
    else
    {
        for(int iatm=0;iatm<natms_lcl;iatm++,xvec+=__dim__,x_dvec+=__dim__,hvec+=__dim__)
            Algebra::V_mul_MLT_add_in(xvec,h.A,x_dvec);
    }
    
        
    
}
/*--------------------------------------------
 init before a run
 --------------------------------------------*/
void MinCG::init()
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
    
    dynamic=new DynamicMD(atoms,ff,chng_box,{},{atoms->x_dof,h.vecs[0],x0.vecs[0],x_d.vecs[0],f0.vecs[0]},{atoms->x_d});
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
 finishing minimization
 --------------------------------------------*/
void MinCG::fin()
{
    if(xprt)
    {
        xprt->fin();
        xprt->atoms=NULL;
    }
    
    dynamic->fin();
    delete dynamic;
    dynamic=NULL;
    
    f0.~VecTens();
    x_d.~VecTens();
    x0.~VecTens();
    h.~VecTens();
    f.~VecTens();
    x.~VecTens();
}
/*--------------------------------------------
 ff_test
 --------------------------------------------*/
#include "random.h"
void MinCG::ff_test(int seed,type0 __max_dx,type0 __max_st,int __n_desc)
{
    
    bool __chng_box=chng_box;
    chng_box=true;
    bool __H_dof[__dim__][__dim__];
    Algebra::V_eq<__dim__*__dim__>(&H_dof[0][0],&__H_dof[0][0]);
    Algebra::DoLT<__dim__>::func([this](int i,int j)
    {H_dof[i][j]=H_dof[j][i]=true;});
    
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
    
    dynamic=new DynamicMD(atoms,ff,chng_box,{},{atoms->x_dof,h.vecs[0],x0.vecs[0],x_d.vecs[0]},{atoms->x_d});
    dynamic->init();
    
    
    
    /* creating random h*/
    Random rand(seed+atoms->comm_rank);
    
    int natms_lcl=atoms->natms_lcl;
    type0* __h=h.vecs[0]->begin();
    for(int i=0;i<__dim__*natms_lcl;i++)
        __h[i]=(2.0*rand.uniform()-1.0)*__max_dx;
    
    Algebra::DoLT<__dim__>::func([this,&rand,&__max_st](int i,int j)
    {h.A[i][j]=(2.0*rand.uniform()-1.0)*__max_st;});
    
    
    x0=x;
    force_calc();
    f_h=f*h;
    
    
    prep();
    type0 alpha=0.0;
    type0 dalpha=1.0/static_cast<type0>(__n_desc);
    
    type0 u0=F(alpha);
    type0 du=0.0;
    type0 du_alpha=0.0;
    
    ThermoDynamics thermo(12,"alpha",alpha,"delta_U",du,"dU*alpha",du_alpha);
    thermo.init();
    for(int i=0;i<__n_desc+1;i++)
    {
        du=F(alpha)-u0;
        du_alpha=-alpha*f_h;
        thermo.print(i);
        alpha+=dalpha;
    }
    thermo.fin();
    
    dynamic->fin();
    delete dynamic;
    dynamic=NULL;

    x_d.~VecTens();
    x0.~VecTens();
    h.~VecTens();
    f.~VecTens();
    x.~VecTens();
    Algebra::V_eq<__dim__*__dim__>(&__H_dof[0][0],&H_dof[0][0]);
    chng_box=__chng_box;
}
/*--------------------------------------------
 min
 --------------------------------------------*/
void MinCG::run(int nsteps)
{
    if(dynamic_cast<LineSearchGoldenSection*>(ls))
        return run(dynamic_cast<LineSearchGoldenSection*>(ls),nsteps);
    
    if(dynamic_cast<LineSearchBrent*>(ls))
        return run(dynamic_cast<LineSearchBrent*>(ls),nsteps);
    
    if(dynamic_cast<LineSearchBackTrack*>(ls))
        return run(dynamic_cast<LineSearchBackTrack*>(ls),nsteps);
}
/*--------------------------------------------
 
 --------------------------------------------*/
type0 MinCG::F(type0 alpha)
{
    x=x0+alpha*x_d;
    if(chng_box)
    atoms->update_H();
    
    dynamic->update(atoms->x);
    return ff->value_timer();
}
/*--------------------------------------------
 inner product of f and h
 --------------------------------------------*/
type0 MinCG::dF(type0 alpha,type0& drev)
{
    x=x0+alpha*x_d;
    if(chng_box)
        atoms->update_H();
    
    dynamic->update(atoms->x);
    force_calc();
    
    drev=-(f*h);
    return atoms->pe;
}
/*--------------------------------------------
 find maximum h
 lets find the next sensible number 
 
 x=x_0+h*alpha
 (x-x0)/alpha=sqrt(eps)/alpha

 --------------------------------------------*/
void MinCG::ls_prep(type0& dfa,type0& h_norm,type0& max_a)
{
    
    h_norm=h*h;
    
    dfa=-f_h;
    
    if(h_norm==0.0)
    {
        max_a=0.0;
        dfa=0.0;
        return;
    }
    
    if(dfa>=0.0)
    {
        max_a=0.0;
        dfa=1.0;
        return;
    }
    
    h_norm=sqrt(h_norm);
    
    type0 max_x_d_lcl=0.0;
    type0 max_x_d;
    type0* x_dvec=x_d.vecs[0]->begin();
    const int n=atoms->natms_lcl*__dim__;
    for(int i=0;i<n;i++)
    max_x_d_lcl=MAX(max_x_d_lcl,fabs(x_dvec[i]));
    
    MPI_Allreduce(&max_x_d_lcl,&max_x_d,1,Vec<type0>::MPI_T,MPI_MAX,atoms->world);
    max_a=fabs(max_dx/max_x_d);
    
}
/*--------------------------------------------
 reset to initial position
 --------------------------------------------*/
void MinCG::F_reset()
{
    x=x0;
    if(chng_box) atoms->update_H();
    dynamic->update(atoms->x);
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
    TypeObject.tp_name="mapp.md.min_cg";
    TypeObject.tp_doc=R"---(
    __init__(e_tol=1.0e-8,H_dof=[[False],[False,False],[False,False,False]],affine=False,max_dx=1.0,ls=mapp.ls_bt())
    
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
    ls : mapp.ls
       Line search method
    
    Notes
    -----
    Cojugate Gradient (CG) algorithm for minimization, see :cite:`press_numerical_2007`.
    
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
    
    int ichk=PyType_Ready(&TypeObject);
    if(ichk<0) return ichk;
    Py_INCREF(&TypeObject);
    /*
     this is a shitty hack since python does not have a slot 
     for __init__, __new__, __call__, and etc. they use
     a wrapper_desriptor with a default doc here I change it
     */
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
    getset_affine(getset[2]);
    getset_max_dx(getset[3]);
    getset_ls(getset[4]);
    getset_ntally(getset[5]);
    getset_export(getset[6]);
}
/*--------------------------------------------*/
PyMethodDef MinCG::methods[]=EmptyPyMethodDef(3);
/*--------------------------------------------*/
void MinCG::setup_tp_methods()
{
    ml_run(methods[0]);
    ml_ff_test(methods[1]);
}
/*--------------------------------------------
 
 --------------------------------------------*/
void MinCG::getset_export(PyGetSetDef& getset)
{
    getset.name=(char*)"export";
    getset.doc=(char*)R"---(
    (mapp.md.export) export object
    
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
    atoms : mapp.md.atoms
        System of interest
    max_nsteps : int
        Maximum number of steps to achieve energy minimization
        
    Returns
    -------
    None

    )---";
}
/*--------------------------------------------
 
 --------------------------------------------*/
void MinCG::ml_ff_test(PyMethodDef& tp_methods)
{
    tp_methods.ml_flags=METH_VARARGS | METH_KEYWORDS;
    tp_methods.ml_name="ff_test";
    tp_methods.ml_meth=(PyCFunction)(PyCFunctionWithKeywords)(
    [](PyObject* self,PyObject* args,PyObject* kwds)->PyObject*
    {
        Object* __self=reinterpret_cast<Object*>(self);
        FuncAPI<OP<AtomsMD>,int,type0,type0,int> f("ff_test",{"atoms","seed","max_dx","max_strain","N"});
        f.logics<1>()[0]=VLogics("ge",0);
        f.logics<2>()[0]=VLogics("ge",0.0);
        f.logics<3>()[0]=VLogics("ge",0.0);
        f.logics<4>()[0]=VLogics("gt",0);
        if(f(args,kwds)) return NULL;
        
        AtomsMD* __atoms=reinterpret_cast<AtomsMD::Object*>(f.val<0>().ob)->atoms;
        ForceFieldMD* __ff=reinterpret_cast<AtomsMD::Object*>(f.val<0>().ob)->ff;
        
        __self->min->atoms=__atoms;
        __self->min->ff=__ff;
        
        __self->min->ff_test(f.val<1>(),f.val<2>(),f.val<3>(),f.val<4>());
        
        
        
        __self->min->ff=NULL;
        __self->min->atoms=NULL;
        
        Py_RETURN_NONE;
    });


    tp_methods.ml_doc=(char*)R"---(
    ff_test,(atoms,seed,max_dx,max_strain,N)
   
    for the purpose of testing
    
    Parameters
    ----------
    atoms : mapp.md.atoms
        System of interest
    max_nsteps : int
        Maximum number of steps to achieve energy minimization
        
    Returns
    -------
    None

    )---";
}
