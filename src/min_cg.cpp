#include "min_cg.h"
#include <stdlib.h>
using namespace MAPP_NS;
/*--------------------------------------------
 constructor
 --------------------------------------------*/
MinCG::MinCG():
atoms(NULL),
ff(NULL),
Min()
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
void MinCG::force_calc()
{
    ff->reset();
    if(chng_box)
    {
        ff->derivative_timer(f.A);
        Algebra::DoLT<__dim__>::func([this](int i,int j)
        {
            f.A[i][j]*=H_dof[i][j];
        });
    }
    else
        ff->derivative_timer();
}
/*--------------------------------------------
 
 --------------------------------------------*/
void MinCG::prepare_affine_h()
{
    const int natms_lcl=atoms->natms_lcl;
    type0 MLT[__dim__][__dim__];
    Algebra::MLT_mul_MLT(atoms->B,f.A,MLT);
    type0* xvec=x0.vecs[0]->begin();
    type0* hvec=h.vecs[0]->begin();
    for(int iatm=0;iatm<natms_lcl;iatm++,xvec+=__dim__,hvec+=__dim__)
        Algebra::V_mul_MLT(xvec,MLT,hvec);
}
/*--------------------------------------------
 
 --------------------------------------------*/
type0 MinCG::calc_ndofs()
{
    return 0.0;
}
/*--------------------------------------------
 init before a run
 --------------------------------------------*/
void MinCG::init()
{
    x.~VecTens();
    new (&x) VecTens<type0,1>(atoms,chng_box,atoms->H,atoms->x);
    f.~VecTens();
    new (&f) VecTens<type0,1>(atoms,chng_box,ff->f);
    h.~VecTens();
    new (&h) VecTens<type0,1>(atoms,chng_box,__dim__);
    x0.~VecTens();
    new (&x0) VecTens<type0,1>(atoms,chng_box,__dim__);
    f0.~VecTens();
    new (&f0) VecTens<type0,1>(atoms,chng_box,__dim__);
    
    dynamic=new DynamicMD(atoms,ff,chng_box,{},{atoms->dof,h.vecs[0],x0.vecs[0],f0.vecs[0]},{atoms->x_d});
    dynamic->init();
}
/*--------------------------------------------
 finishing minimization
 --------------------------------------------*/
void MinCG::fin()
{
    dynamic->fin();
    delete dynamic;
    dynamic=NULL;
    
    f0.~VecTens();
    x0.~VecTens();
    h.~VecTens();
    f.~VecTens();
    x.~VecTens();
}
/*--------------------------------------------
 min
 --------------------------------------------*/
void MinCG::run(int nsteps)
{
    if(ls)
    {
        if(dynamic_cast<LineSearchGoldenSection*>(ls))
            return run(dynamic_cast<LineSearchGoldenSection*>(ls),nsteps);
        
        if(dynamic_cast<LineSearchBrent*>(ls))
            return run(dynamic_cast<LineSearchBrent*>(ls),nsteps);
        
        if(dynamic_cast<LineSearchBackTrack*>(ls))
            return run(dynamic_cast<LineSearchBackTrack*>(ls),nsteps);
    }
    LineSearchBackTrack* __ls=new LineSearchBackTrack();
    run(__ls,nsteps);
    delete __ls;
}

/*--------------------------------------------
 
 --------------------------------------------*/
type0 MinCG::F(type0 alpha)
{
    x=x0+alpha*h;
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
    x=x0+alpha*h;
    if(chng_box)
        atoms->update_H();
    
    dynamic->update(atoms->x);
    force_calc();
    
    drev=-(f*h);
    return ff->nrgy_strss[0];
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
    
    //type0 max_a_lcl=max_dx;
    type0 max_h_lcl=0.0;
    type0 max_h;
    type0* hvec=h.vecs[0]->begin();
    const int n=atoms->natms_lcl*__dim__;
    for(int i=0;i<n;i++)
        max_h_lcl=MAX(max_h_lcl,fabs(hvec[i]));
    
    MPI_Allreduce(&max_h_lcl,&max_h,1,Vec<type0>::MPI_T,MPI_MAX,atoms->world);
    max_a=fabs(max_dx/max_h);
    
}
/*--------------------------------------------
 reset to initial position
 --------------------------------------------*/
void MinCG::F_reset()
{
    x=x0;
    if(chng_box)
        atoms->update_H();
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
    FuncAPI<> f("__init__");
    
    if(f(args,kwds)==-1) return -1;
    Object* __self=reinterpret_cast<Object*>(self);
    __self->min=new MinCG();
    __self->ls=NULL;
    return 0;
}
/*--------------------------------------------
 
 --------------------------------------------*/
PyObject* MinCG::__alloc__(PyTypeObject* type,Py_ssize_t)
{
    Object* __self=new Object;
    __self->ob_type=type;
    __self->ob_refcnt=1;
    __self->min=NULL;
    __self->ls=NULL;
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
    delete __self;
}
/*--------------------------------------------*/
PyTypeObject MinCG::TypeObject={PyObject_HEAD_INIT(NULL)};
/*--------------------------------------------*/
void MinCG::setup_tp()
{
    TypeObject.tp_name="mapp.md.min_cg";
    TypeObject.tp_doc="conjugate gradient minimization";
    
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
PyGetSetDef MinCG::getset[]={[0 ... 6]={NULL,NULL,NULL,NULL,NULL}};
/*--------------------------------------------*/
void MinCG::setup_tp_getset()
{
    getset_max_dx(getset[0]);
    getset_e_tol(getset[1]);
    getset_affine(getset[2]);
    getset_H_dof(getset[3]);
    getset_ntally(getset[4]);
    getset_ls(getset[5]);
}
/*--------------------------------------------*/
PyMethodDef MinCG::methods[]={[0 ... 1]={NULL}};
/*--------------------------------------------*/
void MinCG::setup_tp_methods()
{
    ml_run(methods[0]);
}
/*--------------------------------------------
 
 --------------------------------------------*/
void MinCG::ml_run(PyMethodDef& tp_methods)
{
    tp_methods.ml_flags=METH_VARARGS | METH_KEYWORDS;
    tp_methods.ml_name="run";
    tp_methods.ml_doc="run energy minimization";
    
    tp_methods.ml_meth=(PyCFunction)(PyCFunctionWithKeywords)
    [](PyObject* self,PyObject* args,PyObject* kwds)->PyObject*
    {
        Object* __self=reinterpret_cast<Object*>(self);
        FuncAPI<OP<AtomsMD>,int> f("run",{"atoms","max_nsteps"});
        f.logics<1>()[0]=VLogics("ge",0);
        if(f(args,kwds)) return NULL;
        
        AtomsMD* __atoms=reinterpret_cast<AtomsMD::Object*>(f.val<0>().ob)->atoms;
        ForceFieldMD* __ff=reinterpret_cast<AtomsMD::Object*>(f.val<0>().ob)->ff;
        try
        {
            __self->min->pre_run_chk(__atoms,__ff);
        }
        catch(std::string err_msg)
        {
            PyErr_SetString(PyExc_TypeError,err_msg.c_str());
            return NULL;
        }
        
        __self->min->atoms=__atoms;
        __self->min->ff=__ff;
        
        
        try
        {
            __self->min->init();
        }
        catch(std::string err_msg)
        {
            __self->min->fin();
            __self->min->ff=NULL;
            __self->min->atoms=NULL;
            PyErr_SetString(PyExc_TypeError,err_msg.c_str());
            return NULL;
        }
        
        __self->min->run(f.val<1>());
        
        __self->min->fin();
        __self->min->ff=NULL;
        __self->min->atoms=NULL;
        
        Py_RETURN_NONE;
    };
}
