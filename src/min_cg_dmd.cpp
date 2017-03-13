#include "min_cg_dmd.h"
#include <stdlib.h>
using namespace MAPP_NS;
/*--------------------------------------------
 constructor
 --------------------------------------------*/
MinCGDMD::MinCGDMD():
atoms(NULL),
ff(NULL),
max_dalpha(0.01),
Min()
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
 
 --------------------------------------------*/
void MinCGDMD::force_calc()
{
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
void MinCGDMD::prepare_affine_h()
{
    if(!chng_box) return;
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
type0 MinCGDMD::calc_ndofs()
{
    return 0.0;
}
/*--------------------------------------------
 init before a run
 --------------------------------------------*/
void MinCGDMD::init()
{
    const int c_dim=atoms->c->dim;
    x.~VecTens();
    new (&x) VecTens<type0,2>(atoms,chng_box,atoms->H,atoms->x,atoms->alpha);
    f.~VecTens();
    new (&f) VecTens<type0,2>(atoms,chng_box,ff->f,ff->f_alpha);
    h.~VecTens();
    new (&h) VecTens<type0,2>(atoms,chng_box,__dim__,c_dim);
    x0.~VecTens();
    new (&x0) VecTens<type0,2>(atoms,chng_box,__dim__,c_dim);
    f0.~VecTens();
    new (&f0) VecTens<type0,2>(atoms,chng_box,__dim__,c_dim);
    
    dynamic=new DynamicDMD(atoms,ff,chng_box,{},
    {atoms->dof,h.vecs[0],h.vecs[1],x0.vecs[0],x0.vecs[1],f0.vecs[0],f0.vecs[1]},{});
    
    dynamic->init();
    
    uvecs[0]=atoms->x;
    uvecs[1]=atoms->alpha;
}
/*--------------------------------------------
 finishing minimization
 --------------------------------------------*/
void MinCGDMD::fin()
{
    uvecs[1]=NULL;
    uvecs[0]=NULL;
    
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
void MinCGDMD::run(int nsteps)
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
type0 MinCGDMD::F(type0 alpha)
{
    x=x0+alpha*h;
    type0 max_alpha_lcl=0.0;
    const int n=atoms->natms_lcl*atoms->alpha->dim;
    type0* alpha_vec=atoms->alpha->begin();
    type0* c_vec=atoms->c->begin();
    for(int i=0;i<n;i++)
        if(c_vec[i]>0.0) max_alpha_lcl=MAX(max_alpha_lcl,alpha_vec[i]);
    MPI_Allreduce(&max_alpha_lcl,&atoms->max_alpha,1,Vec<type0>::MPI_T,MPI_MAX,atoms->world);
    
    if(chng_box)
        atoms->update_H();
    
    dynamic->update(uvecs,2);
    return ff->value_timer();
}
/*--------------------------------------------
 inner product of f and h
 --------------------------------------------*/
type0 MinCGDMD::dF(type0 alpha,type0& drev)
{
    x=x0+alpha*h;
    type0 max_alpha_lcl=0.0;
    const int n=atoms->natms_lcl*atoms->alpha->dim;
    type0* alpha_vec=atoms->alpha->begin();
    type0* c_vec=atoms->c->begin();
    for(int i=0;i<n;i++)
        if(c_vec[i]>0.0) max_alpha_lcl=MAX(max_alpha_lcl,alpha_vec[i]);
    MPI_Allreduce(&max_alpha_lcl,&atoms->max_alpha,1,Vec<type0>::MPI_T,MPI_MAX,atoms->world);
    
    if(chng_box)
        atoms->update_H();
    
    dynamic->update(uvecs,2);
    force_calc();
    
    drev=-(f*h);
    return ff->nrgy_strss[0];
}
/*--------------------------------------------
 find maximum h
 --------------------------------------------*/
void MinCGDMD::ls_prep(type0& dfa,type0& h_norm,type0& max_a)
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
    type0* hvec=h.vecs[0]->begin();
    const int natms_lcl=atoms->natms_lcl;
    int n=natms_lcl*__dim__;
    for(int i=0;i<n;i++)
        max_h_lcl=MAX(max_h_lcl,fabs(hvec[i]));

    
    type0 max_alpha_ratio=std::numeric_limits<type0>::infinity();
    type0* h_alphavec=h.vecs[1]->begin();
    type0* alphavec=x.vecs[1]->begin();
    n=natms_lcl*atoms->alpha->dim;
    type0 max_h_alpha_lcl=0.0;
    for(int i=0;i<n;i++)
    {
        if(h_alphavec[i]<0.0)
            max_alpha_ratio=MIN(-alphavec[i]/h_alphavec[i],max_alpha_ratio);
        max_h_alpha_lcl=MAX(max_h_alpha_lcl,fabs(h_alphavec[i]));
    }
    type0 max_a_lcl=MIN(fabs(max_dx/max_h_lcl),max_alpha_ratio);
    max_a_lcl=MIN(max_a_lcl,fabs(max_dalpha/max_h_alpha_lcl));

    MPI_Allreduce(&max_a_lcl,&max_a,1,Vec<type0>::MPI_T,MPI_MIN,atoms->world);
    
}
/*--------------------------------------------
 reset to initial position
 --------------------------------------------*/
void MinCGDMD::F_reset()
{
    x=x0;
    if(chng_box) atoms->update_H();
    dynamic->update(atoms->x);
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
    FuncAPI<> f("__init__");
    
    if(f(args,kwds)==-1) return -1;
    Object* __self=reinterpret_cast<Object*>(self);
    __self->min=new MinCGDMD();
    __self->ls=NULL;
    __self->xprt=NULL;
    return 0;
}
/*--------------------------------------------
 
 --------------------------------------------*/
PyObject* MinCGDMD::__alloc__(PyTypeObject* type,Py_ssize_t)
{
    Object* __self=new Object;
    __self->ob_type=type;
    __self->ob_refcnt=1;
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
void MinCGDMD::setup_tp()
{
    TypeObject.tp_name="mapp.dmd.min_cg";
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
PyGetSetDef MinCGDMD::getset[]={[0 ... 7]={NULL,NULL,NULL,NULL,NULL}};
/*--------------------------------------------*/
void MinCGDMD::setup_tp_getset()
{
    getset_max_dx(getset[0]);
    getset_e_tol(getset[1]);
    getset_affine(getset[2]);
    getset_H_dof(getset[3]);
    getset_ntally(getset[4]);
    getset_ls(getset[5]);
    getset_max_dalpha(getset[6]);
}
/*--------------------------------------------*/
PyMethodDef MinCGDMD::methods[]={[0 ... 1]={NULL}};
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
    getset.doc=(char*)"maximum alpha displacement";
    getset.get=[](PyObject* self,void*)->PyObject*
    {
        return var<type0>::build(reinterpret_cast<Object*>(self)->min->max_dx,NULL);
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
void MinCGDMD::ml_run(PyMethodDef& tp_methods)
{
    tp_methods.ml_flags=METH_VARARGS | METH_KEYWORDS;
    tp_methods.ml_name="run";
    tp_methods.ml_doc="run energy minimization";
    
    tp_methods.ml_meth=(PyCFunction)(PyCFunctionWithKeywords)
    [](PyObject* self,PyObject* args,PyObject* kwds)->PyObject*
    {
        Object* __self=reinterpret_cast<Object*>(self);
        FuncAPI<OP<AtomsDMD>,int> f("run",{"atoms","max_nsteps"});
        f.logics<1>()[0]=VLogics("ge",0);
        if(f(args,kwds)) return NULL;
        
        AtomsDMD* __atoms=reinterpret_cast<AtomsDMD::Object*>(f.val<0>().ob)->atoms;
        ForceFieldDMD* __ff=reinterpret_cast<AtomsDMD::Object*>(f.val<0>().ob)->ff;
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
