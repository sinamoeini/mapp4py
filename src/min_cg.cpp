#include "min_cg.h"
#include <stdlib.h>
#include "MAPP.h"
#include "ff.h"
#include "ff_md.h"
#include "thermo_dynamics.h"
#include "dynamic_md.h"
#include "atoms_styles.h"
using namespace MAPP_NS;
LineSearch* MinCG::ls=NULL;
/*--------------------------------------------
 constructor
 --------------------------------------------*/
MinCG::MinCG(AtomsMD* __atoms,ForceFieldMD* __ff):
atoms(__atoms),
ff(__ff),
world(__atoms->comm.world),
Min(__atoms,__ff)
{
}
/*--------------------------------------------
 destructor
 --------------------------------------------*/
MinCG::~MinCG()
{
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
    const int natms=atoms->natms;
    type0 MLT[__dim__][__dim__];
    Algebra::MLT_mul_MLT(atoms->B,f.A,MLT);
    type0* xvec=x0.vecs[0]->begin();
    type0* hvec=h.vecs[0]->begin();
    for(int iatm=0;iatm<natms;iatm++,xvec+=__dim__,hvec+=__dim__)
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
    
    if(!ls) ls=new LineSearchBrent();
}
/*--------------------------------------------
 finishing minimization
 --------------------------------------------*/
void MinCG::fin()
{
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
    init();
        
    if(atoms->x_d)
        dynamic=new DynamicMD(atoms,ff,chng_box,{atoms->elem},{h.vecs[0],x0.vecs[0],f0.vecs[0]},{atoms->x_d});
    else
        dynamic=new DynamicMD(atoms,ff,chng_box,{atoms->elem},{h.vecs[0],x0.vecs[0],f0.vecs[0]});
    if(atoms->dof)
        dynamic->add_xchng(atoms->dof);
    
    dynamic->init();
    force_calc();
    type0 S[__dim__][__dim__];
    
    ThermoDynamics thermo(6,
    "PE",ff->nrgy_strss[0],
    "S[0][0]",S[0][0],
    "S[1][1]",S[1][1],
    "S[2][2]",S[2][2],
    "S[1][2]",S[2][1],
    "S[2][0]",S[2][0],
    "S[0][1]",S[1][0]);
    
    thermo.init();
    Algebra::DyadicV_2_MLT(&ff->nrgy_strss[1],S);
    thermo.print(0);
    
    type0 e_prev,e_curr=ff->nrgy_strss[0];
    type0 f0_f0,f_f,f_f0;
    type0 ratio,alpha;
    int err=nsteps==0? MIN_F_MAX_ITER:LS_S;
    h=f;
    f0_f0=f*f;
    int istep=0;
    for(;istep<nsteps && err==LS_S;istep++)
    {
        if(f0_f0==0.0)
        {
            err=LS_F_GRAD0;
            continue;
        }
        
        x0=x;
        f0=f;
        
        e_prev=e_curr;
        
        
        f_h=f*h;
        if(f_h<0.0)
        {
            h=f;
            f_h=f0_f0;
        }
        if(affine) prepare_affine_h();
        err=ls->line_min(this,e_curr,alpha,1);
        
        if(err!=LS_S)
            continue;
        
        force_calc();
        
        if(e_prev-e_curr<e_tol)
            err=MIN_S_TOLERANCE;
        
        if(istep+1==nsteps)
            err=MIN_F_MAX_ITER;
        
        if((istep+1)%ntally==0)
        {
            Algebra::DyadicV_2_MLT(&ff->nrgy_strss[1],S);
            thermo.print(istep+1);
        }
        
        if(err) continue;
        
        f_f=f*f;
        f_f0=f*f0;
        
        ratio=(f_f-f_f0)/(f0_f0);
        
        h=ratio*h+f;
        f0_f0=f_f;
    }
    
    if(istep%ntally)
    {
        Algebra::DyadicV_2_MLT(&ff->nrgy_strss[1],S);
        thermo.print(istep);
    }

    thermo.fin();
    dynamic->fin();
    delete dynamic;
    dynamic=NULL;
    fin();
    
    fprintf(MAPP::mapp_out,"%s",err_msgs[err]);
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
    const int n=atoms->natms*__dim__;
    for(int i=0;i<n;i++)
        max_h_lcl=MAX(max_h_lcl,fabs(hvec[i]));
    
    MPI_Allreduce(&max_h_lcl,&max_h,1,Vec<type0>::MPI_T,MPI_MAX,world);
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
    FuncAPI<OP<AtomsMD>> f("__init__",{"atoms"});
    
    if(f(args,kwds)==-1) return -1;
    Object* __self=reinterpret_cast<Object*>(self);
    AtomsMD::Object* atoms=reinterpret_cast<AtomsMD::Object*>(f.val<0>().ob);
    __self->min=new MinCG(atoms->atoms,atoms->ff);
    __self->atoms=atoms;
    Py_INCREF(atoms);
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
    __self->atoms=NULL;
    return reinterpret_cast<PyObject*>(__self);
}
/*--------------------------------------------
 
 --------------------------------------------*/
void MinCG::__dealloc__(PyObject* self)
{
    Object* __self=reinterpret_cast<Object*>(self);
    delete __self->min;
    __self->min=NULL;
    if(__self->atoms) Py_DECREF(__self->atoms);
    __self->atoms=NULL;
    delete __self;
}
/*--------------------------------------------*/
PyMethodDef MinCG::methods[]={[0 ... 1]={NULL}};
/*--------------------------------------------*/
void MinCG::setup_tp_methods()
{
    ml_run(methods[0]);
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
PyGetSetDef MinCG::getset[]={[0 ... 5]={NULL,NULL,NULL,NULL,NULL}};
/*--------------------------------------------*/
void MinCG::setup_tp_getset()
{
    getset_max_dx(getset[0]);
    getset_e_tol(getset[1]);
    getset_affine(getset[2]);
    getset_H_dof(getset[3]);
    getset_ntally(getset[4]);
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
        FuncAPI<int> f("run",{"max_nsteps"});
        f.logics<0>()[0]=VLogics("ge",0);
        if(f(args,kwds)) return NULL;
        __self->min->run(f.val<0>());
        Py_RETURN_NONE;
    };
}
