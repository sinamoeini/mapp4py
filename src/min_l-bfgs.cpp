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
#include <stdlib.h>
#include "min_l-bfgs.h"
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
MinLBFGS::MinLBFGS(AtomsMD*& __atoms,ForceFieldMD*& __ff,int __m):
MinCG(__atoms,__ff),
m(__m),
s(NULL),
y(NULL),
alpha(NULL),
rho(NULL)
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
    MinCG::init();
    
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
        Memory::alloc(rho,m);
        Memory::alloc(alpha,m);
    }
}
/*--------------------------------------------
 run
 --------------------------------------------*/
void MinLBFGS::run(int nsteps)
{
    init();
    
    if(atoms->x_d)
        dynamic=new DynamicMD(atoms,ff,chng_box,{atoms->elem},{h.vecs[0],x0.vecs[0],f0.vecs[0]},{atoms->x_d});
    else
        dynamic=new DynamicMD(atoms,ff,chng_box,{atoms->elem},{h.vecs[0],x0.vecs[0],f0.vecs[0]});
    
    if(atoms->dof)
        dynamic->add_xchng(atoms->dof);
    
    for(int i=0;i<m;i++)
    {
        dynamic->add_xchng(s[i].vecs[0]);
        dynamic->add_xchng(y[i].vecs[0]);
    }
    
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
        
    type0 alpha_m,gamma;
    type0 inner0,inner1;
    
    
    
    int k=0;
    gamma=1.0;
    int err=LS_S;
    

    int istep=0;
    for(;istep<nsteps && err==LS_S;istep++)
    {
        x0=x;
        h=f0=f;
        
        for(int i=0;i<k;i++)
        {
            alpha[i]=-rho[i]*(s[i]*h);
            h+=alpha[i]*y[i];
        }
        
        h*=gamma;
        
        for(int i=k-1;i>-1;i--)
        h+=(-alpha[i]-rho[i]*(y[i]*h))*s[i];
        
        e_prev=e_curr;
        
        
        f_h=f*h;
        if(f_h<0.0)
        {
            h=f;
            k=0;
            f_h=f*f;
        }
        if(affine) prepare_affine_h();
        err=ls->line_min(this,e_curr,alpha_m,0);
        
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
        
        if(m)
        {
            if(k!=m) k++;
            
            s[0].cyclic_shift(k);
            y[0].cyclic_shift(k);
            
            for(int i=m-1;i>0;i--)
                rho[i]=rho[i-1];
            
            s[0]=x-x0;
            y[0]=f0-f;
            
            inner0=s[0]*y[0];
            inner1=y[0]*y[0];
            
            gamma=inner0/inner1;
            rho[0]=1.0/inner0;
        }
        else
            gamma=(x*f0-x*f-x0*f0+x0*f)/(f*f+f0*f0-2.0*(f*f0));
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
 fin after a run
 --------------------------------------------*/
void MinLBFGS::fin()
{
    Memory::dealloc(rho);
    Memory::dealloc(alpha);
    rho=alpha=NULL;
    
    delete [] s;
    delete [] y;
    s=y=NULL;
    MinCG::fin();
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
    FuncAPI<OP<AtomsMD>,int> f("__init__",{"atoms","m"});
    f.noptionals=1;
    f.logics<1>()[0]=VLogics("ge",0);
    f.val<1>()=2;
    
    
    if(f(args,kwds)==-1) return -1;
    Object* __self=reinterpret_cast<Object*>(self);
    AtomsMD::Object* atoms=reinterpret_cast<AtomsMD::Object*>(f.val<0>().ob);
    __self->min=new MinLBFGS(atoms->atoms,atoms->ff,f.val<1>());
    __self->atoms=atoms;
    Py_INCREF(atoms);
    return 0;
}
/*--------------------------------------------
 
 --------------------------------------------*/
PyObject* MinLBFGS::__alloc__(PyTypeObject* type,Py_ssize_t)
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
void MinLBFGS::__dealloc__(PyObject* self)
{
    Object* __self=reinterpret_cast<Object*>(self);
    delete __self->min;
    __self->min=NULL;
    if(__self->atoms) Py_DECREF(__self->atoms);
    __self->atoms=NULL;
    delete __self;
}
/*--------------------------------------------*/
PyTypeObject MinLBFGS::TypeObject ={PyObject_HEAD_INIT(NULL)};
/*--------------------------------------------*/
void MinLBFGS::setup_tp()
{
    TypeObject.tp_name="min_lbfgs";
    TypeObject.tp_doc="l-BFGS minimization";
    
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
}
/*--------------------------------------------*/
PyMethodDef MinLBFGS::methods[]={[0 ... 0]={NULL}};
/*--------------------------------------------*/
void MinLBFGS::setup_tp_methods()
{
}
/*--------------------------------------------*/
PyGetSetDef MinLBFGS::getset[]={[0 ... 1]={NULL,NULL,NULL,NULL,NULL}};
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
    getset.doc=(char*)"number of vectors to store in memory";
    getset.get=[](PyObject* self,void*)->PyObject*
    {
        return var<int>::build(reinterpret_cast<Object*>(self)->min->m,NULL);
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
