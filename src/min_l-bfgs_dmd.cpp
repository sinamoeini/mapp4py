#include "min_l-bfgs_dmd.h"
#include "memory.h"
#include "dynamic_dmd.h"
#include "thermo_dynamics.h"
#include "ff_dmd.h"
#include "MAPP.h"
using namespace MAPP_NS;
/*--------------------------------------------
 constructor
 --------------------------------------------*/
MinLBFGSDMD::MinLBFGSDMD(AtomsDMD* __atoms,ForceFieldDMD* __ff,int __m):
MinCGDMD(__atoms,__ff),
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
MinLBFGSDMD::~MinLBFGSDMD()
{
}
/*--------------------------------------------
 init before a run
 --------------------------------------------*/
void MinLBFGSDMD::init()
{
    MinCGDMD::init();
    const int c_dim=atoms->c->dim;
    if(m)
    {
        s=new VecTens<type0,2>[m];
        y=new VecTens<type0,2>[m];
        for(int i=0;i<m;i++)
        {
            s[i].~VecTens();
            new (s+i) VecTens<type0,2>(atoms,chng_box,__dim__,c_dim);
            y[i].~VecTens();
            new (y+i) VecTens<type0,2>(atoms,chng_box,__dim__,c_dim);
        }
        Memory::alloc(rho,m);
        Memory::alloc(alpha,m);
    }
}
/*--------------------------------------------
 fin after a run
 --------------------------------------------*/
void MinLBFGSDMD::fin()
{
    Memory::dealloc(rho);
    Memory::dealloc(alpha);
    rho=alpha=NULL;
    
    delete [] s;
    delete [] y;
    s=y=NULL;
    MinCGDMD::fin();
}
/*--------------------------------------------
 run
 --------------------------------------------*/
void MinLBFGSDMD::run(int nsteps)
{
    init();
    
    uvecs[0]=atoms->x;
    uvecs[1]=atoms->alpha;
    dynamic=new DynamicDMD(atoms,ff,chng_box,{atoms->elem,atoms->c},
    {h.vecs[0],h.vecs[1],x0.vecs[0],x0.vecs[1],f0.vecs[0],f0.vecs[1]},{});
    
    if(atoms->dof)
        dynamic->add_xchng(atoms->dof);
    
    for(int i=0;i<m;i++)
    {
        dynamic->add_xchng(s[i].vecs[0]);
        dynamic->add_xchng(s[i].vecs[1]);
        dynamic->add_xchng(y[i].vecs[0]);
        dynamic->add_xchng(y[i].vecs[1]);
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
    int err=nsteps==0? MIN_F_MAX_ITER:LS_S;
    

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
/*------------------------------------------------------------------------------------------------------------------------------------
 
 ------------------------------------------------------------------------------------------------------------------------------------*/
PyObject* MinLBFGSDMD::__new__(PyTypeObject* type,PyObject* args,PyObject* kwds)
{
    Object* __self=reinterpret_cast<Object*>(type->tp_alloc(type,0));
    PyObject* self=reinterpret_cast<PyObject*>(__self);
    return self;
}
/*--------------------------------------------
 
 --------------------------------------------*/
int MinLBFGSDMD::__init__(PyObject* self,PyObject* args,PyObject* kwds)
{
    FuncAPI<OP<AtomsDMD>,int> f("__init__",{"atoms","m"});
    f.noptionals=1;
    f.logics<1>()[0]=VLogics("ge",0);
    f.val<1>()=2;
    
    
    if(f(args,kwds)==-1) return -1;
    Object* __self=reinterpret_cast<Object*>(self);
    AtomsDMD::Object* atoms=reinterpret_cast<AtomsDMD::Object*>(f.val<0>().ob);
    __self->min=new MinLBFGSDMD(atoms->atoms,atoms->ff,f.val<1>());
    __self->atoms=atoms;
    Py_INCREF(atoms);
    return 0;
}
/*--------------------------------------------
 
 --------------------------------------------*/
PyObject* MinLBFGSDMD::__alloc__(PyTypeObject* type,Py_ssize_t)
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
void MinLBFGSDMD::__dealloc__(PyObject* self)
{
    Object* __self=reinterpret_cast<Object*>(self);
    delete __self->min;
    __self->min=NULL;
    if(__self->atoms) Py_DECREF(__self->atoms);
    __self->atoms=NULL;
    delete __self;
}
/*--------------------------------------------*/
PyTypeObject MinLBFGSDMD::TypeObject ={PyObject_HEAD_INIT(NULL)};
/*--------------------------------------------*/
void MinLBFGSDMD::setup_tp()
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
    
    TypeObject.tp_base=&MinCGDMD::TypeObject;
}
/*--------------------------------------------*/
PyMethodDef MinLBFGSDMD::methods[]={[0 ... 0]={NULL}};
/*--------------------------------------------*/
void MinLBFGSDMD::setup_tp_methods()
{
}
/*--------------------------------------------*/
PyGetSetDef MinLBFGSDMD::getset[]={[0 ... 1]={NULL,NULL,NULL,NULL,NULL}};
/*--------------------------------------------*/
void MinLBFGSDMD::setup_tp_getset()
{
    getset_m(getset[0]);
}
/*--------------------------------------------
 
 --------------------------------------------*/
void MinLBFGSDMD::getset_m(PyGetSetDef& getset)
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


