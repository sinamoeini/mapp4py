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
MinLBFGSDMD::MinLBFGSDMD(int __m):
MinCGDMD(),
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
    if(dynamic_cast<LineSearchGoldenSection*>(ls))
        return run(dynamic_cast<LineSearchGoldenSection*>(ls),nsteps);
    
    if(dynamic_cast<LineSearchBrent*>(ls))
        return run(dynamic_cast<LineSearchBrent*>(ls),nsteps);
    
    if(dynamic_cast<LineSearchBackTrack*>(ls))
        return run(dynamic_cast<LineSearchBackTrack*>(ls),nsteps);
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
    FuncAPI<int> f("__init__",{"m"});
    f.noptionals=1;
    f.logics<0>()[0]=VLogics("ge",0);
    f.val<0>()=2;
    
    
    if(f(args,kwds)==-1) return -1;
    Object* __self=reinterpret_cast<Object*>(self);
    __self->min=new MinLBFGSDMD(f.val<0>());
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
    return reinterpret_cast<PyObject*>(__self);
}
/*--------------------------------------------
 
 --------------------------------------------*/
void MinLBFGSDMD::__dealloc__(PyObject* self)
{
    Object* __self=reinterpret_cast<Object*>(self);
    delete __self->min;
    __self->min=NULL;
    delete __self;
}
/*--------------------------------------------*/
PyTypeObject MinLBFGSDMD::TypeObject ={PyObject_HEAD_INIT(NULL)};
/*--------------------------------------------*/
void MinLBFGSDMD::setup_tp()
{
    TypeObject.tp_name="mapp.dmd.min_lbfgs";
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


