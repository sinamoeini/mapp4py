#include "min.h"
#include "ls.h"
#include "ff_styles.h"
#include "xmath.h"
using namespace MAPP_NS;
const char* Min::err_msgs[]=
{
    [LS_F_DOWNHILL]="line search failed: not downhill direction\n",
    [LS_F_GRAD0]="line search failed: gradient is zero\n",
    [LS_MIN_ALPHA]="line search failed: minimum alpha reached\n",
    [MIN_S_TOLERANCE]="minimization finished: energy tolerance reached\n",
    [MIN_F_MAX_ITER]="minimization finished: maximum iteration reached\n",
    [B_F_DOWNHILL]="bracketing failed: not downhill direction\n",
    [B_F_MAX_ALPHA]="bracketing failed: maximum alpha reached\n"
};
/*--------------------------------------------
 constructor
 --------------------------------------------*/
Min::Min(Atoms* __atoms,ForceField* __ff):
atoms(__atoms),
ff(__ff),
H_dof{[0 ... __dim__-1]={[0 ... __dim__-1]=false}},
chng_box(false),
e_tol(sqrt(std::numeric_limits<type0>::epsilon())),
affine(false),
max_dx(1.0),
ntally(1000)
{
}
/*--------------------------------------------
 destructor
 --------------------------------------------*/
Min::~Min()
{
}
/*------------------------------------------------------------------------------------------------------------------------------------
 
 ------------------------------------------------------------------------------------------------------------------------------------*/
PyObject* Min::__new__(PyTypeObject* type,PyObject* args,PyObject* kwds)
{
    Object* __self=reinterpret_cast<Object*>(type->tp_alloc(type,0));
    PyObject* self=reinterpret_cast<PyObject*>(__self);
    return self;
}
/*--------------------------------------------
 
 --------------------------------------------*/
int Min::__init__(PyObject* self,PyObject* args,PyObject* kwds)
{
    FuncAPI<OP<Atoms>> f("__init__",{"atoms"});
    
    if(f(args,kwds)==-1) return -1;
    Object* __self=reinterpret_cast<Object*>(self);
    Atoms::Object* atoms=reinterpret_cast<Atoms::Object*>(f.val<0>().ob);
    __self->min=new Min(atoms->atoms,atoms->ff);
    __self->atoms=atoms;
    Py_INCREF(atoms);
    return 0;
}
/*--------------------------------------------
 
 --------------------------------------------*/
PyObject* Min::__alloc__(PyTypeObject* type,Py_ssize_t)
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
void Min::__dealloc__(PyObject* self)
{
    Object* __self=reinterpret_cast<Object*>(self);
    delete __self->min;
    __self->min=NULL;
    if(__self->atoms) Py_DECREF(__self->atoms);
    __self->atoms=NULL;
    delete __self;
}
/*--------------------------------------------*/
PyMethodDef Min::methods[]={[0 ... 0]={NULL}};
/*--------------------------------------------*/
void Min::setup_tp_methods()
{
}
/*--------------------------------------------*/
PyTypeObject Min::TypeObject={PyObject_HEAD_INIT(NULL)};
/*--------------------------------------------*/
void Min::setup_tp()
{
    TypeObject.tp_name="min_cg";
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
PyGetSetDef Min::getset[]={[0 ... 5]={NULL,NULL,NULL,NULL,NULL}};
/*--------------------------------------------*/
void Min::setup_tp_getset()
{
    getset_max_dx(getset[0]);
    getset_e_tol(getset[1]);
    getset_affine(getset[2]);
    getset_H_dof(getset[3]);
    getset_ntally(getset[4]);
}
/*--------------------------------------------
 
 --------------------------------------------*/
void Min::getset_affine(PyGetSetDef& getset)
{
    getset.name=(char*)"affine";
    getset.doc=(char*)"set to true if transformation is affine";
    getset.get=[](PyObject* self,void*)->PyObject*
    {
        return var<bool>::build(reinterpret_cast<Object*>(self)->min->affine,NULL);
    };
    getset.set=[](PyObject* self,PyObject* op,void*)->int
    {
        VarAPI<bool> affine("affine");
        int ichk=affine.set(op);
        if(ichk==-1) return -1;
        reinterpret_cast<Object*>(self)->min->affine=affine.val;
        return 0;
    };
}
/*--------------------------------------------
 
 --------------------------------------------*/
void Min::getset_H_dof(PyGetSetDef& getset)
{
    getset.name=(char*)"H_dof";
    getset.doc=(char*)"unitcell degrees of freedom";
    getset.get=[](PyObject* self,void*)->PyObject*
    {
        return var<symm<bool[__dim__][__dim__]>>::build(reinterpret_cast<Object*>(self)->min->H_dof,NULL);
    };
    getset.set=[](PyObject* self,PyObject* op,void*)->int
    {
        VarAPI<symm<bool[__dim__][__dim__]>> H_dof("H_dof");
        int ichk=H_dof.set(op);
        if(ichk==-1) return -1;
        
        bool (&__H_dof)[__dim__][__dim__]=reinterpret_cast<Object*>(self)->min->H_dof;
        bool chng_box=false;
        Algebra::DoLT<__dim__>::func([&H_dof,&__H_dof,&chng_box](int i,int j)
        {
            __H_dof[i][j]=__H_dof[j][i]=H_dof.val[i][j];
            if(__H_dof[i][j]) chng_box=true;
        });
        reinterpret_cast<Object*>(self)->min->chng_box=chng_box;
        return 0;
    };
}
/*--------------------------------------------
 
 --------------------------------------------*/
void Min::getset_e_tol(PyGetSetDef& getset)
{
    getset.name=(char*)"e_tol";
    getset.doc=(char*)"energy tolerance criterion for stopping minimization";
    getset.get=[](PyObject* self,void*)->PyObject*
    {
        return var<type0>::build(reinterpret_cast<Object*>(self)->min->e_tol,NULL);
    };
    getset.set=[](PyObject* self,PyObject* op,void*)->int
    {
        VarAPI<type0> e_tol("e_tol");
        e_tol.logics[0]=VLogics("ge",0.0);
        int ichk=e_tol.set(op);
        if(ichk==-1) return -1;
        reinterpret_cast<Object*>(self)->min->e_tol=e_tol.val;
        return 0;
    };
}
/*--------------------------------------------
 
 --------------------------------------------*/
void Min::getset_max_dx(PyGetSetDef& getset)
{
    getset.name=(char*)"max_dx";
    getset.doc=(char*)"maximum displacement";
    getset.get=[](PyObject* self,void*)->PyObject*
    {
        return var<type0>::build(reinterpret_cast<Object*>(self)->min->max_dx,NULL);
    };
    getset.set=[](PyObject* self,PyObject* op,void*)->int
    {
        VarAPI<type0> max_dx("max_dx");
        max_dx.logics[0]=VLogics("gt",0.0);
        int ichk=max_dx.set(op);
        if(ichk==-1) return -1;
        reinterpret_cast<Object*>(self)->min->max_dx=max_dx.val;
        return 0;
    };
}
/*--------------------------------------------
 
 --------------------------------------------*/
void Min::getset_ntally(PyGetSetDef& getset)
{
    getset.name=(char*)"ntally";
    getset.doc=(char*)"tally thermodynamic quantities every ntally steps";
    getset.get=[](PyObject* self,void*)->PyObject*
    {
        return var<int>::build(reinterpret_cast<Object*>(self)->min->ntally,NULL);
    };
    getset.set=[](PyObject* self,PyObject* op,void*)->int
    {
        VarAPI<int> ntally("ntally");
        ntally.logics[0]=VLogics("gt",0);
        int ichk=ntally.set(op);
        if(ichk==-1) return -1;
        reinterpret_cast<Object*>(self)->min->ntally=ntally.val;
        return 0;
    };
}


