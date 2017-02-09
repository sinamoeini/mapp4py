#include "atoms_dmd.h"
#include "xmath.h"
/*--------------------------------------------
 
 --------------------------------------------*/
AtomsDMD::AtomsDMD(MPI_Comm& world,int __N):
N(__N),
Atoms(world),
c(NULL),
c_d(NULL),
elem(NULL),
alpha(NULL),
xi(new type0[__N]),
wi(new type0[__N])
{
    XMath::quadrature_hg(N,xi,wi);
}
/*--------------------------------------------
 
 --------------------------------------------*/
AtomsDMD::AtomsDMD(Communication& comm,int __N):
N(__N),
Atoms(comm),
c(NULL),
c_d(NULL),
elem(NULL),
alpha(NULL),
xi(new type0[__N]),
wi(new type0[__N])
{
    XMath::quadrature_hg(N,xi,wi);
}
/*--------------------------------------------
 
 --------------------------------------------*/
AtomsDMD::~AtomsDMD()
{
    delete [] wi;
    delete [] xi;
    delete alpha;
    delete c;
    delete c_d;
    delete elem;
}
/*------------------------------------------------------------------------------------------------------------------------------------
 
 ------------------------------------------------------------------------------------------------------------------------------------*/
#include "ff_styles.h"
PyObject* AtomsDMD::__new__(PyTypeObject* type,PyObject* args,PyObject* kwds)
{
    Object* __self=reinterpret_cast<Object*>(type->tp_alloc(type,0));
    PyObject* self=reinterpret_cast<PyObject*>(__self);
    return self;
}
/*--------------------------------------------
 
 --------------------------------------------*/
int AtomsDMD::__init__(PyObject* self,PyObject* args,PyObject* kwds)
{
    FuncAPI<> func("__init__");
    if(func(args,kwds)==-1) return -1;
    Object* __self=reinterpret_cast<Object*>(self);
    __self->atoms=NULL;
    __self->ff=NULL;
    return 0;
}
/*--------------------------------------------
 
 --------------------------------------------*/
PyObject* AtomsDMD::__alloc__(PyTypeObject* type,Py_ssize_t)
{
    Object* __self=new Object;
    __self->ob_type=type;
    __self->ob_refcnt=1;
    __self->atoms=NULL;
    __self->ff=NULL;
    return reinterpret_cast<PyObject*>(__self);
}
/*--------------------------------------------
 
 --------------------------------------------*/
void AtomsDMD::__dealloc__(PyObject* self)
{
    Object* __self=reinterpret_cast<Object*>(self);
    delete __self->ff;
    delete __self->atoms;
    delete __self;
}
/*--------------------------------------------*/
PyTypeObject AtomsDMD::TypeObject ={PyObject_HEAD_INIT(NULL)};
/*--------------------------------------------*/
void AtomsDMD::setup_tp()
{
    TypeObject.tp_name="atoms_dmd";
    TypeObject.tp_doc="I will add doc here";
    
    TypeObject.tp_flags=Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE;
    TypeObject.tp_basicsize=sizeof(Object);
    
    TypeObject.tp_new=__new__;
    TypeObject.tp_init=__init__;
    TypeObject.tp_alloc=__alloc__;
    TypeObject.tp_dealloc=__dealloc__;
    
    setup_tp_getset();
    TypeObject.tp_getset=getset;
    setup_tp_methods();
    TypeObject.tp_methods=methods;
    
    TypeObject.tp_base=&Atoms::TypeObject;
}
/*--------------------------------------------*/
PyGetSetDef AtomsDMD::getset[]={[0 ... 0]={NULL,NULL,NULL,NULL,NULL}};
/*--------------------------------------------*/
void AtomsDMD::setup_tp_getset()
{
    
}
/*--------------------------------------------*/
PyMethodDef AtomsDMD::methods[]={[0 ... 3]={NULL,NULL,0,NULL}};
/*--------------------------------------------*/
void AtomsDMD::setup_tp_methods()
{
    ForceFieldEAMDMD::ml_new(methods[0],methods[1],methods[2]);
    
}

