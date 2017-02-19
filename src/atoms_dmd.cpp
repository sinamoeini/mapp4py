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
    TypeObject.tp_name="atoms";
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
PyGetSetDef AtomsDMD::getset[]={[0 ... 12]={NULL,NULL,NULL,NULL,NULL}};
/*--------------------------------------------*/
void AtomsDMD::setup_tp_getset()
{
    getset_step(getset[0]);
    getset_h(getset[1]);
    getset_kB(getset[2]);
    getset_H(getset[3]);
    getset_B(getset[4]);
    getset_vol(getset[5]);
    getset_elems(getset[6]);
    getset_skin(getset[7]);
    getset_comm_rank(getset[8]);
    getset_comm_size(getset[9]);
    getset_comm_coords(getset[10]);
    getset_comm_dims(getset[11]);
}
/*--------------------------------------------*/
PyMethodDef AtomsDMD::methods[]={[0 ... 4]={NULL,NULL,0,NULL}};
/*--------------------------------------------*/
void AtomsDMD::setup_tp_methods()
{
    ml_strain(methods[0]);
    ForceFieldEAMDMD::ml_new(methods[1],methods[2],methods[3]);
    
}

