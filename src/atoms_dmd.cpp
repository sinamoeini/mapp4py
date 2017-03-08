#include "atoms_dmd.h"
#include "xmath.h"
/*--------------------------------------------
 
 --------------------------------------------*/
AtomsDMD::AtomsDMD(MPI_Comm& world,int __c_dim,int __N):
c_dim(__c_dim),
N(__N),
Atoms(world),
xi(new type0[__N]),
wi(new type0[__N])
{
    XMath::quadrature_hg(N,xi,wi);
    alpha=new Vec<type0>(this,c_dim);
    c=new Vec<type0>(this,c_dim);
    elem=new Vec<elem_type>(this,c_dim);
    alpha_dof=new Vec<bool>(this,c_dim);
    alpha_dof->empty(true);
    c_dof=new Vec<bool>(this,c_dim);
    c_dof->empty(true);
    c_d=new Vec<type0>(this,c_dim);
    c_d->empty(0.0);
}
/*--------------------------------------------
 
 --------------------------------------------*/
AtomsDMD::~AtomsDMD()
{
    delete [] wi;
    delete [] xi;
    
    delete alpha_dof;
    delete c_dof;
    delete c_d;
    delete elem;
    delete c;
    delete alpha;
    
}
/*--------------------------------------------
 
 --------------------------------------------*/
#include "elements.h"
AtomsDMD& AtomsDMD::operator=(const Atoms& r)
{
    elements=r.elements;
    comm=r.comm;
    natms_lcl=r.natms_lcl;
    natms_ph=r.natms_ph;
    natms=r.natms;
    step=r.step;
    
    max_cut=r.max_cut;
    kB=r.kB;
    h=r.h;
    
    vol=r.vol;
    for(int i=0;i<__dim__;i++)
    {
        depth_inv[i]=r.depth_inv[i];
        for(int j=0;j<__dim__;j++)
        {
            H[i][j]=r.H[i][j];
            B[i][j]=r.B[i][j];
        }
    }
    
    for(int i=0;i<nvecs;i++)
        if(!vecs[i]->is_empty())
            vecs[i]->resize(natms_lcl);
    
    
    memcpy(x->begin(),r.x->begin(),natms_lcl*__dim__*sizeof(type0));
    memcpy(id->begin(),r.id->begin(),natms_lcl*sizeof(unsigned int));
    return* this;
}
/*--------------------------------------------
 
 --------------------------------------------*/
#include "memory.h"
void AtomsDMD::sort_stack(
vec**& xchng,int& nxchng,
vec**& updt,int& nupdt,
vec**& arch,int& narch,
vec**& empty,int& nempty)
{
    
    /*
     go through exchange, update, and archive if they
     are any empty vectors segregate them into empty
     vectors category
     */
    nempty=0;
    delete [] empty;
    int empty_cpcty=nxchng+nupdt+narch;
    Memory::alloc(empty,empty_cpcty);
    
    auto fill_empty_vecs=[&nempty,empty](vec**& vs,int& nvs)->void
    {
        int __nvs=nvs;
        for(int i=0;i<__nvs;)
            if(vs[i]->is_empty())
            {
                empty[nempty++]=vs[i];
                
                
                vs[i]=vs[__nvs-1];
                __nvs--;
                
            }
            else
                i++;
        
        Memory::shrink_to_fit(vs,__nvs,nvs);
        nvs=__nvs;
    };
    
    fill_empty_vecs(xchng,nxchng);
    fill_empty_vecs(updt,nupdt);
    fill_empty_vecs(arch,narch);

    /*
     pop all the empty vectors from vector stack.
     pop remaining archive vectors from vector stack.
     */
    for(int i=0;i<nempty;i++)
        pop(empty[i]);
    for(int i=0;i<narch;i++)
        pop(arch[i]);
    
    /*
     now sort the vector stack as follows
     0. all exchange vectors
     1. all update vectors
     2. remaining vectors (passive)
     ** the VERY FIRST exchange vector: id
     ** the VERY FIRST update vector: x
     ** the SECOND update vector: alpha
     ** the THIRD update vector: c
     ** the FOURTH update vector elem
     */
    vec** __vecs=NULL;
    if(nvecs) __vecs=new vec*[nvecs];
    
    int __nvecs=0;
    __vecs[__nvecs++]=id;
    
    for(int i=0;i<nxchng;i++)
        __vecs[__nvecs++]=xchng[i];
    nxchng++;
    
    __vecs[__nvecs++]=x;
    __vecs[__nvecs++]=alpha;
    __vecs[__nvecs++]=c;
    __vecs[__nvecs++]=elem;
    for(int i=0;i<nupdt;i++)
        __vecs[__nvecs++]=updt[i];
    
    nupdt+=4;
    
    int ___nevcs=__nvecs;
    
    auto is_in=[](const vec* v,vec* const * vs,int nvs)->bool
    {
        for(int i=0;i<nvs;i++)
            if(vs[i]==v) return true;
        return false;
    };
    
    for(int i=0;i<nvecs;i++)
        if(!is_in(vecs[i],__vecs,___nevcs))
            __vecs[__nvecs++]=vecs[i];
}
/*------------------------------------------------------------------------------------------------------------------------------------
 
 ------------------------------------------------------------------------------------------------------------------------------------*/
#include "ff_styles.h"
#include "read_styles.h"
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
    TypeObject.tp_name="mapp.dmd.atoms";
    
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
    
    TypeObject.tp_doc=R"---(
    atoms()
    
    Container class for system configuration and force field.

    
    )---";
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
PyMethodDef AtomsDMD::methods[]={[0 ... 5]={NULL,NULL,0,NULL}};
/*--------------------------------------------*/
void AtomsDMD::setup_tp_methods()
{
    ml_strain(methods[0]);
    ForceFieldEAMDMD::ml_new(methods[1],methods[2],methods[3]);
    ReadCFGDMD::ml_import(methods[4]);
}

