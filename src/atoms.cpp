#include "atoms.h"
#include "comm.h"
#include "elements.h"
#include "xmath.h"
using namespace MAPP_NS;
template<> MPI_Datatype Vec<char>::MPI_T=MPI_CHAR;
template<> MPI_Datatype Vec<short>::MPI_T=MPI_SHORT;
template<> MPI_Datatype Vec<int>::MPI_T=MPI_INT;
template<> MPI_Datatype Vec<long>::MPI_T=MPI_LONG;
template<> MPI_Datatype Vec<long long>::MPI_T=MPI_LONG_LONG;
template<> MPI_Datatype Vec<unsigned char>::MPI_T=MPI_UNSIGNED_CHAR;
template<> MPI_Datatype Vec<unsigned short>::MPI_T=MPI_UNSIGNED_SHORT;
template<> MPI_Datatype Vec<unsigned int>::MPI_T=MPI_UNSIGNED;
template<> MPI_Datatype Vec<unsigned long int>::MPI_T=MPI_UNSIGNED_LONG;
template<> MPI_Datatype Vec<unsigned long long>::MPI_T=MPI_UNSIGNED_LONG_LONG;
template<> MPI_Datatype Vec<float>::MPI_T=MPI_FLOAT;
template<> MPI_Datatype Vec<double>::MPI_T=MPI_DOUBLE;
template<> MPI_Datatype Vec<long double>::MPI_T=MPI_LONG_DOUBLE;
#ifdef MPI_CXX_BOOL
template<> MPI_Datatype Vec<bool>::MPI_T=MPI_CXX_BOOL;
#else
template<> MPI_Datatype Vec<bool>::MPI_T=MPI_BYTE;
#endif
/*---------------------------------------------------------------------------
      ___   _____   _____       ___  ___   _____
     /   | |_   _| /  _  \     /   |/   | /  ___/
    / /| |   | |   | | | |    / /|   /| | | |___
   / / | |   | |   | | | |   / / |__/ | | \___  \
  / /  | |   | |   | |_| |  / /       | |  ___| |
 /_/   |_|   |_|   \_____/ /_/        |_| /_____/
 
 ---------------------------------------------------------------------------*/
Atoms::Atoms(MPI_Comm& __world):
nvecs(0),
vecs(NULL),
tot_natms(0),
natms(0),
natms_ph(0),
comm(__world),
world(comm.world),
comm_rank(Communication::get_rank(__world)),
comm_size(Communication::get_size(__world)),
s_hi(comm.s_hi),
s_lo(comm.s_lo),
elements(new Elements()),
h(1.0),
kB(1.0),
vol(0.0),
H{[0 ... __dim__-1]={[0 ... __dim__-1]=0.0}},
B{[0 ... __dim__-1]={[0 ... __dim__-1]=0.0}},
depth_inv{[0 ... __dim__-1]=0.0},
dof(NULL),
step(0)
{
    x=new Vec<type0>(this,__dim__);
    id= new Vec<unsigned int>(this,1);
}
/*--------------------------------------------
 
 --------------------------------------------*/
Atoms::Atoms(Communication& __comm):
nvecs(0),
vecs(NULL),
tot_natms(0),
natms(0),
natms_ph(0),
comm(__comm),
world(comm.world),
comm_rank(__comm.rank),
comm_size(__comm.size),
s_hi(comm.s_hi),
s_lo(comm.s_lo),
elements(new Elements()),
h(1.0),
kB(1.0),
H{[0 ... __dim__-1]={[0 ... __dim__-1]=0.0}},
B{[0 ... __dim__-1]={[0 ... __dim__-1]=0.0}},
depth_inv{[0 ... __dim__-1]=0.0},
dof(NULL),
step(0)
{
    x=new Vec<type0>(this,__dim__);
    id=new Vec<unsigned int>(this,1);
}
/*--------------------------------------------
 destructor
 --------------------------------------------*/
Atoms::~Atoms()
{
    delete elements;
    while(nvecs) delete vecs[0];
}
/*--------------------------------------------
 use this function after changing H
 --------------------------------------------*/
void Atoms::update_H()
{
    Algebra::MLT_inv(H,B);
    vol=Algebra::MLT_det(H);
    Algebra::MLT_depth(B,depth_inv);
}
/*--------------------------------------------
 add a new vec
 --------------------------------------------*/
void Atoms::push(vec* v)
{
    if(v->name)
        for(int ivec=0;ivec<nvecs;ivec++)
            if(vecs[ivec]->name && !strcmp(v->name,vecs[ivec]->name))
                throw "internal error: vector name already exist";
    
    vec** __vecs=new vec*[nvecs+1];
    memcpy(__vecs,vecs,nvecs*sizeof(vec*));
    delete [] vecs;
    vecs=__vecs;
    vecs[nvecs++]=v;
}
/*--------------------------------------------
 remove a vector
 --------------------------------------------*/
void Atoms::pop(vec* v)
{
    int ivec=nvecs-1;
    for(;vecs[ivec]!=v;ivec--){}
    vecs[ivec]=vecs[nvecs-1];
    nvecs--;

    vec** __vecs=NULL;
    if(nvecs) __vecs=new vec*[nvecs];
    
    memcpy(__vecs,vecs,(nvecs)*sizeof(vec*));
}
/*--------------------------------------------
 add a new vec with name
 --------------------------------------------*/
vec* Atoms::find_vec(const char* name)
{
    for(int ivec=0;ivec<nvecs;ivec++)
        if(vecs[ivec]->name && strcmp(name,vecs[ivec]->name)==0)
            return vecs[ivec];
    return NULL;
}
/*--------------------------------------------
 x2s
 --------------------------------------------*/
void Atoms::x2s(int no)
{
    type0* x_vec=x->begin();
    int x_dim=x->dim;
    for(int i=0;i<no;i++)
        XMatrixVector::x2s<__dim__>(x_vec+i*x_dim,B);
}
/*--------------------------------------------
 s2x
 --------------------------------------------*/
void Atoms::s2x(int no)
{
    type0* x_vec=x->begin();
    int x_dim=x->dim;
    for(int i=0;i<no;i++)
        XMatrixVector::s2x<__dim__>(x_vec+i*x_dim,H);
}
/*--------------------------------------------
 x2s
 --------------------------------------------*/
void Atoms::x2s_lcl()
{
    type0* x_vec=x->begin();
    int x_dim=x->dim;
    for(int i=0;i<natms;i++)
        XMatrixVector::x2s<__dim__>(x_vec+i*x_dim,B);

}
/*--------------------------------------------
 s2x
 --------------------------------------------*/
void Atoms::s2x_lcl()
{
    type0* x_vec=x->begin();
    int x_dim=x->dim;
    for(int i=0;i<natms;i++)
        XMatrixVector::s2x<__dim__>(x_vec+i*x_dim,H);

}
/*--------------------------------------------
 x2s
 --------------------------------------------*/
void Atoms::x2s_all()
{
    type0* x_vec=x->begin();
    int x_dim=x->dim;
    int nall=natms+natms_ph;
    for(int i=0;i<nall;i++)
        XMatrixVector::x2s<__dim__>(x_vec+i*x_dim,B);
}
/*--------------------------------------------
 s2x
 --------------------------------------------*/
void Atoms::s2x_all()
{
    type0* x_vec=x->begin();
    int x_dim=x->dim;
    int nall=natms+natms_ph;
    for(int i=0;i<nall;i++)
        XMatrixVector::s2x<__dim__>(x_vec+i*x_dim,H);
}
/*--------------------------------------------
 insert a number of atoms
 --------------------------------------------*/
void Atoms::insert(byte* buff,vec** vecs_,int nvecs_,int natms_)
{
    int stride=0;
    for(int ivec=0;ivec<nvecs_;ivec++) stride+=vecs_[ivec]->byte_sz;
    

    for(int ivec=0;ivec<nvecs_;ivec++)
        vecs_[ivec]->pst(buff,stride,natms_);
    
    natms+=natms_;
}
/*--------------------------------------------
 make room for some local atoms and phantom 
 atoms; this is used for grand canocical monte 
 carlo, when a successfull insertion trial has
 occured. using this function we make room for 
 the new entries. The new entries are insrerted 
 manually by GCMC
 
 *** we might need a better name for this 
 function
 --------------------------------------------*/
void Atoms::add()
{
    for(int i=0;i<nvecs;i++)
        vecs[i]->add();
    natms++;
}
/*--------------------------------------------
 delete some local atoms and phantom atoms; 
 this is used for grand canocical monte carlo, 
 when a successfull deletion trial has occured.
 it takes a the list of local atoms and phantoms
 
 !! it is assumed that both lists are ascending
 
 *** we might need a better name for this
 function
 --------------------------------------------*/
void Atoms::del(int& del_idx)
{
    for(int i=0;i<nvecs;i++)
        vecs[i]->del(del_idx);
    natms--;
}
/*--------------------------------------------
 restart
 --------------------------------------------*/
void Atoms::restart()
{
    tot_natms=0;
    natms=0;
    natms_ph=0;
    while(nvecs)
        delete vecs[0];
}
/*------------------------------------------------------------------------------------------------------------------------------------
 
 ------------------------------------------------------------------------------------------------------------------------------------*/
#include "ff_styles.h"
PyObject* Atoms::__new__(PyTypeObject* type,PyObject* args,PyObject* kwds)
{
    Object* __self=reinterpret_cast<Object*>(type->tp_alloc(type,0));
    PyObject* self=reinterpret_cast<PyObject*>(__self);
    return self;
}
/*--------------------------------------------
 
 --------------------------------------------*/
int Atoms::__init__(PyObject* self,PyObject* args,PyObject* kwds)
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
PyObject* Atoms::__alloc__(PyTypeObject* type,Py_ssize_t)
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
void Atoms::__dealloc__(PyObject* self)
{
    Object* __self=reinterpret_cast<Object*>(self);
    delete __self->ff;
    delete __self->atoms;
    delete __self;
}
/*--------------------------------------------*/
PyMethodDef Atoms::methods[]={[0 ... 1]={NULL}};
/*--------------------------------------------*/
void Atoms::setup_tp_methods()
{
    ml_strain(methods[0]);
}
/*--------------------------------------------*/
PyTypeObject Atoms::TypeObject ={PyObject_HEAD_INIT(NULL)};
/*--------------------------------------------*/
void Atoms::setup_tp()
{
    TypeObject.tp_name="atoms";
    TypeObject.tp_doc="I will add doc here";
    
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
PyGetSetDef Atoms::getset[]={[0 ... 12]={NULL,NULL,NULL,NULL,NULL}};
/*--------------------------------------------*/
void Atoms::setup_tp_getset()
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
/*--------------------------------------------
 
 --------------------------------------------*/
void Atoms::getset_step(PyGetSetDef& getset)
{
    getset.name=(char*)"step";
    getset.doc=(char*)"step number";
    getset.get=[](PyObject* self,void*)->PyObject*
    {
        return var<int>::build(reinterpret_cast<Object*>(self)->atoms->step,NULL);
    };
    getset.set=[](PyObject* self,PyObject* val,void*)->int
    {
        VarAPI<int> step("step");
        step.logics[0]=VLogics("ge",0);
        if(step.set(val)==-1) return -1;
        reinterpret_cast<Object*>(self)->atoms->step=step.val;
        return 0;
    };
}
/*--------------------------------------------
 
 --------------------------------------------*/
void Atoms::getset_h(PyGetSetDef& getset)
{
    getset.name=(char*)"h";
    getset.doc=(char*)"Planck constant";
    getset.get=[](PyObject* self,void*)->PyObject*
    {
        return var<type0>::build(reinterpret_cast<Object*>(self)->atoms->h,NULL);
    };
    getset.set=[](PyObject* self,PyObject* val,void*)->int
    {
        Object* __self=reinterpret_cast<Object*>(self);
        if(!__self->atoms)
        {
            PyErr_Format(PyExc_TypeError,"cannot set 'h' prior to loading system configuration");
            return -1;
        }
        
        VarAPI<type0> h("h");
        h.logics[0]=VLogics("gt",0.0);
        if(h.set(val)==-1) return -1;
        reinterpret_cast<Object*>(self)->atoms->h=h.val;
        return 0;
    };
}
/*--------------------------------------------
 
 --------------------------------------------*/
void Atoms::getset_kB(PyGetSetDef& getset)
{
    getset.name=(char*)"kB";
    getset.doc=(char*)"Boltzmann constant";
    getset.get=[](PyObject* self,void*)->PyObject*
    {
        return var<type0>::build(reinterpret_cast<Object*>(self)->atoms->h,NULL);
    };
    getset.set=[](PyObject* self,PyObject* val,void*)->int
    {
        Object* __self=reinterpret_cast<Object*>(self);
        if(!__self->atoms)
        {
            PyErr_Format(PyExc_TypeError,"cannot set 'kB' prior to loading system configuration");
            return -1;
        }
        
        VarAPI<type0> kB("kB");
        kB.logics[0]=VLogics("gt",0.0);
        if(kB.set(val)==-1) return -1;
        reinterpret_cast<Object*>(self)->atoms->kB=kB.val;
        return 0;
    };
}
/*--------------------------------------------
 
 --------------------------------------------*/
void Atoms::getset_skin(PyGetSetDef& getset)
{
    getset.name=(char*)"skin";
    getset.doc=(char*)"skin size";
    getset.get=[](PyObject* self,void*)->PyObject*
    {
        return var<type0>::build(reinterpret_cast<Object*>(self)->atoms->comm.skin,NULL);
    };
    getset.set=[](PyObject* self,PyObject* val,void*)->int
    {
        Object* __self=reinterpret_cast<Object*>(self);
        if(!__self->atoms)
        {
            PyErr_Format(PyExc_TypeError,"cannot set 'skin' prior to loading system configuration");
            return -1;
        }
        
        VarAPI<type0> skin("skin");
        skin.logics[0]=VLogics("gt",0.0);
        if(skin.set(val)==-1) return -1;
        reinterpret_cast<Object*>(self)->atoms->comm.skin=skin.val;
        return 0;
    };
}
/*--------------------------------------------
 
 --------------------------------------------*/
void Atoms::getset_H(PyGetSetDef& getset)
{
    getset.name=(char*)"H";
    getset.doc=(char*)"matrix of unitcell vectors";
    getset.get=[](PyObject* self,void*)->PyObject*
    {
        Object* __self=reinterpret_cast<Object*>(self);
        if(!__self->atoms)
        {
            PyErr_Format(PyExc_TypeError,"cannot return 'H' prior to loading system configuration");
            return NULL;
        }

        return var<type0[__dim__][__dim__]>::build(__self->atoms->H,NULL);
    };
    getset.set=[](PyObject*,PyObject*,void*)->int
    {
        PyErr_SetString(PyExc_TypeError,"readonly attribute");
        return -1;
    };
}
/*--------------------------------------------
 
 --------------------------------------------*/
void Atoms::getset_B(PyGetSetDef& getset)
{
    getset.name=(char*)"B";
    getset.doc=(char*)"matrix of reciprocal unitcell vectors";
    getset.get=[](PyObject* self,void*)->PyObject*
    {
        Object* __self=reinterpret_cast<Object*>(self);
        if(!__self->atoms)
        {
            PyErr_Format(PyExc_TypeError,"cannot return 'B' prior to loading system configuration");
            return NULL;
        }

        return var<type0[__dim__][__dim__]>::build(__self->atoms->B,NULL);
    };
    getset.set=[](PyObject*,PyObject*,void*)->int
    {
        PyErr_SetString(PyExc_TypeError,"readonly attribute");
        return -1;
    };
}
/*--------------------------------------------
 
 --------------------------------------------*/
void Atoms::getset_vol(PyGetSetDef& getset)
{
    getset.name=(char*)"vol";
    getset.doc=(char*)"volume of unitcell";
    getset.get=[](PyObject* self,void*)->PyObject*
    {
        Object* __self=reinterpret_cast<Object*>(self);
        if(!__self->atoms)
        {
            PyErr_Format(PyExc_TypeError,"cannot return 'vol' prior to loading system configuration");
            return NULL;
        }

        return var<type0>::build(__self->atoms->vol,NULL);
    };
    getset.set=[](PyObject*,PyObject*,void*)->int
    {
        PyErr_SetString(PyExc_TypeError,"readonly attribute");
        return -1;
    };
}
/*--------------------------------------------
 
 --------------------------------------------*/
void Atoms::getset_elems(PyGetSetDef& getset)
{
    getset.name=(char*)"elems";
    getset.doc=(char*)"dictionary of elements present in the system";
    getset.get=[](PyObject* self,void*)->PyObject*
    {
        Object* __self=reinterpret_cast<Object*>(self);
        if(!__self->atoms)
        {
            PyErr_Format(PyExc_TypeError,"cannot return 'elems' prior to system configuration");
            return NULL;
        }
        return __self->atoms->elements->get_dict();
    };
    getset.set=[](PyObject*,PyObject*,void*)->int
    {
        PyErr_SetString(PyExc_TypeError,"readonly attribute");
        return -1;
    };
}
/*--------------------------------------------
 
 --------------------------------------------*/
void Atoms::getset_comm_rank(PyGetSetDef& getset)
{
    getset.name=(char*)"comm_rank";
    getset.doc=(char*)"rank of the calling process in the communicator";
    getset.get=[](PyObject* self,void*)->PyObject*
    {
        Object* __self=reinterpret_cast<Object*>(self);
        if(!__self->atoms)
        {
            PyErr_Format(PyExc_TypeError,"cannot return 'comm_rank' prior to system configuration");
            return NULL;
        }
        int rank=__self->atoms->comm.rank;
        return var<int>::build(rank,NULL);
    };
    getset.set=[](PyObject*,PyObject*,void*)->int
    {
        PyErr_SetString(PyExc_TypeError,"readonly attribute");
        return -1;
    };
}
/*--------------------------------------------
 
 --------------------------------------------*/
void Atoms::getset_comm_size(PyGetSetDef& getset)
{
    getset.name=(char*)"comm_size";
    getset.doc=(char*)"size of the group associated with the communicator";
    getset.get=[](PyObject* self,void*)->PyObject*
    {
        Object* __self=reinterpret_cast<Object*>(self);
        if(!__self->atoms)
        {
            PyErr_Format(PyExc_TypeError,"cannot return 'comm_size' prior to system configuration");
            return NULL;
        }
        int size=__self->atoms->comm.size;
        return var<int>::build(size,NULL);
    };
    getset.set=[](PyObject*,PyObject*,void*)->int
    {
        PyErr_SetString(PyExc_TypeError,"readonly attribute");
        return -1;
    };
}
/*--------------------------------------------
 
 --------------------------------------------*/
void Atoms::getset_comm_coords(PyGetSetDef& getset)
{
    getset.name=(char*)"comm_coords";
    getset.doc=(char*)"coordinates of calling process in cartesian structure";
    getset.get=[](PyObject* self,void*)->PyObject*
    {
        Object* __self=reinterpret_cast<Object*>(self);
        if(!__self->atoms)
        {
            PyErr_Format(PyExc_TypeError,"cannot return 'comm_coords' prior to system configuration");
            return NULL;
        }
        return var<int[__dim__]>::build(__self->atoms->comm.coords,NULL);
    };
    getset.set=[](PyObject*,PyObject*,void*)->int
    {
        PyErr_SetString(PyExc_TypeError,"readonly attribute");
        return -1;
    };
}
/*--------------------------------------------
 
 --------------------------------------------*/
void Atoms::getset_comm_dims(PyGetSetDef& getset)
{
    getset.name=(char*)"comm_dims";
    getset.doc=(char*)"number of processes for each cartesian dimension";
    getset.get=[](PyObject* self,void*)->PyObject*
    {
        Object* __self=reinterpret_cast<Object*>(self);
        if(!__self->atoms)
        {
            PyErr_Format(PyExc_TypeError,"cannot return 'comm_dims' prior to system configuration");
            return NULL;
        }
        return var<int[__dim__]>::build(__self->atoms->comm.dims,NULL);
    };
    getset.set=[](PyObject* self,PyObject* val,void*)->int
    {
        VarAPI<int[__dim__]> dims("dims");
        dims.logics[0]=VLogics("ge",0);
        int iset=dims.set(val);
        if(iset==-1) return -1;
        int size=1;
        int undet_dims=0;
        for(int i=0;i<__dim__;i++)
        {
            if(dims.val[i]) size*=dims.val[i];
            else undet_dims++;
        }
        Object* __self=reinterpret_cast<Object*>(self);
        if(__self->atoms->comm.size%size!=0)
        {
            PyErr_Format(PyExc_TypeError,"'comm_size' is not divisble by the product of nonzero components of 'comm_dims'");
            return -1;
        }
        
        if(undet_dims==0 && __self->atoms->comm.size!=size)
        {
            PyErr_Format(PyExc_TypeError,"'comm_size' is not equal to the product of nonzero components of 'comm_dims'");
            return -1;
        }
        
        __self->atoms->comm.grid(dims.val);
        return 0;
    };
}
/*--------------------------------------------
 python constructor
 --------------------------------------------*/
void Atoms::ml_strain(PyMethodDef& tp_methods)
{
    tp_methods.ml_flags=METH_VARARGS | METH_KEYWORDS;
    tp_methods.ml_name="strain";
    tp_methods.ml_doc="strain the box";
    
    tp_methods.ml_meth=(PyCFunction)(PyCFunctionWithKeywords)
    [](PyObject* self,PyObject* args,PyObject* kwds)->PyObject*
    {
        FuncAPI<type0[__dim__][__dim__]> f("strain",{"E"});
        if(f(args,kwds)) return NULL;
        
        Atoms::Object* __self=reinterpret_cast<Atoms::Object*>(self);
        
        type0 H[__dim__][__dim__]{[0 ... __dim__-1]={[0 ... __dim__-1]=0.0}};
        type0 B[__dim__][__dim__]{[0 ... __dim__-1]={[0 ... __dim__-1]=0.0}};
        Algebra::DoLT<__dim__>::func([&H,&B,&__self](int i,int j)
        {
            H[i][j]=__self->atoms->H[i][j];
            B[i][j]=__self->atoms->B[i][j];
        });
        
        type0 (&strain)[__dim__][__dim__]=f.val<0>();
        
        Algebra::Do<__dim__>::func([&strain](int i){strain[i][i]++;});
        type0 F[__dim__][__dim__]{[0 ... __dim__-1]={[0 ... __dim__-1]=0.0}};
        type0 __H[__dim__][__dim__]{[0 ... __dim__-1]={[0 ... __dim__-1]=0.0}};
        Algebra::MLT_mul_MSQ(H,strain,__H);
        
        Algebra::MSQ_2_MLT(__H,H);
        Algebra::MLT_mul_MLT(B,H,F);
        if(Algebra::MLT_det(F)<=0.0)
        PyErr_SetString(PyExc_TypeError,"strain should result in transformation tensor with positive determinant");
        Algebra::MLT_inv(H,B);
        
        Algebra::V_eq<__dim__*__dim__>(&H[0][0],&(__self->atoms->H[0][0]));
        __self->atoms->update_H();
        
        
        type0* x=__self->atoms->x->begin();
        int x_dim=__self->atoms->x->dim;
        int natms=__self->atoms->natms;
        for(int i=0;i<natms;i++,x+=x_dim)
            Algebra::V_mul_MLT(x,F,x);
        
        Py_RETURN_NONE;
    };
}
