#include "atoms.h"
#include "comm.h"
#include "elements.h"
#include "xmath.h"
#include <limits>
using namespace MAPP_NS;
#ifdef MPI_CXX_BOOL
template<> MPI_Datatype Vec<bool>::MPI_T=MPI_CXX_BOOL;
#else
template<> MPI_Datatype Vec<bool>::MPI_T=MPI_BYTE;
#endif
template<> MPI_Datatype Vec<char>::MPI_T=MPI_CHAR;
template<> MPI_Datatype Vec<short>::MPI_T=MPI_SHORT;
template<> MPI_Datatype Vec<int>::MPI_T=MPI_INT;
template<> MPI_Datatype Vec<long int>::MPI_T=MPI_LONG;
template<> MPI_Datatype Vec<long long>::MPI_T=MPI_LONG_LONG;
template<> MPI_Datatype Vec<unsigned char>::MPI_T=MPI_UNSIGNED_CHAR;
template<> MPI_Datatype Vec<unsigned short>::MPI_T=MPI_UNSIGNED_SHORT;
template<> MPI_Datatype Vec<unsigned int>::MPI_T=MPI_UNSIGNED;
template<> MPI_Datatype Vec<unsigned long int>::MPI_T=MPI_UNSIGNED_LONG;
template<> MPI_Datatype Vec<unsigned long long>::MPI_T=MPI_UNSIGNED_LONG_LONG;
template<> MPI_Datatype Vec<float>::MPI_T=MPI_FLOAT;
template<> MPI_Datatype Vec<double>::MPI_T=MPI_DOUBLE;
template<> MPI_Datatype Vec<long double>::MPI_T=MPI_LONG_DOUBLE;


template<> const char* Vec<bool>::print_format="%d ";
template<> const char* Vec<char>::print_format="%c ";
template<> const char* Vec<short>::print_format="%d ";
template<> const char* Vec<int>::print_format="%d ";
template<> const char* Vec<long int>::print_format="%ld ";
template<> const char* Vec<long long>::print_format="%lld ";
template<> const char* Vec<unsigned char>::print_format="%c ";
template<> const char* Vec<unsigned short>::print_format="%d ";
template<> const char* Vec<unsigned int>::print_format="%ud ";
template<> const char* Vec<unsigned long int>::print_format="%uld ";
template<> const char* Vec<unsigned long long>::print_format="%ulld ";
template<> const char* Vec<float>::print_format="%e ";
template<> const char* Vec<double>::print_format="%.16lf ";
template<> const char* Vec<long double>::print_format="%Le ";

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
dynamic_vecs(NULL),
ndynamic_vecs(0),
natms(0),
natms_lcl(0),
natms_ph(0),
comm(__world),
world(comm.world),
comm_rank(Communication::get_rank(__world)),
comm_size(Communication::get_size(__world)),
s_hi(comm.s_hi),
s_lo(comm.s_lo),
elements(Elements()),
hP(NAN),
kB(NAN),
vol(0.0),
H{DESIG2(__dim__,__dim__,0.0)},
B{DESIG2(__dim__,__dim__,0.0)},
depth_inv{DESIG(__dim__,0.0)},
step(0)
{
    Algebra::zero<__nvoigt__>(__h);
    Algebra::zero<__nvoigt__>(__b);
    x=new Vec<type0>(this,__dim__,"x");
    id= new Vec<unsigned int>(this,1,"id");
    dof=new Vec<bool>(this,__dim__,"dof");
    dof->empty(true);
}
/*--------------------------------------------
 destructor
 --------------------------------------------*/
Atoms::~Atoms()
{
    delete dof;
    delete id;
    delete x;
    delete [] dynamic_vecs;
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
    Algebra::MLT_2_V(H,__h);
    Algebra::MLT_2_V(B,__b);
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
    
    if(ndynamic_vecs)
    {
        vec** __dyamic_vecs=new vec*[ndynamic_vecs+1];
        memcpy(__dyamic_vecs,dynamic_vecs,ndynamic_vecs*sizeof(vec*));
        delete [] dynamic_vecs;
        dynamic_vecs=__dyamic_vecs;
        dynamic_vecs[ndynamic_vecs++]=v;
    }
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
    
    if(ndynamic_vecs)
    {
        ivec=ndynamic_vecs-1;
        for(;dynamic_vecs[ivec]!=v;ivec--){}
        dynamic_vecs[ivec]=dynamic_vecs[ndynamic_vecs-1];
        ndynamic_vecs--;
        vec** __dyamic_vecs=NULL;
        if(ndynamic_vecs)
            __dyamic_vecs=new vec*[ndynamic_vecs];
        
        
        
        
        memcpy(__dyamic_vecs,dynamic_vecs,ndynamic_vecs*sizeof(vec*));
        delete [] dynamic_vecs;
        dynamic_vecs=__dyamic_vecs;
    }
}
/*--------------------------------------------
 remove a vector
 --------------------------------------------*/
int Atoms::is_dynamic(vec* v)
{
    for(int i=ndynamic_vecs-1;i>-1;i--)
        if(v==dynamic_vecs[i]) return i;
    return -1;
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
    Algebra::X2S<__dim__>(__b,no,x->begin());
}
/*--------------------------------------------
 s2x
 --------------------------------------------*/
void Atoms::s2x(int no)
{
    Algebra::S2X<__dim__>(__h,no,x->begin());
}
/*--------------------------------------------
 x2s
 --------------------------------------------*/
void Atoms::x2s_lcl()
{
    Algebra::X2S<__dim__>(__b,natms_lcl,x->begin());
}
/*--------------------------------------------
 s2x
 --------------------------------------------*/
void Atoms::s2x_lcl()
{
    Algebra::S2X<__dim__>(__h,natms_lcl,x->begin());
}
/*--------------------------------------------
 x2s
 --------------------------------------------*/
void Atoms::x2s_all()
{
    Algebra::X2S<__dim__>(__b,natms_lcl+natms_ph,x->begin());
}
/*--------------------------------------------
 s2x
 --------------------------------------------*/
void Atoms::s2x_all()
{
    Algebra::S2X<__dim__>(__h,natms_lcl+natms_ph,x->begin());
}
/*--------------------------------------------
 x2s
 --------------------------------------------*/
void Atoms::x2s_dump()
{
    Algebra::X2S<__dim__>(__b,natms,x->begin_dump());
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
    
    natms_lcl+=natms_;
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
    natms_lcl++;
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
    natms_lcl--;
}
/*--------------------------------------------
 restart
 --------------------------------------------*/
void Atoms::restart()
{
    natms=0;
    natms_lcl=0;
    natms_ph=0;
    while(nvecs)
        delete vecs[0];
}
/*--------------------------------------------
 reset the domain
 --------------------------------------------*/
#include "exchange.h"
void Atoms::reset_domain()
{
    ndynamic_vecs=0;
    for(int i=0;i<nvecs;i++)
        if(!vecs[i]->is_empty())
            ndynamic_vecs++;
    if(ndynamic_vecs==0) return;
    
    dynamic_vecs=new vec*[ndynamic_vecs];
    ndynamic_vecs=0;
    for(int i=0;i<nvecs;i++)
        if(!vecs[i]->is_empty())
            dynamic_vecs[ndynamic_vecs++]=vecs[i];
    
    
    Exchange xchng(this,ndynamic_vecs);
    x2s_lcl();
    xchng.full_xchng();
    s2x_lcl();
    
    
    delete [] dynamic_vecs;
    dynamic_vecs=NULL;
    ndynamic_vecs=0;
    
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
    Py_TYPE(__self)=type;
    Py_REFCNT(__self)=1;
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
PyTypeObject Atoms::TypeObject ={PyObject_HEAD_INIT(NULL)};
/*--------------------------------------------*/
int Atoms::setup_tp()
{
    TypeObject.tp_name="atoms";
    TypeObject.tp_doc="container class";
    
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
    
    int ichk=PyType_Ready(&TypeObject);
    if(ichk<0) return ichk;
    Py_INCREF(&TypeObject);
    
    return ichk;
}
/*--------------------------------------------*/
PyGetSetDef Atoms::getset[]=EmptyPyGetSetDef(13);
/*--------------------------------------------*/
void Atoms::setup_tp_getset()
{
    getset_step(getset[0]);
    getset_hP(getset[1]);
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
PyMethodDef Atoms::methods[]=EmptyPyMethodDef(2);
/*--------------------------------------------*/
void Atoms::setup_tp_methods()
{
    ml_strain(methods[0]);
}
/*--------------------------------------------
 
 --------------------------------------------*/
void Atoms::getset_step(PyGetSetDef& getset)
{
    getset.name=(char*)"step";
    getset.doc=(char*)R"---(
    (int) step number
    
    Current step number
    )---";
    getset.get=[](PyObject* self,void*)->PyObject*
    {
        return var<int>::build(reinterpret_cast<Object*>(self)->atoms->step);
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
void Atoms::getset_hP(PyGetSetDef& getset)
{
    getset.name=(char*)"hP";
    getset.doc=(char*)R"---(
    (double) Planck constant
    
    Planck constant
    )---";
    getset.get=[](PyObject* self,void*)->PyObject*
    {
        return var<type0>::build(reinterpret_cast<Object*>(self)->atoms->hP);
    };
    getset.set=[](PyObject* self,PyObject* val,void*)->int
    {
        Object* __self=reinterpret_cast<Object*>(self);
        if(!__self->atoms)
        {
            PyErr_Format(PyExc_TypeError,"cannot set 'hP' prior to loading system configuration");
            return -1;
        }
        
        VarAPI<type0> hP("hP");
        hP.logics[0]=VLogics("gt",0.0);
        if(hP.set(val)==-1) return -1;
        reinterpret_cast<Object*>(self)->atoms->hP=hP.val;
        return 0;
    };
}
/*--------------------------------------------
 
 --------------------------------------------*/
void Atoms::getset_kB(PyGetSetDef& getset)
{
    getset.name=(char*)"kB";
    getset.doc=(char*)R"---(
    (double) Boltzmann constant
    
    Boltzmann constant
    )---";
    getset.get=[](PyObject* self,void*)->PyObject*
    {
        return var<type0>::build(reinterpret_cast<Object*>(self)->atoms->kB);
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
    getset.doc=(char*)R"---(
    (double) size of skin
    
    Thickness of neighbor list skin
    )---";
    getset.get=[](PyObject* self,void*)->PyObject*
    {
        return var<type0>::build(reinterpret_cast<Object*>(self)->atoms->comm.skin);
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
    getset.doc=(char*)R"---(
    (symm<double[dim][dim]>) unitcell matrix
    
    Matrix of unitcell vectors, here dim is the dimension of simulation
    )---";
    getset.get=[](PyObject* self,void*)->PyObject*
    {
        Object* __self=reinterpret_cast<Object*>(self);
        if(!__self->atoms)
        {
            PyErr_Format(PyExc_TypeError,"cannot return 'H' prior to loading system configuration");
            return NULL;
        }

        return var<type0[__dim__][__dim__]>::build(__self->atoms->H);
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
    getset.doc=(char*)R"---(
    (symm<double[dim][dim]>) unitcell matrix inverse (:math:`\mathbf{H}^{-1}`)
    
    Inverse of matrix of unitcell vectors, here dim is the dimension of simulation
    )---";
    getset.get=[](PyObject* self,void*)->PyObject*
    {
        Object* __self=reinterpret_cast<Object*>(self);
        if(!__self->atoms)
        {
            PyErr_Format(PyExc_TypeError,"cannot return 'B' prior to loading system configuration");
            return NULL;
        }

        return var<type0[__dim__][__dim__]>::build(__self->atoms->B);
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
    getset.doc=(char*)R"---(
    (double) unitcell volume
    
    Volume of unitcell
    )---";
    getset.get=[](PyObject* self,void*)->PyObject*
    {
        Object* __self=reinterpret_cast<Object*>(self);
        if(!__self->atoms)
        {
            PyErr_Format(PyExc_TypeError,"cannot return 'vol' prior to loading system configuration");
            return NULL;
        }

        return var<type0>::build(__self->atoms->vol);
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
    getset.doc=(char*)R"---(
    (dict{int:string}) elements dictionary
    
    Dictionary of elements present in the system
    )---";
    getset.get=[](PyObject* self,void*)->PyObject*
    {
        Object* __self=reinterpret_cast<Object*>(self);
        if(!__self->atoms)
        {
            PyErr_Format(PyExc_TypeError,"cannot return 'elems' prior to system configuration");
            return NULL;
        }
        return __self->atoms->elements.get_dict();
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
    getset.doc=(char*)R"---(
    (int) communication rank
    
    Rank of the calling process in the communicator
    )---";
    getset.get=[](PyObject* self,void*)->PyObject*
    {
        Object* __self=reinterpret_cast<Object*>(self);
        if(!__self->atoms)
        {
            PyErr_Format(PyExc_TypeError,"cannot return 'comm_rank' prior to system configuration");
            return NULL;
        }
        int rank=__self->atoms->comm.rank;
        return var<int>::build(rank);
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
    getset.doc=(char*)R"---(
    (int) communication size
    
    Size of the group associated with the communicator
    )---";
    getset.get=[](PyObject* self,void*)->PyObject*
    {
        Object* __self=reinterpret_cast<Object*>(self);
        if(!__self->atoms)
        {
            PyErr_Format(PyExc_TypeError,"cannot return 'comm_size' prior to system configuration");
            return NULL;
        }
        int size=__self->atoms->comm.size;
        return var<int>::build(size);
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
    getset.doc=(char*)R"---(
    (int[dim]) communication coordinates
    
    Communication coordinates of the calling process, here dim is the dimension of simulation
    )---";
    getset.get=[](PyObject* self,void*)->PyObject*
    {
        Object* __self=reinterpret_cast<Object*>(self);
        if(!__self->atoms)
        {
            PyErr_Format(PyExc_TypeError,"cannot return 'comm_coords' prior to system configuration");
            return NULL;
        }
        return var<int[__dim__]>::build(__self->atoms->comm.coords);
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
    getset.doc=(char*)R"---(
    (int[dim]) communication dimensions
    
    Number of processes for each cartesian dimension, here dim is the dimension of simulation
    )---";
    getset.get=[](PyObject* self,void*)->PyObject*
    {
        Object* __self=reinterpret_cast<Object*>(self);
        if(!__self->atoms)
        {
            PyErr_Format(PyExc_TypeError,"cannot return 'comm_dims' prior to system configuration");
            return NULL;
        }
        return var<int[__dim__]>::build(__self->atoms->comm.dims);
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
        __self->atoms->reset_domain();
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
    tp_methods.ml_meth=(PyCFunction)(PyCFunctionWithKeywords)(
    [](PyObject* self,PyObject* args,PyObject* kwds)->PyObject*
    {
        FuncAPI<type0[__dim__][__dim__]> f("strain",{"E"});
        if(f(args,kwds)) return NULL;
        
        Atoms::Object* __self=reinterpret_cast<Atoms::Object*>(self);
        
        type0 H[__dim__][__dim__]{DESIG2(__dim__,__dim__,0.0)};
        type0 B[__dim__][__dim__]{DESIG2(__dim__,__dim__,0.0)};
        Algebra::DoLT<__dim__>::func([&H,&B,&__self](int i,int j)
        {
            H[i][j]=__self->atoms->H[i][j];
            B[i][j]=__self->atoms->B[i][j];
        });
        
        type0 (&strain)[__dim__][__dim__]=f.val<0>();
        
        Algebra::Do<__dim__>::func([&strain](int i){strain[i][i]++;});
        type0 F[__dim__][__dim__]{DESIG2(__dim__,__dim__,0.0)};
        type0 __H[__dim__][__dim__]{DESIG2(__dim__,__dim__,0.0)};
        Algebra::MLT_mul_MSQ(H,strain,__H);
        
        Algebra::MSQ_2_MLT(__H,H);
        Algebra::MLT_mul_MLT(B,H,F);
        if(Algebra::MLT_det(F)<=0.0)
        {
            PyErr_SetString(PyExc_TypeError,"strain should result in transformation tensor with positive determinant");
            return NULL;
        }
        
        Algebra::MLT_inv(H,B);
        
        __self->atoms->x2s_lcl();
        Algebra::V_eq<__dim__*__dim__>(&H[0][0],&(__self->atoms->H[0][0]));
        __self->atoms->update_H();
        __self->atoms->s2x_lcl();
        /*
        type0* x=__self->atoms->x->begin();
        int x_dim=__self->atoms->x->dim;
        int natms_lcl=__self->atoms->natms_lcl;
        for(int i=0;i<natms_lcl;i++,x+=x_dim)
            Algebra::V_mul_MLT(x,F,x);
        */
        Py_RETURN_NONE;
    });
    
    tp_methods.ml_doc=R"---(
    strain(E)
    
    Strain the system
        
    Parameters
    ----------
    E : double[dim][dim]
       Strain tensor, here dim is the dimension of simulation
    
    Returns
    -------
    None
   
    )---";
}
/*--------------------------------------------
 python constructor
 --------------------------------------------*/
void Atoms::ml_mul(PyMethodDef& tp_methods)
{
    tp_methods.ml_flags=METH_VARARGS | METH_KEYWORDS;
    tp_methods.ml_name="mul";
    tp_methods.ml_meth=(PyCFunction)(PyCFunctionWithKeywords)(
    [](PyObject* self,PyObject* args,PyObject* kwds)->PyObject*
    {
        FuncAPI<int[__dim__]> f("mul",{"n"});
        
        f.get<0>().logics[0]=VLogics("ge",1);
        if(f(args,kwds)) return NULL;
        type0 n[__dim__];
        Algebra::Do<__dim__>::func([&n,&f](int i){n[i]=static_cast<type0>(f.val<0>()[i]);});
        /*
        
        Atoms::Object* __self=reinterpret_cast<Atoms::Object*>(self);
        
        type0 H[__dim__][__dim__]{DESIG2(__dim__,__dim__,0.0)};
        type0 B[__dim__][__dim__]{DESIG2(__dim__,__dim__,0.0)};
        Algebra::DoLT<__dim__>::func([&H,&B,&__self](int i,int j)
        {
            H[i][j]=__self->atoms->H[i][j];
            B[i][j]=__self->atoms->B[i][j];
        });
        
        type0 (&strain)[__dim__][__dim__]=f.val<0>();
        
        Algebra::Do<__dim__>::func([&strain](int i){strain[i][i]++;});
        type0 F[__dim__][__dim__]{DESIG2(__dim__,__dim__,0.0)};
        type0 __H[__dim__][__dim__]{DESIG2(__dim__,__dim__,0.0)};
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
        int natms_lcl=__self->atoms->natms_lcl;
        for(int i=0;i<natms_lcl;i++,x+=x_dim)
            Algebra::V_mul_MLT(x,F,x);*/
        
        Py_RETURN_NONE;
    });
    
    tp_methods.ml_doc=R"---(
    strain(E)
    
    Strain the system
        
    Parameters
    ----------
    E : double[dim][dim]
       Strain tensor, here dim is the dimension of simulation
    
    Returns
    -------
    None
   
    )---";
}
/*--------------------------------------------
 
 --------------------------------------------*/
void Atoms::ml_do(PyMethodDef& tp_method)
{
    tp_method.ml_flags=METH_VARARGS | METH_KEYWORDS;
    tp_method.ml_name="do";
    tp_method.ml_meth=(PyCFunction)(PyCFunctionWithKeywords)(
    [](PyObject* self,PyObject* args,PyObject* kwds)->PyObject*
    {
        
        FuncAPI<OB<PyFunctionObject,PyFunction_Type>> f("do",{"func"});
        if(f(args,kwds)) return NULL;
        
        PyObject* op=f.val<0>();
        
        Atoms::Object* __self=reinterpret_cast<Atoms::Object*>(self);
        try
        {
            __self->atoms->DO(op);
        }
        catch(std::string& err_msg)
        {
            PyErr_SetString(PyExc_TypeError,err_msg.c_str());
            return NULL;
        }

        Py_RETURN_NONE;
    });
    tp_method.ml_doc=R"---(
    do(func)
    
    acts a given function on every atom
        
    Parameters
    ----------
    func : callable python object
       function
    
    Returns
    -------
    None
   
    )---";
}
