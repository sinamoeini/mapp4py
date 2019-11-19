#ifndef __MAPP__atoms__
#define __MAPP__atoms__
#include "api.h"
#include <typeinfo>
#include <cmath>
#include "comm.h"
#include "elements.h"
#include "macros.h"
#include "global.h"

enum {snd_to_bhnd,snd_to_frnt};
enum {rcv_fm_frnt,rcv_fm_bhnd};
/*-----------------------
 _     _   _____   _____  
| |   / / | ____| /  ___| 
| |  / /  | |__   | |     
| | / /   |  __|  | |     
| |/ /    | |___  | |___  
|___/     |_____| \_____| 
 -----------------------*/
namespace MAPP_NS
{
    class vec
    {
    private:
    protected:
    public:
        
        bool __is_empty__;
        int dim;
        int byte_sz;
        int vec_sz;
        int vec_cpcty;
        class Atoms* atoms;
        const char* name;
        byte* data;
        byte* data_dump;
        static constexpr unsigned int vec_grw=128;
        
        void* begin()const{return data;};
        void* end()const{return data+vec_sz*byte_sz;};

        vec(class Atoms*,int,size_t,const char*);
        vec(class Atoms*,const vec&);
        virtual ~vec();
        
        void change_dim(int);
        
        void reserve(int);
        void resize(int);
        void shrink_to_fit();
        void replicate(int);
        void rearrange(int*,int*,int,int);
        void pst_to(byte*&,int);
        
        void add();
        void del(int&);
        
        void pop_out(byte*&,int);
        void pop_in(byte*&);
        
        void cpy(byte*&,int);
        void pst(byte*&,int,int);
        void cpy_pst(int);
        
        void cpy(byte*&,int*,int);
        void pst(byte*&,int);
        void cpy_pst(int*,int);
        bool is_empty()const{return __is_empty__;};
        
        
        void* begin_dump()const{return data_dump;}
        virtual int ndim_dump()const{return dim;};
        virtual void init_dump();
        virtual void fin_dump();
        virtual void print(FILE*,int)=0;
    };
}
#include <cstring>
using namespace MAPP_NS;
/*--------------------------------------------
 
 --------------------------------------------*/
inline void vec::change_dim(int new_dim)
{
    if(new_dim==dim) return;

    int new_byte_sz=(byte_sz/dim)*new_dim;
    int min_byte_sz=MIN(new_byte_sz,byte_sz);
    byte* __data=new byte[new_byte_sz*vec_cpcty];
    if(!is_empty())
    {
        for(int i=0;i<vec_sz;i++)
            memcpy(__data+i*new_byte_sz,data+i*byte_sz,min_byte_sz);
        delete [] data;
        data=__data;
    }
    byte_sz=new_byte_sz;
    dim=new_dim;
}
/*--------------------------------------------
 
 --------------------------------------------*/
inline void vec::reserve(int xtra_atms)
{
    if(vec_cpcty>=vec_sz+xtra_atms)
        return;
    
    vec_cpcty=vec_sz+xtra_atms+vec_grw;
    byte* __data=new byte[vec_cpcty*byte_sz];
    
    memcpy(__data,data,byte_sz*vec_sz);

    delete [] data;
    data=__data;
    
}
/*--------------------------------------------
 
 --------------------------------------------*/
inline void vec::resize(int n)
{
    if(vec_cpcty>=n)
    {
        vec_sz=n;
        return;
    }
    
    
    vec_cpcty=n+vec_grw;
    byte* __data=NULL;
    if(n) __data=new byte[vec_cpcty*byte_sz];
    memcpy(__data,data,byte_sz*vec_sz);
    delete [] data;
    data=__data;
    vec_sz=n;
}
/*--------------------------------------------
 
 --------------------------------------------*/
inline void vec::shrink_to_fit()
{
    if(vec_cpcty==vec_sz)
        return;
    
    byte* __data=NULL;
    if(vec_sz) __data=new byte[vec_sz*byte_sz];
    memcpy(__data,data,vec_sz*byte_sz);
    
    delete [] data;
    data=__data;
    vec_cpcty=vec_sz;
}
/*--------------------------------------------
 
 --------------------------------------------*/
inline void vec::replicate(int n)
{
    if(__is_empty__) return;
    
    int old_vec_sz=vec_sz;
    resize(n*vec_sz);
    byte* __data=data+old_vec_sz*byte_sz;
    for(int i=1;i<n;i++,__data+=old_vec_sz*byte_sz)
        memcpy(__data,data,old_vec_sz*byte_sz);
        
}
/*--------------------------------------------
 
 --------------------------------------------*/
inline void vec::rearrange(int* old_pos
,int* new_pos,int sz,int new_vec_cpcty)
{
    byte* __data=NULL;
    if(new_vec_cpcty) __data=new byte[new_vec_cpcty*byte_sz];
    
    for(int i=0;i<sz;i++)
        memcpy(__data+new_pos[i]*byte_sz,data+old_pos[i]*byte_sz,byte_sz);
    
    delete [] data;
    data=__data;
    vec_sz=0;
    vec_cpcty=new_vec_cpcty;
}
/*--------------------------------------------
 
 --------------------------------------------*/
inline void vec::pst_to(byte*& buff,int iatm)
{
    memcpy(data+iatm*byte_sz,buff,byte_sz);
    buff+=byte_sz;
}
/*-----------------------------------------------------------------
  _____   _____       ___  ___   _____
 /  ___| /  ___|     /   |/   | /  ___|
 | |     | |        / /|   /| | | |
 | |  _  | |       / / |__/ | | | |
 | |_| | | |___   / /       | | | |___
 \_____/ \_____| /_/        |_| \_____|
 -----------------------------------------------------------------*/
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
inline void vec::add()
{
    if(__is_empty__) return;
    resize(vec_sz+1);
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
inline void vec::del(int &del_idx)
{
    if(__is_empty__) return;
    vec_sz--;
    if(del_idx==vec_sz) return;
    memcpy(data+del_idx*byte_sz,data+vec_sz*byte_sz,byte_sz);
}
/*-----------------------------------------------------------------
 _____  __    __  _____   _   _       ___   __   _   _____   _____  
| ____| \ \  / / /  ___| | | | |     /   | |  \ | | /  ___| | ____| 
| |__    \ \/ /  | |     | |_| |    / /| | |   \| | | |     | |__   
|  __|    }  {   | |     |  _  |   / / | | | |\   | | |  _  |  __|  
| |___   / /\ \  | |___  | | | |  / /  | | | | \  | | |_| | | |___  
|_____| /_/  \_\ \_____| |_| |_| /_/   |_| |_|  \_| \_____/ |_____|
 
 -----------------------------------------------------------------*/
/*--------------------------------------------
 1. removes atom iatm and copy it, using buff
 as destenation
 2. replaces iatm with last atom
 3. moves buff to the end of copied data
 4. decreases the size
 ** used during the EXCHANGE of atoms
 --------------------------------------------*/
inline void vec::pop_out(byte*& buff,int iatm)
{
    memcpy(buff,data+iatm*byte_sz,byte_sz);
    memcpy(data+iatm*byte_sz,data+(vec_sz-1)*byte_sz,byte_sz);
    buff+=byte_sz;
    vec_sz--;
}
/*--------------------------------------------
 1. adds one atom to the bottom of the vector
 using buff as source
 2. moves buff to the end of copied data
 3. increases the size
 !! does NOT check for capicity
  ** used during the EXCHANGE of atoms
 --------------------------------------------*/
inline void vec::pop_in(byte*& buff)
{
    memcpy(data+vec_sz*byte_sz,buff,byte_sz);
    buff+=byte_sz;
    vec_sz++;
}
/*------------------------------------------------------------------------------------------------
 _____   _   _       ___   __   _   _____   _____       ___  ___        _       _   _____   _____
|  _  \ | | | |     /   | |  \ | | |_   _| /  _  \     /   |/   |      | |     | | /  ___/ |_   _| 
| |_| | | |_| |    / /| | |   \| |   | |   | | | |    / /|   /| |      | |     | | | |___    | |   
|  ___/ |  _  |   / / | | | |\   |   | |   | | | |   / / |__/ | |      | |     | | \___  \   | |   
| |     | | | |  / /  | | | | \  |   | |   | |_| |  / /       | |      | |___  | |  ___| |   | |   
|_|     |_| |_| /_/   |_| |_|  \_|   |_|   \_____/ /_/        |_|      |_____| |_| /_____/   |_|   
 
 ------------------------------------------------------------------------------------------------*/
/*--------------------------------------------
 1. copies atom iatm, using buff as 
 destenation
 2. moves buff to the end of copied data
 3. does NOT increase/decrease the size
 ** used during the constructing the PHANTOM
 atoms
 --------------------------------------------*/
inline void vec::cpy(byte*& buff,int iatm)
{
    memcpy(buff,data+iatm*byte_sz,byte_sz);
    buff+=byte_sz;
}
/*--------------------------------------------
 1. adds xtra_natms number of atom to the 
 bottom of the vector using buff as source
 2. moves buff to the end of 1st atom copied
 3. increases the size by xtra_natms
 !! DOES check for capicity
 ** used during the constructing the PHANTOM
 atoms
 --------------------------------------------*/
inline void vec::pst(byte*& buff,int stride,int xtra_natms)
{
    reserve(xtra_natms);
    byte* __data=data+vec_sz*byte_sz;
    for(int i=0;i<xtra_natms;i++)
    {
        memcpy(__data,buff,byte_sz);
        buff+=stride;
        __data+=byte_sz;
    }
    buff+=byte_sz-stride*xtra_natms;
    vec_sz+=xtra_natms;
}
/*--------------------------------------------
 1. copies from atom iatm to the bottom of vector
 2. increases the size by one
 !! DOES check for capicity
 ** used during the constructing the PHANTOM
 atoms
 --------------------------------------------*/
inline void vec::cpy_pst(int iatm)
{
    reserve(1);
    memcpy(data+vec_sz*byte_sz,data+iatm*byte_sz,byte_sz);
    vec_sz++;
}
/*-----------------------------------------------
 _   _   _____   _____       ___   _____   _____
| | | | |  _  \ |  _  \     /   | |_   _| | ____| 
| | | | | |_| | | | | |    / /| |   | |   | |__   
| | | | |  ___/ | | | |   / / | |   | |   |  __|  
| |_| | | |     | |_| |  / /  | |   | |   | |___  
\_____/ |_|     |_____/ /_/   |_|   |_|   |_____| 
 
 -----------------------------------------------*/
/*--------------------------------------------
 1. copies no atoms, given by lst,
 using buff as destenation
 2. moves buff to the end of copied data
 3. does NOT increase/decrease the size
 ** used during the UPDATING the phantom atoms
 --------------------------------------------*/
inline void vec::cpy(byte*& buff,int* lst,int no)
{
    for(int i=0;i<no;i++)
    {
        memcpy(buff,data+lst[i]*byte_sz,byte_sz);
        buff+=byte_sz;
    }
}
/*--------------------------------------------
 1. adds xtra_natms number of atom to the
 bottom of the vector using buff as source
 2. moves buff by specific stride
 3. increases the size by xtra_natms
 !! does NOT check for capicity
 ** used during the UPDATING the phantom atoms
 --------------------------------------------*/
inline void vec::pst(byte*& buff,int xtra_natms)
{
    memcpy(data+vec_sz*byte_sz,buff,byte_sz*xtra_natms);
    buff+=byte_sz*xtra_natms;
    vec_sz+=xtra_natms;
}
/*--------------------------------------------
 1. copies no atoms, given by lst, to the
 bottom of vector
 2. increases the size by no
 !! does NOT check for capicity
 ** used during the UPDATING the phantom atoms
 --------------------------------------------*/
inline void vec::cpy_pst(int* lst,int no)
{
    byte* vec_buff=data+vec_sz*byte_sz;
    for(int i=0;i<no;i++)
    {
        memcpy(vec_buff,data+lst[i]*byte_sz,byte_sz);
        vec_buff+=byte_sz;
    }
    vec_sz+=no;
}
/*-----------------------
 _     _   _____   _____  
| |   / / | ____| /  ___| 
| |  / /  | |__   | |     
| | / /   |  __|  | |     
| |/ /    | |___  | |___  
|___/     |_____| \_____| 
 -----------------------*/
namespace MAPP_NS
{
    
    template<typename T>
    class Vec: public vec
    {
    private:
    protected:
    public:
        T empty_val;
        Vec(Atoms*,int,const char*);
        Vec(Atoms*,int);
        Vec(Atoms*,const Vec&);
        ~Vec();
        unsigned int size()const{return dim*vec_sz;}
        static MPI_Datatype MPI_T;
        static const char* print_format;
        T* begin()const{return reinterpret_cast<T*>(data);};
        T* end()const{return reinterpret_cast<T*>(data+vec_sz*byte_sz);};
        T* begin_dump()const{return reinterpret_cast<T*>(data_dump);};
        void empty(const T);
        void fill();
        void change_dim(int,const T);
        void append(const Vec&,const T);
        void resize(int,const T);
        virtual void print(FILE*,int);
    };
}
/*--------------------------------------------*/
template<> MPI_Datatype Vec<bool>::MPI_T;
template<> MPI_Datatype Vec<char>::MPI_T;
template<> MPI_Datatype Vec<short>::MPI_T;
template<> MPI_Datatype Vec<int>::MPI_T;
template<> MPI_Datatype Vec<long int>::MPI_T;
template<> MPI_Datatype Vec<long long>::MPI_T;
template<> MPI_Datatype Vec<unsigned char>::MPI_T;
template<> MPI_Datatype Vec<unsigned short>::MPI_T;
template<> MPI_Datatype Vec<unsigned int>::MPI_T;
template<> MPI_Datatype Vec<unsigned long int>::MPI_T;
template<> MPI_Datatype Vec<unsigned long long>::MPI_T;
template<> MPI_Datatype Vec<float>::MPI_T;
template<> MPI_Datatype Vec<double>::MPI_T;
template<> MPI_Datatype Vec<long double>::MPI_T;

template<> const char* Vec<bool>::print_format;
template<> const char* Vec<char>::print_format;
template<> const char* Vec<short>::print_format;
template<> const char* Vec<int>::print_format;
template<> const char* Vec<long int>::print_format;
template<> const char* Vec<long long>::print_format;
template<> const char* Vec<unsigned char>::print_format;
template<> const char* Vec<unsigned short>::print_format;
template<> const char* Vec<unsigned int>::print_format;
template<> const char* Vec<unsigned long int>::print_format;
template<> const char* Vec<unsigned long long>::print_format;
template<> const char* Vec<float>::print_format;
template<> const char* Vec<double>::print_format;
template<> const char* Vec<long double>::print_format;
/*--------------------------------------------*/
/*-------------------------------------------------
      ___   _____   _____       ___  ___   _____
     /   | |_   _| /  _  \     /   |/   | /  ___/
    / /| |   | |   | | | |    / /|   /| | | |___
   / / | |   | |   | | | |   / / |__/ | | \___  \
  / /  | |   | |   | |_| |  / /       | |  ___| |
 /_/   |_|   |_|   \_____/ /_/        |_| /_____/
 -------------------------------------------------*/
namespace MAPP_NS
{
    class Atoms
    {
    private:
    protected:
    public:
        
        Elements elements;
        Communication comm;
        MPI_Comm& world;
        type0 (&s_lo)[__dim__];
        type0 (&s_hi)[__dim__];
        
        
        const int comm_size;
        const int comm_rank;
        int natms_lcl;
        int natms_ph;
        int natms;
        int step;
        
        //keeep these
        //dont know what to do with this
        type0 max_cut;
        
        
        // boltzmann constant
        type0 kB;
        //planck constant
        type0 hP;
        type0 vol;
        type0 depth_inv[__dim__];
        type0 H[__dim__][__dim__];
        type0 B[__dim__][__dim__];
        type0 __h[__nvoigt__];
        type0 __b[__nvoigt__];
        
        
        void update_H();
        Vec<type0>* x;
        Vec<id_type>* id;
        Vec<bool>* x_dof;
        vec** vecs;
        int nvecs;
        vec** dynamic_vecs;
        int ndynamic_vecs;
        
        Atoms(MPI_Comm&);
        Atoms(const Atoms&);
        virtual ~Atoms();
        Atoms& operator=(const Atoms&);
        Atoms& operator+(const Atoms&);
        Atoms& operator+=(const Atoms&);
        
        void import_vecs(const Atoms&);
        vec* find_vec(const char*);
        void push(vec*);
        void pop(vec*);
        void x2s(int);
        void s2x(int);
        void x2s_lcl();
        void s2x_lcl();
        void x2s_all();
        void s2x_all();
        void x2s_dump();
        id_type get_max_id();
        
        void insert(byte*,vec**,int,int);
        void add();
        void del(int&);
        void restart();
        int is_dynamic(vec*);
        
        void reset_domain();
        bool xchng_id(unsigned long&);
        virtual void DO(PyObject*){}
        
        
        typedef struct
        {
            PyObject_HEAD
            class Atoms* atoms;
            class ForceField* ff;
        }Object;
        
        static PyTypeObject TypeObject;
        static PyObject* __new__(PyTypeObject*,PyObject*, PyObject*);
        static int __init__(PyObject*, PyObject*,PyObject*);
        static PyObject* __alloc__(PyTypeObject*,Py_ssize_t);
        static void __dealloc__(PyObject*);
        
        
        static PyMethodDef methods[];
        static void setup_tp_methods();
        static void ml_strain(PyMethodDef&);
        static void ml_mul(PyMethodDef&);
        static void ml_cell_change(PyMethodDef&);
        static void ml_autogrid(PyMethodDef&);
        static void ml_do(PyMethodDef&);
        
        
        static PyGetSetDef getset[];
        static void setup_tp_getset();
        static void getset_natms(PyGetSetDef&);
        static void getset_step(PyGetSetDef&);
        static void getset_hP(PyGetSetDef&);
        static void getset_kB(PyGetSetDef&);
        static void getset_H(PyGetSetDef&);
        static void getset_B(PyGetSetDef&);
        static void getset_vol(PyGetSetDef&);
        static void getset_elems(PyGetSetDef&);
        static void getset_skin(PyGetSetDef&);
        static void getset_comm_rank(PyGetSetDef&);
        static void getset_comm_size(PyGetSetDef&);
        static void getset_comm_coords(PyGetSetDef&);
        static void getset_comm_dims(PyGetSetDef&);
        
        
        static int setup_tp();

        
    };
}
/*-----------------------
 _     _   _____   _____  
| |   / / | ____| /  ___| 
| |  / /  | |__   | |     
| | / /   |  __|  | |     
| |/ /    | |___  | |___  
|___/     |_____| \_____| 
 
 -----------------------*/
/*--------------------------------------------
 
 --------------------------------------------*/
inline vec::vec(Atoms* __atoms,int __dim,size_t __t_size,const char* __name):
__is_empty__(false),
dim(__dim),
byte_sz(static_cast<int>(__t_size)*__dim),
vec_sz(0),
vec_cpcty(0),
atoms(__atoms),
name(__name),
data(NULL),
data_dump(NULL)
{
    atoms->push(this);
    reserve(atoms->natms_lcl+atoms->natms_ph);
    resize(atoms->natms_lcl);
}
/*--------------------------------------------
  copy constructor
 --------------------------------------------*/
inline vec::vec(Atoms* __atoms,const vec& other):
__is_empty__(other.__is_empty__),
dim(other.dim),
byte_sz(other.byte_sz),
vec_sz(other.vec_sz),
vec_cpcty(other.vec_cpcty),
atoms(__atoms),
name(other.name),
data(NULL),
data_dump(NULL)
{
    atoms->push(this);
    if(!other.is_empty() && vec_cpcty)
    {
        data=new byte[vec_cpcty*byte_sz];
        memcpy(data,other.data,byte_sz*vec_sz);
    }
}
/*--------------------------------------------
 
 --------------------------------------------*/
inline vec::~vec()
{
    delete [] data;
    atoms->pop(this);
}
/*--------------------------------------------
 
 --------------------------------------------*/
inline void vec::init_dump()
{
    
    if(atoms->comm_rank)
        data_dump=data;
    else
    {
        int natms=atoms->natms;
        if(natms)
            data_dump=new byte[natms*byte_sz];
        memcpy(data_dump,data,byte_sz*atoms->natms_lcl);
    }
}
/*--------------------------------------------
 
 --------------------------------------------*/
inline void vec::fin_dump()
{
    if(atoms->comm_rank==0)
        delete [] data_dump;
    data_dump=NULL;
}
/*--------------------------------------------
 constructor with name;
 --------------------------------------------*/
template<typename T>
Vec<T>::Vec(Atoms* __atoms,int __dim,const char* __name):
vec(__atoms,__dim,sizeof(T),__name)
{}
/*--------------------------------------------
 constructor without name;
 --------------------------------------------*/
template<typename T>
Vec<T>::Vec(Atoms* __atoms,int __dim):
vec(__atoms,__dim,sizeof(T),NULL)
{}
/*--------------------------------------------
 constructor without name;
 --------------------------------------------*/
template<typename T>
Vec<T>::Vec(Atoms* __atoms,const Vec<T>& other):
vec(__atoms,other)
{
    empty_val=other.empty_val;
}
/*--------------------------------------------
 destructor
 --------------------------------------------*/
template<typename T>
inline Vec<T>::~Vec()
{}
/*--------------------------------------------
 filling
 --------------------------------------------*/
template<typename T>
inline void Vec<T>::empty(const T __empty_val)
{
    empty_val=__empty_val;
    if(__is_empty__) return;
    delete [] data;
    data=NULL;
    vec_sz=vec_cpcty=0;
    __is_empty__=true;
}
/*--------------------------------------------
 filling
 --------------------------------------------*/
template<typename T>
inline void Vec<T>::fill()
{
    if(__is_empty__)
    {
        reserve(atoms->natms_lcl+atoms->natms_ph);
        vec::resize(atoms->natms_lcl);
    }
    
    
    
    T* __data=begin();
    for(int i=0;i<dim*vec_sz;i++) __data[i]=empty_val;
    __is_empty__=false;
}
/*--------------------------------------------
 this function changes dimension of the vector
 if __dim is bigger than the original it will
 assign fill_val to the rest of components.
 Otherwise, it will discard extra ones.
 --------------------------------------------*/
template<typename T>
inline void Vec<T>::change_dim(int __dim,const T fill_val)
{
    int old_dim=dim;
    vec::change_dim(__dim);
    if(__dim<=old_dim)
        return;
    
    T* __data=reinterpret_cast<T*>(data);
    for(int i=0;i<vec_sz;i++,__data+=dim)
        for(int j=old_dim;j<dim;j++)
            __data[j]=fill_val;
}
/*--------------------------------------------
 like regular resize, but extra values will
 be filled in
 --------------------------------------------*/
template<typename T>
inline void Vec<T>::resize(int n,const T fill_val)
{
    int n0=vec_sz*dim;
    vec::resize(n);
    T* __data=reinterpret_cast<T*>(data);
    int n1=vec_sz*dim;
    for(int i=n0;i<n1;i++)
        __data[i]=fill_val;
}
/*--------------------------------------------
 append a vec to the bottom of this vec
 
 high level function
 NOT EFFICIENT DO NOT USE FREQUENTLY
 --------------------------------------------*/
template<typename T>
inline void Vec<T>::append(const Vec<T>& other,const T fill_val)
{
    if((is_empty()&& other.is_empty()) && empty_val==other.empty_val)
        return;
    
    if(is_empty()) fill();
    if(other.is_empty())
    {
        resize(vec_sz+other.atoms->natms_lcl,other.empty_val);
        return;
    }
    
    int __vec_sz=vec_sz;
    resize(vec_sz+other.vec_sz,fill_val);
    T* __data=reinterpret_cast<T*>(data)+__vec_sz*dim;
    T* __other_data=reinterpret_cast<T*>(other.data);
    if(dim==other.dim)
    {
        memcpy(__data,__other_data,dim*other.vec_sz*sizeof(T));
        return;
    }

    int min_dim=MIN(dim,other.dim);
    for(int i=0;i<other.vec_sz;i++)
    {
        memcpy(__data,__other_data,min_dim*sizeof(T));
        __data+=dim;
        __other_data+=other.dim;
    }
}
/*--------------------------------------------

 --------------------------------------------*/
template<typename T>
inline void Vec<T>::print(FILE* fp,int i)
{
    T* __data_dump=begin_dump()+dim*i;
    for(int j=0;j<dim;j++)
        fprintf(fp,print_format,__data_dump[j]);
}














/*------------------------------------------------------------------------------------------------
  
 ------------------------------------------------------------------------------------------------*/

namespace MAPP_NS
{
    template<typename T,class F>
    class VecPy
    {
    private:
    protected:
    public:
        size_t ipos;
        bool inc;
        PyObject* op;
        int dim;
        std::remove_pointer<npy_intp>::type npy_dim;
        
        T* data;
        T* arr_data;
        T* head;
        
        F& func;
        Vec<T>* vec;
        
        VecPy(Vec<T>*,F&);
        virtual ~VecPy();
        
        void init();
        void fin();
        virtual void setup_arr_data();
        virtual void pre_iter();
        virtual void post_iter();
        //bool iter();
        
    };
}


/*--------------------------------------------
 
 --------------------------------------------*/
template<typename T,class F>
VecPy<T,F>::VecPy(Vec<T>* __vec,F& __func):
ipos(0),
inc(false),
op(NULL),
dim(__vec->dim),
data(NULL),
arr_data(NULL),
func(__func),
vec(__vec)
{
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<typename T,class F>
void VecPy<T,F>::setup_arr_data()
{
    arr_data=new T[dim];
    npy_dim=static_cast<std::remove_pointer<npy_intp>::type>(dim);
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<typename T,class F>
void VecPy<T,F>::init()
{
    int natms_lcl=vec->atoms->natms_lcl;
    if(natms_lcl) data=new T[natms_lcl*dim];
    if(vec->is_empty())
    {
        T empty_val=vec->empty_val;
        for(int i=0;i<natms_lcl*dim;i++)
        {
            data[i]=empty_val;
        }
    }
    else
    {
        memcpy(data,vec->begin(),natms_lcl*dim*sizeof(T));
    }
    
    setup_arr_data();
//    if(npy_dim==1)
//    {
//        op=PyArray_SimpleNewFromData(0,&npy_dim,cpp_type2type_num<T>::type_num(),arr_data);
//        std::cout<< "ptr before "<< op << " CNT " << Py_REFCNT(op)<<std::endl;
//        op=PyArray_Return((PyArrayObject* )op);
//
//        std::cout<< "ptr after "<< op << " CNT " << Py_REFCNT(op)<<std::endl;
//        std::cout<<"SIZE  SCALAR  " << PyArray_DescrFromScalar(op)->elsize << " VS " <<sizeof(T) <<std::endl;
//    }
//    else
    op=PyArray_SimpleNewFromData(1,&npy_dim,cpp_type2type_num<T>::type_num(),arr_data);
    
    //op=PyArray_SimpleNewFromData(1,&npy_dim,cpp_type2type_num<T>::type_num(),arr_data);
    
    head=data;
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<typename T,class F>
void VecPy<T,F>::pre_iter()
{
    memcpy(arr_data,head,sizeof(T)*dim);
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<typename T,class F>
void VecPy<T,F>::post_iter()
{
//    if(npy_dim==1)
//        PyArray_ScalarAsCtype(op,arr_data);
    
    try
    {
        func(head,arr_data);
    }
    catch(std::string& err_msg)
    {
        throw err_msg;
    }
    
    memcpy(head,arr_data,sizeof(T)*dim);
    head+=dim;
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<typename T,class F>
void VecPy<T,F>::fin()
{
    head=NULL;
    if(op)
    {
        Py_DECREF(op);
        op=NULL;
    }
    
    
    
    delete [] arr_data;
    arr_data=NULL;
    if(vec->is_empty())
    {
        int natms_lcl=vec->atoms->natms_lcl;
        T empty_val=vec->empty_val;
        int i=0;
        while(i<natms_lcl*dim && data[i]==empty_val)
            i++;
        
        int is_empty_lcl=0,is_empty;
        if(i!=natms_lcl*dim) is_empty_lcl=1;
        MPI_Allreduce(&is_empty_lcl,&is_empty,1,MPI_INT,MPI_MAX,vec->atoms->world);
        
        if(is_empty)
        {
            vec->fill();
            memcpy(vec->begin(),data,natms_lcl*dim*sizeof(T));
        }
        
    }
    else
    {
        int natms_lcl=vec->atoms->natms_lcl;
        memcpy(vec->begin(),data,natms_lcl*dim*sizeof(T));
    }
    
    delete [] data;
    data=NULL;
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<typename T,class F>
VecPy<T,F>::~VecPy()
{
    if(op)
    {
        Py_DECREF(op);
        op=NULL;
    }
    
    delete [] arr_data;
    delete [] data;
}

/*------------------------------------------------------------------------------------------------
 
 ------------------------------------------------------------------------------------------------*/
namespace MAPP_NS
{
    class VecPyFunc
    {
    private:
    protected:
    public:
        template <typename T,class F, class... Ss>
        static void find(const size_t i,const char* name,VecPy<T,F>& __v,Ss&... ss)
        {
            
            if(strcmp(name,__v.vec->name)==0)
            {
                __v.inc=true;
                __v.ipos=i;
            }
            else
                return find(i,name,ss...);
        }
        
        
        
        
        
        template <typename T,class F, class... Ss>
        static void init(PyObject* tuple,VecPy<T,F>& __v,Ss&... ss)
        {
            if(__v.inc)
            {
                __v.init();
                PyTuple_SET_ITEM(tuple,__v.ipos,__v.op);
                Py_INCREF(__v.op);
            }
            
            init(tuple,ss...);
        }
        
        
        template <typename T,class F, class... Ss>
        static void fin(VecPy<T,F>& __v,Ss&... ss)
        {
            if(__v.inc)
                __v.fin();
            fin(ss...);
        }
        
        
        
        
        template <typename T,class F, class... Ss>
        static void pre_iter(VecPy<T,F>& __v,Ss&... ss)
        {
            if(__v.inc)
                __v.pre_iter();
            
            pre_iter(ss...);
        }
        
        
        template <typename T,class F, class... Ss>
        static void post_iter(VecPy<T,F>& __v,Ss&... ss)
        {
            try
            {
                if(__v.inc)   __v.post_iter();
            }
            catch(std::string& err_msg)
            {
                throw err_msg;
            }
            post_iter(ss...);
        }
        
        
        
        
        static void find(const size_t,const char* name)
        {throw std::string("no vriable named ")+name+std::string(" exists in the system");}
        static void init(PyObject*){}
        static void fin(){}
        static void pre_iter(){}
        static void post_iter(){}
        
        
        
        
        template <class... Ss>
        static void Do(Atoms*,PyObject*,Ss&... ss);
        

        
    };
}
/*--------------------------------------------
 
 --------------------------------------------*/
template <class... Ss>
void VecPyFunc::Do(Atoms* atoms,PyObject* op,Ss&... ss)
{
    PyCodeObject* co=(PyCodeObject *)PyFunction_GET_CODE(op);
    PyObject* co_varnames=co->co_varnames;
    size_t co_nvars=static_cast<size_t>(co->co_argcount);
    
    
    
    try{
        for(size_t i=0;i<co_nvars;i++)
        {
            const char* var_name=__PyString_AsString(PyTuple_GetItem(co_varnames,i));
            find(i,var_name,ss...);
        }
    }
    catch(std::string& err_msg)
    {
        throw err_msg;
    }
    
    
    PyObject* tuple=PyTuple_New(co_nvars);
    init(tuple,ss...);
    
    int natms_lcl=atoms->natms_lcl;
    int err_lcl=0;
    std::string err_msg;
    PyObject* ans=NULL;
    int* del_idx_lcl=natms_lcl==0 ? NULL:new int[natms_lcl];
    int ndel_idx_lcl=0;
    for(int i=0;i<natms_lcl && err_lcl==0;i++)
    {
        pre_iter(ss...);
        
        // DO THE FUNCTION HERE AND FIND OUT IF THERE IS A PROBLEM
        ans=PyEval_CallObject(op,tuple);
        if(PyErr_Occurred()!=NULL)
        {
            err_lcl=1;
            err_msg=std::string("An error occured in the provided function");
            if(ans) Py_DECREF(ans);
            PyErr_Clear();
            continue;
        }
        
        if(ans!=Py_None && ans!=Py_True && ans!=Py_False)
        {
            Py_DECREF(ans);
            err_lcl=1;
            err_msg=std::string("provided function should only return True, False, or None");
            continue;
        }
        
        if(ans==Py_False) del_idx_lcl[ndel_idx_lcl++]=i;
        
        
        Py_DECREF(ans);


        // NO ERROR CHECK YET
        
        
        
        try
        {
            post_iter(ss...);
        }
        catch(std::string& __err_msg)
        {
            err_lcl=1;
            err_msg=__err_msg;
        }
    }
    
    int err;
    MPI_Allreduce(&err_lcl,&err,1,MPI_INT,MPI_MAX,atoms->world);
    if(err)
    {
        int first_proc;
        int comm_rank=0;
        if(err_lcl)
            comm_rank=atoms->comm_rank;
        else
            comm_rank=atoms->comm_size;
        
        MPI_Allreduce(&comm_rank,&first_proc,1,MPI_INT,MPI_MIN,atoms->world);
        int err_msg_length=0;
        if(atoms->comm_rank==first_proc)
            err_msg_length=static_cast<int>(err_msg.size())+1;
        MPI_Bcast(&err_msg_length,1,MPI_INT,first_proc,atoms->world);
        char* __err_msg=NULL;
        __err_msg=new char[err_msg_length];
        if(atoms->comm_rank==first_proc)
            memcpy(__err_msg,err_msg.c_str(),err_msg_length);
        MPI_Bcast(__err_msg,err_msg_length,MPI_CHAR,first_proc,atoms->world);
        err_msg=std::string(__err_msg);
        delete [] __err_msg;
        
        delete [] del_idx_lcl;
        throw err_msg;
        
    }
    
    
    
    Py_DECREF(tuple);
    fin(ss...);
    
    
    int ndel_idx=0;
    MPI_Allreduce(&ndel_idx_lcl,&ndel_idx,1,Vec<int>::MPI_T,MPI_SUM,atoms->world);
    if(ndel_idx)
    {
        for(int i=ndel_idx_lcl-1;i>-1;i--)
            atoms->del(del_idx_lcl[i]);
        
        
        // redo the ids
        natms_lcl=atoms->natms_lcl;
        int ntams_prev=0;
        MPI_Scan(&natms_lcl,&ntams_prev,1,Vec<int>::MPI_T,MPI_SUM,atoms->world);
        ntams_prev-=natms_lcl;
        id_type st=static_cast<id_type>(ntams_prev);
        id_type* __id=atoms->id->begin();
        for(int i=0;i<natms_lcl;i++,st++)
            __id[i]=st;
        
        //update the number of atoms 
        MPI_Allreduce(&natms_lcl,&(atoms->natms),1,Vec<int>::MPI_T,MPI_SUM,atoms->world);
    }
    
    delete [] del_idx_lcl;
    
    
}

#endif
