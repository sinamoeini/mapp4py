#ifndef __MAPP__atoms__
#define __MAPP__atoms__
#include "api.h"
#include <typeinfo>
#include <cmath>
#include "comm.h"
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
        class Atoms* atoms;
        const char* name;
        byte* data;
        int dim;
        int byte_sz;
        unsigned int vec_sz;
        unsigned int vec_cpcty;
        static constexpr unsigned int vec_grw=128;
        
        void* begin()const{return data;};
        void* end()const{return data+vec_sz*byte_sz;};

        vec(class Atoms*,int,size_t,const char*);
        virtual ~vec();
        
        void change_dim(int);
        
        void reserve(int);
        void resize(int);
        void shrink_to_fit();
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
    
    };
}
#include <cstring>
#include <new>
using namespace MAPP_NS;
/*--------------------------------------------
 
 --------------------------------------------*/
inline void vec::change_dim(int new_dim)
{
    if(new_dim==dim) return;

    int new_byte_sz=(byte_sz/dim)*new_dim;
    int min_byte_sz=MIN(new_byte_sz,byte_sz);
    byte* __data=new byte[new_byte_sz*vec_cpcty];
    for(int i=0;i<vec_sz;i++)
        memcpy(__data+i*new_byte_sz,data+i*byte_sz,min_byte_sz);
    
    delete [] data;
    data=__data;
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
        Vec(Atoms*,int,const char*);
        Vec(Atoms*,int);
        ~Vec();
        unsigned int size()const{return dim*vec_sz;}
        static MPI_Datatype MPI_T;
        T* begin()const{return reinterpret_cast<T*>(data);};
        T* end()const{return reinterpret_cast<T*>(data+vec_sz*byte_sz);};
    };
}
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
        const int comm_size;
        const int comm_rank;
        class Elements* elements;
        Communication comm;
        MPI_Comm& world;
        type0 (&s_lo)[__dim__];
        type0 (&s_hi)[__dim__];
        
        //keeep these
        int natms;
        int natms_ph;
        int tot_natms;
        int step;
        type0 H[__dim__][__dim__];
        type0 B[__dim__][__dim__];
        type0 depth_inv[__dim__];
        type0 vol;
        void update_H();
        Vec<type0>* x;
        Vec<unsigned int>* id;
        Vec<bool>* dof;
        vec** vecs;
        int nvecs;
        
        
        Atoms(Communication&);
        Atoms(MPI_Comm&);
        virtual ~Atoms();
        vec* find_vec(const char*);
        void push(vec*);
        void pop(vec*);
        void x2s(int);
        void s2x(int);
        void x2s_lcl();
        void s2x_lcl();
        void x2s_all();
        void s2x_all();
        
        void insert(byte*,vec**,int,int);
        void add();
        void del(int&);
        void restart();
        
        //dont know what to do with this
        type0 max_cut;
        
        
        // boltzmann constant
        type0 kB;
        //planck constant
        type0 h;
        
        
        
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
        
        
        static PyGetSetDef getset[];
        static void setup_tp_getset();
        static void getset_step(PyGetSetDef&);
        static void getset_h(PyGetSetDef&);
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
        
        
        static void setup_tp();

        
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
dim(__dim),
byte_sz(static_cast<int>(__t_size)*__dim),
vec_sz(0),
vec_cpcty(0),
data(NULL),
name(__name),
atoms(__atoms)
{
    atoms->push(this);
    reserve(atoms->natms+atoms->natms_ph);
    resize(atoms->natms);
}
/*--------------------------------------------
 
 --------------------------------------------*/
inline vec::~vec()
{
    delete [] data;
    atoms->pop(this);
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
 destructor
 --------------------------------------------*/
template<typename T>
inline Vec<T>::~Vec()
{}
#endif
