#ifndef __MAPP__atoms_dmd__
#define __MAPP__atoms_dmd__
#include "atoms.h"


namespace MAPP_NS
{
    template<class>
    class DMDVec;
    class AtomsDMD: public Atoms
    {
    private:
    protected:
    public:
        DMDVec<type0>* alpha;
        DMDVec<type0>* c;
        DMDVec<type0>* c_d;
        DMDVec<bool>* c_dof;
        DMDVec<bool>* alpha_dof;
        Vec<elem_type>* elem;

        AtomsDMD(MPI_Comm&,int,int,int);
        ~AtomsDMD();
        AtomsDMD& operator=(const Atoms&);
        
        type0 max_alpha;
        const int c_dim;
        const int nelems;
        
        const int N;
        type0* xi;
        type0* wi;
        
        typedef struct
        {
            PyObject_HEAD
            class AtomsDMD* atoms;
            class ForceFieldDMD* ff;
        }Object;
        
        static PyTypeObject TypeObject;
        static PyObject* __new__(PyTypeObject*,PyObject*, PyObject*);
        static int __init__(PyObject*, PyObject*,PyObject*);
        static PyObject* __alloc__(PyTypeObject*,Py_ssize_t);
        static void __dealloc__(PyObject*);
        
        static PyGetSetDef getset[];
        static void setup_tp_getset();
        
        static PyMethodDef methods[];
        static void setup_tp_methods();
        
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
namespace MAPP_NS
{
    template<typename T>
    class DMDVec: public Vec<T>
    {
    private:
        class AtomsDMD* atoms_dmd;
    protected:
    public:
        const T def_val;
        const int c_dim;
        const int real_dim;
        DMDVec(class AtomsDMD* __atoms,T __def_val,const char* __name):
        Vec<T>(__atoms,__atoms->c_dim,__name),
        def_val(__def_val),
        atoms_dmd(__atoms),
        c_dim(__atoms->c_dim),
        real_dim(__atoms->nelems)
        {}
        DMDVec(class AtomsDMD* __atoms,T __def_val):
        Vec<T>(__atoms,__atoms->c_dim),
        def_val(__def_val),
        atoms_dmd(__atoms),
        c_dim(__atoms->c_dim),
        real_dim(__atoms->nelems)
        {}
        ~DMDVec(){}
        int ndim_dump()const{return real_dim;};
        void init_dump();
        void fin_dump();
        void print(FILE*,int);
    };
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<typename T>
inline void DMDVec<T>::init_dump()
{
    if(c_dim==real_dim)
        return vec::init_dump();
    
    
    int rank=atoms_dmd->comm_rank;
    int natms_lcl=atoms_dmd->natms_lcl;
    size_t sz;
    if(rank)
        sz=real_dim*sizeof(T)*natms_lcl;
    else
        sz=real_dim*sizeof(T)*atoms_dmd->natms;
    
    if(sz) vec::data_dump=new byte[sz];
    
    T* __data_dump=reinterpret_cast<T*>(vec::data_dump);
    T* __data=reinterpret_cast<T*>(vec::data);
    elem_type* elem=atoms_dmd->elem->begin();
    
    
    for(int i=0;i<natms_lcl;i++)
    {
        for(int j=0;j<real_dim;j++)
            __data_dump[j]=def_val;
        
        for(int j=0;j<c_dim;j++)
            __data_dump[elem[j]]=__data[j];
        
        __data_dump+=real_dim;
        __data+=vec::dim;
        elem+=vec::dim;
    }
    vec::byte_sz=(vec::byte_sz*real_dim)/c_dim;
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<typename T>
inline void DMDVec<T>::fin_dump()
{
    if(c_dim==real_dim)
        return vec::fin_dump();
    
    vec::byte_sz=(vec::byte_sz*c_dim)/real_dim;
    delete [] vec::data_dump;
    vec::data_dump=NULL;
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<typename T>
inline void DMDVec<T>::print(FILE* fp,int i)
{
    T* __data_dump=Vec<T>::begin_dump()+real_dim*i;
    for(int j=0;j<real_dim;j++)
        fprintf(fp,Vec<T>::print_format,__data_dump[j]);
}
#endif
