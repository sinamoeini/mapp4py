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

        AtomsDMD(MPI_Comm&,int,int);
        ~AtomsDMD();
        AtomsDMD& operator=(const Atoms&);
        
        type0 max_alpha;
        const int c_dim;
        
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
        real_dim(static_cast<int>(__atoms->elements.nelems))
        {}
        DMDVec(class AtomsDMD* __atoms,T __def_val):
        Vec<T>(__atoms,__atoms->c_dim),
        def_val(__def_val),
        atoms_dmd(__atoms),
        c_dim(__atoms->c_dim),
        real_dim(static_cast<int>(__atoms->elements.nelems))
        {}
        ~DMDVec(){}
        void gather_dump();
        void print(FILE*,int);
    };
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<typename T>
inline void DMDVec<T>::gather_dump()
{
    if(c_dim!=real_dim)
        return Vec<T>::gather_dump();
    
    int natms_lcl=atoms_dmd->natms_lcl;
    T* dump_data_lcl=NULL;
    if(natms_lcl) dump_data_lcl=new T[natms_lcl*real_dim];
    
    T* __dump_data_lcl=dump_data_lcl;
    T* __data=reinterpret_cast<T*>(vec::data);
    elem_type* elem=atoms_dmd->elem->begin();
    for(int i=0;i<natms_lcl;i++)
    {
        for(int j=0;j<real_dim;j++)
            __dump_data_lcl[j]=def_val;
        
        for(int j=0;j<c_dim;j++)
            __dump_data_lcl[elem[j]]=__data[j];
        
        __dump_data_lcl+=real_dim;
        __data+=vec::dim;
        elem+=vec::dim;
    }
    

    int natms=atoms_dmd->natms;
    if(atoms_dmd->comm_rank==0 && natms)
        Vec<T>::dump_data=new T[natms*real_dim];
    
    MPI_Gather(dump_data_lcl,natms_lcl*real_dim,Vec<T>::MPI_T,Vec<T>::dump_data,natms*real_dim,Vec<T>::MPI_T,0,atoms_dmd->world);

    delete [] dump_data_lcl;
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<typename T>
inline void DMDVec<T>::print(FILE* fp,int i)
{
    for(int j=0;j<real_dim;j++)
        fprintf(fp,Vec<T>::print_format,Vec<T>::dump_data[real_dim*i+j]);
}
#endif
