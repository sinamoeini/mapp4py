#ifndef __MAPP__atoms_dmd__
#define __MAPP__atoms_dmd__
#include "atoms.h"


namespace MAPP_NS
{
    template<class> class DMDVec;
    class AtomsDMD: public Atoms
    {
    private:
    protected:
    public:
        type0 S_fe[__dim__][__dim__];
        type0 fe,pe;
        type0 s;

        type0 max_alpha;
        int c_dim;
        
        const int N;
        type0* xi;
        type0* wi;
        DMDVec<type0>* alpha;
        DMDVec<type0>* c;
        DMDVec<bool>* c_dof;
        DMDVec<bool>* alpha_dof;
        Vec<elem_type>* elem;

        void DO(PyObject*);
        AtomsDMD(MPI_Comm&,int,int);
        AtomsDMD(const AtomsDMD&);
        ~AtomsDMD();
        AtomsDMD& operator=(const AtomsDMD&);
        AtomsDMD& operator+(const AtomsDMD&);
        AtomsDMD& operator+=(const AtomsDMD&);
        void import_vecs(const AtomsDMD&);
        void update_max_alpha();
        AtomsDMD& operator=(const Atoms&);
        
        
        type0 vac_msd();
        type0* ave_comp();
        type0 temp;

        
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
        static void getset_temp(PyGetSetDef&);
        static void getset_S_fe(PyGetSetDef&);
        static void getset_fe(PyGetSetDef&);
        static void getset_pe(PyGetSetDef&);
        static void getset_s(PyGetSetDef&);
        static void getset_ave_mu(PyGetSetDef&);
        static void getset_ave_comp(PyGetSetDef&);
        static void getset_ext_mu(PyGetSetDef&);

        static PyMethodDef methods[];
        static void setup_tp_methods();
        
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
namespace MAPP_NS
{
    template<typename T>
    class DMDVec: public Vec<T>
    {
    private:
        class AtomsDMD* atoms_dmd;
        int dump_dim;
    protected:
    public:
        const T def_val;
        DMDVec(class AtomsDMD*,T,const char*);
        DMDVec(class AtomsDMD*,T);
        DMDVec(class AtomsDMD*,const DMDVec&);
        ~DMDVec(){}
        int ndim_dump()const{return static_cast<int>(atoms_dmd->elements.nelems);};
        void init_dump();
        void fin_dump();
        void init_do();
        void fin_do();
        void print(FILE*,int);
        byte* def_val_buff();
    };
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<typename T>
inline DMDVec<T>::DMDVec(class AtomsDMD* __atoms,T __def_val,const char* __name):
Vec<T>(__atoms,__atoms->c_dim,__name),
atoms_dmd(__atoms),
def_val(__def_val)
{
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<typename T>
inline DMDVec<T>::DMDVec(class AtomsDMD* __atoms,T __def_val):
Vec<T>(__atoms,__atoms->c_dim),
def_val(__def_val),
atoms_dmd(__atoms)
{
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<typename T>
inline DMDVec<T>::DMDVec(class AtomsDMD* __atoms,const DMDVec& other):
Vec<T>(__atoms,other),
atoms_dmd(__atoms),
def_val(other.def_val)
{
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<typename T>
inline void DMDVec<T>::init_dump()
{
    dump_dim=static_cast<int>(atoms_dmd->elements.nelems);
    if(vec::dim==dump_dim)
        return vec::init_dump();
    
    
    int rank=atoms_dmd->comm_rank;
    int natms_lcl=atoms_dmd->natms_lcl;
    size_t sz;
    if(rank)
        sz=dump_dim*sizeof(T)*natms_lcl;
    else
        sz=dump_dim*sizeof(T)*atoms_dmd->natms;
    
    if(sz) vec::data_dump=new byte[sz];
    
    T* __data_dump=reinterpret_cast<T*>(vec::data_dump);
    T* __data=reinterpret_cast<T*>(vec::data);
    elem_type* elem=atoms_dmd->elem->begin();
    
    
    for(int i=0;i<natms_lcl;i++)
    {
        for(int j=0;j<dump_dim;j++)
            __data_dump[j]=def_val;
        
        for(int j=0;j<vec::dim;j++)
            __data_dump[elem[j]]=__data[j];
        
        __data_dump+=dump_dim;
        __data+=vec::dim;
        elem+=vec::dim;
    }
    vec::byte_sz=(vec::byte_sz*dump_dim)/vec::dim;
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<typename T>
inline void DMDVec<T>::print(FILE* fp,int i)
{
    T* __data_dump=Vec<T>::begin_dump()+dump_dim*i;
    for(int j=0;j<dump_dim;j++)
        fprintf(fp,Vec<T>::print_format,__data_dump[j]);
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<typename T>
inline void DMDVec<T>::fin_dump()
{
    if(vec::dim==dump_dim)
        return vec::fin_dump();
    
    vec::byte_sz=(vec::byte_sz*vec::dim)/dump_dim;
    delete [] vec::data_dump;
    vec::data_dump=NULL;
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<typename T>
inline void DMDVec<T>::init_do()
{
    dump_dim=static_cast<int>(atoms_dmd->elements.nelems);
    if(vec::dim==dump_dim)
        return;
    
    int natms_lcl=atoms_dmd->natms_lcl;
    size_t sz=dump_dim*natms_lcl;

    T* __new_data=NULL;
    if(sz) __new_data=new byte[sz];
    T* __data=reinterpret_cast<T*>(vec::data);
    elem_type* elem=atoms_dmd->elem->begin();
    
    
    for(int i=0;i<natms_lcl;i++)
    {
        for(int j=0;j<dump_dim;j++)
            __new_data[j]=def_val;
        
        for(int j=0;j<vec::dim;j++)
            __new_data[elem[j]]=__data[j];
        
        __new_data+=dump_dim;
        __data+=vec::dim;
        elem+=vec::dim;
    }
    vec::byte_sz=(vec::byte_sz*dump_dim)/vec::dim;
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<typename T>
inline void DMDVec<T>::fin_do()
{
    if(vec::dim==dump_dim)
        return vec::fin_dump();
    
    vec::byte_sz=(vec::byte_sz*vec::dim)/dump_dim;
    delete [] vec::data_dump;
    vec::data_dump=NULL;
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<typename T>
inline byte* DMDVec<T>::def_val_buff()
{
    byte* __def_val_buff=new byte[sizeof(T)];
    memcpy(__def_val_buff,&def_val,sizeof(T));
    return __def_val_buff;
}
/*------------------------------------------------------------------------------------------------
 
 ------------------------------------------------------------------------------------------------*/

namespace MAPP_NS
{
    template<typename T,class F>
    class DMDVecPy : public VecPy<T,F>
    {
    private:
    protected:
    public:
        DMDVecPy(DMDVec<T>* __vec,F& __func);
        ~DMDVecPy();
        void setup_arr_data();
        void pre_iter();
        void post_iter();
        int tmp_dim;
        
        elem_type* elem;
        elem_type* elem_head;
        
        int* map_elem;
        T* post_arr_data;
        T def_val;
        //bool iter();
        
    };
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<typename T,class F>
DMDVecPy<T,F>::DMDVecPy(DMDVec<T>* __vec,F& __func):
VecPy<T,F>::VecPy(__vec,__func)
{
    elem=reinterpret_cast<AtomsDMD*>(__vec->atoms)->elem->begin();
    def_val=__vec->def_val;
    tmp_dim=__vec->ndim_dump();
    
    map_elem=NULL;
    post_arr_data=NULL;
    if(tmp_dim)
    {
        map_elem=new int[tmp_dim];
        post_arr_data=new T[tmp_dim];
    }
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<typename T,class F>
DMDVecPy<T,F>::~DMDVecPy()
{
    delete [] post_arr_data;
    delete [] map_elem;
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<typename T,class F>
void DMDVecPy<T,F>::setup_arr_data()
{
    VecPy<T,F>::arr_data=new T[tmp_dim];
    VecPy<T,F>::npy_dim=static_cast<std::remove_pointer<npy_intp>::type>(tmp_dim);
    elem_head=elem;
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<typename T,class F>
void DMDVecPy<T,F>::pre_iter()
{
    for(int i=0;i<tmp_dim;i++)
    {
        map_elem[i]=-1;
        VecPy<T,F>::arr_data[i]=post_arr_data[i]=def_val;
    }
    for(int i=0;i<VecPy<T,F>::dim;i++)
    {
        map_elem[elem_head[i]]=i;
        VecPy<T,F>::arr_data[elem_head[i]]=VecPy<T,F>::head[i];
    }
    
    
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<typename T,class F>
void DMDVecPy<T,F>::post_iter()
{
//    if(VecPy<T,F>::npy_dim==1)
//        PyArray_ScalarAsCtype(VecPy<T,F>::op,VecPy<T,F>::arr_data);
    
    for(int i=0;i<tmp_dim;i++)
    {
        if(map_elem[i]==-1 && VecPy<T,F>::arr_data[i]!=def_val)
        {
            throw std::string("cannot change it");
        }
        if(map_elem[i]!=-1)
            post_arr_data[map_elem[i]]=VecPy<T,F>::arr_data[i];
    }
    
    try
    {
        VecPy<T,F>::func(VecPy<T,F>::head,post_arr_data);
    }
    catch(std::string& err_msg)
    {
        throw err_msg;
    }
    for(int i=0;i<tmp_dim;i++)
        if(map_elem[i]!=-1)
            VecPy<T,F>::head[map_elem[i]]=VecPy<T,F>::arr_data[i];
    VecPy<T,F>::head+=VecPy<T,F>::dim;
    elem_head+=VecPy<T,F>::dim;
}
#endif
