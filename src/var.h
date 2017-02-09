#ifndef __MAPP__var__
#define __MAPP__var__
#include <Python/Python.h>
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#define NO_IMPORT_ARRAY
#define PY_ARRAY_UNIQUE_SYMBOL ARRAY_API
#include <numpy/arrayobject.h>
#include <mpi.h>
#include "print.h"

/*------------------------------------------------------------------------------------------------------------------------------------
 
 ------------------------------------------------------------------------------------------------------------------------------------*/
namespace MAPP_NS
{
    template<const int>class type_num2cpp_type{public:typedef void T;};
    template<>class type_num2cpp_type<NPY_BOOL>{public:typedef bool T;};
    template<>class type_num2cpp_type<NPY_BYTE>{public:typedef char T;};
    template<>class type_num2cpp_type<NPY_UBYTE>{public:typedef unsigned char T;};
    template<>class type_num2cpp_type<NPY_SHORT>{public:typedef short T;};
    template<>class type_num2cpp_type<NPY_USHORT>{public:typedef unsigned short T;};
    template<>class type_num2cpp_type<NPY_INT>{public:typedef int T;};
    template<>class type_num2cpp_type<NPY_UINT>{public:typedef unsigned int T;};
    template<>class type_num2cpp_type<NPY_LONG>{public:typedef long T;};
    template<>class type_num2cpp_type<NPY_ULONG>{public:typedef unsigned long T;};
    template<>class type_num2cpp_type<NPY_LONGLONG>{public:typedef long long T;};
    template<>class type_num2cpp_type<NPY_ULONGLONG>{public:typedef unsigned long long T;};
    template<>class type_num2cpp_type<NPY_FLOAT>{public:typedef float T;};
    template<>class type_num2cpp_type<NPY_DOUBLE>{public:typedef double T;};
    template<>class type_num2cpp_type<NPY_LONGDOUBLE>{public:typedef long double T;};
    template<>class type_num2cpp_type<NPY_STRING>{public:typedef std::string T;};
    
    template<class T>class cpp_type2type_num{public:static constexpr int type_num(){return -1;};};
    template<>class cpp_type2type_num<bool>{public:static constexpr int type_num(){return NPY_BOOL;};};
    template<>class cpp_type2type_num<char>{public:static constexpr int type_num(){return NPY_BYTE;};};
    template<>class cpp_type2type_num<unsigned char>{public:static constexpr int type_num(){return NPY_UBYTE;};};
    template<>class cpp_type2type_num<short>{public:static constexpr int type_num(){return NPY_SHORT;};};
    template<>class cpp_type2type_num<unsigned short>{public:static constexpr int type_num(){return NPY_USHORT;};};
    template<>class cpp_type2type_num<int>{public:static constexpr int type_num(){return NPY_INT;};};
    template<>class cpp_type2type_num<unsigned int>{public:static constexpr int type_num(){return NPY_UINT;};};
    template<>class cpp_type2type_num<long>{public:static constexpr int type_num(){return NPY_LONG;};};
    template<>class cpp_type2type_num<unsigned long>{public:static constexpr int type_num(){return NPY_ULONG;};};
    template<>class cpp_type2type_num<long long>{public:static constexpr int type_num(){return NPY_LONGLONG;};};
    template<>class cpp_type2type_num<unsigned long long>{public:static constexpr int type_num(){return NPY_ULONGLONG;};};
    template<>class cpp_type2type_num<float>{public:static constexpr int type_num(){return NPY_FLOAT;};};
    template<>class cpp_type2type_num<double>{public:static constexpr int type_num(){return NPY_DOUBLE;};};
    template<>class cpp_type2type_num<long double>{public:static constexpr int type_num(){return NPY_LONGDOUBLE;};};
    template<>class cpp_type2type_num<std::string>{public:static constexpr int type_num(){return NPY_STRING;};};

    
    
    template<class T>class symm{};
    
    template<class T>class type_attr{};
    template<>class type_attr<bool>    
    {public: static const std::string name(); static const bool zero;};
    template<>class type_attr<char>
    {public: static const std::string name(); static const char zero;};
    template<>class type_attr<unsigned char>
    {public: static const std::string name(); static const unsigned char zero;};
    template<>class type_attr<short>
    {public: static const std::string name(); static const short zero;};
    template<>class type_attr<unsigned short>
    {public: static const std::string name(); static const unsigned short zero;};
    template<>class type_attr<int>
    {public: static const std::string name(); static const int zero;};
    template<>class type_attr<unsigned int>
    {public: static const std::string name(); static const unsigned int zero;};
    template<>class type_attr<long>
    {public: static const std::string name(); static const long zero;};
    template<>class type_attr<unsigned long>
    {public: static const std::string name(); static const unsigned long zero;};
    template<>class type_attr<long long>
    {public: static const std::string name(); static const long long zero;};
    template<>class type_attr<unsigned long long>
    {public: static const std::string name(); static const unsigned long long zero;};
    template<>class type_attr<float>
    {public: static const std::string name(); static const float zero;};
    template<>class type_attr<double>
    {public: static const std::string name(); static const double zero;};
    template<>class type_attr<long double>
    {public: static const std::string name(); static const long double zero;};
    template<>class type_attr<std::string>
    {public: static const std::string name(); static const std::string zero;};
    
    
    template<class T> class remove_matrix_attr{public: typedef T type;};
    template<class T> class remove_matrix_attr<symm<T>>{public: typedef T type;};
}
using namespace MAPP_NS;
/*------------------------------------------------------------------------------------------------------------------------------------
 
 ------------------------------------------------------------------------------------------------------------------------------------*/
namespace MAPP_NS
{
    template<class T>
    class py_var
    {
        friend class py_var<T*>;
    private:
        py_var(size_t*,void*&);
        py_var(int,long*,PyObject**&);
        py_var(py_var<T>&&);
    public:
        typedef T T_BASE;
        static constexpr int get_rank(){return 0;}
        T val;
        py_var(PyObject*);
        py_var();
        ~py_var();
        
        py_var(const py_var<T>&);
        operator T() const{return val;}
        py_var<T>& operator=(const py_var<T>&);
        py_var<T>& operator=(py_var<T>&&);
        bool operator==(const py_var<T>&);
        static std::string err_msg();
        const static std::string type_name();
    };
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<class T>
py_var<T>::py_var():
val(type_attr<T>::zero)
{
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<class T>
py_var<T>::py_var(PyObject* op):
val(type_attr<T>::zero)
{
    int type_num=cpp_type2type_num<T>::type_num();
    PyArray_Descr* descr=PyArray_DescrFromType(type_num);
    PyObject* __arr__=PyArray_FromAny(op,descr,0,0,0,NULL);
    
    if(__arr__)
    {
        PyArrayObject* arr=(PyArrayObject*)__arr__;
        if(std::is_same<T,bool>::value)
        {
            if(*static_cast<unsigned char*>(PyArray_DATA(arr))==NPY_TRUE)
                val=true;
            else
                val=false;
        }
        else
            val=*static_cast<T*>(PyArray_DATA(arr));
        
        Py_DECREF(__arr__);
        return;
    }
    
    throw 1;
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<>
inline py_var<std::string>::py_var(PyObject* op):
val(type_attr<std::string>::zero)
{
    int type_num=cpp_type2type_num<std::string>::type_num();
    PyArray_Descr* descr=PyArray_DescrFromType(type_num);
    PyObject* __arr__=PyArray_FromAny(op,descr,0,0,0,NULL);
    
    if(__arr__)
    {
        PyArrayObject* arr=(PyArrayObject*)__arr__;
        char* data=static_cast<char*>(PyArray_DATA(arr));
        val=std::string(data,PyArray_DESCR(arr)->elsize);
        Py_DECREF(__arr__);
        return;
    }
    
    throw 1;
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<class T>
py_var<T>::py_var(int depth,long* sz,PyObject**& objs)
{
    this->~py_var<T>();
    try
    {
        new (this) py_var<T>(*objs);
    }
    catch(int)
    {
        throw 1;
    }
    objs++;
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<class T>
py_var<T>::py_var(size_t*,void*& data):
val(type_attr<T>::zero)
{
    if(std::is_same<T,bool>::value)
    {
        unsigned char* d=static_cast<unsigned char*>(data);
        if(*d==NPY_TRUE)
            val=true;
        else
            val=false;
        data=d+1;
    }
    else
    {
        T* d=(T*)data;
        val=*d;
        data=d+1;
    }
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<class T>
py_var<T>::~py_var()
{
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<class T>
py_var<T>::py_var(const py_var<T>& r)
{
    val=r.val;
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<class T>
py_var<T>& py_var<T>::operator=(const py_var<T>& r)
{
    val=r.val;
    return *this;
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<class T>
py_var<T>::py_var(py_var<T>&& r)
{
    val=r.val;
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<class T>
py_var<T>& py_var<T>::operator=(py_var<T>&& r)
{
    val=r.val;
    return *this;
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<class T>
bool py_var<T>::operator==(const py_var<T>& r)
{
    return (this->val==r.val);
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<class T>
std::string py_var<T>::err_msg()
{
    return "fialure to deduce C++ type <"+type_attr<T_BASE>::name+">";
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<class T>
const std::string py_var<T>::type_name()
{
    return type_attr<T_BASE>::name();
}
/*------------------------------------------------------------------------------------------------------------------------------------
 
 ------------------------------------------------------------------------------------------------------------------------------------*/
namespace MAPP_NS
{
    template<class T>
    class py_var<T*>
    {
        friend class py_var<T**>;
    private:
        py_var(size_t*,void*&);
        py_var(int,long*,PyObject**&);
        py_var(py_var<T*>&&);
    public:
        typedef typename py_var<T>::T_BASE T_BASE;
        static constexpr int get_rank(){ return 1+py_var<T>::get_rank();}
        py_var<T>* vars;
        size_t size;
        py_var(PyObject*);
        py_var(PyObject**);
        py_var();
        ~py_var();
        
        py_var(const py_var<T*>&);
        py_var<T*>& operator=(const py_var<T*>&);
        py_var<T*>& operator=(py_var<T*>&&);
        bool operator==(const py_var<T*>&);
        py_var<T>& operator[](const size_t);
        
        bool is_square();
        bool is_square(size_t);
        bool is_triangular();
        bool is_triangular(size_t);
        bool is_symmetric();
        bool is_symmetric(size_t);
        bool is_lower_triangular();
        bool is_lower_triangular(size_t);
        bool is_upper_triangular();
        bool is_upper_triangular(size_t);
        
        void triangular_to_symmetric();
        void triangular_to_lower_triangular();
        void triangular_to_upper_triangular();
        
        static std::string err_msg();
        static const std::string type_name();
        
    };
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<class T>
py_var<T*>::py_var():
size(0),
vars(NULL)
{
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<class T>
py_var<T*>::py_var(PyObject* op):
size(0),
vars(NULL)
{
    int type_num=cpp_type2type_num<T_BASE>::type_num();
    int rank=get_rank();
    PyObject* __arr__;
    PyArray_Descr* descr;
    
    if(std::is_same<T_BASE,char>::value)
    {
        descr=PyArray_DescrFromType(NPY_STRING);
        __arr__=PyArray_FromAny(op,descr,rank-1,rank-1,0,NULL);
        if(__arr__)
        {
            PyArrayObject* arr=(PyArrayObject*)__arr__;
            size_t sz[rank];
            npy_intp* dims=PyArray_DIMS(arr);
            for(int i=0;i<rank-1;i++)
                sz[i]=dims[i];
            sz[rank-1]=PyArray_DESCR(arr)->elsize;
            void* data=PyArray_DATA(arr);
            this->~py_var<T*>();
            new (this) py_var<T*>(sz,data);
            Py_DECREF(__arr__);
            return;
        }
        PyErr_Clear();
    }
    
    descr=PyArray_DescrFromType(type_num);
    __arr__=PyArray_FromAny(op,descr,rank,rank,NPY_ARRAY_CARRAY_RO,NULL);
    if(__arr__)
    {
        PyArrayObject* arr=(PyArrayObject*)__arr__;
        size_t sz[rank];
        npy_intp* dims=PyArray_DIMS(arr);
        for(int i=0;i<rank;i++)
            sz[i]=dims[i];
        void* data=PyArray_DATA(arr);
        this->~py_var<T*>();
        new (this) py_var<T*>(sz,data);
        Py_DECREF(__arr__);
        return;
    }
    PyErr_Clear();
    
    if(rank>1)
    {
        descr=PyArray_DescrFromType(NPY_OBJECT);
        __arr__=PyArray_FromAny(op,descr,1,rank-1,NPY_ARRAY_CARRAY_RO,NULL);
        PyErr_Clear();
        
        if(__arr__)
        {
            PyArrayObject* arr=(PyArrayObject*)__arr__;
            int depth=PyArray_NDIM(arr);
            PyObject** objs=static_cast<PyObject**>(PyArray_DATA(arr));
            long* sz=PyArray_DIMS(arr);
            this->~py_var<T*>();
            try
            {
                new (this) py_var<T*>(depth,sz,objs);
            }
            catch(int)
            {
                Py_DECREF(__arr__);
                throw 1;
            }
            Py_DECREF(__arr__);
            return;
        }
        PyErr_Clear();
    }
    
    throw 1;
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<class T>
py_var<T*>::py_var(PyObject** op_ptr):
size(0),
vars(NULL)
{
    size=1;
    vars=new py_var<T>[1];
    try
    {
        vars->~py_var<T>();
        new (vars) py_var<T>(*op_ptr);
    }
    catch(int)
    {
        delete [] vars;
        throw 1;
    }
    
    size_t dim=0;
    while(dim*(dim+1)/2<vars[0].size)
        dim++;
    
    if(dim*(dim+1)/2!=vars[0].size)
    {
        delete [] vars;
        throw 1;
    }
    size=dim;
    
    typedef typename std::remove_pointer<T>::type S;
    py_var<S*>* __vars=NULL;
    if(size) __vars=new py_var<S*>[size];
    
    
    for(size_t i=0;i<size;i++)
    {
        __vars[i].size=size;
        __vars[i].vars=new py_var<S>[size];
    }
    size_t it=0;
    
    for(size_t i=0;i<size;i++)
    {
        __vars[i][i]=std::move(vars[0][it++]);
    }
    
    int dir=-1;
    for(size_t i=size-1;i>0;i--)
    {
        if(dir==-1)
            for(size_t j=i;j>0;j--)
            {
                __vars[i][j-1]=std::move(vars[0][it++]);
                __vars[j-1][i]=__vars[i][j-1];
            }
        else
            for(size_t j=0;j<i;j++)
            {
                __vars[i][j]=std::move(vars[0][it++]);
                __vars[j][i]=__vars[i][j];
            }
        
        dir*=-1;
    }
    
    delete [] vars;
    vars=__vars;
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<class T>
py_var<T*>::py_var(size_t* sz,void*& data):
size(*sz),
vars(NULL)
{
    if(std::is_same<T,char>::value)
    {
        char* d=static_cast<char*>(data);
        while(d[size-1]=='\0') size--;
        size++;
        if(size) vars=new py_var<T>[size];
        for(size_t i=0;i<size-1;i++)
        {
            (vars+i)-> ~py_var<T>();
            new (vars+i) py_var<T>(sz+1,data);
        }
        
        char null_char='\0';
        void* _data=&null_char;
        (vars+size-1)-> ~py_var<T>();
        new (vars+size-1) py_var<T>(sz+1,_data);
    }
    else
    {
        if(size) vars=new py_var<T>[size];
        for(size_t i=0;i<size;i++)
        {
            (vars+i)-> ~py_var<T>();
            new (vars+i) py_var<T>(sz+1,data);
        }
    }
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<class T>
py_var<T*>::py_var(int depth,long* sz,PyObject**& objs):
size(0),
vars(NULL)
{
    if(depth==0)
    {
        this->~py_var<T*>();
        try
        {
            new (this) py_var<T*>(*objs);
        }
        catch(int)
        {
            delete [] vars;
            throw 1;
        }
        objs++;
        return;
    }
    else
    {
        size=*sz;
        if(size) vars=new py_var<T>[*sz];
        for(size_t i=0;i<size;i++)
        {
            try
            {
                //vars[i]=std::move(py_var<T>(depth-1,sz+1,objs));
                (vars+i)->~py_var<T>();
                new (vars+i) py_var<T>(depth-1,sz+1,objs);
            }
            catch(int)
            {
                delete [] vars;
                throw 1;
            }
        }
    }
    
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<class T>
py_var<T*>::~py_var()
{
    delete [] vars;
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<class T>
py_var<T*>::py_var(const py_var<T*>& r):
vars(NULL)
{
    size=r.size;
    if(size) vars=new py_var<T>[size];
    for(size_t i=0;i<size;i++)
        vars[i]=r.vars[i];
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<class T>
py_var<T*>& py_var<T*>::operator=(const py_var<T*>& r)
{
    delete [] this->vars;
    this->vars=NULL;
    this->size=r.size;
    if(this->size) this->vars=new py_var<T>[size];
    for(size_t i=0;i<this->size;i++)
        this->vars[i]=r.vars[i];
    return *this;
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<class T>
py_var<T*>::py_var(py_var<T*>&& r)
{
    size=r.size;
    r.size=0;
    vars=r.vars;
    r.vars=NULL;
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<class T>
py_var<T*>& py_var<T*>::operator=(py_var<T*>&& r)
{
    delete [] this->vars;
    this->vars=NULL;
    this->size=r.size;
    r.size=0;
    this->vars=r.vars;
    r.vars=NULL;
    return *this;
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<class T>
bool py_var<T*>::operator==(const py_var<T*>& r)
{
    if(this->size!=r.size) return false;
    for(size_t i=0;i<size;i++)
        if(!(this->vars[i]==r.vars[i])) return false;
    return true;
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<class T>
py_var<T>& py_var<T*>::operator[](const size_t i)
{
    return vars[i];
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<class T>
bool py_var<T*>::is_square()
{
    if(get_rank()<2) return false;
    for(size_t i=0;i<size;i++)
        if(vars[i].size!=size) return false;
    return true;
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<class T>
bool py_var<T*>::is_square(size_t __size)
{
    if(size!=__size) return false;
    return is_square();
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<class T>
bool py_var<T*>::is_triangular()
{
    if(get_rank()<2) return false;
    for(size_t i=0;i<size;i++)
        if(vars[i].size!=i+1) return false;
    return true;
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<class T>
bool py_var<T*>::is_triangular(size_t __size)
{
    if(size!=__size) return false;
    return is_triangular();
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<class T>
bool py_var<T*>::is_symmetric()
{
    if(is_triangular())
    {
        triangular_to_symmetric();
        return true;
    }
    
    if(!is_square()) return false;
    for(size_t i=0;i<size;i++)
        for(size_t j=0;j<i;j++)
            if(!(vars[i][j]==vars[j][i]))
                return false;
    return true;
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<class T>
bool py_var<T*>::is_symmetric(size_t __size)
{
    if(size!=__size) return false;
    return is_symmetric();
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<class T>
bool py_var<T*>::is_lower_triangular()
{
    if(is_triangular())
    {
        triangular_to_lower_triangular();
        return true;
    }
    
    if(!is_square()) return false;
    typedef typename std::remove_pointer<T>::type S;
    py_var<S> null_var;
    for(size_t i=0;i<size;i++)
        for(size_t j=0;j<i;j++)
            if(!(vars[j][i]==null_var))
                return false;
    return true;
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<class T>
bool py_var<T*>::is_lower_triangular(size_t __size)
{
    if(size!=__size) return false;
    return is_lower_triangular();
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<class T>
bool py_var<T*>::is_upper_triangular()
{
    if(is_triangular())
    {
        triangular_to_upper_triangular();
        return true;
    }
    
    if(!is_square()) return false;
    typedef typename std::remove_pointer<T>::type S;
    py_var<S> null_var;
    for(size_t i=0;i<size;i++)
        for(size_t j=0;j<i;j++)
            if(!(vars[i][j]==null_var))
                return false;
    return true;
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<class T>
bool py_var<T*>::is_upper_triangular(size_t __size)
{
    if(size!=__size) return false;
    return is_upper_triangular();
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<class T>
void py_var<T*>::triangular_to_symmetric()
{
    typedef typename std::remove_pointer<T>::type S;
    py_var<S*>* __vars=NULL;
    if(size) __vars=new py_var<S*>[size];
    
    
    for(size_t i=0;i<size;i++)
    {
        __vars[i].size=size;
        __vars[i].vars=new py_var<S>[size];
        for(size_t j=0;j<=i;j++)
            __vars[i][j]=std::move(vars[i][j]);
    }
    
    for(size_t i=0;i<size;i++)
        for(size_t j=i+1;j<size;j++)
            __vars[i][j]=__vars[j][i];
    
    delete [] vars;
    vars=__vars;
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<class T>
void py_var<T*>::triangular_to_lower_triangular()
{
    typedef typename std::remove_pointer<T>::type S;
    py_var<S*>* __vars=NULL;
    if(size) __vars=new py_var<S*>[size];
    
    
    for(size_t i=0;i<size;i++)
    {
        __vars[i].size=size;
        __vars[i].vars=new py_var<S>[size];
        for(size_t j=0;j<=i;j++)
            __vars[i][j]=std::move(vars[i][j]);
    }
    
    delete [] vars;
    vars=__vars;
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<class T>
void py_var<T*>::triangular_to_upper_triangular()
{
    typedef typename std::remove_pointer<T>::type S;
    py_var<S*>* __vars=NULL;
    if(size) __vars=new py_var<S*>[size];
    
    
    for(size_t i=0;i<size;i++)
    {
        __vars[i].size=size;
        __vars[i].vars=new py_var<S>[size];
        for(size_t j=0;j<=i;j++)
            __vars[j][i]=std::move(vars[i][j]);
    }
    
    delete [] vars;
    vars=__vars;
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<class T>
std::string py_var<T*>::err_msg()
{
    constexpr int rank=get_rank();
    std::string err="fialure to deduce C++ type <"+type_attr<T_BASE>::name+"[]";
    
    for(int i=1;i<rank;i++) err+="[]";
    err+=">";
    return err;
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<class T>
const std::string py_var<T*>::type_name()
{
    constexpr int rank=get_rank();
    std::string err=type_attr<T_BASE>::name()+"*";
    for(int i=1;i<rank;i++) err+="*";
    return err;
}
/*------------------------------------------------------------------------------------------------------------------------------------
 
 ------------------------------------------------------------------------------------------------------------------------------------*/
namespace MAPP_NS
{
    template<class> class var;
    class Var
    {
    private:
        bool pushed;
    protected:
        virtual const void* get_ptr()=0;
        bool root;
        size_t** size_list();
        void size_list(size_t*);
        void size_list(size_t**);
    public:
        static Var** var_stack;
        static size_t stack_size;
        static size_t stack_capacity;
        static size_t stack_grow;
        static void push(Var*);
        static void pop(Var*);
        static void replace(Var*,Var*);
        template<class T>
        static Var* find_var(T&);
        
        template<class T,class...Ts>
        static void assign_dynamic_size(T**,size_t&,Ts&&...);
        template<class T,class...Ts>
        static void assign_dynamic_size(T**,size_t&&,Ts&&...);
        template<class T>
        static void assign_dynamic_size(T**){}
                                        
        bool obj_set;
        const int rank;
        const size_t hash_code;
        size_t size;
        std::string name;
        
        
        Var(int,size_t,const char*);
        Var(int,size_t,std::string&&);
        Var(int,size_t);
        virtual ~Var();
        Var(const Var&);
        Var(Var&&);
        virtual Var* clone()=0;
        
        virtual void set(PyObject*)=0;
        virtual bool set_nothrow(PyObject*)=0;
        virtual PyObject* get()=0;
        virtual Var& operator [] (const size_t)=0;
        
        size_t* is_array(const int);
        size_t is_square();
        size_t is_triangular();
        
        bool is_set();
        bool belongs_to(Var&);
        virtual bool operator == (Var&);
        virtual bool operator != (Var&);
        virtual bool operator > (Var&);
        virtual bool operator >= (Var&);
        virtual bool operator <= (Var&);
        virtual bool operator < (Var&);

        
        
        
        template<class Logics>
        void log_chek(Logics* logics);
    };
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<class Logics>
void Var::log_chek(Logics* logics)
{
    std::string err_msg=(*logics)(this);
    if(!err_msg.empty()) throw err_msg;
    if(rank)
        for(size_t i=0;i<size;i++)
            (*this)[i].log_chek(logics+1);
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<class T>
Var* Var::find_var(T& val)
{
    if(!stack_size) return NULL;
    int var_rank=var<T>::get_rank();
    const void* ptr=&val;
    size_t ivar=stack_size-1;
    for(;ivar!=0 && (var_stack[ivar]->get_ptr()!=ptr || var_stack[ivar]->rank!=var_rank) ;ivar--);
    if(ivar) return var_stack[ivar];
    else if(var_stack[0]->get_ptr()==ptr && var_stack[ivar]->rank==var_rank)
        return var_stack[0];
    return NULL;
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<class T,class...Ts>
void Var::assign_dynamic_size(T** __dsizes__,size_t& sz,Ts&&...szs)
{
    *__dsizes__=reinterpret_cast<T*>(find_var(sz));
    assign_dynamic_size(__dsizes__+1,szs...);
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<class T,class...Ts>
void Var::assign_dynamic_size(T** __dsizes__,size_t&& sz,Ts&&...szs)
{
    *__dsizes__=NULL;
    assign_dynamic_size(__dsizes__+1,szs...);
}
/*------------------------------------------------------------------------------------------------------------------------------------
 
 ------------------------------------------------------------------------------------------------------------------------------------*/
namespace MAPP_NS
{
    template<class T>
    class var: public Var
    {
    private:
        T val;
    protected:
        const void* get_ptr();
    public:
        typedef typename std::remove_const<T>::type T_BASE;
        typedef T_BASE T_EQUIV;
        static constexpr int get_rank();
        static size_t base_hash_code();
        static void allocate(void**,T&,size_t*);
        static void deallocate(T&);
        static void accum_size(py_var<T_EQUIV>&,size_t*){};
        static bool is_size_compatible_nothrow(py_var<T_EQUIV>&,var<size_t>**);
        static void is_size_compatible(py_var<T_EQUIV>&,const std::string&,var<size_t>**);
        static std::string type_name(var<size_t>**);
        static PyObject* build(T&,size_t**);
        
        T* ptr;
        
        var(T&,const char*);
        
        var(T&,std::string&&,py_var<T_EQUIV>&,void**);
        var();
        ~var();
        var(T&&);
        var(const var<T>&);
        var(var<T>&&);
        var<T>& operator=(const var<T>&);
        var<T>& operator=(var<T>&&);
        var<T>& operator=(const T&);
        Var* clone();
        
        operator T() const;
        var<T>& operator [] (const size_t);
        void set(PyObject*);
        bool set_nothrow(PyObject*);
        PyObject* get();
        
        std::string type_name(){return type_attr<T_BASE>::name;}
        
        bool operator==(Var&);
        bool operator!=(Var&);
        bool operator> (Var&);
        bool operator< (Var&);
        bool operator>=(Var&);
        bool operator<=(Var&);
    };
    
}
/*--------------------------------------------
 
 --------------------------------------------*/
template <class T>
constexpr int var<T>::get_rank()
{
    return 0;
}
/*--------------------------------------------
 
 --------------------------------------------*/
template <class T>
size_t var<T>::base_hash_code()
{
    return typeid(T_BASE).hash_code();
}
/*--------------------------------------------
 
 --------------------------------------------*/
template <class T>
void var<T>::allocate(void**,T&,size_t*)
{}
/*--------------------------------------------
 
 --------------------------------------------*/
template <class T>
void var<T>::deallocate(T&)
{}
/*--------------------------------------------
 
 --------------------------------------------*/
template <class T>
bool var<T>::is_size_compatible_nothrow(py_var<T_EQUIV>&,var<size_t>**)
{
    return true;
}
/*--------------------------------------------
 
 --------------------------------------------*/
template <class T>
void var<T>::is_size_compatible(py_var<T_EQUIV>&,const std::string&,var<size_t>**)
{
}
/*--------------------------------------------
 
 --------------------------------------------*/
template <class T>
std::string var<T>::type_name(var<size_t>**)
{
    return  std::string();
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<class T>
PyObject* var<T>::build(T& v,size_t**)
{
    
    if(std::is_integral<T_BASE>::value)
    {
        long i=0;
        if(std::is_same<T_BASE,bool>::value)
        {
            if(v) i=1;
            return PyBool_FromLong(i);
        }
        else
        {
            i=static_cast<long>(v);
            return PyInt_FromLong(i);
        }
    }
    else
    {
        
        double i=static_cast<double>(v);
        return PyFloat_FromDouble(i);
    }
    Py_RETURN_NONE;
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<>
inline PyObject* var<std::string>::build(std::string& v,size_t**)
{return PyString_FromString(v.c_str());}
/*--------------------------------------------
 constructor:
 usage: direct
 --------------------------------------------*/
template<class T>
var<T>::var(T& v,const char* name_):
val(type_attr<T>::zero),
ptr(&v),
Var(get_rank(),base_hash_code(),name_)
{
    size=1;
}
/*--------------------------------------------
 constructor:
 usage: indirect by same template class with
 one upper rank
 --------------------------------------------*/
template<class T>
var<T>::var(T& v,std::string&& name_,py_var<T_EQUIV>& pv,void**):
val(type_attr<T>::zero),
ptr(&v),
Var(get_rank(),base_hash_code(),std::move(name_))
{
    size=1;
    v=pv;
}
/*--------------------------------------------
 default constructor:
 usage: indirect by same template class with
 one upper rank
 --------------------------------------------*/
template<class T>
var<T>::var():
ptr(NULL),
Var(get_rank(),base_hash_code()),
val(type_attr<T>::zero)
{
    size=0;
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<class T>
var<T>::~var()
{
}
/*--------------------------------------------
 constructor:
 usage: indirect by Logics class
 it is only used for all literals except for
 string literals
 --------------------------------------------*/
template<class T>
var<T>::var(T&& v):
val(v),
ptr(&val),
Var(get_rank(),base_hash_code())
{
    size=1;
    name=std::to_string(val);
}
/*--------------------------------------------
 copy constructor:
 --------------------------------------------*/
template<class T>
var<T>::var(const var<T>& r):
Var(static_cast<const Var&>(r)),
val(r.val),
ptr(r.ptr==&(r.val)?&val:r.ptr)
{
    
}
/*--------------------------------------------
 move constructor:
 --------------------------------------------*/
template<class T>
var<T>::var(var<T>&& r):
Var(static_cast<Var&&>(r)),
val(r.val),
ptr(r.ptr==&(r.val)?&val:r.ptr)
{
    r.ptr=NULL;
}
/*--------------------------------------------
 copy assignment
 --------------------------------------------*/
template<class T>
var<T>& var<T>::operator=(const var<T>& r)
{
    if(this==&r) return *this;
    this->~var<T>();
    new (this) var<T>(r);
    return *this;
}
/*--------------------------------------------
 move assignment
 --------------------------------------------*/
template<class T>
var<T>& var<T>::operator=(var<T>&& r)
{
    if(this==&r) return *this;
    this->~var<T>();
    new (this) var<T>(std::move(r));
    return *this;
}
/*--------------------------------------------
 copy assignment from variable
 --------------------------------------------*/
template<class T>
var<T>& var<T>::operator = (const T& v)
{
    typedef typename std::remove_const<T>::type S;
    *const_cast<S*>(ptr)=v;
    return *this;
}
/*--------------------------------------------
 clone used by Logics class
 --------------------------------------------*/
template<class T>
Var* var<T>::clone()
{
    return new var<T>(*this);
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<class T>
const void* var<T>::get_ptr()
{
    return ptr;
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<class T>
var<T>::operator T() const
{
    return *ptr;
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<class T>
var<T>& var<T>::operator[](const size_t)
{
    return *this;
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<class T>
void var<T>::set(PyObject* py_obj)
{
    py_var<T_EQUIV> pv;
    try
    {
        pv=py_var<T_EQUIV>(py_obj);
    }
    catch(int)
    {
        throw "fialure to deduce C++ type <"+py_var<T_EQUIV>::type_name()+">";
    }
    *ptr=pv;
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<class T>
bool var<T>::set_nothrow(PyObject* py_obj)
{
    py_var<T_EQUIV> pv;
    try
    {
        pv=py_var<T_EQUIV>(py_obj);
    }
    catch(int)
    {
        return false;
    }
    *ptr=pv;
    return true;
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<class T>
PyObject* var<T>::get()
{
    if(std::is_integral<T_BASE>::value)
    {
        long i=0;
        if(std::is_same<T_BASE,bool>::value)
        {
            if(*ptr) i=1;
            return PyBool_FromLong(i);
        }
        else
        {
            i=static_cast<long>(*ptr);
            return PyInt_FromLong(i);
        }
    }
    else
    {
        double i=static_cast<double>(*ptr);
        return PyFloat_FromDouble(i);
    }
    Py_RETURN_NONE;
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<>
inline PyObject* var<std::string>::get()
{return PyString_FromString(ptr->c_str());}
/*--------------------------------------------
 
 --------------------------------------------*/
template<class T>
bool var<T>::operator == (Var& r)
{
    var<T>* var_ptr0=dynamic_cast<var<T>*>(&r);
    if(var_ptr0)
        return (*ptr==*(var_ptr0->ptr));
    
    if(std::is_const<T>::value)
    {
        typedef typename std::remove_const<T>::type S;
        var<S>* var_ptr0=dynamic_cast<var<S>*>(&r);
        if(var_ptr0)
            return (*ptr==*(var_ptr0->ptr));
    }
    else
    {
        var<const T>* var_ptr0=dynamic_cast<var<const T>*>(&r);
        if(var_ptr0)
            return (*ptr==*(var_ptr0->ptr));
    }
    return false;
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<class T>
bool var<T>::operator!=(Var& r)
{
    var<T>* var_ptr0=dynamic_cast<var<T>*>(&r);
    if(var_ptr0)
        return (*ptr!=*(var_ptr0->ptr));
    
    if(std::is_const<T>::value)
    {
        typedef typename std::remove_const<T>::type S;
        var<S>* var_ptr0=dynamic_cast<var<S>*>(&r);
        if(var_ptr0)
            return (*ptr!=*(var_ptr0->ptr));
    }
    else
    {
        var<const T>* var_ptr0=dynamic_cast<var<const T>*>(&r);
        if(var_ptr0)
            return (*ptr!=*(var_ptr0->ptr));
    }
    return false;
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<class T>
bool var<T>::operator>(Var& r)
{
    var<T>* var_ptr0=dynamic_cast<var<T>*>(&r);
    if(var_ptr0)
        return (*ptr>*(var_ptr0->ptr));
    
    if(std::is_const<T>::value)
    {
        typedef typename std::remove_const<T>::type S;
        var<S>* var_ptr0=dynamic_cast<var<S>*>(&r);
        if(var_ptr0)
            return (*ptr>*(var_ptr0->ptr));
    }
    else
    {
        var<const T>* var_ptr0=dynamic_cast<var<const T>*>(&r);
        if(var_ptr0)
            return (*ptr>*(var_ptr0->ptr));
    }
    return false;
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<class T>
bool var<T>::operator<(Var& r)
{
    var<T>* var_ptr0=dynamic_cast<var<T>*>(&r);
    if(var_ptr0)
        return (*ptr<*(var_ptr0->ptr));
    
    if(std::is_const<T>::value)
    {
        typedef typename std::remove_const<T>::type S;
        var<S>* var_ptr0=dynamic_cast<var<S>*>(&r);
        if(var_ptr0)
            return (*ptr<*(var_ptr0->ptr));
    }
    else
    {
        var<const T>* var_ptr0=dynamic_cast<var<const T>*>(&r);
        if(var_ptr0)
            return (*ptr<*(var_ptr0->ptr));
    }
    return false;
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<class T>
bool var<T>::operator>=(Var& r)
{
    var<T>* var_ptr0=dynamic_cast<var<T>*>(&r);
    if(var_ptr0)
        return (*ptr>=*(var_ptr0->ptr));
    
    if(std::is_const<T>::value)
    {
        typedef typename std::remove_const<T>::type S;
        var<S>* var_ptr0=dynamic_cast<var<S>*>(&r);
        if(var_ptr0)
            return (*ptr>=*(var_ptr0->ptr));
    }
    else
    {
        var<const T>* var_ptr0=dynamic_cast<var<const T>*>(&r);
        if(var_ptr0)
            return (*ptr>=*(var_ptr0->ptr));
    }
    return false;
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<class T>
bool var<T>::operator<=(Var& r)
{
    var<T>* var_ptr0=dynamic_cast<var<T>*>(&r);
    if(var_ptr0)
        return (*ptr<=*(var_ptr0->ptr));
    
    if(std::is_const<T>::value)
    {
        typedef typename std::remove_const<T>::type S;
        var<S>* var_ptr0=dynamic_cast<var<S>*>(&r);
        if(var_ptr0)
            return (*ptr<=*(var_ptr0->ptr));
    }
    else
    {
        var<const T>* var_ptr0=dynamic_cast<var<const T>*>(&r);
        if(var_ptr0)
            return (*ptr<=*(var_ptr0->ptr));
    }
    return false;
}
/*------------------------------------------------------------------------------------------------------------------------------------
 
 ------------------------------------------------------------------------------------------------------------------------------------*/
namespace MAPP_NS
{
    template<class T,size_t N>
    class var<T[N]>:public Var
    {
    private:
    protected:
        var<size_t>** __dsizes__;
        const void* get_ptr();
    public:
        typedef typename var<T>::T_BASE T_BASE;
        typedef typename std::add_pointer<typename var<T>::T_EQUIV>::type T_EQUIV;
        static constexpr int get_rank();
        static size_t base_hash_code();
        static void allocate(void**,T(&)[N],size_t*);
        static void deallocate(T(&)[N]);
        static void accum_size(py_var<T_EQUIV>&,size_t*);
        static bool is_size_compatible_nothrow(py_var<T_EQUIV>&,var<size_t>**);
        static void is_size_compatible(py_var<T_EQUIV>&,const std::string&,var<size_t>**);
        static std::string type_name(var<size_t>**);
        static PyObject* build(T(&)[N],size_t**);
        
        var<T> vars[N];
        T (*ptr)[N];
        
        var(T(&)[N],const char*);
        template<class... Ts>
        var(T(&)[N],const char*,Ts&&...);
        var(T(&)[N],std::string&&,py_var<T_EQUIV>&,void**);
        var();
        virtual ~var();
        var(T(&)[N]);
        var(const var<T[N]>&);
        var(var<T[N]>&&);
        var<T[N]>& operator = (const var<T[N]>&);
        var<T[N]>& operator = (var<T[N]>&&);
        var<T[N]>& operator = (const T (&)[N]);
        Var* clone();
        
        var<T>& operator [] (const size_t);
        void set(PyObject*);
        bool set_nothrow(PyObject*);
        void set(py_var<T_EQUIV>&);
        PyObject* get();
        std::string type_name()
        {return type_attr<T_BASE>::name()+" "+type_name(__dsizes__);}
        template<class ... Ts>
        void dynamic_size(Ts&&...);
    };
}
/*--------------------------------------------
 
 --------------------------------------------*/
template <class T,size_t N>
constexpr int var<T[N]>::get_rank()
{
    return 1+var<T>::get_rank();
}
/*--------------------------------------------
 
 --------------------------------------------*/
template <class T,size_t N>
size_t var<T[N]>::base_hash_code()
{
    return typeid(T_BASE).hash_code();
}
/*--------------------------------------------
 
 --------------------------------------------*/
template <class T,size_t N>
void var<T[N]>::allocate(void** ptr,T (&v)[N],size_t* sz)
{
    *ptr=NULL;
    if(get_rank()==1) return;
    var<T>::allocate(ptr+1,v[0],sz+1);
}
/*--------------------------------------------
 
 --------------------------------------------*/
template <class T,size_t N>
void var<T[N]>::deallocate(T (&v)[N])
{
    var<T>::deallocate(v[0]);
}
/*--------------------------------------------
 
 --------------------------------------------*/
template <class T,size_t N>
void var<T[N]>::accum_size(py_var<T_EQUIV>& pv,size_t* sz)
{
    *sz+=N;
    if(get_rank()==1) return;
    for(size_t i=0;i<N;i++)
        var<T>::accum_size(pv[i],sz+1);
}
/*--------------------------------------------
 
 --------------------------------------------*/
template <class T,size_t N>
bool var<T[N]>::is_size_compatible_nothrow(py_var<T_EQUIV>& pv,var<size_t>** sz)
{
    if(N!=pv.size)
        return false;
    
    if(sz) sz+=1;

    for(size_t i=0;i<pv.size;i++)
        if(!var<T>::is_size_compatible_nothrow(pv[i],sz))
            return false;

    return true;
}
/*--------------------------------------------
 
 --------------------------------------------*/
template <class T,size_t N>
void var<T[N]>::is_size_compatible(py_var<T_EQUIV>& pv,const std::string& name,var<size_t>** sz)
{
    if(N!=pv.size)
        throw "expected size "+std::to_string(N)+" for argument '"+name+"'";

    if(sz)  sz+=1;
    
    for(size_t i=0;i<pv.size;i++)
    {
        try
        {
            var<T>::is_size_compatible(pv[i],name,sz);
        }
        
        catch(std::string& err_msg)
        {
            throw err_msg+"["+std::to_string(i)+"]";
        }
    }
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<class T,size_t N>
std::string var<T[N]>::type_name(var<size_t>** sz)
{
    if(sz)
        return std::string("[")+std::to_string(N)
        +std::string("]")+var<T>::type_name(sz+1);
    return std::string("[")+std::to_string(N)
    +std::string("]")+var<T>::type_name(sz);
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<class T,size_t N>
PyObject* var<T[N]>::build(T(&v)[N],size_t** sz)
{
    
    PyObject* py_obj=PyList_New(N);
    if(sz)
    {
        for(size_t i=0;i<N;i++)
            PyList_SET_ITEM(py_obj,i,var<T>::build(v[i],sz+1));
    }
    else
        for(size_t i=0;i<N;i++)
            PyList_SET_ITEM(py_obj,i,var<T>::build(v[i],sz));
        
    return py_obj;
}
/*--------------------------------------------
 constructor:
 usage: direct
 --------------------------------------------*/
template<class T,size_t N>
var<T[N]>::var(T (&v)[N],const char* name_):
ptr(&v),
__dsizes__(NULL),
Var(get_rank(),base_hash_code(),name_)
{
    size=0;
}
/*--------------------------------------------
 constructor:
 usage: direct
 --------------------------------------------*/
template<class T,size_t N>template<class... Ts>
var<T[N]>::var(T (&v)[N],const char* name_,Ts&&...____dsizes__):
Var(get_rank(),base_hash_code(),name_),
__dsizes__(NULL),
ptr(&v)
{
    assert(sizeof...(____dsizes__)==get_rank());
    size=0;
    __dsizes__=new var<size_t>*[rank];
    for(int i=0;i<rank;i++) __dsizes__[i]=NULL;
    Var::assign_dynamic_size(__dsizes__,____dsizes__...);
}
/*--------------------------------------------
 constructor:
 usage: indirect by same template class with
 one upper rank
 --------------------------------------------*/
template<class T,size_t N>
var<T[N]>::var(T(&v)[N],std::string&& name_,py_var<T_EQUIV>& pv,void** data_ptr):
ptr(&v),
__dsizes__(NULL),
Var(get_rank(),base_hash_code(),std::move(name_))
{
    size=N;
    for(size_t i=0;i<size;i++)
    {
        (vars+i)->~var<T>();
        new (vars+i) var<T>(v[i],name+"["+std::to_string(i)+"]",pv[i],data_ptr+1);
    }
}
/*--------------------------------------------
 default constructor:
 usage: indirect by same template class with 
 one upper rank
 --------------------------------------------*/
template<class T,size_t N>
var<T[N]>::var():
ptr(NULL),
__dsizes__(NULL),
Var(get_rank(),base_hash_code())
{
    size=N;
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<class T,size_t N>
var<T[N]>::~var()
{
    delete [] __dsizes__;
}
/*--------------------------------------------
 constructor:
 usage: indirect by Logics class
 it is only used for string literals
 --------------------------------------------*/
template<class T,size_t N>
var<T[N]>::var(T (&v)[N]):
ptr(&v),
__dsizes__(NULL),
Var(get_rank(),base_hash_code())
{
    size=N;
    if(std::is_same<const char,T>::value)
    {
        const char* v_str=reinterpret_cast<const char*>(v);
        name=std::string(v_str);
        char* name_=new char[2];
        name_[1]='\0';
        for(size_t i=0;i<size;i++)
        {
            name_[0]=v[i];
            //vars[i]=std::move(var<T>(v[i],name_));
            (vars+i)->~var<T>();
            new (vars+i) var<T>(v[i],std::string(v[i]));
        }
    }
    else
        throw "this is not allowed";
}
/*--------------------------------------------
 copy constructor:
 --------------------------------------------*/
template<class T,size_t N>
var<T[N]>::var(const var<T[N]>& r):
Var(static_cast<const Var&>(r)),
ptr(r.ptr),
__dsizes__(NULL)
{
    if(r.__dsizes__)
    {
        __dsizes__=new var<size_t>*[rank];
        memcpy(__dsizes__,r.__dsizes__,rank*sizeof(var<size_t>*));
    }
    for(size_t i=0;i<size;i++)
        vars[i]=r.vars[i];
}
/*--------------------------------------------
 move constructor:
 --------------------------------------------*/
template<class T,size_t N>
var<T[N]>::var(var<T[N]>&& r):
Var(static_cast<Var&&>(r)),
ptr(r.ptr),
__dsizes__(r.__dsizes__)
{
    r.ptr=NULL;
    r.__dsizes__=NULL;
    for(size_t i=0;i<size;i++)
        vars[i]=std::move(r.vars[i]);
}
/*--------------------------------------------
 copy assignment
 --------------------------------------------*/
template<class T,size_t N>
var<T[N]>& var<T[N]>::operator=(const var<T[N]>& r)
{
    if(this==&r) return *this;
    this-> ~var<T[N]>();
    new (this) var<T[N]>(r);
    return *this;
}
/*--------------------------------------------
 move assignment
 --------------------------------------------*/
template<class T,size_t N>
var<T[N]>& var<T[N]>::operator=(var<T[N]>&& r)
{
    if(this==&r) return *this;
    this-> ~var<T[N]>();
    new (this) var<T[N]>(std::move(r));
    return *this;
}
/*--------------------------------------------
 copy assignment from variable
 --------------------------------------------*/
template<class T,size_t N>
var<T[N]>& var<T[N]>::operator=(const T (&v)[N])
{
    size=N;
    for(size_t i=0;i<size;i++)
        vars[i]=v[i];
    for(size_t i=0;i<this->size;i++)
    {
        (this->vars+i)->~var<T>();
        new (this->vars+i) var<T>(v[i],this->name+"["+std::to_string(i)+"]");
    }
    
    return *this;
}
/*--------------------------------------------
 clone used by Logics class
 --------------------------------------------*/
template<class T,size_t N>
Var* var<T[N]>::clone()
{
    return new var<T[N]>(*this);
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<class T,size_t N>
const void* var<T[N]>::get_ptr()
{
    return ptr;
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<class T,size_t N>
var<T>& var<T[N]>::operator[](const size_t i)
{
    return vars[i];
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<class T,size_t N>
void var<T[N]>::set(PyObject* py_obj)
{
    py_var<T_EQUIV> pv;
    try
    {
        pv=py_var<T_EQUIV>(py_obj);
    }
    catch(int)
    {
        throw "fialure to deduce C++ type <"+py_var<T_EQUIV>::type_name()+">";
    }
    
    try
    {
        is_size_compatible(pv,name,__dsizes__);
    }
    catch (std::string& err_msg)
    {
        throw "fialure to deduce type <"+type_name()+">: "+err_msg;
    }
    
    set(pv);
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<class T,size_t N>
bool var<T[N]>::set_nothrow(PyObject* py_obj)
{
    py_var<T_EQUIV> pv;
    try
    {
        pv=py_var<T_EQUIV>(py_obj);
    }
    catch(int)
    {
        return false;
    }
    
    if(!is_size_compatible_nothrow(pv,__dsizes__))
        return false;
    set(pv);
    return true;
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<class T,size_t N>
void var<T[N]>::set(py_var<T_EQUIV>& pv)
{
    size=pv.size;
    
    deallocate(*ptr);
    size_t sz[rank];
    for(int i=0;i<rank;i++) sz[i]=0;
    accum_size(pv,sz);
    void* data_ptr[rank];
    for(int i=0;i<rank;i++) data_ptr[i]=NULL;
    allocate(data_ptr,*ptr,sz);
    for(size_t i=0;i<size;i++)
    {
        //vars[i]=std::move(var<T>((*ptr)[i],name+"["+std::to_string(i)+"]",pv[i],data_ptr+1));
        (vars+i)->~var<T>();
        new (vars+i) var<T>((*ptr)[i],name+"["+std::to_string(i)+"]",pv[i],data_ptr+1);
    }
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<class T,size_t N>
PyObject* var<T[N]>::get()
{
    if(std::is_same<T,char>::value || std::is_same<T,const char>::value)
        return PyString_FromString(reinterpret_cast<char*>(*ptr));
    
    PyObject* py_obj=PyList_New(size);
    for(size_t i=0;i<size;i++)
        PyList_SET_ITEM(py_obj,i,vars[i].get());
    return py_obj;
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<class T,size_t N>template<class ... Ts>
void  var<T[N]>::dynamic_size(Ts&&... ____dsizes__)
{
    assert(root);
    assert(sizeof...(____dsizes__)==get_rank());
    delete [] __dsizes__;
    __dsizes__=new var<size_t>*[rank];
    for(int i=0;i<rank;i++) __dsizes__[i]=NULL;
    Var::assign_dynamic_size(__dsizes__,____dsizes__...);
}
/*------------------------------------------------------------------------------------------------------------------------------------
 
 ------------------------------------------------------------------------------------------------------------------------------------*/
namespace MAPP_NS
{
    template<class T>
    class var<T*>: public Var
    {
    private:
    protected:
        var<size_t>** __dsizes__;
        const void* get_ptr();
    public:
        typedef typename var<T>::T_BASE T_BASE;
        typedef typename std::add_pointer<typename var<T>::T_EQUIV>::type T_EQUIV;
        constexpr static int get_rank();
        static size_t base_hash_code();
        static void allocate(void**,T*&,size_t*);
        static void deallocate(T*&);
        static void accum_size(py_var<T_EQUIV>&,size_t*);
        static bool is_size_compatible_nothrow(py_var<T_EQUIV>&,var<size_t>**);
        static void is_size_compatible(py_var<T_EQUIV>&,const std::string&,var<size_t>**);
        static std::string type_name(var<size_t>**);
        static PyObject* build(T*&,size_t**);
        var<T>* vars;
        T** ptr;
        
        var(T*&,const char*);
        template<class... Ts>
        var(T*&,const char*,Ts&&...);
        var(T*&,std::string&&,py_var<T_EQUIV>&,void**);
        var();
        virtual ~var();
        var(T*&);
        var(const var<T*>&);
        var(var<T*>&&);
        var<T*>& operator = (const var<T*>&);
        var<T*>& operator = (var<T*>&&);
        template<size_t N>
        var<T*>& operator = (const T(&)[N]);
        Var* clone();
        
        var<T>& operator [] (const size_t);
        void set(PyObject*);
        bool set_nothrow(PyObject*);
        void set(py_var<T_EQUIV>&);
        PyObject* get();
        std::string type_name(){return type_attr<T_BASE>::name()+" "+type_name(__dsizes__);}
        template<class ... Ts>
        void dynamic_size(Ts&&...);
    };
}
/*--------------------------------------------
 
 --------------------------------------------*/
template <class T>
constexpr int var<T*>::get_rank()
{
    return 1+var<T>::get_rank();
}
/*--------------------------------------------
 
 --------------------------------------------*/
template <class T>
size_t var<T*>::base_hash_code()
{
    return typeid(T_BASE).hash_code();
}
/*--------------------------------------------
 
 --------------------------------------------*/
template <class T>
void var<T*>::allocate(void** ptr,T*& v,size_t* sz)
{
    v=NULL;
    if(*sz) v=new T[*sz];
    *ptr=v;
    if(get_rank()==1) return;
    var<T>::allocate(ptr+1,v[0],sz+1);
}
/*--------------------------------------------
 
 --------------------------------------------*/
template <class T>
void var<T*>::deallocate(T*& v)
{
    if(v) var<T>::deallocate(v[0]);
    delete [] v;
    v=NULL;
}
/*--------------------------------------------
 
 --------------------------------------------*/
template <class T>
void var<T*>::accum_size(py_var<T_EQUIV>& pv,size_t* sz)
{
    *sz+=pv.size;
    if(get_rank()==1) return;
    
    for(size_t i=0;i<pv.size;i++)
        var<T>::accum_size(pv[i],sz+1);
}
/*--------------------------------------------
 
 --------------------------------------------*/
template <class T>
bool var<T*>::is_size_compatible_nothrow(py_var<T_EQUIV>& pv,var<size_t>** sz)
{
    if(sz)
    {
        if(*sz && **sz!=pv.size)
            return false;
        sz+=1;
    }
    
    for(size_t i=0;i<pv.size;i++)
        if(!var<T>::is_size_compatible_nothrow(pv[i],sz))
            return false;
    return true;
}
/*--------------------------------------------
 
 --------------------------------------------*/
template <class T>
void var<T*>::is_size_compatible(py_var<T_EQUIV>& pv,const std::string& name,var<size_t>** sz)
{
    if(sz)
    {
        if(*sz && **sz!=pv.size)
        {
            throw "expected size "+std::to_string(*((*sz)->ptr))+" ("+(*sz)->name+") for argument '"+name+"'";
        }
        sz+=1;
    }
    
    for(size_t i=0;i<pv.size;i++)
    {
        try
        {
            var<T>::is_size_compatible(pv[i],name,sz);
        }
        
        catch(std::string& err_msg)
        {
            throw err_msg+"["+std::to_string(i)+"]";
        }
    }
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<class T>
std::string var<T*>::type_name(var<size_t>** sz)
{
    if(sz)
    {
        if(*sz)
            return std::string("[")+std::string((*sz)->name)+std::string("]")+var<T>::type_name(sz+1);
        else
            return std::string("[]")+var<T>::type_name(sz+1);
    }
    return std::string("[]")+var<T>::type_name(sz);
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<class T>
PyObject* var<T*>::build(T*& v,size_t** sz)
{
    
    PyObject* py_obj=PyList_New(**sz);
    for(size_t i=0;i<**sz;i++)
        PyList_SET_ITEM(py_obj,i,var<T>::build(v[i],sz+1));
    return py_obj;
}
/*--------------------------------------------
 constructor:
 usage: direct
 --------------------------------------------*/
template<class T>
var<T*>::var(T*& v,const char* name_):
vars(NULL),
ptr(&v),
__dsizes__(NULL),
Var(get_rank(),base_hash_code(),name_)
{
    size=0;
    v=NULL;
}
/*--------------------------------------------
 constructor:
 usage: direct
 --------------------------------------------*/
template<class T>template<class... Ts>
var<T*>::var(T*& v,const char* name_,Ts&&...____dsizes__):
Var(get_rank(),base_hash_code(),name_),
vars(NULL),
__dsizes__(NULL),
ptr(&v)
{
    assert(sizeof...(____dsizes__)==get_rank());
    size=0;
    v=NULL;
    __dsizes__=new var<size_t>*[rank];
    for(int i=0;i<rank;i++) __dsizes__[i]=NULL;
    Var::assign_dynamic_size(__dsizes__,____dsizes__...);
        
}
/*--------------------------------------------
 constructor:
 usage: indirect by same template class with
 one upper rank
 --------------------------------------------*/
template<class T>
var<T*>::var(T*& v,std::string&& name_,py_var<T_EQUIV>& pv,void** data_ptr):
Var(get_rank(),base_hash_code(),std::move(name_)),
vars(NULL),
__dsizes__(NULL),
ptr(&v)
{
    v=static_cast<T*>(*data_ptr);
    size=pv.size;
    vars=new var<T>[size];
    for(size_t i=0;i<size;i++)
    {
        //vars[i]=std::move(var<T>((*ptr)[i],name+"["+std::to_string(i)+"]",pv[i],data_ptr+1));
        (vars+i)->~var<T>();
        new (vars+i) var<T>((*ptr)[i],name+"["+std::to_string(i)+"]",pv[i],data_ptr+1);
    }
    *data_ptr=v+size;
}
/*--------------------------------------------
 default constructor:
 usage: indirect by same template class with
 one upper rank
 --------------------------------------------*/
template<class T>
var<T*>::var():
Var(get_rank(),base_hash_code()),
vars(NULL),
__dsizes__(NULL),
ptr(NULL)
{
    size=0;
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<class T>
var<T*>::~var()
{
    delete [] __dsizes__;
    delete [] vars;
    if(root && ptr) var<T*>::deallocate(*ptr);
}
/*--------------------------------------------
 constructor:
 usage: indirect by Logics class
 it is only used for string literals
 --------------------------------------------*/
template<class T>
var<T*>::var(T*& v):
vars(NULL),
ptr(&v),
__dsizes__(NULL),
Var(get_rank(),base_hash_code())
{
    size=0;
    if(std::is_same<const char,T>::value)
    {
        const char* v_str=reinterpret_cast<const char*>(v);
        name=std::string(v_str);
        size=strlen(v_str)+1;
        char* name_=new char[2];
        name_[1]='\0';
        for(size_t i=0;i<size;i++)
        {
            name_[0]=v_str[i];
            //vars[i]=std::move(var<T>(v[i],name_));
            (vars+i)->~var<T>();
            new (vars+i) var<T>(v[i],std::string(v[i]));
        }
    }
    else
        throw "this is not allowed";
}
/*--------------------------------------------
 copy contructor:
 --------------------------------------------*/
template<class T>
var<T*>::var(const var<T*>& r):
Var(static_cast<const Var&>(r)),
ptr(r.ptr),
__dsizes__(NULL)
{
    if(r.__dsizes__)
    {
        __dsizes__=new var<size_t>*[rank];
        memcpy(__dsizes__,r.__dsizes__,rank*sizeof(var<size_t>*));
    }
    vars=new var<T>[r.size];
    for(size_t i=0;i<size;i++)
        vars[i]=r.vars[i];
}
/*--------------------------------------------
 move contructor:
 --------------------------------------------*/
template<class T>
var<T*>::var(var<T*>&& r):
Var(static_cast<Var&&>(r)),
ptr(r.ptr),
__dsizes__(r.__dsizes__)
{
    r.ptr=NULL;
    r.__dsizes__=NULL;
    vars=r.vars;
    r.vars=NULL;
}
/*--------------------------------------------
 copy assignment
 --------------------------------------------*/
template<class T>
var<T*>& var<T*>::operator=(const var<T*>& r)
{
    if(this==&r) return *this;
    this->~var<T*>();
    new (this) var<T*>(r);
    return *this;
}
/*--------------------------------------------
 move assignment
 --------------------------------------------*/
template<class T>
var<T*>& var<T*>::operator=(var<T*>&& r)
{
    if(this==&r) return *this;
    this->~var<T*>();
    new (this) var<T*>(std::move(r));
    return *this;
}
/*--------------------------------------------
 copy assignment from variable
 --------------------------------------------*/
template<class T>template<size_t N>
var<T*>& var<T*>::operator=(const T (&v_arr)[N])
{
    assert(root);
    delete [] this->vars;
    if(ptr) delete [] *(this->ptr);
    this->size=N;
    this->vars=new var<T>[this->size];
    T* v=new T[this->size];
    *(this->ptr)=v;

    for(size_t i=0;i<this->size;i++)
    {
        (this->vars+i)->~var<T>();
        new (this->vars+i) var<T>(v[i],this->name+"["+std::to_string(i)+"]");
    }
    return *this;
}
/*--------------------------------------------
 clone used by Logics class
 --------------------------------------------*/
template<class T>
Var* var<T*>::clone()
{
    return new var<T*>(*this);
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<class T>
const void* var<T*>::get_ptr()
{
    return ptr;
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<class T>
var<T>& var<T*>::operator[](const size_t i)
{
    return vars[i];
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<class T>
void var<T*>::set(PyObject* py_obj)
{
    py_var<T_EQUIV> pv;
    try
    {
        pv=py_var<T_EQUIV>(py_obj);
    }
    catch(int)
    {
        throw "fialure to deduce C++ type <"+py_var<T_EQUIV>::type_name()+">";
    }
    
    try
    {
        is_size_compatible(pv,name,__dsizes__);
    }
    catch (std::string& err_msg)
    {
        throw "fialure to deduce type <"+type_name()+">: "+err_msg;
    }
    
    set(pv);
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<class T>
bool var<T*>::set_nothrow(PyObject* py_obj)
{
    py_var<T_EQUIV> pv;
    try
    {
        pv=py_var<T_EQUIV>(py_obj);
    }
    catch(int)
    {
        return false;
    }
    
    if(!is_size_compatible_nothrow(pv,__dsizes__))
        return false;
    
    delete [] vars;
    vars=NULL;
    size=pv.size;
    if(size) vars=new var<T>[size];
    
    set(pv);
    return true;
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<class T>
void var<T*>::set(py_var<T_EQUIV>& pv)
{
    delete [] vars;
    vars=NULL;
    size=pv.size;
    if(size) vars=new var<T>[size];
    
    deallocate(*ptr);
    size_t sz[rank];
    for(int i=0;i<rank;i++) sz[i]=0;
    accum_size(pv,sz);
    void* data_ptr[rank];
    for(int i=0;i<rank;i++) data_ptr[i]=NULL;
    allocate(data_ptr,*ptr,sz);
    for(size_t i=0;i<size;i++)
    {
        (vars+i)->~var<T>();
        new (vars+i) var<T>((*ptr)[i],name+"["+std::to_string(i)+"]",pv[i],data_ptr+1);
    }
    
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<class T>
PyObject* var<T*>::get()
{
    if(std::is_same<T,char>::value || std::is_same<T,const char>::value)
        return PyString_FromString(reinterpret_cast<char*>(*ptr));
    
    PyObject* py_obj=PyList_New(size);
    for(size_t i=0;i<size;i++)
        PyList_SET_ITEM(py_obj,i,vars[i].get());
    return py_obj;
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<class T>template<class ... Ts>
void  var<T*>::dynamic_size(Ts&&... ____dsizes__)
{
    assert(root);
    assert(sizeof...(____dsizes__)==get_rank());
    delete [] __dsizes__;
    __dsizes__=new var<size_t>*[rank];
    for(int i=0;i<rank;i++) __dsizes__[i]=NULL;
    Var::assign_dynamic_size(__dsizes__,____dsizes__...);
}
/*------------------------------------------------------------------------------------------------------------------------------------
 
 ------------------------------------------------------------------------------------------------------------------------------------*/
namespace MAPP_NS
{
    template<class T>
    class var<symm<T>>:public var<T>
    {
    public:
        typedef typename var<T>::T_BASE T_BASE;
        typedef typename var<T>::T_EQUIV T_EQUIV;
        var(T&,const char*);
        template<class... Ts>
        var(T&,const char*,Ts&&...);
        ~var();
        void set(PyObject*);
        bool set_nothrow(PyObject*);
        std::string type_name();
    };
    
}
/*--------------------------------------------
 desstructor:
 --------------------------------------------*/
template<class T>
var<symm<T>>::~var()
{
}
/*--------------------------------------------
 constructor:
 usage: direct
 --------------------------------------------*/
template<class T>
var<symm<T>>::var(T& v,const char* name_):
var<T>(v,name_)
{
}
/*--------------------------------------------
 constructor:
 usage: direct
 --------------------------------------------*/
template<class T>template<class... Ts>
var<symm<T>>::var(T& v,const char* name_,Ts&&...____dsizes__):
var<T>(v,name_,____dsizes__...)
{
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<class T>
void var<symm<T>>::set(PyObject* py_obj)
{
    py_var<T_EQUIV> pv;
    try
    {
        pv=py_var<T_EQUIV>(py_obj);
    }
    catch(int)
    {
        try
        {
            pv=py_var<T_EQUIV>(&py_obj);
        }
        catch(int)
        {
            throw "fialure to deduce C++ type <"+py_var<T_EQUIV>::type_name()+"> or <"
            +py_var<typename std::remove_pointer<T_EQUIV>::type>::type_name()+">";
        }
    }
    
    if(!pv.is_symmetric())
        throw "fialure to deduce type <"+type_name()+">: "+
        var<T>::name+" must be represented in symmetric, triangular or voigt form";
    
    try
    {
        var<T>::is_size_compatible(pv,var<T>::name,var<T>::__dsizes__);
    }
    catch (std::string& err_msg)
    {
        throw "fialure to deduce type <"+type_name()+">: "+err_msg;
    }
    
    var<T>::set(pv);
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<class T>
bool var<symm<T>>::set_nothrow(PyObject* py_obj)
{
    py_var<T_EQUIV> pv;
    try
    {
        pv=py_var<T_EQUIV>(py_obj);
    }
    catch(int)
    {
        try
        {
            pv=py_var<T_EQUIV>(&py_obj);
        }
        catch(int)
        {
            return false;
        }
    }
    
    if(!pv.is_symmetric())
        return false;
    
    if(!var<T>::is_size_compatible_nothrow(pv,var<T>::__dsizes__))
        return false;
    
    var<T>::set(pv);
    return true;
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<class T>
std::string var<symm<T>>::type_name()
{
    return std::string("symmetric ")+var<T>::type_name();
}



#endif
