#ifndef __MAPP__var_op__
#define __MAPP__var_op__
#include "var.h"
/*------------------------------------------------------------------------------------------------------------------------------------
 
 ------------------------------------------------------------------------------------------------------------------------------------*/
namespace MAPP_NS
{
    template<class O,PyTypeObject& T>
    class OB
    {
    private:
    protected:
    public:
        PyObject* ob;
        OB();
        ~OB();
        OB(O* __ob);
        OB(PyObject* __ob);
        OB(const OB& r);
        OB(OB&& r);
        OB& operator=(const OB& r);
        OB& operator=(OB&& r);
        OB& operator=(PyObject* r);
        operator PyObject*&();
        
    };
}
using namespace MAPP_NS;
/*--------------------------------------------
 
 --------------------------------------------*/
template<class O,PyTypeObject& T>
OB<O,T>::OB()
{
    ob=NULL;
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<class O,PyTypeObject& T>
OB<O,T>::~OB()
{
    if(ob) Py_DECREF(ob);
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<class O,PyTypeObject& T>
OB<O,T>::OB(O* __ob):
ob(reinterpret_cast<PyObject*>(__ob))
{
    if(ob) Py_INCREF(ob);
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<class O,PyTypeObject& T>
OB<O,T>::OB(PyObject* __ob):
ob(__ob)
{
    if(ob) Py_INCREF(ob);
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<class O,PyTypeObject& T>
OB<O,T>::OB(const OB& r):
ob(r.ob)
{
    if(ob) Py_INCREF(ob);
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<class O,PyTypeObject& T>
OB<O,T>::OB(OB&& r):
ob(r.ob)
{
    r.ob=NULL;
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<class O,PyTypeObject& T>
OB<O,T>& OB<O,T>::operator=(const OB& r)
{
    if(ob) Py_DECREF(ob);
    ob=r.ob;
    if(ob) Py_INCREF(ob);
    return *this;
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<class O,PyTypeObject& T>
OB<O,T>& OB<O,T>::operator=(OB&& r)
{
    if(ob) Py_DECREF(ob);
    ob=r.ob;
    r.ob=NULL;
    return *this;
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<class O,PyTypeObject& T>
OB<O,T>& OB<O,T>::operator=(PyObject* r)
{
    if(r)
    {
        
        if(r->ob_type!=&T && r->ob_type->tp_base!=&T) throw 0;
        if(ob) Py_DECREF(ob);
        ob=r;
        if(ob) Py_INCREF(ob);
        return *this;
    }
    
    if(ob) Py_DECREF(ob);
    ob=NULL;
    return *this;
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<class O,PyTypeObject& T>
OB<O,T>::operator PyObject*&()
{
    return ob;
}
/*------------------------------------------------------------------------------------------------------------------------------------
 
 ------------------------------------------------------------------------------------------------------------------------------------*/
namespace MAPP_NS
{
    template<class C>
    using OP=OB<typename C::Object,C::TypeObject>;
    
    template<class O,PyTypeObject& T>class cpp_type2type_num<OB<O,T>>
    {public:static constexpr int type_num(){return NPY_OBJECT;};};
    template<class O,PyTypeObject& T>class type_attr<OB<O,T>>
    {public: static const std::string name(); static const OB<O,T> zero;};
}
/*------------------------------------------------------------------------------------------------------------------------------------
 
 ------------------------------------------------------------------------------------------------------------------------------------*/
template<class O,PyTypeObject& T>
const std::string type_attr<OB<O,T>>::name(){return std::string(T.tp_name);};
template<class O,PyTypeObject& T>
const OB<O,T> type_attr<OB<O,T>>::zero=OB<O,T>();
/*------------------------------------------------------------------------------------------------------------------------------------
 
 ------------------------------------------------------------------------------------------------------------------------------------*/
namespace MAPP_NS
{
    template<class O,PyTypeObject& T>
    class py_var<OB<O,T>>
    {
        friend class py_var<OB<O,T>*>;
    private:
        py_var(size_t*,void*&);
        py_var(int,long*,PyObject**&);
        py_var(py_var<OB<O,T>>&&);
    public:
        typedef OB<O,T> T_BASE;
        static constexpr int get_rank(){return 0;}
        OB<O,T> val;
        py_var(PyObject*);
        py_var();
        ~py_var();
        
        py_var(const py_var<OB<O,T>>&);
        operator OB<O,T>() const{return val;}
        py_var<OB<O,T>>& operator=(const py_var<OB<O,T>>&);
        py_var<OB<O,T>>& operator=(py_var<OB<O,T>>&&);
        bool operator==(const py_var<OB<O,T>>&);
        static std::string err_msg();
        const static std::string type_name();
    };
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<class O,PyTypeObject& T>
py_var<OB<O,T>>::py_var():
val(type_attr<OB<O,T>>::zero)
{
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<class O,PyTypeObject& T>
py_var<OB<O,T>>::py_var(PyObject* op):
val(type_attr<OB<O,T>>::zero)
{
    try
    {
        val=op;
    }
    catch(int)
    {
        throw 1;
    }
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<class O,PyTypeObject& T>
py_var<OB<O,T>>::py_var(int depth,long* sz,PyObject**& objs)
{
    this->~py_var<OB<O,T>>();
    try
    {
        new (this) py_var<OB<O,T>>(*objs);
    }
    catch(int)
    {
        throw 1;
    }
    objs++;
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<class O,PyTypeObject& T>
py_var<OB<O,T>>::py_var(size_t*,void*& data)
{
    PyObject** d=reinterpret_cast<PyObject**>(data);
    try
    {
        val=*d;
    }
    catch(int)
    {
        throw 1;
    }
    data=d+1;
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<class O,PyTypeObject& T>
py_var<OB<O,T>>::~py_var()
{
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<class O,PyTypeObject& T>
py_var<OB<O,T>>::py_var(const py_var<OB<O,T>>& r):
val(r.val)
{
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<class O,PyTypeObject& T>
py_var<OB<O,T>>& py_var<OB<O,T>>::operator=(const py_var<OB<O,T>>& r)
{
    this->~py_var<OB<O,T>>();
    new (this) py_var<OB<O,T>>(r);
    return *this;
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<class O,PyTypeObject& T>
py_var<OB<O,T>>::py_var(py_var<OB<O,T>>&& r):
val(std::move(r.val))
{
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<class O,PyTypeObject& T>
py_var<OB<O,T>>& py_var<OB<O,T>>::operator=(py_var<OB<O,T>>&& r)
{
    this->~py_var<OB<O,T>>();
    new (this) py_var<OB<O,T>>(std::move(r));
    return *this;
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<class O,PyTypeObject& T>
bool py_var<OB<O,T>>::operator==(const py_var<OB<O,T>>& r)
{
    return (this->val==r.val);
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<class O,PyTypeObject& T>
std::string py_var<OB<O,T>>::err_msg()
{
    return "fialure to deduce C++ type <"+std::string(T.tp_name)+">";
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<class O,PyTypeObject& T>
const std::string py_var<OB<O,T>>::type_name()
{
    return std::string(T.tp_name);
}
/*------------------------------------------------------------------------------------------------------------------------------------
 
 ------------------------------------------------------------------------------------------------------------------------------------*/
namespace MAPP_NS
{
    template<class O,PyTypeObject& T>
    class var<OB<O,T>>:public Var
    {
    private:
    protected:
        const void* get_ptr();
    public:
        
        typedef OB<O,T> T_BASE;
        typedef T_BASE T_EQUIV;
        static constexpr int get_rank();
        static size_t base_hash_code();
        static void allocate(void**,OB<O,T>&,size_t*);
        static void deallocate(OB<O,T>&);
        static void accum_size(py_var<T_EQUIV>&,size_t*){};
        static bool is_size_compatible_nothrow(py_var<T_EQUIV>&,var<size_t>**);
        static void is_size_compatible(py_var<T_EQUIV>&,const std::string&,var<size_t>**);
        static std::string type_name(var<size_t>**);
        static PyObject* build(O*&,size_t**);
        
        OB<O,T>* ptr;
        
        var(OB<O,T>&,const char*);
        
        var(OB<O,T>&,std::string&&,py_var<T_EQUIV>&,void**);
        var();
        ~var();
        var(OB<O,T>&&);
        var(const var<OB<O,T>>&);
        var(var<OB<O,T>>&&);
        var<OB<O,T>>& operator=(const var<OB<O,T>>&);
        var<OB<O,T>>& operator= (var<OB<O,T>>&&);
        var<OB<O,T>>& operator=(const OB<O,T>&);
        Var* clone();
        
        operator PyObject*() const;
        var<OB<O,T>>& operator [] (const size_t);
        void set(PyObject*);
        bool set_nothrow(PyObject*);
        PyObject* get();
        
        std::string type_name(){return std::string(T.tp_name);}
        
        bool operator==(Var&){return false;};
        bool operator!=(Var&){return false;};
        bool operator> (Var&){return false;};
        bool operator< (Var&){return false;};
        bool operator>=(Var&){return false;};
        bool operator<=(Var&){return false;};
    };
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<class O,PyTypeObject& T>
constexpr int var<OB<O,T>>::get_rank()
{
    return 0;
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<class O,PyTypeObject& T>
size_t var<OB<O,T>>::base_hash_code()
{
    return typeid(O).hash_code();
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<class O,PyTypeObject& T>
void var<OB<O,T>>::allocate(void**,OB<O,T>&,size_t*)
{}
/*--------------------------------------------
 
 --------------------------------------------*/
template<class O,PyTypeObject& T>
void var<OB<O,T>>::deallocate(OB<O,T>&)
{}
/*--------------------------------------------
 
 --------------------------------------------*/
template<class O,PyTypeObject& T>
bool var<OB<O,T>>::is_size_compatible_nothrow(py_var<T_EQUIV>&,var<size_t>**)
{
    return true;
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<class O,PyTypeObject& T>
void var<OB<O,T>>::is_size_compatible(py_var<T_EQUIV>&,const std::string&,var<size_t>**)
{
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<class O,PyTypeObject& T>
std::string var<OB<O,T>>::type_name(var<size_t>**)
{
    return std::string();
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<class O,PyTypeObject& T>
PyObject* var<OB<O,T>>::build(O*& v,size_t**)
{
    return reinterpret_cast<PyObject*>(v);
}
/*--------------------------------------------
 constructor:
 usage: direct
 --------------------------------------------*/
template<class O,PyTypeObject& T>
var<OB<O,T>>::var(OB<O,T>& v,const char* name_):
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
template<class O,PyTypeObject& T>
var<OB<O,T>>::var(OB<O,T>& v,std::string&& name_,py_var<T_EQUIV>& pv,void**):
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
template<class O,PyTypeObject& T>
var<OB<O,T>>::var():
ptr(NULL),
Var(get_rank(),base_hash_code())
{
    size=0;
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<class O,PyTypeObject& T>
var<OB<O,T>>::~var()
{}
/*--------------------------------------------
 copy constructor:
 --------------------------------------------*/
template<class O,PyTypeObject& T>
var<OB<O,T>>::var(const var<OB<O,T>>& r):
Var(static_cast<const Var&>(r)),
ptr(r.ptr)
{}
/*--------------------------------------------
 move constructor:
 --------------------------------------------*/
template<class O,PyTypeObject& T>
var<OB<O,T>>::var(var<OB<O,T>>&& r):
Var(static_cast<Var&&>(r)),
ptr(r.ptr)
{}
/*--------------------------------------------
 copy assignment
 --------------------------------------------*/
template<class O,PyTypeObject& T>
var<OB<O,T>>& var<OB<O,T>>::operator=(const var<OB<O,T>>& r)
{
    this-> ~var<OB<O,T>>();
    new (this) var<OB<O,T>>(r);
    return *this;
}
/*--------------------------------------------
 move assignment
 --------------------------------------------*/
template<class O,PyTypeObject& T>
var<OB<O,T>>& var<OB<O,T>>::operator=(var<OB<O,T>>&& r)
{
    this-> ~var<OB<O,T>>();
    new (this) var<OB<O,T>>(std::move(r));
    return *this;
}
/*--------------------------------------------
 copy assignment from variable
 --------------------------------------------*/
template<class O,PyTypeObject& T>
var<OB<O,T>>& var<OB<O,T>>::operator=(const OB<O,T>& r)
{
    ptr=&r;
    return *this;
}
/*--------------------------------------------
 clone used by Logics class
 --------------------------------------------*/
template<class O,PyTypeObject& T>
Var* var<OB<O,T>>::clone()
{
    return new var<OB<O,T>>(*this);
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<class O,PyTypeObject& T>
const void* var<OB<O,T>>::get_ptr()
{
    return ptr;
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<class O,PyTypeObject& T>
var<OB<O,T>>::operator PyObject*() const
{
    return *ptr;
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<class O,PyTypeObject& T>
var<OB<O,T>>& var<OB<O,T>>::operator[](const size_t)
{
    return *this;
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<class O,PyTypeObject& T>
void var<OB<O,T>>::set(PyObject* py_obj)
{
    try
    {
        *ptr=py_obj;
    }
    catch(int)
    {
        throw "fialure to deduce type <"+type_name()+">";
    }
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<class O,PyTypeObject& T>
bool var<OB<O,T>>::set_nothrow(PyObject* py_obj)
{
    try
    {
        *ptr=py_obj;
    }
    catch(int)
    {
        return false;
    }
    return true;
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<class O,PyTypeObject& T>
PyObject* var<OB<O,T>>::get()
{
    return ptr->ob;
}

/*------------------------------------------------------------------------------------------------------------------------------------
 
 ------------------------------------------------------------------------------------------------------------------------------------*/
namespace MAPP_NS
{
    template<class T0,class T1>
    class py_var<std::map<T0,T1>>
    {
        friend class py_var<std::map<T0,T1>*>;
    private:
        py_var(size_t*,void*&);
        py_var(int,long*,PyObject**&);
        py_var(py_var<std::map<T0,T1>>&&);
    public:
        typedef std::map<T0,T1> T_BASE;
        static constexpr int get_rank(){return 0;}
        std::map<T0,T1> val;
        py_var(PyObject*);
        py_var();
        ~py_var();
        
        py_var(const py_var<std::map<T0,T1>>&);
        operator std::map<T0,T1>() const{return val;}
        py_var<std::map<T0,T1>>& operator=(const py_var<std::map<T0,T1>>&);
        py_var<std::map<T0,T1>>& operator=(py_var<std::map<T0,T1>>&&);
        bool operator==(const py_var<std::map<T0,T1>>&);
        static std::string err_msg();
        const static std::string type_name();
    };
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<class T0,class T1>
py_var<std::map<T0,T1>>::py_var():
val(type_attr<std::map<T0,T1>>::zero)
{
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<class T0,class T1>
py_var<std::map<T0,T1>>::py_var(PyObject* op):
val(type_attr<std::map<T0,T1>>::zero)
{
    if(!PyDict_Check(op))
        throw 1;
    PyObject* keys=PyDict_Keys(op);
    PyObject* vals=PyDict_Values(op);
    py_var<T0*> __keys;
    py_var<T1*> __vals;
    try
    {
        __keys=py_var<T0*>(keys);
        __vals=py_var<T1*>(vals);
        
    }
    catch(int)
    {
        Py_DECREF(keys);
        Py_DECREF(vals);
        throw 1;
    }
    
    val.clear();
    for(size_t i=0;i<Py_SIZE(keys);i++)
        val[__keys[i].val]=__vals[i].val;
    
    Py_DECREF(keys);
    Py_DECREF(vals);
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<class T0,class T1>
py_var<std::map<T0,T1>>::py_var(int depth,long* sz,PyObject**& objs)
{
    this->~py_var<std::map<T0,T1>>();
    try
    {
        new (this) py_var<std::map<T0,T1>>(*objs);
    }
    catch(int)
    {
        throw 1;
    }
    objs++;
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<class T0,class T1>
py_var<std::map<T0,T1>>::py_var(size_t*,void*& data)
{
    PyObject** d=reinterpret_cast<PyObject**>(data);
    if(!PyDict_Check(*d))
        throw 1;
    PyObject* keys=PyDict_Keys(*d);
    PyObject* vals=PyDict_Values(*d);
    py_var<T0*> __keys;
    py_var<T1*> __vals;
    try
    {
        __keys=py_var<T0*>(keys);
        __vals=py_var<T1*>(vals);
        
    }
    catch(int)
    {
        Py_DECREF(keys);
        Py_DECREF(vals);
        throw 1;
    }
    
    val.clear();
    for(size_t i=0;i<Py_SIZE(keys);i++)
        val[__keys[i].val]=__vals[i].val;
    
    Py_DECREF(keys);
    Py_DECREF(vals);
    data=d+1;
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<class T0,class T1>
py_var<std::map<T0,T1>>::~py_var()
{
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<class T0,class T1>
py_var<std::map<T0,T1>>::py_var(const py_var<std::map<T0,T1>>& r):
val(r.val)
{
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<class T0,class T1>
py_var<std::map<T0,T1>>& py_var<std::map<T0,T1>>::operator=(const py_var<std::map<T0,T1>>& r)
{
    this->~py_var<std::map<T0,T1>>();
    new (this) py_var<std::map<T0,T1>>(r);
    return *this;
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<class T0,class T1>
py_var<std::map<T0,T1>>::py_var(py_var<std::map<T0,T1>>&& r):
val(std::move(r.val))
{
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<class T0,class T1>
py_var<std::map<T0,T1>>& py_var<std::map<T0,T1>>::operator=(py_var<std::map<T0,T1>>&& r)
{
    this->~py_var<std::map<T0,T1>>();
    new (this) py_var<std::map<T0,T1>>(std::move(r));
    return *this;
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<class T0,class T1>
bool py_var<std::map<T0,T1>>::operator==(const py_var<std::map<T0,T1>>& r)
{
    return (this->val==r.val);
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<class T0,class T1>
std::string py_var<std::map<T0,T1>>::err_msg()
{
    return "fialure to deduce C++ type <map "+type_attr<T0>::name()+" "+type_attr<T1>::name()+">";
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<class T0,class T1>
const std::string py_var<std::map<T0,T1>>::type_name()
{
    return "map "+type_attr<T0>::name()+" "+type_attr<T1>::name();
}

#endif
