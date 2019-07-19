#ifndef __MAPP__api__
#define __MAPP__api__

#include "logics.h"

using namespace MAPP_NS;
/*------------------------------------------------------------------------------------------------------------------------------------
 
 ------------------------------------------------------------------------------------------------------------------------------------*/
namespace MAPP_NS
{
               
    
    template<class T>
    class VarAPI
    {
    private:
    protected:
    public:
        typedef typename cpp_type<T>::type S;
        S val;
        var<T> __var__;
        Logics logics[var<T>::get_rank()+1];
        
        VarAPI(const char*);
        template<class ... Ts>
        VarAPI(const char*,Ts&&...);
        ~VarAPI();
    
        
        PyObject* get();
        int set(PyObject*);
        bool set_nothrow(PyObject*);
        const char* name() const;
        void* get_ptr(){return &val;}
    };
    

}
/*--------------------------------------------
 
 --------------------------------------------*/
template<class T>
VarAPI<T>::VarAPI(const char* name):
__var__(val,name)
{
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<class T>template<class ... Ts>
VarAPI<T>::VarAPI(const char* name,Ts&&...dsizes):
__var__(val,name,dsizes...)
{
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<class T>
VarAPI<T>::~VarAPI()
{}
/*--------------------------------------------
 
 --------------------------------------------*/
template<class T>
PyObject* VarAPI<T>::get()
{
    if(__var__.obj_set)
        return __var__.get();
    
    PyErr_Format(PyExc_TypeError,"%s has not been set yet",__var__.name.c_str());
    return NULL;
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<class T>
int VarAPI<T>::set(PyObject* op)
{
    try
    {
        __var__.set(op);
    }
    catch(std::string& err_msg)
    {
        PyErr_SetString(PyExc_TypeError,err_msg.c_str());
        return -1;
    }
    try
    {
        Logics::var_log(&__var__,logics+var<T>::get_rank());
    }
    catch(std::string& err_msg)
    {
        PyErr_SetString(PyExc_ValueError,err_msg.c_str());
        return -1;
    }
    __var__.obj_set=true;
    return 0;
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<class T>
bool VarAPI<T>::set_nothrow(PyObject* op)
{
    if(!__var__.set_nothrow(op))
        return false;
    
    if(!Logics::var_log_nothrow(&__var__,logics+var<T>::get_rank()))
        return false;
    
    return true;
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<class T>
const char* VarAPI<T>::name() const
{
    return __var__.name.c_str();
}
/*------------------------------------------------------------------------------------------------------------------------------------
 
 ------------------------------------------------------------------------------------------------------------------------------------*/
namespace MAPP_NS
{
    template <class... Ts>
    class FuncTuple
    {};
    
    template <class T, class... Ts>
    class FuncTuple<T, Ts...>
    {
    private:
    protected:
    public:
        typedef T type;
        
        FuncTuple<Ts...> next_func;
        VarAPI<T> __var_api__;
        PyObject* op;
        
        
        void match_kwd(const char*,PyObject*);
        void match(PyObject*,Py_ssize_t);
        void reset();
        int set();
        void expected(size_t);
        
        template <class... Ss>
        FuncTuple(const char* s, Ss... ss):
        next_func(ss...),
         __var_api__(s),
        op(NULL)
        {}
        
    
        FuncTuple(const char*const* s):
        next_func(s+1),
         __var_api__(*s),
        op(NULL)
        {}
        
        
        
        template<class S,class... TTs>
        PyObject* act(S s,PyObject* cls,TTs ... vars)
        {
            return next_func.act(s,cls,vars..., __var_api__.val);
        }
        
        template<class S>
        PyObject* operator()(S s,PyObject* cls)
        {
            return next_func.act(s,cls,__var_api__.val);
        }
        
    };
    

    template <class T>
    class FuncTuple<T>
    {
    private:
    protected:
    public:
        typedef T type;
        VarAPI<T>  __var_api__;
        PyObject* op;
        
        void match_kwd(const char*,PyObject*);
        void match(PyObject*,Py_ssize_t);
        void reset();
        int set();
        void expected(size_t);
        
        FuncTuple(const char* s):
        __var_api__(s),
        op(NULL)
        {}
        FuncTuple(const char*const* s):
        __var_api__(*s),
        op(NULL)
        {}
        
        template<class S,class... TTs>
        PyObject* act(S s,PyObject* cls,TTs ... vars)
        {
            return s(cls,vars...,__var_api__.val);
        }
        template<class S>
        PyObject* operator()(S s,PyObject* cls)
        {
            return s(cls,__var_api__.val);
        }
    };
}
/*--------------------------------------------
 
 --------------------------------------------*/
template <class T, class... Ts>
void FuncTuple<T, Ts...>::match_kwd(const char* kwd,PyObject* val)
{
    bool match=!strcmp(__var_api__.name(),kwd);
    if(match && op)
        throw 0; /*multiple values*/
    if(match)
    {
        op=val;
        return;
    }
    return next_func.match_kwd(kwd,val);
}
/*--------------------------------------------
 
 --------------------------------------------*/
template <class T>
void FuncTuple<T>::match_kwd(const char* kwd,PyObject* val)
{
    bool match=!strcmp(__var_api__.name(),kwd);
    if(match && op)
        throw 0; /*multiple values*/
    if(match)
    {
        op=val;
        return;
    }
    throw 1; /*unexpected keyword*/
}
/*--------------------------------------------
 
 --------------------------------------------*/
template <class T, class... Ts>
void FuncTuple<T, Ts...>::match(PyObject* tuple,Py_ssize_t i)
{
    op=PySequence_Fast_GET_ITEM(tuple,i);
    if(i==PyObject_Size(tuple)-1)
        return;
    return next_func.match(tuple,i+1);
}
/*--------------------------------------------
 
 --------------------------------------------*/
template <class T>
void FuncTuple<T>::match(PyObject* tuple,Py_ssize_t i)
{
    op=PySequence_Fast_GET_ITEM(tuple,i);
}
/*--------------------------------------------
 
 --------------------------------------------*/
template <class T, class... Ts>
void FuncTuple<T, Ts...>::reset()
{
    op=NULL;
    next_func.reset();
}
/*--------------------------------------------
 
 --------------------------------------------*/
template <class T>
void FuncTuple<T>::reset()
{
    op=NULL;
}
/*--------------------------------------------
 
 --------------------------------------------*/
template <class T, class... Ts>
int FuncTuple<T, Ts...>::set()
{
    if(op && __var_api__.set(op)==-1)
        return -1;
    return next_func.set();
}
/*--------------------------------------------
 
 --------------------------------------------*/
template <class T>
int FuncTuple<T>::set()
{
    if(op && __var_api__.set(op)==-1)
        return -1;
    return 0;
}
/*--------------------------------------------
 
 --------------------------------------------*/
template <class T, class... Ts>
void FuncTuple<T, Ts...>::expected(size_t i)
{
    if(i==0) return;

    if(op==NULL)
        throw __var_api__.name();
    
    return next_func.expected(i-1);
}
/*--------------------------------------------
 
 --------------------------------------------*/
template <class T>
void FuncTuple<T>::expected(size_t i)
{
    if(i==0) return;
    
    if(op==NULL)
        throw __var_api__.name();
}
/*------------------------------------------------------------------------------------------------------------------------------------
 
 ------------------------------------------------------------------------------------------------------------------------------------*/
namespace MAPP_NS
{
    
    
    template <size_t I, class T>
    class  FuncElement;
    
    template <size_t I,class T,class ... Ts>
    class  FuncElement<I,FuncTuple<T,Ts...>>
    {
    public:
        typedef typename FuncElement<I-1,FuncTuple<Ts...>>::type type;
        typedef typename FuncElement<I-1,FuncTuple<Ts...>>::val_type val_type;
        
        static VarAPI<type>& var_api(FuncTuple<T,Ts...>& f)
        {
            return FuncElement<I-1,FuncTuple<Ts...>>::var_api(f.next_func);
        }
        
        static ::var<type>& var(FuncTuple<T,Ts...>& f)
        {
            return FuncElement<I-1,FuncTuple<Ts...>>::var(f.next_func);
        }
        
        static val_type& val(FuncTuple<T,Ts...>& f)
        {
            return FuncElement<I-1,FuncTuple<Ts...>>::val(f.next_func);
        }
    };
    template <class T,class ... Ts>
    class  FuncElement<0,FuncTuple<T,Ts...>>
    {
    public:
        typedef T type;
        typedef typename VarAPI<T>::S val_type;
        static VarAPI<type>& var_api(FuncTuple<T,Ts...>& f)
        {
            return f.__var_api__;
        }
        
        static ::var<type>& var(FuncTuple<T,Ts...>& f)
        {
            return f.__var_api__.__var__;
        }
        static val_type& val(FuncTuple<T,Ts...>& f)
        {
            return f.__var_api__.val;
        }
    };

    
    template <class... Ts>
    class FuncAPI;

    template <class T,class... Ts>
    class FuncAPI<T,Ts ...>
    {
    private:
    protected:
    public:
        FuncTuple<T,Ts...> func;
        size_t nvars;
        size_t noptionals;
        std::string name;
        
        template<class ... Ss>
        FuncAPI(const char* __name,Ss ... vs):
        name(__name+std::string("()")),
        func(vs...),
        nvars(sizeof...(vs)),
        noptionals(0){}
        
        FuncAPI(const char* __name,std::initializer_list<const char*> __names):
        func(__names.begin()),
        nvars(__names.size()),
        noptionals(0),
        name(__name+std::string("()")){}
        
        int operator()(PyObject*,PyObject*);
        template<class F>
        PyObject* operator()(F,PyObject*,PyObject*,PyObject*);

        
        template <size_t I>
        VarAPI<typename FuncElement<I, FuncTuple<T,Ts...> >::type>& get()
        {return FuncElement<I, FuncTuple<T,Ts...> >::var_api(func);}
        
        template <size_t I>
        var<typename FuncElement<I, FuncTuple<T,Ts...> >::type>& v()
        {return get<I>().__var__;}
        
        template <size_t I>
        typename FuncElement<I, FuncTuple<T,Ts...> >::val_type& val()
        {return get<I>().val;}
        template <size_t I>
        typename FuncElement<I, FuncTuple<T,Ts...> >::val_type&& mov()
        {return std::move(get<I>().val);}
        
        template <size_t I>
        Logics* logics()
        {return get<I>().logics;}
        
        /*
        template<size_t I,size_t I0,size_t ... Is>
        void remap(const char*,const std::string*,const size_t);
        */
        
        template<size_t I,size_t I0,size_t ... Is>
        void ___remap(const char* name,const std::string* elems,const size_t nelems)
        {
            
        }
        
        template<size_t,size_t,size_t ... Is,class ... TTs>
        int remap(const char*,const std::string*,const size_t nelems,TTs& ...);
        template<size_t,class ... TTs>
        int remap(const char*,const std::string*,const size_t,TTs& ...);
    };
}
/*--------------------------------------------
 
 --------------------------------------------*/
template <class T,class... Ts>
int FuncAPI<T,Ts...>::operator()(PyObject* args,PyObject* kwds)
{
    size_t max_nargs=nvars;
    size_t min_nargs=nvars-noptionals;
    
    size_t nargs=PyObject_Size(args);
    size_t nkwds=0;
    if(kwds) nkwds=PyObject_Size(kwds);
    nargs+=nkwds;
    
    if(nargs<min_nargs)
    {
        if(min_nargs==max_nargs)
        {
            PyErr_Format(PyExc_TypeError,"%s takes exactly %zu %s (%zu given)",
            name.c_str(),min_nargs,min_nargs==1?"argument":"arguments",nargs);
            return -1;
        }
        else
        {
            PyErr_Format(PyExc_TypeError,"%s takes at least %zu %s (%zu given)",
            name.c_str(),min_nargs,min_nargs==1?"argument":"arguments",nargs);
            return -1;
        }
    }
    if(nargs>max_nargs)
    {
        PyErr_Format(PyExc_TypeError,"%s takes at most %zu %s (%zu given)",
        name.c_str(),max_nargs,max_nargs==1?"argument":"arguments",nargs);
        return -1;
    }
        
    if(PyObject_Size(args))
        func.match(args,0);

    if(nkwds)
    {
        PyObject* keys=PyDict_Keys(kwds);
        PyObject* vals=PyDict_Values(kwds);
        for(size_t i=0;i<nkwds;i++)
        {
            try
            {
                func.match_kwd(__PyString_AsString(PyList_GetItem(keys,i)),PyList_GetItem(vals,i));
            }
            catch (int __err)
            {
                if(__err)
                    PyErr_Format(PyExc_TypeError,"%s got an unexpected keyword "
                    "argument keyword argument '%s'",name.c_str(),
                    __PyString_AsString(PyList_GetItem(keys,i)));
                else
                    PyErr_Format(PyExc_TypeError,"%s got multiple values for "
                    "keyword argument '%s'",name.c_str(),
                    __PyString_AsString(PyList_GetItem(keys,i)));
                func.reset();
                Py_DECREF(keys);
                Py_DECREF(vals);
                return -1;
            }
        }
        
        Py_DECREF(keys);
        Py_DECREF(vals);
    }
    
    try
    {
        func.expected(min_nargs);
    }
    catch(const char* var_name)
    {
        PyErr_Format(PyExc_TypeError,"%s expects "
        "value for argument '%s'",name.c_str(),var_name);
        func.reset();
        return -1;
    }
    
    if(func.set()==-1)
    {
        func.reset();
        return -1;
    }
    return 0;
}
/*--------------------------------------------
 
 --------------------------------------------*/
template <class T,class... Ts>template<class F>
PyObject* FuncAPI<T,Ts...>::operator()(F f,PyObject* cls,PyObject* args,PyObject* kwds)
{
    if(this->operator()(args,kwds)==-1)
        return NULL;

    PyObject* ans=func(f,cls);
    func.reset();
    return ans;
}
/*------------------------------------------------------------------------------------------------------------------------------------
 
 ------------------------------------------------------------------------------------------------------------------------------------*/
namespace MAPP_NS
{
    template <>
    class FuncAPI<>
    {
    private:
        template<class T,class ... Ts>
        static void same_size(const std::string&,const size_t,T&,Ts&...);
        template<class T>
        static void same_size(const std::string&,const size_t,T&);
        static void same_size(const std::string&,const size_t){}
        template<class T,class ... Ts>
        static void reset_n_remap(const size_t*,const size_t,T&,Ts& ...);
        template<class T>
        static void reset_n_remap(const size_t*,const size_t,T&);
        static void reset_n_remap(const size_t*){};
        static size_t* mapping(const std::string*,const size_t,const VarAPI<std::string*>&);
    protected:
    public:
        std::string name;
        FuncAPI(const char* __name):
        name(__name){}
        int operator()(PyObject*,PyObject*);
        template<class F>
        PyObject* operator()(F,PyObject*,PyObject*,PyObject*);
        
        template<class ... Ts>
        static void __remap(const char*,const std::string*,const size_t,const VarAPI<std::string*>&,Ts&...);
        
    };
}
/*--------------------------------------------
 
 --------------------------------------------*/
inline int FuncAPI<>::operator()(PyObject* args,PyObject* kwds)
{
    size_t nargs=PyObject_Size(args);
    size_t nkwds=0;
    if(kwds) nkwds=PyObject_Size(kwds);
    nargs+=nkwds;
    
    if(nargs!=0)
    {
        PyErr_Format(PyExc_TypeError,"%s takes no arguments (%zu given)",
        name.c_str(),nargs);
        return -1;
    }
    
    return 0;
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<class F>
PyObject* FuncAPI<>::operator()(F f,PyObject* cls,PyObject* args,PyObject* kwds)
{
    if(this->operator()(args,kwds)==-1)
        return NULL;
    return f();
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<class T,class ... Ts>
void FuncAPI<>::same_size(const std::string& name,const size_t size,T& v,Ts& ... vs)
{
    if(v.__var__.rank>0 && size!=v.__var__.size)
        throw "'"+v.__var__.name+"' is expected to have the same size as '"+name+"'";
    
    same_size(name,size,vs...);
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<class T>
void FuncAPI<>::same_size(const std::string& name,const size_t size,T& v)
{
    if(v.__var__.rank>0 && size!=v.__var__.size)
        throw "'"+v.__var__.name+"' is expected to have the same size as '"+name+"'";
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<class T,class ... Ts>
void FuncAPI<>::reset_n_remap(const size_t* map,const size_t map_sz,T& v,Ts& ... vs)
{
    if(v.__var__.rank>0)
        v.__var__.reset_n_remap(map,map_sz);
    reset_n_remap(map,map_sz,vs...);
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<class T>
void FuncAPI<>::reset_n_remap(const size_t* map,const size_t map_sz,T& v)
{
    if(v.__var__.rank>0)
        v.__var__.reset_n_remap(map,map_sz);
}
/*--------------------------------------------
 
 --------------------------------------------*/
inline size_t* FuncAPI<>::mapping(const std::string* elems,const size_t nelems,const VarAPI<std::string*>& v)
{
    
    for(size_t i=0;i<v.__var__.size;i++)
        for(size_t j=i+1;j<v.__var__.size;j++)
            if(!strcmp(v.val[i].c_str(),v.val[j].c_str()))
                throw "duplicate element '"+v.val[i]+"' in '"+v.__var__.name+"'";
    
    
    size_t* map=new size_t[nelems];
    
    for(size_t i=0;i<nelems;i++)
    {
        bool found=false;
        for(size_t j=0;j<v.__var__.size && !found;j++)
            if(!strcmp(v.val[j].c_str(),elems[i].c_str()))
            {
                map[i]=j;
                found=true;
            }
        
        if(!found)
        {
            delete [] map;
            throw "unable to find elment '"+elems[i]+"' in '"+v.__var__.name+"'";
        }
    }
    
    return map;
    
    
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<class ... Ts>
void FuncAPI<>::__remap(const char* __elems_name,const std::string* __elems,const size_t __nelems,const VarAPI<std::string*>& __elem_api,Ts&... vs)
{
    if(__elem_api.val)
    {
        try
        {
            same_size(__elem_api.__var__.name,__elem_api.__var__.size,vs...);
        }
        catch(std::string& err_msg)
        {
            throw err_msg;
        }
        
        size_t* map=NULL;
        try
        {
            map=mapping(__elems,__nelems,__elem_api);
        }
        catch(std::string& err_msg)
        {
            throw err_msg;
        }
        
        reset_n_remap(map,__nelems,vs ...);
        delete [] map;
    }
    else
    {
        try
        {
            same_size(std::string(__elems_name),__nelems,vs...);
        }
        catch(std::string& err_msg)
        {
            throw err_msg;
        }
    }
    
}
/*--------------------------------------------
 
 --------------------------------------------*/
template <class T,class... Ts> template<size_t M,size_t I,size_t ... Is,class ... TTs>
int FuncAPI<T,Ts...>::remap(const char* name,const std::string* elems,const size_t nelems,TTs& ... vs)
{
    return remap<M,Is...>(name,elems,nelems,vs ...,get<I>());
}
/*--------------------------------------------
 
 --------------------------------------------*/
template <class T,class... Ts> template<size_t M,class ... TTs>
int FuncAPI<T,Ts...>::remap(const char* name,const std::string* elems,const size_t nelems,TTs& ... vs)
{
    try
    {
        FuncAPI<>::__remap(name,elems,nelems,get<M>(),vs ...);
    }
    catch(std::string& err_msg)
    {
        PyErr_SetString(PyExc_TypeError,err_msg.c_str());
        return -1;
    }
    return 0;
}


#endif
