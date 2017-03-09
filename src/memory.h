#ifndef __MAPP__memory__
#define __MAPP__memory__
#include <exception>
#include <new>
#include <cstring>
#include "macros.h"
namespace MAPP_NS
{

    template<class T>
    class Mem
    {
    public:
        static void alloc(T& data)
        {}
        static void dealloc(T& data)
        {}
    };
    
    template<class T,size_t N>
    class Mem<T[N]>
    {
    public:
        template<class...Ts>
        static void alloc(T(&data)[N],size_t sz,Ts ... szs)
        {
            alloc(data[0],N*sz,szs...);
            if(std::is_pointer<T>::value)
            {
                T* __data=data;
                for(size_t i=1;i<sz;i++) __data[i]=__data[i-1]+N;
            }
        }
        static void alloc(T(&data)[N])
        {
        }
        static void dealloc(T(&data)[N])
        {
            Mem<T>::dealloc(data[0]);
        }
    };
    template<class T>
    class Mem<T*>
    {
    public:
        template<class...Ts>
        static void alloc(T*& data,size_t sz0,size_t sz1,Ts ... szs)
        {
            if(sz0==0)
            {
                data=NULL;
                return;
            }
            data=new T[sz0];
            Mem<T>::alloc(*data,sz0*sz1,szs...);
            for(size_t i=1;i<sz0;i++) data[i]=data[i-1]+sz1;
        }
        static void alloc(T*& data, size_t sz0)
        {
            if(sz0==0)
            {
                data=NULL;
                return;
            }
            data=new T[sz0];
        }
        
        static void dealloc(T*& data)
        {
            if(data) Mem<T>::dealloc(*data);
            delete [] data;
            data=NULL;
        }
    };
    
    class Memory
    {
    private:
    protected:
    public:
        
        
        template<class T,class...Ts>
        static inline void alloc(T& data,Ts ... szs)
        {
            Mem<T>::alloc(data,szs...);
        }
        template<class T>
        static inline void dealloc(T& data)
        {
            Mem<T>::dealloc(data);
        }
        
        
        template<class T>
        static inline void shrink_to_fit(T*& data,size_t size,size_t capacity)
        {
            if(capacity==size) return;
            
            T* __data=NULL;
            if(size) __data=new T[size];
            memcpy(__data,data,size*sizeof(T));
            delete [] data;
            data=__data;
        }
        
        template <typename TYPE>
        static TYPE* create(TYPE*&,long,const char*
        ,int,const char*,const char*);

        template <typename TYPE>
        static TYPE** create(TYPE**&,long,long,const char*
        ,int,const char*,const char*);
        
        template <typename TYPE>
        static TYPE** create_2d(TYPE**&,long,long,const char*
        ,int,const char*,const char*);
        
        template <typename TYPE>
        static void del_2d(TYPE**&);
        
        template <typename TYPE>
        static TYPE* grow(TYPE*&,long,long,const char*
        ,int,const char*,const char*);
        template <typename T>
        static void grow(T*&,size_t,size_t);
    };
}
using namespace MAPP_NS;
/*--------------------------------------------
 create 1d vector
 --------------------------------------------*/
template<typename TYPE>
TYPE* Memory::create(TYPE*& array,long d0
,const char* name,int line_no,const char* file
,const char* function)
{
    try
    {
        array = new TYPE [d0];
    }
    catch(std::bad_alloc&)
    {
        /*-------temp_remove-------
        Error::abort("memory allocation failure "
        "in file  %s, function %s, line: %d for "
        "variable: %s",file,function,line_no,name);
         */
    }
    return array;
}
/*--------------------------------------------
 create 2d vector
 --------------------------------------------*/
template <typename TYPE>
TYPE** Memory::create(TYPE**& array,long d0
,long d1,const char* name,int line_no,
const char* file,const char* function)
{
    create(array,d0,name,line_no,file,function);
    for(int i=0;i<d0;i++)
        create(array[i],d1,name,line_no,file,function);
    return array;
}
/*--------------------------------------------
 create 2d vector
 --------------------------------------------*/
template <typename TYPE>
TYPE** Memory::create_2d(TYPE**& array,long d0
,long d1,const char* name,int line_no,
const char* file,const char* function)
{
    create(array,d0,name,line_no,file,function);
    create(*array,d0*d1,name,line_no,file,function);
    
    for(int i=1;i<d0;i++)
        array[i]=array[i-1]+d1;
    return array;
}
/*--------------------------------------------
 create 2d vector
 --------------------------------------------*/
template <typename TYPE>
void Memory::del_2d(TYPE**& array)
{
    if(array==NULL)
    {
        delete [] array;
        return;
    }
    delete [] *array;
    delete [] array;    
}
/*--------------------------------------------
 grow 1d vector
 --------------------------------------------*/
template <typename TYPE>
TYPE* Memory::grow(TYPE*& array,long oldsize,
long newsize,const char* name,int line_no,
const char* file,const char* function)
{
    if(oldsize==0)
    {
        return create(array,newsize,name,line_no,file,function);
    }
    else if(oldsize==newsize)
    {
        return array;
    }
    else
    {
        
        TYPE* newarray=array;
        try
        {
            long size1=newsize;
            long size2=oldsize;
            long size=MIN(size1,size2);
            newarray = new TYPE[newsize];
            memcpy(newarray,array,size*sizeof(TYPE));
            delete [] array;
            array=newarray;
        }
        catch (std::bad_alloc&)
        {
            /*-------temp_remove-------
            Error::abort("reallocation "
            "failure in file  %s, function %s, "
            "line: %d for variable: %s",file,function,line_no,name);
             */
        }
        return array;
    }
}
/*--------------------------------------------
 grow 1d vector
 --------------------------------------------*/
template <typename T>
void Memory::grow(T*& v,size_t old_sz,size_t new_sz)
{
    T* v_=new T[new_sz];
    memcpy(v_,v,old_sz*sizeof(T));
    delete [] v;
    v=v_;
}
#endif
