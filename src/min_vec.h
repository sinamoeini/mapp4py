#ifndef __MAPP__min_vec__
#define __MAPP__min_vec__
#include "atoms.h"
#include "xmath.h"
#include "memory.h"
namespace MAPP_NS
{
    
    
    template<class T,const int N,class E>
    class VecTensExpr
    {
    private:
    protected:
    public:
        int size(const int i) const
        {return static_cast<const E&>(*this).size(i);}
        template<const int I>
        T operator[](int i) const {return static_cast<const E&>(*this).operator[]<I>(i);}
        T operator()(int i,int j)const{return static_cast<const E&>(*this)(i,j);}
        
        T get(int I,int i) const {return static_cast<const E&>(*this).get(I,i);}
    };
    
    template<class T,const int N>
    class VecTens:public VecTensExpr<T,N,VecTens<T,N>>
    {
    private:
    protected:
    public:
        Vec<T>* vecs[N];
        bool vecs_alloc[N];
        
        T (*A)[__dim__];
        bool A_alloc;
        
        bool box_chng;
        VecTens();
        ~VecTens();
        template<class ...Ts>
        VecTens(Atoms*,bool,Ts...);
        template<class ...Ts>
        VecTens(Atoms*,bool,T(&)[__dim__][__dim__],Ts...);
        
        
        VecTens(const VecTens&);
        VecTens(VecTens&&);
        VecTens& operator=(VecTens&&);
        
        template<class... Ts>
        void assign(Atoms* atoms,int i,Vec<T>* v,Ts... vs )
        {
            assign(atoms,i,v);
            assign(atoms,i+1,vs...);
        }
        void assign(Atoms* atoms,int i,Vec<T>* v)
        {
            vecs[i]=v;
            vecs_alloc[i]=false;
        }
        
        template<class... Ts>
        void assign(Atoms* atoms,int i,int d,Ts... vs)
        {
            assign(atoms,i,d);
            assign(atoms,i+1,vs...);
        }
        void assign(Atoms* atoms,int i,int d)
        {
            vecs[i]=new Vec<T>(atoms,d);
            vecs_alloc[i]=true;
        }
        
        
        int size(const int i) const
        { return vecs[i]->atoms->natms_lcl*vecs[i]->dim;}
        
        template<const int I>
        T& operator[](int i){return vecs[I]->begin()[i];}
        template<const int I>
        T operator[](int i)const {return vecs[I]->begin()[i];}
        T operator()(int i,int j)const{return A[i][j];}
        T& operator()(int i,int j){return A[i][j];}
        
        T get(int I,int i)const{return vecs[I]->begin()[i];}
        T& get(int I,int i){return vecs[I]->begin()[i];}
        
        VecTens& operator*=(const T&);
        void cyclic_shift(int);
        
        T operator*(const VecTens&);
        
        template<class E>
        VecTens& operator+=(const VecTensExpr<T,N,E>&);
        VecTens& operator+=(const VecTens&);
        template<class E>
        VecTens& operator-=(const VecTensExpr<T,N,E>&);
        VecTens& operator-=(const VecTens&);
        template<class E>
        VecTens& operator=(const VecTensExpr<T,N,E>&);
        VecTens& operator=(const VecTens&);
        
    };
    
    template <class T,const int N,class E0,class E1>
    class VecTensSum:public VecTensExpr<T,N,VecTensSum<T,N,E0,E1>>
    {
        E0 const& v0;
        E1 const& v1;
    public:
        VecTensSum(E0 const& __v0, E1 const& __v1):v0(__v0),v1(__v1){}
        template<const int I>
        T operator[](int i)const{return v0.operator[]<I>(i)+v1.operator[]<I>(i);}
        T get(int I,int i) const {return v0.get(I,i)+v1.get(I,i);}
        T operator()(int i,int j)const{return v0(i,j)+v1(i,j);}
        int size(const int i) const
        {return v0.size(i);}
    };
    
    template <class T,const int N,class E0,class E1>
    class VecTensSub:public VecTensExpr<T,N,VecTensSub<T,N,E0,E1>>
    {
        E0 const& v0;
        E1 const& v1;
    public:
        VecTensSub(E0 const& __v0, E1 const& __v1):v0(__v0),v1(__v1){}
        template<const int I>
        T operator[](int i)const{return v0.operator[]<I>(i)-v1.operator[]<I>(i);}
        T get(int I,int i) const {return v0.get(I,i)-v1.get(I,i);}
        T operator()(int i,int j)const{return v0(i,j)-v1(i,j);}
        int size(const int i) const
        {return v0.size(i);}
    };
    
    template <class T,const int N,class E>
    class VecTensMulScl:public VecTensExpr<T,N,VecTensMulScl<T,N,E>>
    {
        E const& v;
        T const& scl;
    public:
        VecTensMulScl(T const& __scl, E const& __v):scl(__scl),v(__v){}
        template<const int I>
        T operator[](int i)const{return v.operator[]<I>(i)*scl;}
        T get(int I,int i) const {return v.get(I,i)*scl;}
        T operator()(int i,int j)const{return v(i,j)*scl;}
        int size(const int i) const
        {return v.size(i);}
    };
    
    
    template<class T,const int N,class ...Ts,template<class,const int,class...> class E>
    VecTensMulScl<T,N,E<T,N,Ts...>> operator*(T const& u,const E<T,N,Ts...>& v)
    {
        return VecTensMulScl<T,N,E<T,N,Ts...>>(u, v);
    }
    
    template<class T,const int N,class ...Ts0,class ...Ts1,template<class,const int,class...> class E0,template<class,const int,class...> class E1>
    VecTensSum<T,N,E0<T,N,Ts0...>,E1<T,N,Ts1...>> operator+(const E0<T,N,Ts0...>& v0,const E1<T,N,Ts1...>& v1)
    {
        return VecTensSum<T,N,E0<T,N,Ts0...>,E1<T,N,Ts1...>>(v0,v1);
    }
    
    template<class T,const int N,class ...Ts0,class ...Ts1,template<class,const int,class...> class E0,template<class,const int,class...> class E1>
    VecTensSub<T,N,E0<T,N,Ts0...>,E1<T,N,Ts1...>> operator-(const E0<T,N,Ts0...>& v0,const E1<T,N,Ts1...>& v1)
    {
        return VecTensSub<T,N,E0<T,N,Ts0...>,E1<T,N,Ts1...>>(v0,v1);
    }
    
    
    template<const int N,const int I=N-1>
    class VecOp
    {
    public:
        template<class T,class E,class F>
        static inline void func(VecTens<T,N>& v0,const VecTensExpr<T,N,E>& v1,F f)
        {
            int size=v0.size(I);
            T* __v0=v0.vecs[I]->begin();
            for(int i=0;i<size;i++)
                f(__v0[i],v1.get(I,i));
          
            VecOp<N,I-1>::func(v0,v1,f);
        }
    };
    template<const int N>
    class VecOp<N,0>
    {
    public:
        template<class T,class E,class F>
        static inline void func(VecTens<T,N>& v0,const VecTensExpr<T,N,E>& v1,F f)
        {
            int size=v0.size(0);
            T* __v0=v0.vecs[0]->begin();
            for(int i=0;i<size;i++)
            {
                f(__v0[i],v1.get(0,i));
            }
        }
    };
}

using namespace MAPP_NS;
/*--------------------------------------------
 
 --------------------------------------------*/
template<typename T,const int N>
VecTens<T,N>::VecTens():
A(NULL),
A_alloc(false)
{
    Algebra::Do<N>::func([this](int i)
    {vecs[i]=NULL; vecs_alloc[i]=false;});
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<typename T,const int N>
VecTens<T,N>::~VecTens()
{
    Algebra::Do<N>::func([this](int i)
    {
        if(vecs_alloc[i]) delete vecs[i];
        vecs[i]=NULL;
    });
    if(A_alloc) delete [] A;
    A_alloc=false;
    A=NULL;
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<typename T,const int N>template<class... Ts>
VecTens<T,N>::VecTens(Atoms* atoms,bool __box_chng,Ts... vs):
box_chng(__box_chng)
{
    static_assert(sizeof...(Ts)==N,"not enough vecs");
    A_alloc=false;
    if(__box_chng)
    {
        A=new T[__dim__][__dim__];
        A_alloc=true;
    }
    assign(atoms,0,vs...);
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<typename T,const int N>template<class... Ts>
VecTens<T,N>::VecTens(Atoms* atoms,bool __box_chng,T(&__A)[__dim__][__dim__],Ts... vs):
box_chng(__box_chng),
A_alloc(false),
A(__A)
{
    static_assert(sizeof...(Ts)==N,"not enough vecs");
    assign(atoms,0,vs...);
}
    
/*--------------------------------------------
 
 --------------------------------------------*/
template<typename T,const int N>
VecTens<T,N>::VecTens(const VecTens<T,N>& r):
box_chng(r.box_chng)
{
    Algebra::Do<N>::func([this,&r](int i)
    {
        vecs_alloc[i]=true;
        vecs[i]=new Vec<T>(r.vecs[i]->atoms,r.vecs[i]->dim);
        memcpy(vecs[i]->begin(),r.vecs[i]->begin(),vecs[i]->dim*r.vecs[i]->atoms->natms_lcl*sizeof(T));
    });
    A_alloc=false;
    if(box_chng)
    {
        A_alloc=true;
        A=new T[__dim__][__dim__];
        memcpy(&A[0][0],&r.A[0][0],__dim__*__dim__*sizeof(T));
    }
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<typename T,const int N>
VecTens<T,N>::VecTens(VecTens<T,N>&& r):
box_chng(r.box_chng)
{
    Algebra::Do<N>::func([this,&r](int i)
    {
        vecs_alloc[i]=r.vecs_alloc[i];
        r.vecs_alloc[i]=false;
        vecs[i]=r.vecs[i];
        r.vecs[i]=NULL;
    });
    
    A_alloc=r.A_alloc;
    r.A_alloc=false;
    A=r.A;
    r.A=NULL;
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<typename T,const int N>
VecTens<T,N>& VecTens<T,N>::operator=(VecTens<T,N>&& r)
{
    box_chng=r.box_chng;
    Algebra::Do<N>::func([this,&r](int i)
    {
        vecs_alloc[i]=r.vecs_alloc[i];
        r.vecs_alloc[i]=false;
        vecs[i]=r.vecs[i];
        r.vecs[i]=NULL;
    });
    
    A_alloc=r.A_alloc;
    r.A_alloc=false;
    A=r.A;
    r.A=NULL;
    return *this;
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<typename T,const int N>
T VecTens<T,N>::operator*(const VecTens<T,N>& r)
{
    T ans_lcl=0.0;
    const int natms_lcl=vecs[0]->atoms->natms_lcl;
    Algebra::Do<N>::func([this,&ans_lcl,&natms_lcl,&r](int j)
    {
        const int dim=vecs[j]->dim;
        T* v0=vecs[j]->begin();
        T* v1=r.vecs[j]->begin();
        for(int i=0;i<natms_lcl*dim;i++)
            ans_lcl+=v0[i]*v1[i];
    });
    
    T ans;
    MPI_Allreduce(&ans_lcl,&ans,1,Vec<T>::MPI_T,MPI_SUM,vecs[0]->atoms->world);
    if(box_chng)
        Algebra::DoLT<__dim__>::func([this,&r,&ans](int i,int j){ans+=this->A[i][j]*r.A[i][j];});
    return ans;
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<typename T,const int N>
VecTens<T,N>& VecTens<T,N>::operator*=(const T& scl)
{
    
    const int natms_lcl=vecs[0]->atoms->natms_lcl;
    Algebra::Do<N>::func([this,&scl,&natms_lcl](int j)
    {
        const int dim=vecs[j]->dim;
        T* v0=vecs[j]->begin();
        for(int i=0;i<natms_lcl*dim;i++)
            v0[i]*=scl;
    });
    

    if(box_chng)
        Algebra::DoLT<__dim__>::func([this,&scl](int i,int j){this->A[i][j]*=scl;});
    return *this;
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<typename T,const int N>template<class E>
VecTens<T,N>& VecTens<T,N>::operator+=(const VecTensExpr<T,N,E>& expr)
{
    VecOp<N>::func(*this,expr,[](T& l,const T& r){l+=r;});
    if(box_chng)
        Algebra::DoLT<__dim__>::func([this,&expr](int i,int j){this->A[i][j]+=expr(i,j);});
    return *this;
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<typename T,const int N>
VecTens<T,N>& VecTens<T,N>::operator+=(const VecTens<T,N>& r)
{
    const int natms_lcl=vecs[0]->atoms->natms_lcl;
    Algebra::Do<N>::func([this,&natms_lcl,&r](int j)
    {
        const int dim=vecs[j]->dim;
        T* v0=vecs[j]->begin();
        const T* v1=r.vecs[j]->begin();
        for(int i=0;i<natms_lcl*dim;i++)
            v0[i]+=v1[i];
    });
    if(box_chng)
        Algebra::DoLT<__dim__>::func([this,&r](int i,int j){this->A[i][j]+=r.A[i][j];});
    return *this;
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<typename T,const int N>template<class E>
VecTens<T,N>& VecTens<T,N>::operator-=(const VecTensExpr<T,N,E>& expr)
{
    VecOp<N>::func(*this,expr,[](T& l,const T& r){l-=r;});
    if(box_chng)
        Algebra::DoLT<__dim__>::func([this,&expr](int i,int j){this->A[i][j]-=expr(i,j);});
    return *this;
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<typename T,const int N>
VecTens<T,N>& VecTens<T,N>::operator-=(const VecTens<T,N>& r)
{
    const int natms_lcl=vecs[0]->atoms->natms_lcl;
    Algebra::Do<N>::func([this,&natms_lcl,&r](int j)
    {
        const int dim=vecs[j]->dim;
        T* v0=vecs[j]->begin();
        const T* v1=r.vecs[j]->begin();
        for(int i=0;i<natms_lcl*dim;i++)
            v0[i]-=v1[i];
    });
    if(box_chng)
        Algebra::DoLT<__dim__>::func([this,&r](int i,int j){this->A[i][j]-=r.A[i][j];});
    return *this;
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<typename T,const int N>template<class E>
VecTens<T,N>& VecTens<T,N>::operator=(const VecTensExpr<T,N,E>& expr)
{
    VecOp<N>::func(*this,expr,[](T& l,const T& r){l=r;});
    if(box_chng)
        Algebra::DoLT<__dim__>::func([this,&expr](int i,int j){this->A[i][j]=expr(i,j);});
    return *this;
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<typename T,const int N>
VecTens<T,N>& VecTens<T,N>::operator=(const VecTens<T,N>& r)
{
    const int natms_lcl=vecs[0]->atoms->natms_lcl;
    Algebra::Do<N>::func([this,&natms_lcl,&r](int j)
    {
        const int dim=vecs[j]->dim;
        T* v0=vecs[j]->begin();
        const T* v1=r.vecs[j]->begin();
        for(int i=0;i<natms_lcl*dim;i++)
            v0[i]=v1[i];
    });
    if(box_chng)
        Algebra::DoLT<__dim__>::func([this,&r](int i,int j){this->A[i][j]=r.A[i][j];});
    return *this;
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<typename T,const int N>
void VecTens<T,N>::cyclic_shift(int n)
{
    Algebra::Do<N>::func([this,&n](int j)
    {
        Vec<T>* __v=(this+n-1)->vecs[j];
        for(int i=n-1;i>0;i--)
            (this+i)->vecs[j]=(this+i-1)->vecs[j];
        
        this->vecs[j]=__v;
    });
    if(box_chng)
    {
        T(*__A)[__dim__]=(this+n-1)->A;
        for(int i=n-1;i>0;i--)
            (this+i)->A=(this+i-1)->A;
        this->A=__A;
    }
}
/*--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
     
 --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------*/
namespace MAPP_NS
{
    template<class T,class E>
    class __VecTensExpr
    {
    private:
    protected:
    public:
        E& operator()(){return static_cast<E&>(*this); }
        const E& operator()() const{ return static_cast<const E&>(*this);}
        T operator[](int i) const {return static_cast<const E&>(*this)[i];}
        T operator()(int i,int j)const{return static_cast<const E&>(*this)(i,j);}
    };
    
    template<class T>
    class __VecTens: public __VecTensExpr<T,__VecTens<T>>
    {
    private:
    protected:
    public:
        Vec<T>* vec;
        T (*A)[__dim__];
        bool box_chng;
        unsigned int alloc_flag;
        
        __VecTens();
        __VecTens(Atoms*,bool);
        __VecTens(Atoms*,Vec<T>*&,bool);
        __VecTens(Atoms*,Vec<T>*& ,T(&)[__dim__][__dim__],bool);
        __VecTens(__VecTens&&);
        __VecTens(const __VecTens&);
        ~__VecTens();
    
        Vec<T>* operator()();
        T operator*(const __VecTens&);
        __VecTens& operator*= (T&);
        __VecTens& operator=(__VecTens&&);
        __VecTens& operator=(const __VecTens&);
        __VecTens& operator-=(const __VecTens&);
        __VecTens& operator+=(const __VecTens&);
        template <class E>
        __VecTens& operator=(const __VecTensExpr<T,E>&&);
        template <class E>
        __VecTens& operator-=(const __VecTensExpr<T,E>&&);
        template <class E>
        __VecTens& operator+=(const __VecTensExpr<T,E>&&);
        T operator[](int i) const {return vec->begin()[i];}
        T& operator[](int i){return vec->begin()[i];}
        T operator()(int i,int j)const{return A[i][j];}
        T& operator()(int i,int j){return A[i][j];}
    };
    
    
    
    
    
    
    template <class T,class E0,class E1>
    class __VecTensSum:public __VecTensExpr<T,__VecTensSum<T,E0,E1>>
    {
        E0 const& v0;
        E1 const& v1;
    public:
        __VecTensSum(E0 const& __v0, E1 const& __v1):v0(__v0),v1(__v1)
        {}
        T operator[](int i)const{return v0[i]+v1[i];}
        T operator()(int i,int j)const{return v0(i,j)+v1(i,j);}
    };
    template <class T,class E0,class E1>
    class __VecTensSub:public __VecTensExpr<T,__VecTensSub<T,E0,E1>>
    {
        E0 const& v0;
        E1 const& v1;
    public:
        __VecTensSub(E0 const& __v0, E1 const& __v1):v0(__v0),v1(__v1)
        {}
        T operator[](int i)const{return v0[i]-v1[i];}
        T operator()(int i,int j)const{return v0(i,j)-v1(i,j);}
    };
    
    template <class T,class E>
    class __VecTensMulScl : public __VecTensExpr<T,__VecTensMulScl<T,E>>
    {
        E const& v;
        T const& scl;
    public:
        __VecTensMulScl(T const& __scl, E const& __v):scl(__scl),v(__v)
        {}
        T operator[](int i)const{return scl*v[i];}
        T operator()(int i,int j)const{return scl*v(i,j);}
    };
    
    template<class T,class ...Ts,template<class,class...> class E>
    __VecTensMulScl<T,E<T,Ts...>> operator*(T const& u,const E<T,Ts...>& v)
    {
        return __VecTensMulScl<T,E<T,Ts...>>(u, v);
    }
    
    template<class T,class ...Ts0,class ...Ts1,template<class,class...> class E0,template<class,class...> class E1>
    __VecTensSum<T,E0<T,Ts0...>,E1<T,Ts1...>> operator+(const E0<T,Ts0...>& v0,const E1<T,Ts1...>& v1)
    {
        return __VecTensSum<T,E0<T,Ts0...>,E1<T,Ts1...>>(v0,v1);
    }
    
    template<class T,class ...Ts0,class ...Ts1,template<class,class...> class E0,template<class,class...> class E1>
    __VecTensSub<T,E0<T,Ts0...>,E1<T,Ts1...>> operator-(const E0<T,Ts0...>& v0,const E1<T,Ts1...>& v1)
    {
        return __VecTensSub<T,E0<T,Ts0...>,E1<T,Ts1...>>(v0,v1);
    }
    
}
using namespace MAPP_NS;
/*--------------------------------------------
 
 --------------------------------------------*/
template<typename T>
__VecTens<T>::__VecTens():
vec(NULL),
A(NULL),
box_chng(false),
alloc_flag(0)
{
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<typename T>
__VecTens<T>::~__VecTens()
{
    if(alloc_flag & 1)
        delete vec;
    if(alloc_flag & 2)
        delete [] A;
    
    vec=NULL;
    A=NULL;
    alloc_flag=0;
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<typename T>
__VecTens<T>::__VecTens(Atoms* __atoms,Vec<T>*& __vec,T(&__A)[__dim__][__dim__],bool __box_chng):
vec(__vec),
A(__A),
box_chng(__box_chng),
alloc_flag(0)
{
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<typename T>
__VecTens<T>::__VecTens(Atoms* __atoms,Vec<T>*& __vec,bool __box_chng):
vec(__vec),
A(__box_chng?new T[__dim__][__dim__]:NULL),
box_chng(__box_chng),
alloc_flag(__box_chng?2:0)
{
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<typename T>
__VecTens<T>::__VecTens(Atoms* __atoms,bool __box_chng):
vec(new Vec<T>(__atoms,__atoms->x->dim)),
A(__box_chng ? new T[__dim__][__dim__]:NULL),
box_chng(__box_chng),
alloc_flag(__box_chng?3:1)
{
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<typename T>
__VecTens<T>::__VecTens(__VecTens&& other):
vec(other.vec),
box_chng(other.box_chng),
A(other.A),
alloc_flag(other.alloc_flag)
{
    other.A=NULL;
    other.vec=NULL;
    other.alloc_flag=0;
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<typename T>
__VecTens<T>::__VecTens(const __VecTens& other):
vec(other.vec),
box_chng(other.box_chng),
A(other.A),
alloc_flag(other.alloc_flag)
{
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<typename T>
Vec<T>* __VecTens<T>::operator()()
{
    return vec;
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<typename T>
T __VecTens<T>::operator*(const __VecTens<T>& rhs)
{
    T ans_lcl=0.0,ans;
    T* vec0=this->vec->begin();
    T* vec1=rhs.vec->begin();
    const int n=this->vec->atoms->natms_lcl*this->vec->dim;
    for(int i=0;i<n;i++)
        ans_lcl+=vec0[i]*vec1[i];
    
    MPI_Allreduce(&ans_lcl,&ans,1,Vec<T>::MPI_T,MPI_SUM,this->vec->atoms->world);
    if(!box_chng) return ans;
    Algebra::DoLT<__dim__>::func([this,&rhs,&ans](int i,int j)
    {ans+=this->A[i][j]*rhs.A[i][j];});
    return ans;
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<typename T>
__VecTens<T>& __VecTens<T>::operator*=(T& a)
{
    T* vec0=this->vec->begin();
    const int n=this->vec->atoms->natms_lcl*this->vec->dim;
    for(int i=0;i<n;i++)
        vec0[i]*=a;
    if(!box_chng) return *this;
    Algebra::DoLT<__dim__>::func([this,&a](int i,int j)
    {this->A[i][j]*=a;});
    return *this;
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<typename T>
__VecTens<T>& __VecTens<T>::operator=(__VecTens&& other)
{
    this->vec=other.vec;
    other.vec=NULL;
    this->A=other.A;
    other.A=NULL;
    this->alloc_flag=other.alloc_flag;
    other.alloc_flag=0;
    this->box_chng=other.box_chng;
    other.box_chng=false;
    return *this;
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<typename T>
__VecTens<T>& __VecTens<T>::operator=(const __VecTens& other)
{
    memcpy(this->vec->begin(),other.vec->begin(),this->vec->atoms->natms_lcl*this->vec->dim*sizeof(T));
    if(!box_chng) return *this;
    Algebra::DoLT<__dim__>::func([this,&other](int i,int j)
    {this->A[i][j]=other.A[i][j];});
    return *this;
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<typename T>
__VecTens<T>& __VecTens<T>::operator-=(const __VecTens& other)
{
    T* vec0=this->vec->begin();
    T* vec1=other.vec->begin();
    const int n=this->vec->atoms->natms_lcl*this->vec->dim;
    for(int i=0;i<n;i++) vec0[i]-=vec1[i];
    
    if(!box_chng) return *this;
    T (&__A)[__dim__][__dim__]=other.vec_tens.A;
    Algebra::DoLT<__dim__>::func([this,__A](int i,int j)
    {this->A[i][j]-=__A[i][j];});
    return *this;
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<typename T>
__VecTens<T>& __VecTens<T>::operator+=(const __VecTens& other)
{
    T* vec0=this->vec->begin();
    T* vec1=other.vec->begin();
    const int n=this->vec->atoms->natms_lcl*this->vec->dim;
    for(int i=0;i<n;i++) vec0[i]+=vec1[i];
    
    if(!box_chng) return *this;
    T (&__A)[__dim__][__dim__]=other.vec_tens.A;
    Algebra::DoLT<__dim__>::func([this,__A](int i,int j)
    {this->A[i][j]+=__A[i][j];});
    return *this;
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<typename T>template<class E>
__VecTens<T>& __VecTens<T>::operator=(const __VecTensExpr<T,E>&& other)
{
    T* vec0=this->vec->begin();
    const int n=this->vec->atoms->natms_lcl*this->vec->dim;
    for(int i=0;i<n;i++)
        vec0[i]=other[i];
    if(!box_chng) return *this;
    Algebra::DoLT<__dim__>::func([this,&other](int i,int j)
    {this->A[i][j]=other(i,j);});
    return *this;
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<typename T>template<class E>
__VecTens<T>& __VecTens<T>::operator-=(const __VecTensExpr<T,E>&& other)
{
    T* vec0=this->vec->begin();
    const int n=this->vec->atoms->natms_lcl*this->vec->dim;
    for(int i=0;i<n;i++)
        vec0[i]-=other[i];
    if(!box_chng) return *this;
    Algebra::DoLT<__dim__>::func([this,&other](int i,int j)
    {this->A[i][j]-=other(i,j);});
    return *this;
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<typename T>template<class E>
__VecTens<T>& __VecTens<T>::operator+=(const __VecTensExpr<T,E>&& other)
{
    T* vec0=this->vec->begin();
    const int n=this->vec->atoms->natms_lcl*this->vec->dim;
    for(int i=0;i<n;i++)
        vec0[i]+=other[i];
    if(!box_chng) return *this;
    Algebra::DoLT<__dim__>::func([this,&other](int i,int j)
    {this->A[i][j]+=other(i,j);});
    return *this;
}



/*--------------------------------------------
 
 --------------------------------------------*/
namespace MAPP_NS
{
    template<class V>
    class __GMRES__
    {
    private:
        const int m;
        V* Q;
        type0** A_hat;
        type0* Ax_hat;
        type0(* cos_sin)[2];
        
        
        type0* x_hat;
        
        
        type0 calc(V&,V&);
        type0 calc(int,V&,V&);
        type0 solve_y(int,V&);
        
    protected:
    public:
        template<class ... Cs>
        __GMRES__(int,Cs ...);
        ~__GMRES__();
        template<class KERNEL>
        bool solve(KERNEL&,V&,type0,type0&,V&);
        
    };
}

/*--------------------------------------------
 
 --------------------------------------------*/
template<class V>template<class ... Cs>
__GMRES__<V>::__GMRES__(int __m,Cs ... cs):
m(__m)
{
    
    
    A_hat=new type0*[m];
    *A_hat=new type0[m*(m+1)/2];
    for(int i=1;i<m;i++)
        A_hat[i]=A_hat[i-1]+i;
    
    Memory::alloc(Ax_hat,m+1);
    Memory::alloc(cos_sin,m+1);
    Memory::alloc(x_hat,m+1);
    
    Q=new V[m+1];
    for(int i=0;i<m+1;i++)
    {
        Q[i].~V();
        new (Q+i) V(cs...);
    }
    
    
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<class V>
__GMRES__<V>::~__GMRES__()
{
    
    delete [] Q;

    Memory::dealloc(x_hat);
    Memory::dealloc(cos_sin);
    Memory::dealloc(Ax_hat);
    delete [] *A_hat;
    delete [] A_hat;
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<class V>
type0 __GMRES__<V>::calc(V& Ax,V& x)
{
    type0 norm=sqrt(Ax*Ax);
    type0 norm_inv=1.0/norm;

    Q[0]=x=norm_inv*Ax;
    
    Ax_hat[0]=norm;
    //printf("\t\t %e\n",fabs(Ax_hat[0]));
    return norm;
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<class V>
type0 __GMRES__<V>::calc(int iter,V& Ax,V& x)
{
    type0* h=A_hat[iter];
    int ivec=iter+1;
    
    
    for(int i=0;i<ivec;i++)
        h[i]=Ax*Q[i];
    
    Q[ivec]=Ax;
    
    for(int i=0;i<ivec;i++)
        Q[ivec]-=h[i]*Q[i];
    
    type0 norm=sqrt(Q[ivec]*Q[ivec]);
    type0 norm_inv=1.0/norm;
    
    Q[ivec]*=norm_inv;
    x=Q[ivec];
    
    
    type0 tmp;
    for(int i=0;i<iter;i++)
    {
        tmp=cos_sin[i][0]*h[i]-cos_sin[i][1]*h[i+1];
        h[i+1]=cos_sin[i][1]*h[i]+cos_sin[i][0]*h[i+1];
        h[i]=tmp;
    }
    
    tmp=sqrt(h[iter]*h[iter]+norm*norm);
    cos_sin[iter][0]=h[iter]/tmp;
    cos_sin[iter][1]=-norm/tmp;
   
    h[iter]=cos_sin[iter][0]*h[iter]-cos_sin[iter][1]*norm;
    Ax_hat[iter+1]=cos_sin[iter][1]*Ax_hat[iter];
    Ax_hat[iter]*=cos_sin[iter][0];
    //printf("\t\t %e\n",fabs(Ax_hat[iter+1]));
    return fabs(Ax_hat[iter+1]);
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<class V>
type0 __GMRES__<V>::solve_y(int nvecs,V& x)
{
    for(int i=nvecs-1;i>-1;i--)
    {
        x_hat[i]=Ax_hat[i];
        for(int j=i+1;j<nvecs;j++)
            x_hat[i]-=A_hat[j][i]*x_hat[j];
        x_hat[i]/=A_hat[i][i];
    }

    x=x_hat[0]*Q[0];
    for(int ivec=1;ivec<nvecs;ivec++)
        x+=x_hat[ivec]*Q[ivec];

    type0 norm=0.0;
    for(int ivec=0;ivec<nvecs;ivec++)
        norm+=x_hat[ivec]*x_hat[ivec];
    
    return sqrt(norm);
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<class V>template<class KERNEL>
bool __GMRES__<V>::solve(KERNEL& A,V& Ax,type0 tol,type0& norm,V& x)
{
    calc(Ax,x);
    for(int i=0;i<m;i++)
    {
        A(x,Ax);
        if(calc(i,Ax,x)<tol)
        {
            norm=solve_y(i+1,x);
            return true;
        }
    }
    
    norm=solve_y(m,x);
    return false;
}


#endif
