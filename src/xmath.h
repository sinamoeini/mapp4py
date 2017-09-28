#ifndef __MAPP__xmath__
#define __MAPP__xmath__
#include <utility>
#include <math.h>
#include <limits>
#include <iostream>
#include "global.h"
using namespace MAPP_NS;
namespace MAPP_NS
{
    namespace XMath
    {
        // recursive form of above functions that feeds into it
        void fac_list_rec(int,int,int,int*&,int*&,int&);
        template<typename,class>
        class XVec;
        // inverse square matrix
        void invert(type0**,type0**,int);
        // inverse lower triangle square matrix
        template<const int d>
        //void invert_lower_triangle(type0 (&)[d][d],type0 (&)[d][d]);
        void invert_lower_triangle(type0**&,type0**&,int&);
        // return the list of all possible groups of integers that their products are equal to specific number
        void fac_list(int,int,int*&,int&);
        
        void square2lo_tri(type0(&)[3][3],type0(&)[3][3]);
        void quadrature_lg(int,type0*,type0*);
        void quadrature_hg(int,type0*,type0*);
        template<const int N>
        void quad_hg(type0(&)[N],type0(&)[N]);
        template<const int N>
        void quad_lg(type0(&)[N],type0(&)[N]);
        bool Msqrt(type0(&)[3][3],type0(&)[3][3]);
        bool Msqrt(type0(&)[2][2],type0(&)[2][2]);
        
        template<typename T0,class COMP,class SWAP>
        void quicksort(T0,T0,COMP,SWAP);
        template<typename T0,class C0>
        void srch_lst_lst(T0*,int,C0*,T0*,int,C0*);
    };
}
/*------------------------------------------------------------------------------------------------------------------------------------
 
 ------------------------------------------------------------------------------------------------------------------------------------*/
template<const int N>
void XMath::quad_hg(type0 (&x)[N],type0 (&w)[N])
{
    constexpr int M=N/2+1;
    const int max_iter=50;
    constexpr type0 tol=std::numeric_limits<type0>::epsilon();
    
    
    
    
    
    type0 p_1[M]{[0 ... M-1]=0.0};
    type0 p_2[M]{[0 ... M-1]=0.0};
    type0 p_coef[M]{[0 ... M-1]=0.0};
    type0 dp_coef[M]{[0 ... M-1]=0.0};
    
    
    p_1[0]=1.0;
    p_coef[0]=1.0;
    type0 ii;
    for(int m,i=0;i<N;i++)
    {
        ii=static_cast<type0>(i+1);
        m=(i+3)/2;
        if(i%2)
        {
            for(int j=0;j<m;j++)
                p_coef[j]=-2.0*(ii-1.0)*p_2[j];
            for(int j=1;j<m;j++)
                p_coef[j]+=2.0*p_1[j-1];
        }
        else
        {
            for(int j=0;j<m;j++)
                p_coef[j]=2.0*p_1[j]-2.0*(ii-1.0)*p_2[j];
        }
        
        for(int j=0;j<m;j++)
        {
            p_2[j]=p_1[j];
            p_1[j]=p_coef[j];
        }
    }
    
    
    for(int i=0;i<M-1;i++)
        dp_coef[i]=p_coef[i+1]*static_cast<type0>(i+1);
    
    int iter,ord,icurs;
    type0 a,u0,inv_u0,f,df,tmp0,tmp1,del_u0;
    
    ord=M;
    icurs=N/2;
    if(N%2==1)
    {
        x[icurs]=0.0;
        icurs++;
    }
    
    a=p_coef[M-1];
    
    for(int i=0;i<M;i++)
        p_coef[i]/=a;
    for(int i=0;i<M-1;i++)
        dp_coef[i]/=a;
    
    
    
    while(ord>1)
    {
        u0=0.0;
        f=1.0;
        del_u0=0.0;
        iter=max_iter;
        while(fabs(f)>tol && iter)
        {
            u0+=del_u0;
            df=f=0.0;
            tmp0=1.0;
            for(int j=0;j<ord-1;j++)
            {
                f+=p_coef[j]*tmp0;
                df+=dp_coef[j]*tmp0;
                tmp0*=u0;
            }
            f+=p_coef[ord-1]*tmp0;
            del_u0=-f/df;
            iter--;
        }
        
        x[icurs]=sqrt(u0);
        
        tmp0=1.0;
        tmp1=0.0;
        for(int i=0;i<M;i++)
        {
            tmp1+=p_2[i]*tmp0;
            tmp0*=u0;
        }
        
        if(N%2==0)
            tmp1*=sqrt(u0);
        
        w[icurs]=1.0/(tmp1*tmp1);
        
        
        inv_u0=1.0/u0;
        p_coef[0]*=-inv_u0;
        for(int i=1;i<ord-1;i++)
        {
            p_coef[i]*=-inv_u0;
            p_coef[i]+=inv_u0*p_coef[i-1];
            dp_coef[i-1]=p_coef[i]*static_cast<type0>(i);
        }
        
        ord--;
        icurs++;
    }
    
    tmp0=sqrt(M_PI)/static_cast<type0>(N);
    for(int i=0;i<N-1;i++)
        tmp0*=static_cast<type0>(2*(i+1));
    
    if(N%2==1)
        w[N/2]=tmp0/(p_2[0]*p_2[0]);
    
    
    for(int i=0;i<N/2;i++)
    {
        icurs--;
        w[icurs]*=tmp0;
        w[i]=w[icurs];
        x[i]=-x[icurs];
    }
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<const int N>
void XMath::quad_lg(type0 (&x)[N],type0 (&w)[N])
{
    constexpr int M=N/2+1;
    constexpr int max_iter=50;
    constexpr type0 tol=std::numeric_limits<type0>::epsilon();
    
    
    int iter,ord,icurs,jcurs;
    type0 a,u0,inv_u0,f,up,df,tmp0;
    type0 ii,jj,del_u0;
    type0 p_coef[M]{[0 ... M-1]=0.0};
    type0 dp_coef[M]{[0 ... M-1]=0.0};
    
    up=0.0;
    
    p_coef[0]=1.0;
    
    
    
    for(int __m,i=0;i<N;i++)
    {
        ii=static_cast<type0>(i+1);
        if(i%2)
        {
            __m=(i+3)/2;
            p_coef[__m-1]=(2.0-1.0/ii)*p_coef[__m-2];
            for(int j=__m-2;j>0;j--)
            {
                jj=static_cast<type0>(j);
                p_coef[j]*=-(2.0*jj+1)/ii;
                p_coef[j]+=(1.0+(2*jj-1.0)/ii)*p_coef[j-1];
            }
            p_coef[0]*=-1.0/ii;
        }
        else
        {
            __m=i/2+1;
            for(int j=0;j<__m-1;j++)
            {
                jj=static_cast<type0>(j);
                p_coef[j]*=1.0+2.0*jj/ii;
                p_coef[j]-=2.0*(jj+1.0)*p_coef[j+1]/ii;
            }
            //2.0-1.0/ii=(1+2*(m-1)/i)
            p_coef[__m-1]*=2.0-1.0/ii;
        }
    }
    
    
    
    for(int i=0;i<M-1;i++)
        dp_coef[i]=p_coef[i+1]*static_cast<type0>(i+1);
    
    ord=M;
    icurs=N-1;
    a=p_coef[M-1];
    
    for(int i=0;i<M;i++)
        p_coef[i]/=a;
    for(int i=0;i<M-1;i++)
        dp_coef[i]/=a;
    
    while(ord>1)
    {
        u0=1.0;
        f=1.0;
        del_u0=0.0;
        iter=max_iter;
        while(fabs(f)>tol && iter)
        {
            u0+=del_u0;
            df=f=0.0;
            tmp0=1.0;
            for(int j=0;j<ord-1;j++)
            {
                f+=p_coef[j]*tmp0;
                df+=dp_coef[j]*tmp0;
                tmp0*=u0;
            }
            f+=p_coef[ord-1]*tmp0;
            del_u0=-f/df;
            iter--;
        }
        
        x[icurs]=sqrt(u0);
        
        inv_u0=1.0/u0;
        p_coef[0]*=-inv_u0;
        for(int i=1;i<ord-1;i++)
        {
            p_coef[i]*=-inv_u0;
            p_coef[i]+=inv_u0*p_coef[i-1];
            dp_coef[i-1]=p_coef[i]*static_cast<type0>(i);
        }
        
        ord--;
        icurs--;
    }
    
    
    if(N%2==0)
    {
        icurs++;
        jcurs=icurs-1;
        for(int i=icurs;i<N;i++)
        {
            tmp0=a;
            for(int j=icurs;j<N;j++)
            {
                if(i!=j)
                    tmp0*=x[i]*x[i]-x[j]*x[j];
                else
                    tmp0*=x[i]+x[j];
            }
            w[i]=2.0/(tmp0*tmp0*(1.0-x[i]*x[i]));
            w[jcurs]=w[i];
            x[jcurs]=-x[i];
            jcurs--;
        }
    }
    else
    {
        x[icurs]=0.0;
        icurs++;
        tmp0=a;
        for(int i=icurs;i<N;i++)
            tmp0*=-x[i]*x[i];
        w[icurs-1]=2.0/(tmp0*tmp0);
        
        jcurs=icurs-2;
        for(int i=icurs;i<N;i++)
        {
            tmp0=a*x[i];
            for(int j=icurs;j<N;j++)
            {
                if(i!=j)
                    tmp0*=x[i]*x[i]-x[j]*x[j];
                else
                    tmp0*=x[i]+x[j];
            }
            w[i]=2.0/(tmp0*tmp0*(1.0-x[i]*x[i]));
            w[jcurs]=w[i];
            x[jcurs]=-x[i];
            jcurs--;
        }
    }
}
/*------------------------------------------------------------------------------------------------------------------------------------
 
 ------------------------------------------------------------------------------------------------------------------------------------*/
namespace MAPP_NS
{
    template<typename T0,class C0>
    class XMath::XVec
    {
    private:
        T0* data;
        int* key;
        int size;
        C0& act;
    protected:
    public:
        
        class iterator
        {
        private:
            XVec* xvec;
            int* ikey;
        protected:
        public:
            iterator(XVec* _xvec,int* _ikey)
            {
                xvec=_xvec;
                ikey=_ikey;
            }
            
            iterator(const iterator& other)
            {
                this->xvec=other.xvec;
                this->ikey=other.ikey;
            }
            
            iterator(iterator&& other)
            {
                this->xvec=other.xvec;
                this->ikey=other.ikey;
            }
            
            iterator& operator =(const iterator& other)
            {
                this->xvec=other.xvec;
                this->ikey=other.ikey;
                return *this;
            }
            
            iterator& operator =(iterator&& other)
            {
                this->xvec=other.xvec;
                this->ikey=other.ikey;
                return *this;
            }
            
            iterator& operator +=(int i)
            {
                ikey+=i;
                return *this;
            }
            
            iterator& operator ++()
            {
                ikey++;
                return *this;
            }
            
            iterator operator ++(int)
            {
                iterator old=*this;
                ikey++;
                return old;
            }
            
            iterator& operator --()
            {
                ikey--;
                return *this;
            }
            
            iterator operator --(int)
            {
                iterator old=*this;
                ikey--;
                return old;
            }
            
            T0& operator *()
            {
                return xvec->data[*ikey];
            }
            
            T0& operator [](int i)
            {
                return xvec->data[*(ikey+i)];
            }
            
            bool operator!=(const iterator& rhs)
            {
                return (this->ikey!=rhs.ikey);
            }
            
            bool operator==(const iterator& rhs)
            {
                return (this->ikey==rhs.ikey);
            }
            void neq()
            {
                xvec->act.neq(*ikey);
            }
            void eq()
            {
                xvec->act.eq(*ikey);
            }
        };
        
        XVec(T0* _data
             ,int _size,C0& _act):
        data(_data),
        size(_size),
        act(_act)
        {
            key=NULL;
            if(!size) return;
            
            key=new int[size];
            for(int i=0;i<size;i++) key[i]=i;
            quicksort(key,key+size,
            [this](int* ikey,int* jkey){return (data[*ikey]<data[*jkey]);},
            [](int* ikey,int* jkey){std::swap(*ikey,*jkey);});
        }
        
        XVec(XVec&& other):
        act(other.act)
        {
            this->data=other.data;
            this->key=other.key;
            this->size=other.size;
            other.key=NULL;
        }
        ~XVec()
        {
            delete [] key;
        }
        
        XVec& operator=(XVec&& other)
        {
            this->act=other.act;
            this->data=other.data;
            this->key=other.key;
            this->size=other.size;
            other.key=NULL;
            return *this;
        }
        
        iterator begin()
        {
            return iterator(this,key);
        }
        
        iterator end()
        {
            return iterator(this,key+size);
        }
        
        iterator rbegin()
        {
            return iterator(this,key+size-1);
        }
        
        iterator rend()
        {
            return iterator(this,key-1);
        }
    };
}
/*------------------------------------------------------------------------------------------------------------------------------------
 
 ------------------------------------------------------------------------------------------------------------------------------------*/
namespace MAPP_NS
{
    namespace XMatrixVector
    {
        const float precision_f(sqrt(2.0*std::numeric_limits<float>::epsilon()));
        const double precision_d(sqrt(2.0*std::numeric_limits<double>::epsilon()));
        const long double precision_ld(sqrt(2.0*std::numeric_limits<long double>::epsilon()));
        
        template<class T>
        class precision
        {};
        
        template<>
        class precision<float>
        {public: static float val(){return precision_f;};};
        template<>
        class precision<double>
        {public: static float val(){return precision_d;};};
        template<>
        class precision<long double>
        {public: static long double val(){return precision_ld;};};
        
        
        template<const int dim>
        class _zero_
        {
        public:
            template<typename T>
            static inline void func(T* v)
            {
                *v=0;
                _zero_<dim-1>::func(v+1);
            }
        };
        template<>
        class _zero_<0>
        {
        public:
            template<typename T>
            static inline void func(T*)
            {
            }
        };
        
        template<const int i>
        class _rsq_
        {
        public:
            template<typename T>
            static inline T func(T* V0,T* V1)
            {
                return _rsq_<i-1>::func(V0+1,V1+1)+(*V0-*V1)*(*V0-*V1);
            }
        };
        
        template<>
        class _rsq_<1>
        {
        public:
            template<typename T>
            static inline T func(T* V0,T* V1)
            {
                return (*V0-*V1)*(*V0-*V1);
            }
        };
        
        template<const int i>
        class _V_V_
        {
        public:
            template<typename T>
            static inline T func(T* V0,T* V1)
            {
                return _V_V_<i-1>::func(V0+1,V1+1)+(*V0)*(*V1);
            }
        };
        
        template<>
        class _V_V_<1>
        {
        public:
            template<typename T>
            static inline T func(T* V0,T* V1)
            {
                return (*V0)*(*V1);
            }
        };
        template<>
        class _V_V_<0>
        {
        public:
            template<typename T>
            static inline T func(T*,T*)
            {
                return 0;
            }
        };
        
        template<const int stride,const int i>
        class _V_V_str_
        {
        public:
            template<typename T>
            static inline T func(T* V0,T* V1)
            {
                return _V_V_str_<stride,i-1>::func(V0+stride,V1+1)+(*V0)*(*V1);
            }
        };
        
        template<const int stride>
        class _V_V_str_<stride,1>
        {
        public:
            template<typename T>
            static inline T func(T* V0,T* V1)
            {
                return (*V0)*(*V1);
            }
        };
        
        template<const int stride,const int i>
        class _V_str_V_str_
        {
        public:
            template<typename T>
            static inline T func(T* V0,T* V1)
            {
                return _V_str_V_str_<stride,i-1>::func(V0+stride,V1+stride)+(*V0)*(*V1);
            }
        };
        
        template<const int stride>
        class _V_str_V_str_<stride,1>
        {
        public:
            template<typename T>
            static inline T func(T* V0,T* V1)
            {
                return (*V0)*(*V1);
            }
        };
        
        
        template<const int dim,const int i>
        class _V_Mlt_
        {
        public:
            template<typename T>
            static inline void func(T* Mlt,T* V)
            {
                *V=_V_V_str_<dim,i>::func(Mlt,V);
                _V_Mlt_<dim,i-1>::func(Mlt+dim+1,V+1);
            }
        };
        
        template<const int dim>
        class _V_Mlt_<dim,1>
        {
        public:
            template<typename T>
            static inline void func(T* Mlt,T* V)
            {
                *V*=*Mlt;
            }
        };
        
        template<const int dim,const int i>
        class _V_V_Mlt_
        {
        public:
            template<typename T>
            static inline void func(T* V1,T* Mlt,T* V0)
            {
                *V1=_V_V_str_<dim,i>::func(Mlt,V0);
                _V_V_Mlt_<dim,i-1>::func(V1+1,Mlt+dim+1,V0+1);
            }
        };
        
        template<const int dim>
        class _V_V_Mlt_<dim,1>
        {
        public:
            template<typename T>
            static inline void func(T* V1, T* Mlt,T* V0)
            {
                *V1*=*V0 * *Mlt;
            }
        };
        
        
        template<const int dim,const int i>
        class _x2s_
        {
        public:
            template<typename T>
            static inline void func(T* Mlt,T* V)
            {
                *V=_V_V_str_<dim,i>::func(Mlt,V);
                while(*V<0.0)
                    (*V)++;
                while(*V>=1.0)
                    (*V)--;
                _x2s_<dim,i-1>::func(Mlt+dim+1,V+1);
            }
        };
        
        template<const int dim>
        class _x2s_<dim,1>
        {
        public:
            template<typename T>
            static inline void func(T* Mlt,T* V)
            {
                *V*=*Mlt;
                while(*V<0.0)
                    (*V)++;
                while(*V>=1.0)
                    (*V)--;
            }
        };
        
        
        template<const int dim,const int i,const int j>
        class _Mlt_inv_
        {
        public:
            template<typename T>
            static inline void func(T* M,T* Minv)
            {
                *Minv=-_V_V_str_<dim,i-j>::func(Minv-(i-j)*dim,M)*(*(Minv+i-j));
                _Mlt_inv_<dim,i,j-1>::func(M-1,Minv-1);
            }
        };
        
        template<const int dim,const int i>
        class _Mlt_inv_<dim,i,i>
        {
        public:
            template<typename T>
            static inline void func(T* M,T* Minv)
            {
                *Minv=1.0/(*M);
                _Mlt_inv_<dim,i,i-1>::func(M-1,Minv-1);
            }
        };
        template<const int dim,const int i>
        class _Mlt_inv_<dim,i,0>
        {
        public:
            template<typename T>
            static inline void func(T* M,T* Minv)
            {
                *Minv=-_V_V_str_<dim,i>::func(Minv-i*dim,M)*(*(Minv+i));
                _Mlt_inv_<dim,i+1,i+1>::func(M+dim+i+1,Minv+dim+i+1);
            }
        };
        
        template<const int dim>
        class _Mlt_inv_<dim,dim,dim>
        {
        public:
            template<typename T>
            static inline void func(T*,T*)
            {}
        };

        template<const int dim,const int i,const int j>
        class _M_2_Mlt_
        {
        public:
            template<typename T>
            static inline void func(T* M_old,T* M_new)
            {
                
                M_new[i*dim+j]=(_V_V_<dim>::func(M_old+j*dim,M_old+i*dim)-
                _V_V_<j>::func(M_new+j*dim,M_new+i*dim))/M_new[j*(dim+1)];
                _M_2_Mlt_<dim,i,j+1>::func(M_old,M_new);
            }
        };
        
        template<const int dim,const int i>
        class _M_2_Mlt_<dim,i,i>
        {
        public:
            template<typename T>
            static inline void func(T* M_old,T* M_new)
            {
                M_new[i*dim+i]=sqrt(_V_V_<dim>::func(M_old+i*dim,M_old+i*dim)-
                    _V_V_<i>::func(M_new+i*dim,M_new+i*dim));
                
                _zero_<dim-i-1>::func(M_new+i*(dim+1)+1);
                _M_2_Mlt_<dim,i+1,0>::func(M_old,M_new);
            }
        };
        
        template<const int dim>
        class _M_2_Mlt_<dim,dim,0>
        {
        public:
            template<typename T>
            static inline void func(T*,T*)
            {}
        };
        
        template<const int dim,const int i,const int j>
        class _Mlt_Mlt_
        {
        public:
            template<typename T>
            static inline void func(T* M0,T* M1,T* B)
            {
                *B=_V_V_str_<dim,i-j+1>::func(M1-(i-j)*dim,M0);
                _Mlt_Mlt_<dim,i,j+1>::func(M0+1,M1+1,B+1);
            }
        };
        
        template<const int dim,const int i>
        class _Mlt_Mlt_<dim,i,i>
        {
        public:
            template<typename T>
            static inline void func(T* M0,T* M1,T* B)
            {
                *B=(*M0)*(*M1);
                _Mlt_Mlt_<dim,i+1,0>::func(M0+dim-i,M1+dim-i,B+dim-i);
            }
        };
        template<const int dim>
        class _Mlt_Mlt_<dim,dim,0>
        {
        public:
            template<typename T>
            static inline void func(T*,T*,T*)
            {}
        };
        
        
        template<int dim,int i>
        class __depth__
        {
        public:
            template<typename T>
            static inline void func(T* A,T* d)
            {
                *d=sqrt(_V_str_V_str_<dim,i>::func(A,A));
                __depth__<dim,i-1>::func(A+dim+1,d+1);
            }
        };
        
        template<int dim>
        class __depth__<dim,1>
        {
        public:
            template<typename T>
            static inline void func(T* A,T* d)
            {
                *d=abs(*A);
            }
        };
        
        template <int i,int dim>
        class __trace__
        {
        public:
            template<typename T>
            static inline T func(T* A)
            {
                return *A+__trace__<i-1,dim>::func(A+dim+1);
            }
        };
        template <int dim>
        class __trace__<1,dim>
        {
        public:
            template<typename T>
            static inline T func(T* A)
            {
                return *A;
            }
        };
        
        
    
        template<int i,int dim>
        class __dyadic
        {
        public:
            template<typename T>
            static inline void func(T& scl,T* x,T* v)
            {
                *v=*x*x[i]*scl;
                __dyadic<i+1,dim>::func(scl,x,v+1);
            }
        };
        
        template<int dim>
        class __dyadic<dim,dim>
        {
        public:
            template<typename T>
            static inline void func(T& scl,T* x,T* v)
            {
                __dyadic<0,dim-1>::func(scl,x+1,v);
            }
        };
        
        template<>
        class __dyadic<0,0>
        {
        public:
            template<typename T>
            static inline void func(T&,T*,T*)
            {
            }
        };
        
        template<int i,int j,int dim>
        class __SymmVec_2_Mlt
        {
        public:
            template<typename T>
            static inline void func(T* v,T* A)
            {
                *A=*v;
                __SymmVec_2_Mlt<i+1,j,dim>::func(v+1,A+dim+1);
            }
        };
        

        template<int j,int dim>
        class __SymmVec_2_Mlt<dim,j,dim>
        {
        public:
            template<typename T>
            static inline void func(T* v,T* A)
            {
                *A=*v;
                __SymmVec_2_Mlt<j+1,j+1,dim>::func(v+1,A+(j-dim+1)*(dim+1)+1);
            }
        };
        
        template<int dim>
        class __SymmVec_2_Mlt<dim,dim,dim>
        {
        public:
            template<typename T>
            static inline void func(T* v,T* A)
            {
                *A=*v;
            }
        };
        
        
        
        /*-------------------------------------------------------------------*/
        template <const int i>
        class UnrolledLoop
        {
        public:
            template < typename FuncType >
            static inline void Do(FuncType func)
            {
                UnrolledLoop <i-1>::Do( func);
                func(i-1);
            }
        };
        template <>
        class UnrolledLoop<1>
        {
        public:
            template < typename FuncType >
            static inline void Do(FuncType func){func(0);}
        };
        
        template<const int dim,typename T>
        inline T rsq(T* V0,T* V1)
        {
            return _rsq_<dim>::func(V0,V1);
        }
        
        template<const int dim,typename T>
        inline T V_V(T* V0,T* V1)
        {
            return _V_V_<dim>(V0,V1);
        }
        
        template<const int dim,typename T>
        inline void V_Mlt(T* V,T** Mlt)
        {
            _V_Mlt_<dim,dim>::func(*Mlt,V);
        }
        template<const int dim,typename T>
        static inline void V_Mlt(T* V,T (&Mlt)[dim][dim])
        {
            _V_Mlt_<dim,dim>::func((T*)Mlt,V);
        }
        
        template<const int dim,typename T>
        inline void V_Mlt(T* V,T** Mlt,T* ans)
        {
            _V_V_Mlt_<dim,dim>::func(ans,*Mlt,V);
        }
        template<const int dim,typename T>
        static inline void V_Mlt(T* V,T (&Mlt)[dim][dim],T* ans)
        {
            _V_V_Mlt_<dim,dim>::func(ans,(T*)Mlt,V);
        }
        
        template<const int dim,typename T>
        static inline T Tr(T(&A)[dim][dim])
        {
            return __trace__<dim,dim>::func((T*)A);
        }
        
        template<const int dim,typename T>
        inline void s2x(T* s,T (&H)[dim][dim])
        {
            _V_Mlt_<dim,dim>::func(&H[0][0],s);
        }
        
        template<const int dim,typename T>
        inline void s2x(T* s,T** H)
        {
            _V_Mlt_<dim,dim>::func(*H,s);
        }
        
        template<const int dim,typename T>
        inline void x2s(T* x,T (&B)[dim][dim])
        {
            _x2s_<dim,dim>::func(&B[0][0],x);
        }
        
        template<const int dim,typename T>
        inline void x2s(T* x,T** B)
        {
            _x2s_<dim,dim>::func(*B,x);
        }
        
        template<const int dim,typename T>
        inline void Mlt_inv(T** A,T** Ainv)
        {
            **Ainv=1.0/(**A);
            _Mlt_inv_<dim,1,1>::func(*A+dim+1,*Ainv+dim+1);
        }
        
        template<const int dim,typename T>
        inline void Mlt_inv(T (&A)[dim][dim],T (&Ainv)[dim][dim])
        {
            **Ainv=1.0/(**A);
            _Mlt_inv_<dim,1,1>::func((T*)A+dim+1,(T*)Ainv+dim+1);
        }
        
        template<const int dim,typename T>
        inline void Mlt_Mlt(T** M0,T** M1,T** B)
        {
            _Mlt_Mlt_<dim,0,0>::func(*M0,*M1,*B);
        }
        
        template<const int dim,typename T>
        inline void Mlt_Mlt(T (&M0)[dim][dim],T (&M1)[dim][dim],T (&B)[dim][dim])
        {
            _Mlt_Mlt_<dim,0,0>::func((T*)M0,(T*)M1,(T*)B);
        }
        
        template<const int dim,typename T>
        inline void M_2_Mlt(T** &M0,T** M1)
        {
            _M_2_Mlt_<dim,0,0>::func(*M0,*M1);
        }
        
        template<const int dim,typename T>
        inline void M_2_Mlt(T (&M0)[dim][dim],T (&M1)[dim][dim])
        {
            _M_2_Mlt_<dim,0,0>::func((T*)M0,(T*)M1);
        }
        
        
        template<const int dim,typename T>
        inline void depth(T (&A)[dim][dim],T (&d)[dim])
        {
            __depth__<dim,dim>::func((T*)A,d);
        }

        template<const int dim,typename T>
        inline void dyadic(T& scl,T* x,T (&v)[dim*(dim+1)/2])
        {
            __dyadic<0,dim>::func(scl,x,v);
        }
        
        template<const int dim,typename T>
        inline void SymmVec_2_Mlt(T (&v)[dim*(dim+1)/2],T(&A)[dim][dim])
        {
            __SymmVec_2_Mlt<0,0,dim-1>::func(v,&A[0][0]);
        }
        
        template<const int dim,typename T>
        inline void zero(T(&v)[dim])
        {
            _zero_<dim>::func(v);
        }

        template<const int dim,const int idim>
        class _opt_comm_grid_
        {
        public:
            static void func(const type0 (&h)[dim],int no,int(&curr_grid)[dim],int(&opt_grid)[dim],type0& opt_size)
            {
                for(int i=1;i<=no;i++)
                    if(no%i==0)
                    {
                        curr_grid[idim-1]=i;
                        _opt_comm_grid_<dim,idim-1>::func(h,no/i,curr_grid,opt_grid,opt_size);
                    }
            }
        };
        
        template<const int dim>
        class _opt_comm_grid_<dim,1>
        {
        public:
            static void func(const type0 (&h)[dim],int no,int(&curr_grid)[dim],int(&opt_grid)[dim],type0& opt_size)
            {
                curr_grid[0]=no;
                type0 curr_size=0.0;
                for(int i=0;i<dim;i++) curr_size+=static_cast<type0>(curr_grid[i])/h[i];
                if(opt_size==0.0 || curr_size<opt_size)
                {
                    opt_size=curr_size;
                    memcpy(opt_grid,curr_grid,dim*sizeof(int));
                }
            }
        };
        
        
        class __opt_comm_grid__
        {
        public:
            static void func(const int idim,const int dim,const type0* h,int no,int* curr_grid,int* opt_grid,type0& opt_size)
            {
                if(idim==1)
                {
                    curr_grid[0]=no;
                    type0 curr_size=0.0;
                    for(int i=0;i<dim;i++) curr_size+=static_cast<type0>(curr_grid[i])/h[i];
                    if(opt_size==0.0 || curr_size<opt_size)
                    {
                        opt_size=curr_size;
                        for(int i=0;i<dim;i++)
                            opt_grid[i]=curr_grid[i];
                    }
                    return;
                }
                
                for(int i=1;i<=no;i++)
                    if(no%i==0)
                    {
                        curr_grid[idim-1]=i;
                        __opt_comm_grid__::func(idim-1,dim,h,no/i,curr_grid,opt_grid,opt_size);
                    }
            }
        };
        
        template<const int dim>
        inline void opt_comm_grid(const type0 (&h)[dim],int no,int (&opt_grid)[dim])
        {
            if(no==1)
            {
                for(int i=0;i<dim;i++)
                    opt_grid[i]=1;
                return;
            }
            
            int curr_grid[dim];
            type0 opt_size=0.0;
            _opt_comm_grid_<dim,dim>::func(h,no,curr_grid,opt_grid,opt_size);
            
        }
        
        inline void opt_comm_grid(const type0 (&)[1],int no,int (&opt_grid)[1])
        {
            opt_grid[0]=no;
        }
        

        inline void opt_comm_grid(const int dim,const type0* h,int no,int* opt_grid)
        {
            if(dim==1)
            {
                opt_grid[0]=no;
                return;
            }
            
            if(no==1)
            {
                for(int i=0;i<dim;i++)
                    opt_grid[i]=1;
                return;
            }
            
            int* curr_grid =new int[dim];
            type0 opt_size=0.0;
            __opt_comm_grid__::func(dim,dim,h,no,curr_grid,opt_grid,opt_size);
            delete [] curr_grid;
        }
    }
}
/*--------------------------------------------
 quick sort algorithm
 
 example to test:
 
 int list_sz=1000;
 int* list=new int[list_sz];
 
 Random* rand=new Random(this,58245647);
 for(int i=0;i<list_sz;++i)
     list[i]=static_cast<int>(rand->uniform()*10000.0);
 
 XMath* xmath=new XMath();
 
 auto comp=
 [] (int* i,int* j)->bool
 {
     return (*i<*j);
 };
 auto swap=
 [] (int* i,int* j)->void
 {
     if(i==j)
         return;
     int k=*i;
     *i=*j;
     *j=k;
 };
 
 xmath->quicksort(list,list+list_sz,comp,swap);
 
 for(int i=1;i<list_sz;i++)
     if(list[i]<list[i-1])
         cout<<"error"<<endl;
 
 delete xmath;
 delete rand;
 delete [] list;
 --------------------------------------------*/
template<typename T0,class COMP,class SWAP>
void XMath::quicksort(T0 start,T0 end,COMP comp,SWAP swap)
{
    if(start+1>=end)
        return;
    T0 pindex=start;
    for(T0 i=start;i!=end-1;i++)
        if(comp(i,end-1))
        {
            swap(i,pindex);
            pindex++;
        }
    swap(end-1,pindex);
    quicksort(start,pindex,comp,swap);
    quicksort(pindex+1,end,comp,swap);
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<typename T0,class C0>
void XMath::srch_lst_lst(T0* ilst,int isize,C0* iact
,T0* jlst,int jsize,C0* jact)
{
    XVec<T0,C0> _ilst(ilst,isize,*iact);
    XVec<T0,C0> _jlst(jlst,jsize,*jact);
    
    if(isize==0 || jsize==0)
    {
        if(isize)
        {
            auto ipos=_ilst.begin();
            auto iend=_ilst.end();
            while(ipos!=iend)
            {
                ipos.neq();
                ++ipos;
            }
        }
        
        if(jsize)
        {
            auto jpos=_jlst.begin();
            auto jend=_jlst.end();
            while(jpos!=jend)
            {
                jpos.neq();
                ++jpos;
            }
        }
        
        return;
    }
    
    auto ipos=_ilst.begin();
    auto jpos=_jlst.begin();
    auto iend=_ilst.end();
    auto jend=_jlst.end();
    
    auto exit=
    [&]()->void
    {
        while(ipos!=iend)
        {
            ipos.neq();
            ++ipos;
        }
        
        while(jpos!=jend)
        {
            jpos.neq();
            ++jpos;
        }
        return;
    };
    
    while(1)
    {
        while(*ipos!=*jpos)
        {
            while(*ipos<*jpos)
            {
                ipos.neq();
                ++ipos;
                if(ipos==iend)
                    return exit();
            }
            std::swap(ipos,jpos);
            std::swap(iend,jend);
        }
        while(*ipos==*jpos)
        {
            ipos.eq();
            jpos.eq();
            ++ipos;
            ++jpos;
            if(ipos==iend || jpos==jend)
                return exit();
        }
    }    
}
/*--------------------------------------------
 
 --------------------------------------------*/
/*
template<const int dim>
void XMath::invert_lower_triangle(type0 (&A)[dim][dim],type0 (&Ainv)[dim][dim])
{
    type0 ATA[dim][dim];
    type0 c[dim];
    type0 x[dim];
    type0 g[dim];
    type0 g0[dim];
    type0 h[dim];
    type0 a0,a1,alpha;
    type0 g0g0,gg,gg0,ratio;
    
    for(int i=0;i<dim;i++)
        for(int j=0;j<dim;j++)
            Ainv[i][j]=ATA[i][j]=0.0;
    
    for(int i=0;i<dim;i++)
        for(int j=0;j<dim;j++)
            for(int k=0;k<dim;k++)
                ATA[i][j]+=A[k][i]*A[k][j];
    
    for(int itry=0;itry<dim;itry++)
    {
        for(int i=0;i<dim;i++)
        {
            c[i]=A[itry][i];
            x[i]=c[i];
        }
        
        
        g0g0=0.0;
        for(int i=0;i<dim;i++)
        {
            h[i]=2.0*c[i];
            for(int j=0;j<dim;j++)
                h[i]-=2.0*ATA[i][j]*x[j];
            g[i]=h[i];
            g0g0+=h[i]*h[i];
        }
        
        int jtry=0;
        double error=1.0;
        while(jtry<dim+1 && error!=0.0)
        {
            
            if(g0g0==0.0)
            {
                error=0.0;
                continue;
            }
            
            
            a0=0.0;
            a1=0.0;
            for(int i=0;i<dim;i++)
            {
                a0+=h[i]*g[i];
                for(int j=0;j<dim;j++)
                    a1+=h[i]*ATA[i][j]*h[j];
            }
            if(a1==0.0)
            {
                error=0.0;
                continue;
            }
            alpha=0.5*a0/a1;
            
            for(int i=0;i<dim;i++)
                x[i]+=alpha*h[i];
            
            //cout << "chk 3" << endl;
            
            gg=0.0;
            gg0=0.0;
            for(int i=0;i<dim;i++)
            {
                g[i]=2.0*c[i];
                for(int j=0;j<dim;j++)
                    g[i]-=2.0*ATA[i][j]*x[j];
                gg+=g[i]*g[i];
                gg0+=g0[i]*g[i];
            }
            
            //cout << "chk 4" << endl;
            ratio=(gg-gg0)/g0g0;
            g0g0=gg;
            
            
            for(int i=0;i<dim;i++)
            {
                h[i]=ratio*h[i]+g[i];
                g0[i]=g[i];
            }
            
            
            error=0.0;
            for(int i=0;i<dim;i++)
            {
                for(int j=0;j<dim;j++)
                    error+=x[i]*ATA[i][j]*x[i];
                error-=2*c[i]*x[i];
            }
            error++;
            
            jtry++;
        }
        
        for(int i=0;i<dim;i++)
            Ainv[i][itry]=x[i];
    }
}
*/
/*------------------------------------------------------------------------------------------------------------------------------------
 
 ------------------------------------------------------------------------------------------------------------------------------------*/
namespace MAPP_NS
{
    namespace Algebra
    {
        
        template<const int dim>
        class __V_eq
        {
        public:
            template<typename T>
            static inline void func(T const * src,T* dst)
            {
                *dst=*src;
                __V_eq<dim-1>::func(src+1,dst+1);
            }
        };
        
        template<>
        class __V_eq<1>
        {
        public:
            template<typename T>
            static inline void func(T const * src,T* dst)
            {
                *dst=*src;
            }
        };
        
        template<const int dim>
        class __V_add
        {
        public:
            template<typename T>
            static inline void func(T const * src,T* dst)
            {
                *dst+=*src;
                __V_add<dim-1>::func(src+1,dst+1);
            }
        };
        
        template<>
        class __V_add<1>
        {
        public:
            template<typename T>
            static inline void func(T const * src,T* dst)
            {
                *dst+=*src;
            }
        };
        
        template<const int dim>
        class __V_add_x_mul_V
        {
        public:
            template<typename T>
            static inline void func(const T& x,T const * src,T* dst)
            {
                *dst+=x**src;
                __V_add_x_mul_V<dim-1>::func(x,src+1,dst+1);
            }
        };
        
        template<>
        class __V_add_x_mul_V<1>
        {
        public:
            template<typename T>
            static inline void func(const T& x,T const * src,T* dst)
            {
                *dst+=x**src;
            }
        };
        
        template<const int dim>
        class __V_eq_x_mul_V
        {
        public:
            template<typename T>
            static inline void func(const T& x,T const * src,T* dst)
            {
                *dst=x**src;
                __V_eq_x_mul_V<dim-1>::func(x,src+1,dst+1);
            }
        };
        
        template<>
        class __V_eq_x_mul_V<1>
        {
        public:
            template<typename T>
            static inline void func(const T& x,T const * src,T* dst)
            {
                *dst=x**src;
            }
        };
        
        
        template<const int dim>
        class __V_zero
        {
        public:
            template<typename T>
            static inline void func(T* v)
            {
                *v=0;
                __V_zero<dim-1>::func(v+1);
            }
        };
        template<>
        class __V_zero<1>
        {
        public:
            template<typename T>
            static inline void func(T* v)
            {
                *v=0;
            }
        };
        
        template<>
        class __V_zero<0>
        {
        public:
            template<typename T>
            static inline void func(T*)
            {
            }
        };
        
        
        
        template<const int dim,const int strd>
        class __V_strd_zero
        {
        public:
            template<typename T>
            static inline void func(T* v)
            {
                *v=0;
                __V_zero<dim-1>::func(v+strd);
            }
        };
        template<const int strd>
        class __V_strd_zero<1,strd>
        {
        public:
            template<typename T>
            static inline void func(T* v)
            {
                *v=0;
            }
        };
        template<const int strd>
        class __V_strd_zero<0,strd>
        {
        public:
            template<typename T>
            static inline void func(T*)
            {
            }
        };
        
        
        
        /* dot product for aligned vectors */
        template<const int i>
        class __V_mul_SCL
        {
        public:
            template<typename T>
            static inline void func(T* vec,T scl)
            {
                *vec*=scl;
                return __V_mul_SCL::func(vec+1,scl);
            }
        };
        
        template<>
        class __V_mul_SCL<1>
        {
        public:
            template<typename T>
            static inline void func(T* vec,T scl)
            {
                *vec*=scl;
            }
        };
        
        
        
        
        /* dot product for aligned vectors */
        template<const int i>
        class __V_mul_V
        {
        public:
            template<typename T>
            static inline T func(T* vec0,T* vec1)
            {
                return *vec0**vec1+__V_mul_V<i-1>::func(vec0+1,vec1+1);
            }
        };
        
        template<>
        class __V_mul_V<1>
        {
        public:
            template<typename T>
            static inline T func(T* vec0,T* vec1)
            {
                return *vec0**vec1;
            }
        };
        
        template<>
        class __V_mul_V<0>
        {
        public:
            template<typename T>
            static inline T func(T*,T*)
            {
                return 0;
            }
        };
        
        
        
        /* dot product for vectors with one of them with stride */
        template<const int i,const int strd>
        class __V_strd_mul_V
        {
        public:
            template<typename T>
            static inline T func(T* vec0,T* vec1)
            {
                return *vec0**vec1+__V_strd_mul_V<i-1,strd>::func(vec0+strd,vec1+1);
            }
        };
        
        template<const int strd>
        class __V_strd_mul_V<1,strd>
        {
        public:
            template<typename T>
            static inline T func(T* vec0,T* vec1)
            {
                return *vec0**vec1;
            }
        };
        
        
        template<const int strd>
        class __V_strd_mul_V<0,strd>
        {
        public:
            template<typename T>
            static inline T func(T*,T*)
            {
                return 0;
            }
        };
        
        
        /* dot product for vectors with both of them with strides */
        template<const int i,const int strd0,const int strd1>
        class __V_strd_mul_V_strd
        {
        public:
            template<typename T>
            static inline T func(T* vec0,T* vec1)
            {
                return *vec0**vec1+__V_strd_mul_V_strd<i-1,strd0,strd1>::func(vec0+strd0,vec1+strd1);
            }
        };
        
        template<const int strd0,const int strd1>
        class __V_strd_mul_V_strd<1,strd0,strd1>
        {
        public:
            template<typename T>
            static inline T func(T* vec0,T* vec1)
            {
                return *vec0**vec1;
            }
        };
        
        
        template<const int strd0,const int strd1>
        class __V_strd_mul_V_strd<0,strd0,strd1>
        {
        public:
            template<typename T>
            static inline T func(T*,T*)
            {
                return 0;
            }
        };
        
        
        
        
        
        
        
        
        
        
        
        /* dot product for vector by lower triangular matrix VMLT=V*MLT */
        template<const int i,const int dim>
        class __V_mul_MLT
        {
        public:
            template<typename T>
            static inline void func(T* V,T* MLT,T* VMLT)
            {
                *VMLT=__V_strd_mul_V<i,dim>::func(MLT,V);
                __V_mul_MLT<i-1,dim>::func(V+1,MLT+dim+1,VMLT+1);
            }
            template<typename T>
            static inline void add_in(T* V,T* MLT,T* VMLT)
            {
                *VMLT+=__V_strd_mul_V<i,dim>::func(MLT,V);
                __V_mul_MLT<i-1,dim>::add_in(V+1,MLT+dim+1,VMLT+1);
            }
        };
        
        template<const int dim>
        class __V_mul_MLT<1,dim>
        {
        public:
            template<typename T>
            static inline void func(T* V,T* MLT,T* VMLT)
            {
                *VMLT=*V**MLT;
            }
            template<typename T>
            static inline void add_in(T* V,T* MLT,T* VMLT)
            {
                *VMLT+=*V**MLT;
            }
        };
        
        /* dot product for vector by lower triangular matrix VMLT=MLT*V */
        /* WARNING: MLT should start at the last row MLT+(dim-1)*dim
                and VMLT should start at the last component VMLT+ dim-1
         */
        template<const int i,const int dim>
        class __MLT_mul_V
        {
        public:
            template<typename T>
            static inline void func(T* MLT,T* V,T* VMLT)
            {
                *VMLT=__V_mul_V<i>::func(MLT,V);
                __MLT_mul_V<i-1,dim>::func(MLT-dim,V,VMLT-1);
            }
        };
        
        template<const int dim>
        class __MLT_mul_V<1,dim>
        {
        public:
            template<typename T>
            static inline void func(T* MLT,T* V,T* VMLT)
            {
                *VMLT=*V**MLT;
            }
        };
        
        
        
        
        /* dot product for vector by upper triangular matrix VMLT=V*MLT */
        /* WARNING: MUT should start at the last column MUT+dim-1
                and VMUT should start at the last component VMUT+ dim-1
         */
        template<const int i,const int dim>
        class __V_mul_MUT
        {
        public:
            template<typename T>
            static inline void func(T* V,T* MUT,T* VMUT)
            {
                *VMUT=__V_strd_mul_V<i,dim>::func(MUT,V);
                __V_mul_MUT<i-1,dim>::func(V,MUT-1,VMUT-1);
            }
        };
        
        template<const int dim>
        class __V_mul_MUT<1,dim>
        {
        public:
            template<typename T>
            static inline void func(T* V,T* MUT,T* VMUT)
            {
                *VMUT=*V**MUT;
            }
        };
        
        
        
        
        
        /* dot product for vector by lower triangular matrix VMUT=MUT*V */
        template<const int i,const int dim>
        class __MUT_mul_V
        {
        public:
            template<typename T>
            static inline void func(T* MUT,T* V,T* VMUT)
            {
                *VMUT=__V_mul_V<i>::func(MUT,V);
                __MUT_mul_V<i-1,dim>::func(MUT+dim+1,V+1,VMUT+1);
            }
        };
        
        template<const int dim>
        class __MUT_mul_V<1,dim>
        {
        public:
            template<typename T>
            static inline void func(T* MUT,T* V,T* VMUT)
            {
                *VMUT=*V**MUT;
            }
        };
        
        
        
        
        
        
        template<const int row,const int diff,const int dim>
        class __MLT_inv
        {
        public:
            template<typename T>
            static inline void func(T* MLT,T* MLT_inv)
            {
                *MLT_inv=-MLT_inv[diff]*__V_strd_mul_V<diff,dim>::func(MLT_inv-diff*dim,MLT);
                __MLT_inv<row,diff+1,dim>::func(MLT-1,MLT_inv-1);
            }
        };
        
        
        template<const int diff,const int dim>
        class __MLT_inv<diff,diff,dim>
        {
        public:
            template<typename T>
            static inline void func(T* MLT,T* MLT_inv)
            {
                *MLT_inv=-MLT_inv[diff]*__V_strd_mul_V<diff,dim>::func(MLT_inv-diff*dim,MLT);
                __MLT_inv<diff+1,0,dim>::func(MLT+dim+diff+1,MLT_inv+dim+diff+1);
            }
        };
        
        template<const int row,const int dim>
        class __MLT_inv<row,0,dim>
        {
        public:
            template<typename T>
            static inline void func(T* MLT,T* MLT_inv)
            {
                *MLT_inv=1.0/(*MLT);
                __MLT_inv<row,1,dim>::func(MLT-1,MLT_inv-1);
            }
        };
        
        
        template<const int dim>
        class __MLT_inv<0,0,dim>
        {
        public:
            template<typename T>
            static inline void func(T* MLT,T* MLT_inv)
            {
                *MLT_inv=1.0/(*MLT);
                __MLT_inv<1,0,dim>::func(MLT+dim+1,MLT_inv+dim+1);
            }
        };
        
        
        template<const int dim>
        class __MLT_inv<dim,0,dim>
        {
        public:
            template<typename T>
            static inline void func(T*,T*)
            {}
        };
        
        
        template<const int i,const int dim>
        class __MLT_depth
        {
        public:
            template<typename T>
            static inline void func(T const * MLT,T* depth)
            {
                *depth=sqrt(__V_strd_mul_V_strd<i,dim,dim>::func(MLT,MLT));
                __MLT_depth<i-1,dim>::func(MLT+dim+1,depth+1);
            }
        };
        
        template<const int dim>
        class __MLT_depth<1,dim>
        {
        public:
            template<typename T>
            static inline void func(T const * MLT,T* depth)
            {
                *depth=fabs(*MLT);
            }
        };
        
        
        
        
        
        template<const int i,const int j,const int dim>
        class __MLT_mul_MSQ
        {
        public:
            template<typename T>
            static inline void func(T* MLT,T* MSQ,T* MLT_MSQ)
            {
                *MLT_MSQ=__V_strd_mul_V<i,dim>::func(MSQ,MLT);
                __MLT_mul_MSQ<i,j-1,dim>::func(MLT,MSQ-1,MLT_MSQ-1);
            }
        };
        
        template<const int i,const int dim>
        class __MLT_mul_MSQ<i,1,dim>
        {
        public:
            template<typename T>
            static inline void func(T* MLT,T* MSQ,T* MLT_MSQ)
            {
                *MLT_MSQ=__V_strd_mul_V<i,dim>::func(MSQ,MLT);
                __MLT_mul_MSQ<i-1,dim,dim>::func(MLT-dim,MSQ+dim-1,MLT_MSQ-1);
            }
        };
        
        template<const int dim>
        class __MLT_mul_MSQ<1,1,dim>
        {
        public:
            template<typename T>
            static inline void func(T* MLT,T* MSQ,T* MLT_MSQ)
            {
                *MLT_MSQ=*MSQ**MLT;
            }
        };
        
        
        
        template<const int i,const int j,const int dim>
        class __MSQ_mul_MLT
        {
        public:
            template<typename T>
            static inline void func(T* MSQ,T* MLT,T* MSQ_MLT)
            {
                *MSQ_MLT=__V_strd_mul_V<j,dim>::func(MLT,MSQ);
                __MSQ_mul_MLT<i,j-1,dim>::func(MSQ+1,MLT+dim+1,MSQ_MLT+1);
            }
        };
        
        template<const int i,const int dim>
        class __MSQ_mul_MLT<i,1,dim>
        {
        public:
            template<typename T>
            static inline void func(T* MSQ,T* MLT,T* MSQ_MLT)
            {
                *MSQ_MLT=*MSQ**MLT;
                __MSQ_mul_MLT<i-1,dim,dim>::func(MSQ+1,MLT-dim*dim+1,MSQ_MLT+1);
            }
        };
        
        template<const int dim>
        class __MSQ_mul_MLT<1,1,dim>
        {
        public:
            template<typename T>
            static inline void func(T* MSQ,T* MLT,T* MSQ_MLT)
            {
                *MSQ_MLT=*MSQ**MLT;
            }
        };
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        template<const int i,const int j,const int dim>
        class ___MSQ_mul_MLT
        {
        public:
            template<typename T>
            static inline void func(T* MSQ,T* MLT,T* MSQ_MLT)
            {
                *MSQ_MLT=__V_strd_mul_V<j,dim>::func(MLT,MSQ);
                ___MSQ_mul_MLT<i,j-1,dim>::func(MSQ+1,MLT+dim+1,MSQ_MLT+1);
            }
        };
        
        template<const int i,const int dim>
        class ___MSQ_mul_MLT<i,i,dim>
        {
        public:
            template<typename T>
            static inline void func(T* MSQ,T* MLT,T* MSQ_MLT)
            {
                *MSQ_MLT=__V_strd_mul_V<i,dim>::func(MLT,MSQ);
                ___MSQ_mul_MLT<i-1,dim,dim>::func(MSQ+i,MLT-(dim-i)*(dim+1),MSQ_MLT+i);
            }
        };
        
        
        template<const int dim>
        class ___MSQ_mul_MLT<1,1,dim>
        {
        public:
            template<typename T>
            static inline void func(T* MSQ,T* MLT,T* MSQ_MLT)
            {
                *MSQ_MLT=*MSQ**MLT;
            }
        };
        
        
        
        
        
        
        
        
        
        
        
        
        
        template<const int i,const int j,const int dim>
        class ___MUT_mul_MSQ_mul_MLT__
        {
        public:
            template<typename T>
            static inline void func(T* MLT,T* MSQ,T* ANS)
            {
                *ANS=__V_strd_mul_V_strd<i,dim,dim>::func(MLT,MSQ);
                ___MUT_mul_MSQ_mul_MLT__<i-1,j,dim>::func(MLT+dim+1,MSQ+dim,ANS+dim);
            }
        };
        
        template<const int j,const int dim>
        class ___MUT_mul_MSQ_mul_MLT__<1,j,dim>
        {
        public:
            template<typename T>
            static inline void func(T* MLT,T* MSQ,T* ANS)
            {
                *ANS=*MLT**MSQ;
                ___MUT_mul_MSQ_mul_MLT__<j-1,j-1,dim>::func(MLT-(j-2)*(dim+1),MSQ-(j-2)*dim+1,ANS-(j-2)*dim+1);
            }
        };
        
        
        template<const int dim>
        class ___MUT_mul_MSQ_mul_MLT__<1,1,dim>
        {
        public:
            template<typename T>
            static inline void func(T* MLT,T* MSQ,T* ANS)
            {
                *ANS=*MLT**MSQ;
            }
        };
        
        
        template<const int i,const int j,const int dim >
        class __MLT_T
        {
        public:
            template<typename T>
            static inline void func(T* RESTRICT MLT,T* RESTRICT MUT)
            {
                *MUT=*MLT;
                __MLT_T<i-1,j,dim>::func(MLT+dim,MUT+1);
            }
        };
        
        
        template<const int j,const int dim >
        class __MLT_T<1,j,dim>
        {
            public:
            template<typename T>
            static inline void func(T* RESTRICT MLT,T* RESTRICT MUT)
            {
                *MUT=*MLT;
                __MLT_T<j-1,j-1,dim>::func(MLT+(2-j)*dim+1,MUT-j+dim+2);
            }
        };
        
        template<const int dim >
        class __MLT_T<1,1,dim>
        {
            public:
            template<typename T>
            static inline void func(T* RESTRICT MLT,T* RESTRICT MUT)
            {
                *MUT=*MLT;
            }
        };
        
        
         template<const int i,const int j,const int dim >
        class __MLT_2_V
        {
        public:
            template<typename T>
            static inline void func(T* RESTRICT MLT,T* RESTRICT mlt)
            {
                *mlt=*MLT;
                __MLT_2_V<i-1,j,dim>::func(MLT+dim,mlt+1);
            }
        };
        
        
        template<const int j,const int dim >
        class __MLT_2_V<1,j,dim>
        {
            public:
            template<typename T>
            static inline void func(T* RESTRICT MLT,T* RESTRICT mlt)
            {
                *mlt=*MLT;
                __MLT_2_V<j-1,j-1,dim>::func(MLT+(2-j)*dim+1,mlt+1);
            }
        };
        
        template<const int dim >
        class __MLT_2_V<1,1,dim>
        {
            public:
            template<typename T>
            static inline void func(T* RESTRICT MLT,T* RESTRICT mlt)
            {
                *mlt=*MLT;
            }
        };

        
        
        
        template<const int column,const int diff,const int dim>
        class __MUT_inv
        {
        public:
            template<typename T>
            static inline void func(T* MUT,T* MUT_inv)
            {
                *MUT_inv=-MUT_inv[diff*dim]*__V_strd_mul_V<diff,dim>::func(MUT,MUT_inv-diff);
                __MUT_inv<column,diff+1,dim>::func(MUT-dim,MUT_inv-dim);
            }
        };
        
        
        template<const int diff,const int dim>
        class __MUT_inv<diff,diff,dim>
        {
        public:
            template<typename T>
            static inline void func(T* MUT,T* MUT_inv)
            {
                *MUT_inv=-MUT_inv[diff*dim]*__V_strd_mul_V<diff,dim>::func(MUT,MUT_inv-diff);
                __MUT_inv<diff+1,0,dim>::func(MUT+(diff+1)*dim+1,MUT_inv+(diff+1)*dim+1);
            }
        };
        
        
        template<const int column,const int dim>
        class __MUT_inv<column,0,dim>
        {
        public:
            template<typename T>
            static inline void func(T* MUT,T* MUT_inv)
            {
                *MUT_inv=1.0/(*MUT);
                __MUT_inv<column,1,dim>::func(MUT-dim,MUT_inv-dim);
            }
        };
        
        
        template<const int dim>
        class __MUT_inv<0,0,dim>
        {
        public:
            template<typename T>
            static inline void func(T* MUT,T* MUT_inv)
            {
                *MUT_inv=1.0/(*MUT);
                __MUT_inv<1,0,dim>::func(MUT+dim+1,MUT_inv+dim+1);
            }
        };
        
        
        
        template<const int dim>
        class __MUT_inv<dim,0,dim>
        {
        public:
            template<typename T>
            static inline void func(T*,T*)
            {}
        };
        
        
        
    
        
        
        
        
        template<const int i,const int j,const int dim>
        class __MSQ_2_MLT
        {
        public:
            template<typename T>
            static inline void func(T* MSQ,T* MLT)
            {
                
                MLT[i*dim+j]=(__V_mul_V<dim>::func(MSQ+j*dim,MSQ+i*dim)-
                __V_mul_V<j>::func(MLT+j*dim,MLT+i*dim))/MLT[j*(dim+1)];
                __MSQ_2_MLT<i,j+1,dim>::func(MSQ,MLT);
            }
        };
        
        template<const int i,const int dim>
        class __MSQ_2_MLT<i,i,dim>
        {
        public:
            template<typename T>
            static inline void func(T* MSQ,T* MLT)
            {
                MLT[i*dim+i]=sqrt(__V_mul_V<dim>::func(MSQ+i*dim,MSQ+i*dim)-
                __V_mul_V<i>::func(MLT+i*dim,MLT+i*dim));
                
                __V_zero<dim-i-1>::func(MLT+i*(dim+1)+1);
                __MSQ_2_MLT<i+1,0,dim>::func(MSQ,MLT);
            }
        };
        
        template<const int dim>
        class __MSQ_2_MLT<dim,0,dim>
        {
        public:
            template<typename T>
            static inline void func(T*,T*)
            {}
        };
        
        
        
        
        
        
        
        
        
        
        template<const int i,const int j,const int dim>
        class __MSQ_2_MUT
        {
        public:
            template<typename T>
            static inline void func(T* M_old,T* M_new)
            {
                
                M_new[i*dim+j]=(__V_strd_mul_V_strd<dim,dim,dim>::func(M_old+j,M_old+i)-
                __V_strd_mul_V_strd<i,dim,dim>::func(M_new+j,M_new+i))/M_new[i*(dim+1)];
                __MSQ_2_MUT<i+1,j,dim>::func(M_old,M_new);
            }
        };
        
        template<const int i,const int dim>
        class __MSQ_2_MUT<i,i,dim>
        {
        public:
            template<typename T>
            static inline void func(T* M_old,T* M_new)
            {
                
                
                M_new[i*dim+i]=sqrt(__V_strd_mul_V_strd<dim,dim,dim>::func(M_old+i,M_old+i)-
                __V_strd_mul_V_strd<i,dim,dim>::func(M_new+i,M_new+i));
                
                __V_strd_zero<dim-i-1,dim>::func(M_new+i*(dim+1)+dim);
                __MSQ_2_MUT<0,i+1,dim>::func(M_old,M_new);
            }
        };
        
        template<const int dim>
        class __MSQ_2_MUT<0,dim,dim>
        {
        public:
            template<typename T>
            static inline void func(T*,T*)
            {}
        };
        
        
        
        
        
        template<const int row,const int diff,const int dim>
        class __MLT_mul_MLT
        {
        public:
            template<typename T>
            static inline void func(T* MLTL,T* MLTR,T* MLT)
            {
                *MLT=__V_strd_mul_V<diff+1,dim>::func(MLTR,MLTL);
                __MLT_mul_MLT<row,diff-1,dim>::func(MLTL+1,MLTR+dim+1,MLT+1);
            }
        };
        
        template<const int row,const int dim>
        class __MLT_mul_MLT<row,0,dim>
        {
        public:
            template<typename T>
            static inline void func(T* MLTL,T* MLTR,T* MLT)
            {
                *MLT=*MLTL**MLTR;
                __MLT_mul_MLT<row+1,row+1,dim>::func(MLTL-row+dim,MLTR-(dim+1)*row,MLT-row+dim);
            }
        };
        
        template<const int dim>
        class __MLT_mul_MLT<dim,dim,dim>
        {
        public:
            template<typename T>
            static inline void func(T*,T*,T*)
            {}
        };
        
        
        
        
        
        
        template<const int column,const int diff,const int dim>
        class __MUT_mul_MUT
        {
        public:
            template<typename T>
            static inline void func(T* MUTL,T* MUTR,T* MUT)
            {
                *MUT=__V_strd_mul_V<diff+1,dim>::func(MUTR,MUTL);
                __MUT_mul_MUT<column,diff-1,dim>::func(MUTL+1+dim,MUTR+dim,MUT+dim);
            }
        };
        template<const int column,const int dim>
        class __MUT_mul_MUT<column,0,dim>
        {
        public:
            template<typename T>
            static inline void func(T* MUTL,T* MUTR,T* MUT)
            {
                *MUT=*MUTL**MUTR;
                __MUT_mul_MUT<column+1,column+1,dim>::func(MUTL-column*(dim+1),MUTR-column*dim+1,MUT-column*dim+1);
            }
        };
        
        template<const int dim>
        class __MUT_mul_MUT<dim,dim,dim>
        {
        public:
            template<typename T>
            static inline void func(T*,T*,T*)
            {}
        };
        
        
        
        
        
        
        
        
        
        template<const int i,const int j,const int dim>
        class __DyadicV_mul_MLT
        {
        public:
            template<typename T>
            static inline void func(T* DyadicV,T* MLT, T* M)
            {
                *M=__V_strd_mul_V<j,dim>::func(MLT,DyadicV);
                __DyadicV_mul_MLT<i,j-1,dim>::func(DyadicV+1,MLT+dim+1,M+1);
                
            }
        };
        
        template<const int i,const int dim>
        class __DyadicV_mul_MLT<i,1,dim>
        {
        public:
            template<typename T>
            static inline void func(T* DyadicV,T* MLT, T* M)
            {
                *M=*MLT**DyadicV;
                __DyadicV_mul_MLT<i-1,i-1,dim>::func(DyadicV+1,MLT-dim*dim+1,M+1);
                
            }
        };
        template<const int dim>
        class __DyadicV_mul_MLT<1,1,dim>
        {
        public:
            template<typename T>
            static inline void func(T* DyadicV,T* MLT, T* M)
            {
                *M=*MLT**DyadicV;
            }
        };
        
        
        
        
        
        
        
        
        
        
        
        
        template<const int i,const int dim>
        class __DyadicV
        {
        public:
            template<typename T>
            static inline void func(T& scl,T* x,T* dyad)
            {
                *dyad+=*x*x[dim-i]*scl;
                __DyadicV<i-1,dim>::func(scl,x,dyad+1);
            }
            template<typename T>
            static inline void func(T* x,T* y,T* dyad)
            {
                *dyad+=*x**y;
                __DyadicV<i-1,dim>::func(x,y+1,dyad+1);
            }
            
            template<typename T>
            static inline void func(T scl,T* x,T* y,T* dyad)
            {
                *dyad+=scl*(*x*y[dim-i]+*y*x[dim-i]);
                __DyadicV<i-1,dim>::func(scl,x,y,dyad+1);
            }
        };
        
        template<const int dim>
        class __DyadicV<1,dim>
        {
        public:
            template<typename T>
            static inline void func(T& scl,T* x,T* dyad)
            {
                *dyad+=*x*x[dim-1]*scl;
                __DyadicV<dim-1,dim-1>::func(scl,x+1,dyad+1);
            }
            template<typename T>
            static inline void func(T* x,T* y,T* dyad)
            {
                *dyad+=*x**y;
                __DyadicV<dim-1,dim-1>::func(x+1,y-dim+2,dyad+1);
            }
            template<typename T>
            static inline void func(T scl,T* x,T* y,T* dyad)
            {
                *dyad+=scl*(*x*y[dim-1]+*y*x[dim-1]);
                __DyadicV<dim-1,dim-1>::func(scl,x+1,y+1,dyad+1);
            }
        };
        
        template<>
        class __DyadicV<1,1>
        {
        public:
            template<typename T>
            static inline void func(T& scl,T* x,T* dyad)
            {
                *dyad+=*x**x*scl;
            }
            template<typename T>
            static inline void func(T* x,T* y,T* dyad)
            {
                *dyad+=*x**y;
            }
            template<typename T>
            static inline void func(T scl,T* x,T* y,T* dyad)
            {
                *dyad+=2.0*scl**x**y;
                
            }
        };
        
        
        
        template<const int i,const int j,const int dim>
        class __DyadicV_2_MLT
        {
        public:
            template<typename T>
            static inline void func(T* dyad,T* MLT)
            {
                *MLT=*dyad;
                __DyadicV_2_MLT<i-1,j,dim>::func(dyad+1,MLT+dim);
            }
        };
        
        
        template<const int j,const int dim>
        class __DyadicV_2_MLT<1,j,dim>
        {
        public:
            template<typename T>
            static inline void func(T* dyad,T* MLT)
            {
                *MLT=*dyad;
                __DyadicV_2_MLT<j-1,j-1,dim>::func(dyad+1,MLT-dim*(j-2)+1);
            }
        };
        
        template<const int dim>
        class __DyadicV_2_MLT<1,1,dim>
        {
        public:
            template<typename T>
            static inline void func(T* dyad,T* MLT)
            {
                *MLT=*dyad;
            }
        };
        
        
        
        template<const int i,const int j,const int dim>
        class __NONAME_DyadicV_mul_MLT
        {
        public:
            template<typename T>
            static inline void func(T* dyad,T* MLT,T* ANS)
            {
                *ANS=__V_strd_mul_V<i,dim>::func(MLT,dyad);
                __NONAME_DyadicV_mul_MLT<i-1,j,dim>::func(dyad+1,MLT+dim+1,ANS+dim);
            }
        };
        
        template<const int j,const int dim>
        class __NONAME_DyadicV_mul_MLT<1,j,dim>
        {
        public:
            template<typename T>
            static inline void func(T* dyad,T* MLT,T* ANS)
            {
                *ANS=*MLT**dyad;
                __NONAME_DyadicV_mul_MLT<j-1,j-1,dim>::func(dyad+1,MLT-(j-2)*(dim+1),ANS-(j-2)*dim+1);
            }
        };
        
        template<const int dim>
        class __NONAME_DyadicV_mul_MLT<1,1,dim>
        {
        public:
            template<typename T>
            static inline void func(T* dyad,T* MLT,T* ANS)
            {
                *ANS=*MLT**dyad;
            }
        };
        
        
        template<int i,int j,int dim>
        class __NONAME_DyadicV_2_MLT
        {
        public:
            template<class T>
            static void func(T* x)
            {
                x[(i-j)*(dim-1)+j*(j-1)/2]=*x;
                __NONAME_DyadicV_2_MLT<i-1,j,dim>::func(x-1);
            };
        };
        
        template<int i,int dim>
        class __NONAME_DyadicV_2_MLT<i,i,dim>
        {
        public:
            template<class T>
            static void func(T* x)
            {
                x[i*(i-1)/2]=*x;
                __NONAME_DyadicV_2_MLT<dim,i-1,dim>::func(x-1);
            };
        };
        
        template<int dim>
        class __NONAME_DyadicV_2_MLT<1,1,dim>
        {
        public:
            template<class T>
            static void func(T*)
            {
            };
        };
        
        template<const int i,const int j,const int dim>
        class __DyadicV_2_MUT
        {
        public:
            template<typename T>
            static inline void func(T* dyad,T* MUT)
            {
                *MUT=*dyad;
                __DyadicV_2_MUT<i-1,j,dim>::func(dyad+1,MUT+1);
            }
        };
        
        template<const int j,const int dim>
        class __DyadicV_2_MUT<1,j,dim>
        {
        public:
            template<typename T>
            static inline void func(T* dyad,T* MUT)
            {
                *MUT=*dyad;
                __DyadicV_2_MUT<j-1,j-1,dim>::func(dyad+1,MUT+dim+2-j);
            }
        };
        
        template<const int dim>
        class __DyadicV_2_MUT<1,1,dim>
        {
        public:
            template<typename T>
            static inline void func(T* dyad,T* MUT)
            {
                *MUT=*dyad;
            }
        };
        
        

        template<const int i,const int j,const int dim>
        class __DyadicV_2_MSY
        {
        public:
            template<typename T>
            static inline void func(T* dyad,T* MSY)
            {
                *MSY=MSY[(j-i)*(dim-1)]=*dyad;
                __DyadicV_2_MSY<i-1,j,dim>::func(dyad+1,MSY+1);
            }
        };
        
        
        template<const int i,const int dim>
        class __DyadicV_2_MSY<i,i,dim>
        {
        public:
            template<typename T>
            static inline void func(T* dyad,T* MSY)
            {
                *MSY=*dyad;
                __DyadicV_2_MSY<i-1,i,dim>::func(dyad+1,MSY+1);
            }
        };
        
        template<const int j,const int dim>
        class __DyadicV_2_MSY<1,j,dim>
        {
        public:
            template<typename T>
            static inline void func(T* dyad,T* MSY)
            {
                *MSY=MSY[(j-1)*(dim-1)]=*dyad;
                __DyadicV_2_MSY<j-1,j-1,dim>::func(dyad+1,MSY+dim+2-j);
            }
        };
        
        template<const int dim>
        class __DyadicV_2_MSY<1,1,dim>
        {
        public:
            template<typename T>
            static inline void func(T* dyad,T* MSY)
            {
                *MSY=*dyad;
            }
        };
        
        
        template <const int i,const int dim>
        class __MLT_det
        {
        public:
            template<typename T>
            static inline T func(T* M)
            {
                return *M*__MLT_det<i-1,dim>::func(M+dim+1);
            }
        };
        template <const int dim>
        class __MLT_det<1,dim>
        {
        public:
            template<typename T>
            static inline T func(T* M)
            {
                return *M;
            }
        };
        
        
        template <const int i,const int dim>
        class __Tr
        {
        public:
            template<typename T>
            static inline T func(T* M)
            {
                return *M+__Tr<i-1,dim>::func(M+dim+1);
            }
        };
        template <const int dim>
        class __Tr<1,dim>
        {
        public:
            template<typename T>
            static inline T func(T* M)
            {
                return *M;
            }
        };
        
        template <const int i,const int dim>
        class __Tr_DyadicV
        {
        public:
            template<typename T>
            static inline T func(T* M)
            {
                return *M+__Tr_DyadicV<i-1,dim>::func(M+i);
            }
        };
        template <const int dim>
        class __Tr_DyadicV<1,dim>
        {
        public:
            template<typename T>
            static inline T func(T* M)
            {
                return *M;
            }
        };
        
        
        template<const int i>
        class __DX_RSQ
        {
        public:
            template<class T>
            static inline T func(T const* xi,T const* xj ,T* dxij)
            {
                *dxij=*xi-*xj;
                return *dxij**dxij+__DX_RSQ<i-1>::func(xi+1,xj+1,dxij+1);
            }
            
            template<class T>
            static inline T func(T const * xi,T const * xj)
            {
                return (*xi-*xj)*(*xi-*xj)+__DX_RSQ<i-1>::func(xi+1,xj+1);
            }
            
            template<class T>
            static inline void __func(T const* xi,T const* xj ,T* dxij)
            {
                *dxij=*xi-*xj;
                __DX_RSQ<i-1>::__func(xi+1,xj+1,dxij+1);
            }
        };
        
        template<>
        class __DX_RSQ<1>
        {
        public:
            template<class T>
            static inline T func(T const* xi,T const* xj ,T* dxij)
            {
                *dxij=*xi-*xj;
                return *dxij**dxij;
            }
            
            template<class T>
            static inline T func(T const * xi,T const * xj)
            {
                return (*xi-*xj)*(*xi-*xj);
            }
            
            template<class T>
            static inline void __func(T const* xi,T const* xj ,T* dxij)
            {
                *dxij=*xi-*xj;
            }
        };
        
        
        
        template<const int dim,const int i>
        class __DX_HAT_R
        {
        public:
            template<class T>
            static inline T func(T const* xi,T const* xj ,T* dxij)
            {
                *dxij=*xi-*xj;
                return *dxij**dxij+__DX_HAT_R<dim,i-1>::func(xi+1,xj+1,dxij+1);
            }
        };
        
        template<const int dim>
        class __DX_HAT_R<dim,1>
        {
        public:
            template<class T>
            static inline T func(T const* xi,T const* xj ,T* dxij)
            {
                *dxij=*xi-*xj;
                return *dxij**dxij;
            }
        };
        
        template<const int dim>
        class __DX_HAT_R<dim,dim>
        {
        public:
            template<class T>
            static inline T func(T const* xi,T const* xj ,T* dxij)
            {
                *dxij=*xi-*xj;
                T r=sqrt(*dxij**dxij+__DX_HAT_R<dim,dim-1>::func(xi+1,xj+1,dxij+1));
                __V_mul_SCL<dim>::func(dxij,1.0/r);
                return r;
            }
        };
        
        
        
        template<const int N,const int M>
        class __pow
        {
        public:
            static constexpr int value=N*__pow<N,M-1>::value;
        };
        
        template<const int N>
        class __pow<N,1>
        {
        public:
            static constexpr int value=N;
        };
        
        
        
        
        template <const int i>
        class Do
        {
        public:
            template<typename FuncType>
            static inline void func(FuncType f)
            {
                Do<i-1>::func(f);
                f(i-1);
            }
        };
        template <>
        class Do<1>
        {
        public:
            template<typename FuncType>
            static inline void func(FuncType f){f(0);}
        };
        
        
        template <const int i,const int j=i>
        class DoLT
        {
        public:
            template<typename FuncType>
            static inline void func(FuncType f)
            {
                DoLT<i,j-1>::func(f);
                f(i-1,j-1);
            }
        };
        
        template <const int i>
        class DoLT<i,1>
        {
        public:
            template<typename FuncType>
            static inline void func(FuncType f)
            {
                DoLT<i-1,i-1>::func(f);
                f(i-1,0);
            }
        };
        
        template <>
        class DoLT<1,1>
        {
        public:
            template<typename FuncType>
            static inline void func(FuncType f)
            {f(0,0);}
        };
        
        

        
        
        template<const int dim>
        class __S2X
        {
            public:
            template<typename T>
            static inline void func(T* RESTRICT h,T* RESTRICT s)
            {
                *s=__V_mul_V<dim>::func(h,s);
                
                __S2X<dim-1>::func(h+dim,s+1);
            }
        };
        
        template<>
        class __S2X<1>
        {
            public:
            template<typename T>
            static inline void func(T* RESTRICT h,T* RESTRICT s)
            {
                *s*=*h;
            }
        };
        
        template<const int dim>
        class __X2S
        {
            public:
            template<typename T>
            static inline void func(T* RESTRICT b,T* RESTRICT x)
            {
                *x=__V_mul_V<dim>::func(b,x);
                while(*x<0.0) ++*x;
                while(*x>=1.0) --*x;
                __X2S<dim-1>::func(b+dim,x+1);
            }
        };
        
        template<>
        class __X2S<1>
        {
            public:
            template<typename T>
            static inline void func(T* RESTRICT b,T* RESTRICT x)
            {
                *x*=*b;
                while(*x<0.0) ++*x;
                while(*x>=1.0) --*x;
            }
        };
        
        
        
        class __opt_comm_grid
        {
            public:
            static void func(const int idim,const int dim,const type0* h,int no,int* curr_grid,int* opt_grid,type0& opt_size)
            {
                if(idim==1)
                {
                    curr_grid[0]=no;
                    type0 curr_size=0.0;
                    for(int i=0;i<dim;i++) curr_size+=static_cast<type0>(curr_grid[i])/h[i];
                    if(opt_size==0.0 || curr_size<opt_size)
                    {
                        opt_size=curr_size;
                        for(int i=0;i<dim;i++)
                        opt_grid[i]=curr_grid[i];
                    }
                    return;
                }
                
                for(int i=1;i<=no;i++)
                if(no%i==0)
                {
                    curr_grid[idim-1]=i;
                    __opt_comm_grid::func(idim-1,dim,h,no/i,curr_grid,opt_grid,opt_size);
                }
            }
        };
        
        

        
        
        /*==========================================================================*/
        template<const int dim,typename T>
        void zero(T* V)
        {
            __V_zero<dim>::func(V);
        }
        /*==========================================================================*/
        template<const int dim,typename T>
        void MSQ_2_MLT(T(&MSQ)[dim][dim],T(&MLT)[dim][dim])
        {
            __MSQ_2_MLT<0,0,dim>::func(&MSQ[0][0],&MLT[0][0]);
        }
        /*==========================================================================*/
        template<const int dim,typename T>
        void MLT_inv(T(&MLT)[dim][dim],T(&MLT_inv)[dim][dim])
        {
            __MLT_inv<0,0,dim>::func(&MLT[0][0],&MLT_inv[0][0]);
        }
        /*==========================================================================*/
        template<const int dim,typename T>
        void MUT_inv(T(&MUT)[dim][dim],T(&MUT_inv)[dim][dim])
        {
            __MUT_inv<0,0,dim>::func(&MUT[0][0],&MUT_inv[0][0]);
        }
        /*==========================================================================*/
        template<const int dim,typename T>
        T V_mul_V(const T* V0,const T* V1)
        {
            return __V_mul_V<dim>::func(V0,V1);
        }
        /*==========================================================================*/
        template<const int dim,typename T>
        void MLT_mul_V(T(&MLT)[dim][dim],T*&V,T*&MLTV)
        {
            __MLT_mul_V<dim,dim> ::func(&MLT[dim-1][0],V,&MLTV[dim-1]);
        }
        template<const int dim,typename T>
        void MLT_mul_V(T(&MLT)[dim][dim],T*&V,T(&MLTV)[dim])
        {
            __MLT_mul_V<dim,dim> ::func(&MLT[dim-1][0],V,&MLTV[dim-1]);
        }
        /*==========================================================================*/
        template<typename T,const int dim>
        void V_mul_MLT(T*& V,T(&MLT)[dim][dim],T*& VMLT)
        {
            __V_mul_MLT<dim,dim> ::func(V,&MLT[0][0],VMLT);
        }
        template<typename T,const int dim>
        void V_mul_MLT(T*& V,T(&MLT)[dim][dim],T(&VMLT)[dim])
        {
            __V_mul_MLT<dim,dim> ::func(V,&MLT[0][0],VMLT);
        }
        template<typename T,const int dim>
        void V_mul_MLT(T (&V)[dim],T(&MLT)[dim][dim],T(&VMLT)[dim])
        {
            __V_mul_MLT<dim,dim> ::func(V,&MLT[0][0],VMLT);
        }

        /*==========================================================================*/
        template<typename T,const int dim>
        void V_mul_MLT_add_in(T*& V,T(&MLT)[dim][dim],T*& VMLT)
        {
            __V_mul_MLT<dim,dim>::add_in(V,&MLT[0][0],VMLT);
        }
        template<typename T,const int dim>
        void V_mul_MLT_add_in(T*& V,T(&MLT)[dim][dim],T(&VMLT)[dim])
        {
            __V_mul_MLT<dim,dim>::add_in(V,&MLT[0][0],VMLT);
        }
        /*==========================================================================*/
        template<typename T,const int dim>
        void MLT_depth(T(&MLT)[dim][dim],T(&depth)[dim])
        {
            __MLT_depth<dim,dim>::func(&(MLT[0][0]),depth);
        }
        /*==========================================================================*/
        template<const int dim,typename T>
        void MUT_mul_V(T(&MUT)[dim][dim],T*&V,T*&MUTV)
        {
            __MUT_mul_V<dim,dim> ::func(&MUT[0][0],V,MUTV);
        }
        template<const int dim,typename T>
        void MUT_mul_V(T(&MUT)[dim][dim],T*&V,T(&MUTV)[dim])
        {
            __MUT_mul_V<dim,dim> ::func(&MUT[0][0],V,MUTV);
        }
        /*==========================================================================*/
        template<const int dim,typename T>
        void V_mul_MUT(T*&V,T(&MUT)[dim][dim],T*&VMUT)
        {
            __V_mul_MUT<dim,dim> ::func(V,&MUT[0][dim-1],VMUT[dim-1]);
        }
        template<const int dim,typename T>
        void V_mul_MUT(T*&V,T(&MUT)[dim][dim],T(&VMUT)[dim])
        {
            __V_mul_MUT<dim,dim> ::func(V,&MUT[0][dim-1],VMUT[dim-1]);
        }
        /*==========================================================================*/
        template<const int dim,typename T>
        void MLT_mul_MLT(T(&MLTL)[dim][dim],T(&MLTR)[dim][dim],T(&MLT)[dim][dim])
        {
            __MLT_mul_MLT<0,0,dim>::func(&MLTL[0][0],&MLTR[0][0],&MLT[0][0]);
        }
        template<const int dim,typename T>
        void MLT_mul_MLT(T(&MLTL)[dim][dim],T(*&MLTR)[dim],T(&MLT)[dim][dim])
        {
            __MLT_mul_MLT<0,0,dim>::func(&MLTL[0][0],&MLTR[0][0],&MLT[0][0]);
        }
        template<const int dim,typename T>
        void MLT_mul_MLT(T(*&MLTL)[dim],T(&MLTR)[dim][dim],T(&MLT)[dim][dim])
        {
            __MLT_mul_MLT<0,0,dim>::func(&MLTL[0][0],&MLTR[0][0],&MLT[0][0]);
        }
        /*==========================================================================*/
        template<const int dim,typename T>
        void MUT_mul_MUT(T(&MUTL)[dim][dim],T(&MUTR)[dim][dim],T(&MUT)[dim][dim])
        {
            __MUT_mul_MUT<0,0,dim>::func(&MUTL[0][0],&MUTR[0][0],&MUT[0][0]);
        }
        /*==========================================================================*/
        template<const int dim,typename T>
        void MSQ_mul_MLT(T(&MSQ)[dim][dim],T(&MLT)[dim][dim],T(&MSQ_MLT)[dim][dim])
        {
            __MSQ_mul_MLT<dim,dim,dim>::func(&MSQ[0][0],&MLT[0][0],&MSQ_MLT[0][0]);
        }
        template<const int dim,typename T>
        void MSQ_mul_MLT_LT(T(&MSQ)[dim][dim],T(&MLT)[dim][dim],T(&MSQ_MLT)[dim][dim])
        {
            ___MSQ_mul_MLT<dim,dim,dim>::func(&MSQ[0][0],&MLT[0][0],&MSQ_MLT[0][0]);
        }
        
        template<const int dim,typename T>
        void MUT_mul_MSQ_mul_MLT(T(&MSQ)[dim][dim],T(&MLT)[dim][dim],T(&MSQ_MLT)[dim][dim])
        {
            
            ___MUT_mul_MSQ_mul_MLT__<dim,dim,dim>::func(&MSQ[0][0],&MLT[0][0],&MSQ_MLT[0][0]);
        }
        template<const int dim,typename T>
        void MUT_mul_MSY_mul_MLT(T(&MSY)[dim][dim],T(&MLT)[dim][dim],T(&__MLT)[dim][dim],T(&ANS)[dim][dim])
        {
            ___MSQ_mul_MLT<dim,dim,dim>::func(&MSY[0][0],&MLT[0][0],&__MLT[0][0]);
            ___MUT_mul_MSQ_mul_MLT__<dim,dim,dim>::func(&MLT[0][0],&__MLT[0][0],&ANS[0][0]);
        }
        /*==========================================================================*/
        template<const int dim,typename T>
        void MLT_mul_MSQ(T(&MLT)[dim][dim],T(&MSQ)[dim][dim],T(&MLT_MSQ)[dim][dim])
        {
            __MLT_mul_MSQ<dim,dim,dim>::func(&MLT[dim-1][0],&MSQ[0][dim-1],&MLT_MSQ[dim-1][dim-1]);
        }
        /*==========================================================================*/
        template<const int dim,typename T>
        void MLT_T(T(&MLT)[dim][dim],T(&MUT)[dim][dim])
        {
            __MLT_T<dim,dim,dim>::func(&MLT[0][0],&MUT[0][0]);
        }
        /*==========================================================================*/
        template<const int dim,typename T>
        void MLT_2_V(T(&MLT)[dim][dim],T* mlt)
        {
            __MLT_2_V<dim,dim,dim>::func(&MLT[0][0],mlt);
        }
        /*==========================================================================*/
        template<const int dim,typename T>
        void DyadicV(T scl,T* x,T* dyad)
        {
            __DyadicV<dim,dim>::func(scl,x,dyad);
        }
        /*
        template<const int dim,typename T>
        void DyadicV(T& scl,T* x,T (&dyad)[dim*(dim+1)/2])
        {
            __DyadicV<dim,dim>::func(scl,x,dyad);
        }
        template<const int dim,typename T>
        void DyadicV(T scl,T (&x)[dim],T* dyad)
        {
            __DyadicV<dim,dim>::func(scl,x,dyad);
        }
        template<const int dim,typename T>
        void DyadicV(T& scl,T(&x)[dim],T (&dyad)[dim*(dim+1)/2])
        {
            __DyadicV<dim,dim>::func(scl,x,dyad);
        }*/
        template<const int dim,typename T>
        void DyadicV(T*&x,T*& y,T* dyad)
        {
            __DyadicV<dim,dim>::func(x,y,dyad);
        }
        /*==========================================================================*/
        template<const int dim,typename T>
        void DyadicV(T scl,T* x,T* y,T* dyad)
        {
            __DyadicV<dim,dim>::func(scl,x,y,dyad);
        }
        /*==========================================================================*/
        /*
        template<const int dim,typename T>
        void DyadicV_2_MLT(T* dyad,T (&MLT)[dim][dim])
        {
            __DyadicV_2_MLT<dim,dim,dim>::func(dyad,&MLT[0][0]);
        }*/
        template<const int dim,typename T>
        void DyadicV_2_MLT(T* dyad,T (*MLT)[dim])
        {
            __DyadicV_2_MLT<dim,dim,dim>::func(dyad,MLT[0]);
        }
        template<const int dim,typename T>
        void DyadicV_2_MSY(T* dyad,T (*MLT)[dim])
        {
            __DyadicV_2_MSY<dim,dim,dim>::func(dyad,MLT[0]);
        }
        /*==========================================================================*/
        template<const int dim,typename T>
        void DyadicV_2_MUT(T(&dyad)[dim*(dim+1)/2],T (&MUT)[dim][dim])
        {
            __DyadicV_2_MUT<dim,dim,dim>::func(dyad,&MUT[0][0]);
        }
        /*==========================================================================*/
        template<const int dim,typename T>
        T Tr(T(&M)[dim][dim])
        {
            return __Tr<dim,dim>::func(&M[0][0]);
        }
        template<const int dim,typename T>
        T Tr_DyadicV(T* DyadicV)
        {
            return __Tr_DyadicV<dim,dim>::func(DyadicV);
        }
        /*==========================================================================*/
        template<const int dim,typename T>
        T MLT_det(T(&M)[dim][dim])
        {
            return __MLT_det<dim,dim>::func(&M[0][0]);
        }
        /*==========================================================================*/
        template<const int dim,typename T>
        T MLT_max(T(&MLT)[dim][dim])
        {
            T max=MLT[0][0];
            Algebra::DoLT<dim>::func([&max,&MLT](int i,int j){max=max<MLT[i][j]?MLT[i][j]:max;});
            return max;
        }
        /*==========================================================================*/
        template<const int dim,typename T>
        T MLT_min(T(&MLT)[dim][dim])
        {
            T min=MLT[0][0];
            Algebra::DoLT<dim>::func([&min,&MLT](int i,int j){min=min>MLT[i][j]?MLT[i][j]:min;});
            return min;
        }
        /*==========================================================================*/
        template<const int dim,typename T>
        void NONAME_DyadicV_mul_MLT(T* DyadicV,T(&MLT)[dim][dim],T(&ANS)[dim][dim])
        {
            __NONAME_DyadicV_mul_MLT<dim,dim,dim>::func(DyadicV,&MLT[0][0],&ANS[0][0]);
        }
        template<const int dim,typename T>
        void NONAME_DyadicV_mul_MLT(T* DyadicV,T(&MLT)[dim][dim],T(*&ANS)[dim])
        {
            __NONAME_DyadicV_mul_MLT<dim,dim,dim>::func(DyadicV,&MLT[0][0],&ANS[0][0]);
        }
        /*==========================================================================*/
        template<int dim,class T>
        void NONAME_DyadicV_2_MLT(T(&MLT)[dim][dim])
        {
            __NONAME_DyadicV_2_MLT<dim,dim,dim>::func(&MLT[0][0]+dim*(dim+1)/2-1);
        }
        /*==========================================================================*/
        template<const int dim,typename T>
        T DX_RSQ(T const * xi,T const * xj,T (&dxij)[dim])
        {
            return __DX_RSQ<dim>::func(xi,xj,dxij);
        }
        /*==========================================================================*/
        template<const int dim,typename T>
        T DX_HAT_R(T const * xi,T const * xj,T (&dxij)[dim])
        {
            return __DX_HAT_R<dim,dim>::func(xi,xj,dxij);
        }
        /*==========================================================================*/
        /*
        template<const int dim,typename T>
        T RSQ(T const *& xi,T const *& xj)
        {
            return __DX_RSQ<dim>::func(xi,xj);
        }*/
        template<const int dim,typename T>
        T RSQ(T const * xi,T const * xj)
        {
            return __DX_RSQ<dim>::func(xi,xj);
        }
        /*==========================================================================*/
        template<const int dim,typename T>
        void DX(T const * xi,T const * xj,T* dxij)
        {
            __DX_RSQ<dim>::__func(xi,xj,dxij);
        }
        /*==========================================================================*/
        template<const int dim,typename T>
        void V_eq(T const * src,T* dst)
        {
            __V_eq<dim>::func(src,dst);
        }
        /*==========================================================================*/
        template<const int dim,typename T>
        void V_add(T const * src,T* dst)
        {
            __V_add<dim>::func(src,dst);
        }
        /*==========================================================================*/
        template<const int dim,typename T>
        void V_add_x_mul_V(const T x,T const * src,T* dst)
        {
            __V_add_x_mul_V<__dim__>::func(x,src,dst);
        }
        /*==========================================================================*/
        template<const int dim,typename T>
        void V_eq_x_mul_V(const T x,T const * src,T* dst)
        {
            __V_eq_x_mul_V<__dim__>::func(x,src,dst);
        }
        /*==========================================================================*/
        template<const int N,const int M>
        constexpr int pow()
        {
            return __pow<N,M>::value;
        }
        /*==========================================================================*/
        template<const int dim,typename T>
        void X2S(T* RESTRICT b,int N,T* RESTRICT x)
        {
            for(int i=0;i<N;i++)
            {
                __X2S<dim>::func(b,x);
                x+=dim;
            }
        }
        /*==========================================================================*/
        template<const int dim,typename T>
        void S2X(T* RESTRICT h,int N,T* RESTRICT s)
        {
            for(int i=0;i<N;i++)
            {
                __S2X<dim>::func(h,s);
                s+=dim;
            }
        }
        
        /*==========================================================================*/
        template<const int dim,typename T>
        void X2S(T* RESTRICT b,T* RESTRICT x)
        {
            __X2S<dim>::func(b,x);
        }
        /*==========================================================================*/
        template<const int dim,typename T>
        void S2X(T* RESTRICT h,T* RESTRICT s)
        {
            __S2X<dim>::func(h,s);
        }
        /*==========================================================================*/
        inline void opt_comm_grid(const int dim,const type0* h,int no,int* opt_grid)
        {
            if(dim==1)
            {
                opt_grid[0]=no;
                return;
            }
            
            if(no==1)
            {
                for(int i=0;i<dim;i++)
                    opt_grid[i]=1;
                return;
            }
            
            int* curr_grid =new int[dim];
            type0 opt_size=0.0;
            __opt_comm_grid::func(dim,dim,h,no,curr_grid,opt_grid,opt_size);
            delete [] curr_grid;
        }
        
        
        
        
    }
}
#endif
