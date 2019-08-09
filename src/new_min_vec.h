#ifndef __MAPP__new_min_vec__
#define __MAPP__new_min_vec__
#include "atoms.h"
#include "xmath.h"
#include "memory.h"
/*--------------------------------------------
 
 --------------------------------------------*/
namespace MAPP_NS
{

    template<class E>
    class VTExpr
    {
    private:
    protected:
    public:
        template<const int I>
        void reset_ptr(){static_cast<E&>(*this).template reset_ptr<I>();}
        type0 operator[](const int& i) const {return static_cast<const E&>(*this).operator[](i);}
        template<const int I,const int J>
        type0 get() const {return static_cast<const E&>(*this).template get<I,J>();}
    };
    
    template<class E>
    class VTMulScl : public VTExpr<VTMulScl<E>>
    {
    private:
    protected:
    public:
        E&  vt;
        type0& scl;
        VTMulScl(type0& __scl,E& __vt):vt(__vt),scl(__scl){}
        template<const int I>
        void reset_ptr() {vt.template reset_ptr<I>();}
        type0 operator[](const int& i) const {return vt.operator[](i)*scl;}
        template<const int I,const int J>
        type0 get() const {return vt.template get<I,J>()*scl;}
    };
    
    template<class El,class Er>
    class VTSum : public VTExpr<VTSum<El,Er>>
    {
    private:
    protected:
    public:
        El& vt_l;
        Er& vt_r;
        VTSum(El& __vt_l,Er& __vt_r):vt_l(__vt_l),vt_r(__vt_r){}
        template<const int I>
        void reset_ptr() {vt_l.template reset_ptr<I>();vt_r.template reset_ptr<I>();}
        type0 operator[](const int& i) const {return vt_l.operator[](i)+vt_r.operator[](i);}
        template<const int I,const int J>
        type0 get() const {return vt_l.template get<I,J>()+vt_r.template get<I,J>();}
    };
    
    template<class El,class Er>
    class VTSub : public VTExpr<VTSub<El,Er>>
    {
    private:
    protected:
    public:
        El& vt_l;
        Er& vt_r;
        VTSub(El& __vt_l,Er& __vt_r):vt_l(__vt_l),vt_r(__vt_r){}
        template<const int I>
        void reset_ptr() {vt_l.template reset_ptr<I>();vt_r.template reset_ptr<I>();}
        type0 operator[](const int& i) const {return vt_l.operator[](i)-vt_r.operator[](i);}
        template<const int I,const int J>
        type0 get() const {return vt_l.template get<I,J>()-vt_r.template get<I,J>();}
    };
    
    template<bool B,bool... Bs>
    class VTN
    {
    public:
        static constexpr int N=1+VTN<Bs...>::N;
        
        template<class ...Vs>
        void assign_vecs(Atoms* __atoms,Vec<type0>** __vs_arr,bool* __vs_arr_alloc,int* __vs_dim,Vec<type0>* __v,Vs... __vs)
        {
            *__vs_arr=__v;
            *__vs_arr_alloc=false;
            *__vs_dim=__v->dim;
            VTN<Bs...>::assign_vecs(__atoms,__vs_arr+1,__vs_arr_alloc+1,__vs_dim+1,__vs...);
        }
        template<class ...Vs>
        void assign_vecs(Atoms* __atoms,Vec<type0>** __vs_arr,bool* __vs_arr_alloc,int* __vs_dim,int __dim,Vs... __vs)
        {
            *__vs_arr=new Vec<type0>(__atoms,__dim);
            *__vs_arr_alloc=true;
            *__vs_dim=__dim;
            VTN<Bs...>::assign_vecs(__atoms,__vs_arr+1,__vs_arr_alloc+1,__vs_dim+1,__vs...);
        }
    };
    
    template<bool... Bs>
    class VTN<false,Bs...>
    {
    public:
        static constexpr int N=VTN<Bs...>::N;
        template<class ...Vs>
        void assign_vecs(Atoms* __atoms,Vec<type0>** __vs_arr,bool* __vs_arr_alloc,int* __vs_dim,Vec<type0>* __v,Vs... __vs)
        {
            VTN<Bs...>::assign_vecs(__atoms,__vs_arr,__vs_arr_alloc,__vs_dim,__vs...);
        }
        template<class ...Vs>
        void assign_vecs(Atoms* __atoms,Vec<type0>** __vs_arr,bool* __vs_arr_alloc,int* __vs_dim,int __dim,Vs... __vs)
        {
            VTN<Bs...>::assign_vecs(__atoms,__vs_arr,__vs_arr_alloc,__vs_dim,__vs...);
        }
    };
    
    template<>
    class VTN<true>
    {
    public:
        static constexpr int N=1;
        void assign_vecs(Atoms* __atoms,Vec<type0>** __vs_arr,bool* __vs_arr_alloc,int* __vs_dim,Vec<type0>* __v)
        {
            *__vs_arr=__v;
            *__vs_arr_alloc=false;
            *__vs_dim=__v->dim;
        }
        
        void assign_vecs(Atoms* __atoms,Vec<type0>** __vs_arr,bool* __vs_arr_alloc,int* __vs_dim,int __dim)
        {
            *__vs_arr=new Vec<type0>(__atoms,__dim);
            *__vs_arr_alloc=true;
            *__vs_dim=__dim;
        }
    };
    template<>
    class VTN<false>
    {
    public:
        static constexpr int N=0;
        void assign_vecs(Atoms*,Vec<type0>**,bool*,int* __vs_dim,Vec<type0>*)
        {
        }
        void assign_vecs(Atoms*,Vec<type0>**,bool*,int* __vs_dim,int)
        {
        }
    };
    
    template<const int I>
    class VOP
    {
    public:
        template<class E0,class E1,class F>
        static void func_expr(const int& natms_lcl,E0& vt_l,E1& vt_r,F&& f)
        {
            VOP<I-1>::func_expr(natms_lcl,vt_l,vt_r,f);
            vt_r.template reset_ptr<I-1>();
            const int __v_size=natms_lcl*vt_l.vecs_dims[I-1];
            type0* __v0=vt_l.vecs[I-1]->begin();
            for(int j=0;j<__v_size;j++) f(__v0[j],vt_r[j]);
        }
        
        template<class E,class F>
        static void func(const int& natms_lcl,E& vt_l,const E& vt_r,F&& f)
        {
            VOP<I-1>::func(natms_lcl,vt_l,vt_r,f);
            const int __v_size=natms_lcl*vt_l.vecs_dims[I-1];
            type0* __v0=vt_l.vecs[I-1]->begin();
            const type0* __v1=vt_r.vecs[I-1]->begin();
            for(int j=0;j<__v_size;j++) f(__v0[j],__v1[j]);
        }
        
        template<class E>
        static void eq(const int& natms_lcl,E& vt_l,const E& vt_r)
        {
            VOP<I-1>::eq(natms_lcl,vt_l,vt_r);
            memcpy(vt_l.vecs[I-1]->begin(),vt_r.vecs[I-1]->begin(),natms_lcl*vt_l.vecs_dims[I-1]*sizeof(type0));
        }
        
        template<class E,class F>
        static void func_scl(const int& natms_lcl,E& vt,const type0& scl,F&& f)
        {
            VOP<I-1>::func_scl(natms_lcl,vt,scl,f);
            const int __v_size=natms_lcl*vt.vecs_dims[I-1];
            type0* __v0=vt.vecs[I-1]->begin();
            for(int j=0;j<__v_size;j++) f(__v0[j],scl);
        }
    };
    
    template<>
    class VOP<1>
    {
    public:
        template<class E0,class E1,class F>
        static void func_expr(const int& natms_lcl,E0& vt_l,E1& vt_r,F&& f)
        {
            vt_r.template reset_ptr<0>();
            const int __v_size=natms_lcl*vt_l.vecs_dims[0];
            type0* __v0=vt_l.vecs[0]->begin();
            for(int j=0;j<__v_size;j++) f(__v0[j],vt_r[j]);
        }
        
        template<class E,class F>
        static void func(const int& natms_lcl,E& vt_l,const E& vt_r,F&& f)
        {
            const int __v_size=natms_lcl*vt_l.vecs_dims[0];
            type0* __v0=vt_l.vecs[0]->begin();
            const type0* __v1=vt_r.vecs[0]->begin();
            for(int j=0;j<__v_size;j++) f(__v0[j],__v1[j]);
        }
        
        template<class E>
        static void eq(const int& natms_lcl,E& vt_l,const E& vt_r)
        {
            memcpy(vt_l.vecs[0]->begin(),vt_r.vecs[0]->begin(),natms_lcl*vt_l.vecs_dims[0]*sizeof(type0));
        }
        
        template<class E,class F>
        static void func_scl(const int& natms_lcl,E& vt,const type0& scl,F&& f)
        {
            const int __v_size=natms_lcl*vt.vecs_dims[0];
            type0* __v0=vt.vecs[0]->begin();
            for(int j=0;j<__v_size;j++) f(__v0[j],scl);
        }
    };
    
    
    template<const int I,const int J=I>
    class TOP
    {
    public:
        template<class E0,class E1,class F>
        static void func_expr(E0& vt_l,const E1& vt_r,F&& f)
        {
            TOP<I,J-1>::func_expr(vt_l,vt_r,f);
            f(vt_l.A[I-1][J-1],vt_r.template get<I-1,J-1>());
        }
        
        template<class E,class F>
        static void func(E& vt_l,const E& vt_r,F&& f)
        {
            TOP<I,J-1>::func(vt_l,vt_r,f);
            f(vt_l.A[I-1][J-1],vt_r.A[I-1][J-1]);
        }
        
        template<class E>
        static void eq(E& vt_l,const E& vt_r)
        {
            TOP<I,J-1>::eq(vt_l,vt_r);
            vt_l.A[I-1][J-1]=vt_r.A[I-1][J-1];
        }
        
        
        template<class E,class F>
        static void func_scl(E& vt,const type0& scl,F&& f)
        {
            TOP<I,J-1>::func_scl(vt,scl,f);
            f(vt.A[I-1][J-1],scl);
        }
    };
    
    template<const int I>
    class TOP<I,1>
    {
    public:
        template<class E0,class E1,class F>
        static void func_expr(E0& vt_l,const E1& vt_r,F&& f)
        {
            TOP<I-1,I-1>::func_expr(vt_l,vt_r,f);
            f(vt_l.A[I-1][0],vt_r.template get<I-1,0>());
        }
        
        template<class E,class F>
        static void func(E& vt_l,const E& vt_r,F&& f)
        {
            TOP<I-1,I-1>::func(vt_l,vt_r,f);
            f(vt_l.A[I-1][0],vt_r.A[I-1][0]);
        }
        
        template<class E>
        static void eq(E& vt_l,const E& vt_r)
        {
            TOP<I-1,I-1>::eq(vt_l,vt_r);
            vt_l.A[I-1][0]=vt_r.A[I-1][0];
        }
        
        template<class E,class F>
        static void func_scl(E& vt,const type0& scl,F&& f)
        {
            TOP<I-1,I-1>::func_scl(vt,scl,f);
            f(vt.A[I-1][0],scl);
        }
        
    };
    
    
    template<>
    class TOP<1,1>
    {
    public:
        template<class E0,class E1,class F>
        static void func_expr(E0& vt_l,const E1& vt_r,F&& f)
        {
            f(vt_l.A[0][0],vt_r.template get<0,0>());
        }
        
        template<class E,class F>
        static void func(E& vt_l,const E& vt_r,F&& f)
        {
            f(vt_l.A[0][0],vt_r.A[0][0]);
        }
        
        template<class E>
        static void eq(E& vt_l,const E& vt_r)
        {
            vt_l.A[0][0]=vt_r.A[0][0];
        }
        
        template<class E,class F>
        static void func_scl(E& vt,const type0& scl,F&& f)
        {
            f(vt.A[0][0],scl);
        }
    };
    
    template<bool>
    class TOPHelper
    {
    public:
        template<class E0,class E1,class F>
        static void func_expr(E0& vt_l,const E1& vt_r,F&& f)
        {
            TOP<__dim__>::func_expr(vt_l,vt_r,f);
        }
        
        template<class E,class F>
        static void func(E& vt_l,const E& vt_r,F&& f)
        {
            TOP<__dim__>::func(vt_l,vt_r,f);
        }
        
        template<class E>
        static void eq(E& vt_l,const E& vt_r)
        {
            TOP<__dim__>::eq(vt_l,vt_r);
        }
        
        template<class E,class F>
        static void func_scl(E& vt,const type0& scl,F&& f)
        {
            TOP<__dim__>::func_scl(vt,scl,f);
        }
    };
    
    template<>
    class TOPHelper<false>
    {
    public:
        template<class E0,class E1,class F>
        static void func_expr(E0&,const E1&,F&&)
        {
        }
        
        template<class E,class F>
        static void func(E&,const E&,F&&)
        {
        }
        
        template<class E>
        static void eq(E&,const E&)
        {
        }
        
        template<class E,class F>
        static void func_scl(E&,const type0&,F&&)
        {
        }
    };
    
    
    template<bool HDOF,const int N>
    class VT : public VTExpr<VT<HDOF,N>>
    {
    private:
        
        void assign_vecs(Atoms* __atoms,Vec<type0>** __vs_arr,bool* __vs_arr_alloc,int* __vs_dim,Vec<type0>* __v,bool __v_dof)
        {
            if(__v_dof)
            {
                *__vs_arr=__v;
                *__vs_arr_alloc=false;
                *__vs_dim=__v->dim;
            }
        }
        
        template<class ...Vs>
        void assign_vecs(Atoms* __atoms,Vec<type0>** __vs_arr,bool* __vs_arr_alloc,int* __vs_dim,Vec<type0>* __v,bool __v_dof,Vs... __vs)
        {
            assign_vecs(__atoms,__vs_arr,__vs_arr_alloc,__vs_dim,__v,__v_dof);
            if(__v_dof)
                assign_vecs(__atoms,__vs_arr+1,__vs_arr_alloc+1,__vs_dim+1,__vs...);
            else
                assign_vecs(__atoms,__vs_arr,__vs_arr_alloc,__vs_dim,__vs...);
            
            
        }
        
        void assign_vecs(Atoms* __atoms,Vec<type0>** __vs_arr,bool* __vs_arr_alloc,int* __vs_dim,int __dim,bool __v_dof)
        {
            if(__v_dof)
            {
                *__vs_arr=new Vec<type0>(__atoms,__dim);
                *__vs_arr_alloc=true;
                *__vs_dim=__dim;
            }
        }
        
        template<class ...Vs>
        void assign_vecs(Atoms* __atoms,Vec<type0>** __vs_arr,bool* __vs_arr_alloc,int* __vs_dim,int __dim,bool __v_dof,Vs... __vs)
        {
            assign_vecs(__atoms,__vs_arr,__vs_arr_alloc,__vs_dim,__dim,__v_dof);
            if(__v_dof)
            
                assign_vecs(__atoms,__vs_arr+1,__vs_arr_alloc+1,__vs_dim+1,__vs...);
            else
                assign_vecs(__atoms,__vs_arr,__vs_arr_alloc,__vs_dim,__vs...);
        }
        
        
        
        void assign_vecs(Atoms* __atoms,Vec<type0>** __vs_arr,bool* __vs_arr_alloc,int* __vs_dim,Vec<type0>** __vs_arr_src,int __dim,bool __v_dof,bool __v_src_dof,bool is_same)
        {
            if(__v_dof)
            {
                if(is_same)
                {
                    // __v_src_dof is definitely true
                    *__vs_arr_alloc=false;
                    *__vs_arr=*__vs_arr_src;
                    *__vs_dim=__dim;
                }
                else
                {
                    // __v_src_dof can be true or false
                    *__vs_arr_alloc=true;
                    *__vs_arr=new Vec<type0>(__atoms,__dim);
                    *__vs_dim=__dim;
                }
            }
        }
        template<class ...Vs>
        void assign_vecs(Atoms* __atoms,Vec<type0>** __vs_arr,bool* __vs_arr_alloc,int* __vs_dim,Vec<type0>** __vs_arr_src,int __dim,bool __v_dof,bool __v_src_dof,bool is_same,Vs... __vs)
        {
            assign_vecs(__atoms,__vs_arr,__vs_arr_alloc,__vs_dim,__vs_arr_src,__dim,__v_dof,__v_src_dof,is_same);
            
            if(__v_dof)
            {
                ++__vs_arr_alloc;
                ++__vs_arr;
                ++__vs_dim;
            }
            
            if(__v_src_dof)
            {
                ++__vs_arr_src;
            }
            
            assign_vecs(__atoms,__vs_arr,__vs_arr_alloc,__vs_dim,__vs_arr_src,__vs...);
            
        }
        
        
        
    protected:
    public:
        type0 A_arr[__dim__][__dim__];
        type0 (*A)[__dim__];
        
        Atoms* atoms;
        Vec<type0>* vecs[N];
        bool vecs_alloc[N];
        int vecs_dims[N];
        
        type0* vec_ptr;
        
        
        
        template<const int I>
        void reset_ptr(){vec_ptr=vecs[I]->begin();}
        type0 operator[](const int& i)const {return vec_ptr[i];}
        template<const int I,const int J>
        type0 get() const{return A[I][J];}
        
        VT()
        {
            Algebra::Do<N>::func([this](int i)
            {
                vecs_alloc[i]=false;
                vecs[i]=NULL;
                vecs_dims[i]=0;
            });
            
            atoms=NULL;
            vec_ptr=NULL;
            A=NULL;
            Algebra::zero<__dim__*__dim__>(A_arr[0]);
        }
        ~VT()
        {
            Algebra::Do<N>::func([this](int i)
            {
                if(vecs_alloc[i]) delete vecs[i];
                vecs_alloc[i]=false;
                vecs[i]=NULL;
                vecs_dims[i]=0;
            });
            
            atoms=NULL;
            vec_ptr=NULL;
            A=NULL;
            Algebra::zero<__dim__*__dim__>(A_arr[0]);
        }

        
        template<class ...Vs>
        VT(Atoms* __atoms,Vs... __vs):
        atoms(__atoms),vec_ptr(NULL)
        {
            
            assign_vecs(__atoms,vecs,vecs_alloc,vecs_dims,__vs...);
            A=A_arr;
            Algebra::zero<__dim__*__dim__>(A_arr[0]);
        }
        
        template<class ...Vs>
        VT(Atoms* __atoms,type0 (&__A)[__dim__][__dim__],Vs... __vs):
        atoms(__atoms),vec_ptr(NULL)
        {
            assign_vecs(__atoms,vecs,vecs_alloc,vecs_dims,__vs...);
            A=__A;
        }
        
        
        template<bool __HDOF,const int __N,class ...Vs>
        VT(Atoms* __atoms,VT<__HDOF,__N>& __src,Vs... __vs):
        atoms(__atoms),vec_ptr(NULL)
        {
            assign_vecs(__atoms,vecs,vecs_alloc,vecs_dims,__src.vecs,__vs...);
            A=A_arr;
            Algebra::zero<__dim__*__dim__>(A_arr[0]);
        }
        
        VT(const VT& r)
        {
            this->atoms=r.atoms;
            const int natms_lcl=this->atoms->natms_lcl;
            Algebra::Do<N>::func([this,&r,&natms_lcl](int i)
            {
                this->vecs_dims[i]=r.vecs_dims[i];
                this->vecs_alloc[i]=true;
                this->vecs[i]=new Vec<type0>(this->atoms,this->vecs_dims[i]);
                memcpy(r.vecs[i]->begin(),this->vecs[i]->begin(),natms_lcl*(this->vecs_dims[i])*sizeof(type0));
            });
            
            this->A=this->A_arr;
            Algebra::V_eq<__dim__*__dim__>(r.A[0],this->A[0]);
        }
        VT(VT&& r)
        {
            this->atoms=r.atoms;
            
            Algebra::Do<N>::func([this,&r](int i)
            {
                this->vecs_dims[i]=r.vecs_dims[i];
                this->vecs_alloc[i]=r.vecs_alloc[i];
                this->vecs[i]=r.vecs[i];
            });
            if(r.A==r.A_arr)
            {
                this->A=this->A_arr;
                Algebra::V_eq<__dim__*__dim__>(r.A_arr[0],this->A_arr[0]);
            }
            else
            {
                Algebra::zero<__dim__*__dim__>(A_arr[0]);
                this->A=r.A;
            }
            
            
            r.~VT();
        }
        VT& operator=(VT&& r)
        {
            this->~VT();
            new (this) VT(std::move(r));
            return *this;
        }
        
        type0 operator*(const VT& r)
        {
            type0 ans_lcl=0.0;
            const int natms_lcl=this->atoms->natms_lcl;
            VOP<N>::func(natms_lcl,*this,r,[&ans_lcl](type0& l,const type0& r){ans_lcl+=l*r;});
            type0 ans;
            MPI_Allreduce(&ans_lcl,&ans,1,Vec<type0>::MPI_T,MPI_SUM,atoms->world);
            TOPHelper<HDOF>::func(*this,r,[&ans](type0& l,const type0& r){ans+=l*r;});
            return ans;
        }
        

        

        VT& operator*=(const type0& scl)
        {
            const int natms_lcl=this->atoms->natms_lcl;
            VOP<N>::func_scl(natms_lcl,*this,scl,[](type0& l,const type0& r){l*=r;});
            TOPHelper<HDOF>::func_scl(*this,scl,[](type0& l,const type0& r){l*=r;});
            return *this;
        }
        VT& operator/=(const type0& scl)
        {
            const int natms_lcl=this->atoms->natms_lcl;
            VOP<N>::func_scl(natms_lcl,*this,scl,[](type0& l,const type0& r){l/=r;});
            TOPHelper<HDOF>::func_scl(*this,scl,[](type0& l,const type0& r){l/=r;});
            return *this;
        }
        
        
        VT& operator=(const VT& r)
        {
            const int natms_lcl=this->atoms->natms_lcl;
            VOP<N>::eq(natms_lcl,*this,r);
            TOPHelper<HDOF>::eq(*this,r);
            return *this;
        }
        
        VT& operator+=(const VT& r)
        {
            const int natms_lcl=this->atoms->natms_lcl;
            VOP<N>::func(natms_lcl,*this,r,[](type0& l,const type0& r){l+=r;});
            TOPHelper<HDOF>::func(*this,r,[](type0& l,const type0& r){l+=r;});
            return *this;
        }
        
        VT& operator-=(const VT& r)
        {
            const int natms_lcl=this->atoms->natms_lcl;
            VOP<N>::func(natms_lcl,*this,r,[](type0& l,const type0& r){l-=r;});
            TOPHelper<HDOF>::func(*this,r,[](type0& l,const type0& r){l-=r;});
            return *this;
        }
        
        
        
        template<class E>
        VT& operator=(VTExpr<E>&& r)
        {
            const int natms_lcl=this->atoms->natms_lcl;
            VOP<N>::func_expr(natms_lcl,*this,r,[](type0& l,const type0& r){l=r;});
            TOPHelper<HDOF>::func_expr(*this,r,[](type0& l,const type0& r){l=r;});
            return *this;
        }
        
        template<class E>
        VT& operator+=(VTExpr<E>&& r)
        {
            const int natms_lcl=this->atoms->natms_lcl;
            VOP<N>::func_expr(natms_lcl,*this,r,[](type0& l,const type0& r){l+=r;});
            TOPHelper<HDOF>::func_expr(*this,r,[](type0& l,const type0& r){l+=r;});
            return *this;
        }
        
        template<class E>
        VT& operator-=(VTExpr<E>&& r)
        {
            const int natms_lcl=this->atoms->natms_lcl;
            VOP<N>::func_expr(natms_lcl,*this,r,[](type0& l,const type0& r){l-=r;});
            TOPHelper<HDOF>::func_expr(*this,r,[](type0& l,const type0& r){l-=r;});
            return *this;
        }
    };
    
    
    
    
    
    template<bool HDOF>
    class VT <HDOF,0>: public VTExpr<VT<HDOF,0>>
    {
    private:
        
    protected:
    public:
        type0 A_arr[__dim__][__dim__];
        type0 (*A)[__dim__];
        
        Atoms* atoms;
        Vec<type0>** vecs;
        
        
        template<const int I>
        void reset_ptr(){}
        type0 operator[](const int& i)const {return 0.0;}
        template<const int I,const int J>
        type0 get() const{return A[I][J];}
        
        VT():
        atoms(NULL),
        vecs(NULL)
        {

            A=NULL;
            Algebra::zero<__dim__*__dim__>(A_arr[0]);
        }
        ~VT()
        {
            A=NULL;
            Algebra::zero<__dim__*__dim__>(A_arr[0]);
        }

        
        template<class ...Vs>
        VT(Atoms* __atoms,Vs... __vs):
        atoms(NULL),
        vecs(NULL)
        {
            A=A_arr;
        }
        
        template<class ...Vs>
        VT(Atoms* __atoms,type0 (&__A)[__dim__][__dim__],Vs... __vs):
        atoms(NULL),
        vecs(NULL)
        {
            A=__A;
            Algebra::zero<__dim__*__dim__>(A_arr[0]);
        }
        
        
        VT(const VT& r):
        atoms(NULL),
        vecs(NULL)
        {
            this->A=this->A_arr;
            Algebra::V_eq<__dim__*__dim__>(r.A[0],this->A[0]);
            
        }
        VT(VT&& r):
        atoms(NULL),
        vecs(NULL)
        {
            if(r.A==r.A_arr)
            {
                this->A=this->A_arr;
                Algebra::V_eq<__dim__*__dim__>(r.A_arr[0],this->A_arr[0]);
            }
            else
            {
                Algebra::zero<__dim__*__dim__>(A_arr[0]);
                this->A=r.A;
            }
            
            
            r.~VT();
        }
        VT& operator=(VT&& r)
        {
            this->~VT();
            new (this) VT(std::move(r));
            return *this;
        }
        
        type0 operator*(const VT& r)
        {
            type0 ans=0.0;
            TOPHelper<HDOF>::func(*this,r,[&ans](type0& l,const type0& r){ans+=l*r;});
            return ans;
        }
        
        VT& operator*=(const type0& scl)
        {
            TOPHelper<HDOF>::func_scl(*this,scl,[](type0& l,const type0& r){l*=r;});
            return *this;
        }
        VT& operator/=(const type0& scl)
        {
            TOPHelper<HDOF>::func_scl(*this,scl,[](type0& l,const type0& r){l/=r;});
            return *this;
        }
        
        
        VT& operator=(const VT& r)
        {
            TOPHelper<HDOF>::eq(*this,r);
            return *this;
        }
        
        VT& operator+=(const VT& r)
        {
            TOPHelper<HDOF>::func(*this,r,[](type0& l,const type0& r){l+=r;});
            return *this;
        }
        
        VT& operator-=(const VT& r)
        {
            TOPHelper<HDOF>::func(*this,r,[](type0& l,const type0& r){l-=r;});
            return *this;
        }
        
        
        
        template<class E>
        VT& operator=(const VTExpr<E>& r)
        {
            TOPHelper<HDOF>::func_expr(*this,r,[](type0& l,const type0& r){l=r;});
            return *this;
        }
        
        template<class E>
        VT& operator+=(const VTExpr<E>& r)
        {
            TOPHelper<HDOF>::func_expr(*this,r,[](type0& l,const type0& r){l+=r;});
            return *this;
        }
        
        template<class E>
        VT& operator-=(const VTExpr<E>& r)
        {
            TOPHelper<HDOF>::func_expr(*this,r,[](type0& l,const type0& r){l-=r;});
            return *this;
        }
    };
    
    
    template<class E_l,template<class>class VTExpr_l,class E_r,template<class>class VTExpr_r>
    VTSum<VTExpr_l<E_l>,VTExpr_r<E_r>> operator+(VTExpr_l<E_l>&& vt_l,VTExpr_r<E_r>&& vt_r)
    {
        return VTSum<VTExpr_l<E_l>,VTExpr_r<E_r>>(vt_l,vt_r);
    }
    
    template<class E_l,template<class>class VTExpr_l,bool HDOF_r,const int N_r>
    VTSum<VTExpr_l<E_l>,VT<HDOF_r,N_r>> operator+(VTExpr_l<E_l>&& vt_l,VT<HDOF_r,N_r>& vt_r)
    {
        return VTSum<VTExpr_l<E_l>,VT<HDOF_r,N_r>>(vt_l,vt_r);
    }
    
    template<bool HDOF_l,const int N_l,class E_r,template<class>class VTExpr_r>
    VTSum<VT<HDOF_l,N_l>,VTExpr_r<E_r>> operator+(VT<HDOF_l,N_l>& vt_l,VTExpr_r<E_r>&& vt_r)
    {
        return VTSum<VT<HDOF_l,N_l>,VTExpr_r<E_r>>(vt_l,vt_r);
    }
    
    template<bool HDOF_l,const int N_l,bool HDOF_r,const int N_r>
    VTSum<VT<HDOF_l,N_l>,VT<HDOF_r,N_r>> operator+(VT<HDOF_l,N_l>& vt_l,VT<HDOF_r,N_r>& vt_r)
    {
        return VTSum<VT<HDOF_l,N_l>,VT<HDOF_r,N_r>>(vt_l,vt_r);
    }
    
    template<class E_l,template<class>class VTExpr_l,class E_r,template<class>class VTExpr_r>
    VTSub<VTExpr_l<E_l>,VTExpr_r<E_r>> operator-(VTExpr_l<E_l>&& vt_l,VTExpr_r<E_r>&& vt_r)
    {
        return VTSub<VTExpr_l<E_l>,VTExpr_r<E_r>>(vt_l,vt_r);
    }
    
    template<class E_l,template<class>class VTExpr_l,bool HDOF_r,const int N_r>
    VTSub<VTExpr_l<E_l>,VT<HDOF_r,N_r>> operator-(VTExpr_l<E_l>&& vt_l,VT<HDOF_r,N_r>& vt_r)
    {
        return VTSub<VTExpr_l<E_l>,VT<HDOF_r,N_r>>(vt_l,vt_r);
    }
    
    template<bool HDOF_l,const int N_l,class E_r,template<class>class VTExpr_r>
    VTSub<VT<HDOF_l,N_l>,VTExpr_r<E_r>> operator-(VT<HDOF_l,N_l>& vt_l,VTExpr_r<E_r>&& vt_r)
    {
        return VTSub<VT<HDOF_l,N_l>,VTExpr_r<E_r>>(vt_l,vt_r);
    }
    
    template<bool HDOF_l,const int N_l,bool HDOF_r,const int N_r>
    VTSub<VT<HDOF_l,N_l>,VT<HDOF_r,N_r>> operator-(VT<HDOF_l,N_l>& vt_l,VT<HDOF_r,N_r>& vt_r)
    {
        return VTSub<VT<HDOF_l,N_l>,VT<HDOF_r,N_r>>(vt_l,vt_r);
    }
    
    template<class E,template<class>class VExpr>
    VTMulScl<VExpr<E>> operator*(type0& scl,VExpr<E>& vt)
    {
        return VTMulScl<VExpr<E>>(scl,vt);
    }
    
    template<bool H_DOF,int N>
    VTMulScl<VT<H_DOF,N>> operator*(type0& scl,VT<H_DOF,N>& vt)
    {
        return VTMulScl<VT<H_DOF,N>>(scl,vt);
    }
    
    template<bool H_DOF,int N>
    void cyclic_shift(VT<H_DOF,N>* vs,int n)
    {
        Vec<type0>* __vecs[N];
        Algebra::V_eq<N>(vs[n-1].vecs,__vecs);
        for(int i=n-1;i>0;i--)
            Algebra::V_eq<N>(vs[i-1].vecs,vs[i].vecs);
        Algebra::V_eq<N>(__vecs,vs[0].vecs);
    }
    
    template<int N>
    void cyclic_shift(VT<true,N>* vs,int n)
    {
        Vec<type0>* __vecs[N];
        type0 __A[__dim__][__dim__];
        Algebra::V_eq<N>(vs[n-1].vecs,__vecs);
        Algebra::V_eq<__dim__*__dim__>(vs[n-1].A[0],__A[0]);
        
        for(int i=n-1;i>0;i--)
        {
            Algebra::V_eq<N>(vs[i-1].vecs,vs[i].vecs);
            Algebra::V_eq<__dim__*__dim__>(vs[i-1].A[0],vs[i].A[0]);
        }
        Algebra::V_eq<N>(__vecs,vs[0].vecs);
        Algebra::V_eq<__dim__*__dim__>(__A[0],vs[0].A[0]);
    }
    
    template<bool H_DOF>
    void cyclic_shift(VT<H_DOF,0>* ,int)
    {}
    
    inline void cyclic_shift(VT<true,0>* vs,int n)
    {
        
        type0 __A[__dim__][__dim__];
        Algebra::V_eq<__dim__*__dim__>(vs[n-1].A[0],__A[0]);
        
        for(int i=n-1;i>0;i--)
        {
            Algebra::V_eq<__dim__*__dim__>(vs[i-1].A[0],vs[i].A[0]);
        }
        Algebra::V_eq<__dim__*__dim__>(__A[0],vs[0].A[0]);
    }
    
    
}


/*--------------------------------------------
   _____   __  __   _____    ______    _____
  / ____| |  \/  | |  __ \  |  ____|  / ____|
 | |  __  | \  / | | |__) | | |__    | (___
 | | |_ | | |\/| | |  _  /  |  __|    \___ \
 | |__| | | |  | | | | \ \  | |____   ____) |
  \_____| |_|  |_| |_|  \_\ |______| |_____/
 --------------------------------------------*/
namespace MAPP_NS
{
    template<class V>
    class GMRES
    {
    private:
        const int m;
        int N;
        V* Q;
        type0** A_hat;
        type0* Ax_hat;
        type0(* cos_sin)[2];
        
        type0* x_hat;
        
        
        type0 calc(V&,V&);
        type0 calc(int,V&,V&);
        type0 solve_y(int,V&);
        
        V RES;
        V dx;
        
    protected:
    public:
        template<class ... Cs>
        GMRES(int,Cs ...);
        template<class ... Cs>
        GMRES(int,int,Cs ...);
        ~GMRES();
        template<class KERNEL>
        bool solve(KERNEL&,V&,type0,type0&,V&);
        template<class KERNEL>
        bool solve_restart(KERNEL&,V&,type0,type0&,V&);
        int iter;
        type0 res;
        
        void res_calc(V&);
        
    };
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<class V>template<class ... Cs>
GMRES<V>::GMRES(int __m,Cs ... cs):
m(__m),
N(0),
RES(cs...),
dx(cs...),
iter(0),
res(0.0)
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
template<class V>template<class ... Cs>
GMRES<V>::GMRES(int __N,int __m,Cs ... cs):
GMRES<V>(__m,cs ...)
{
    N=__N;
    
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<class V>
GMRES<V>::~GMRES()
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
type0 GMRES<V>::calc(V& Ax,V& x)
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
type0 GMRES<V>::calc(int iter,V& Ax,V& x)
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
type0 GMRES<V>::solve_y(int nvecs,V& x)
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
bool GMRES<V>::solve(KERNEL& A,V& b,type0 tol,type0& norm,V& x)
{
    res=calc(b,x);
    for(iter=0;iter<m;iter++)
    {
        A(x,b);
        res=calc(iter,b,x);
        if(res<tol)
        {
            norm=solve_y(iter+1,x);
            iter++;
            return true;
        }
    }
    
    norm=solve_y(m,x);
    return false;
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<class V>template<class KERNEL>
bool GMRES<V>::solve_restart(KERNEL& A,V& b,type0 tol,type0& norm,V& x)
{
    bool conv=solve(A,b,tol,norm,x);
    if(conv) return true;
    if(N==0) return false;
    while(!conv)
    {
        res_calc(RES);
        conv=solve(A,RES,tol,norm,dx);
        x+=dx;
    }
    
    
    
    return true;
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<class V>
void GMRES<V>::res_calc(V& RES)
{
    
    type0 res_norm=Ax_hat[iter];
    
    RES=(res_norm*cos_sin[iter-1][0])*Q[iter];
    res_norm*=cos_sin[iter-1][1];
    for(int i=iter-1;i>0;i--)
    {
        RES+=(res_norm*cos_sin[i-1][0])*Q[i];
        res_norm*=cos_sin[i-1][1];
    }
    
    RES+=res_norm*Q[0];
    
  
    

}

#endif
