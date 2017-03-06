#ifndef __MAPP__gmres__
#define __MAPP__gmres__
#include "atoms.h"
namespace MAPP_NS
{
    template<class T,class C0>
    class __GMRES
    {
    private:
        const int m;
        const int dim;
        int n;
        C0& kernel;
        Vec<T>** vecs;
        T** Q;
        T** H;
        T* b_hat;
        Atoms* atoms;
        T* cos;
        T* sin;
        
        T* y;
        T* ans_lcl;
        
        T calc(int);
        T solve_y(int,type0*);
        void restart();
        void start(T*);
        MPI_Comm& world;
    protected:
    public:
        __GMRES(Atoms*,int,int,C0&);
        ~__GMRES();
        void refresh();
        bool solve(T,Vec<T>*,Vec<T>*,T&,T&);
        bool solve(Vec<T>*,T,T&,Vec<T>*);
        
    };
}
using namespace MAPP_NS;
/*--------------------------------------------
 
 --------------------------------------------*/
template<class T,class C0>
__GMRES<T,C0>::__GMRES(Atoms* __atoms,int __m,int __dim,C0& __kernel):
atoms(__atoms),
world(__atoms->world),
m(__m),
dim(__dim),
kernel(__kernel)
{
    
    Q=new T*[m+1];
    vecs=new Vec<T>*[m+1];
    for(int ivec=0;ivec<m+1;ivec++)
        vecs[ivec]=new Vec<T>(atoms,dim);
    
    H=new T*[m+1];
    *H=new T[m*(m+1)/2+1];
    for(int i=1;i<m+1;i++)
        H[i]=H[i-1]+i;
    
    refresh();
    
    b_hat=new T[m+1];
    cos=new T[m];
    sin=new T[m];
    y=new T[m];
    ans_lcl=new T[m+1];
    
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<class T,class C0>
__GMRES<T,C0>::~__GMRES()
{
    for(int ivec=0;ivec<m+1;ivec++)
        delete vecs[ivec];
    delete [] vecs;
    delete [] Q;
    
    delete [] *H;
    delete [] H;
    delete [] b_hat;
    delete [] cos;
    delete [] sin;
    delete [] y;
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<class T,class C0>
void __GMRES<T,C0>::refresh()
{
    for(int ivec=0;ivec<m+1;ivec++)
        Q[ivec]=vecs[ivec]->begin();
    n=atoms->natms_lcl*dim;
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<class T,class C0>
bool __GMRES<T,C0>::solve(Vec<T>* b_ptr,T tol,T& norm,Vec<T>* x_ptr)
{
    type0* x=x_ptr->begin();
    type0* b=b_ptr->begin();
    for(int i=0;i<n;i++) x[i]=0.0;
    T last_h=0.0,tmp_0,b_norm,res_norm;
    start(b);
    res_norm=b_norm=b_hat[0];
    for(int i=0;i<m;i++)
    {
        kernel(vecs[i],vecs[i+1]);
        
        last_h=calc(i);
        
        for(int j=0;j<i;j++)
        {
            tmp_0=cos[j]*H[i][j]-sin[j]*H[i][j+1];
            H[i][j+1]=sin[j]*H[i][j]+cos[j]*H[i][j+1];
            H[i][j]=tmp_0;
        }
        
        tmp_0=sqrt(H[i][i]*H[i][i]+last_h*last_h);
        
        cos[i]=H[i][i]/tmp_0;
        sin[i]=-last_h/tmp_0;
        
        H[i][i]=cos[i]*H[i][i]-sin[i]*last_h;
        b_hat[i+1]=sin[i]*b_hat[i];
        b_hat[i]*=cos[i];
        
        res_norm=fabs(b_hat[i+1]);
        if(res_norm<tol)
        {
            norm=solve_y(i+1,x);
            return true;
        }
    }
    
    norm=solve_y(m,x);
    return false;
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<class T,class C0>
bool __GMRES<T,C0>::solve(T tol,Vec<T>* b,Vec<T>* x,T& res_norm,T& norm)
{
    for(int i=0;i<n;i++) x[i]=0.0;
    T last_h=0.0,tmp_0,b_norm;
    int max_iter=1;
    start(b);
    
    res_norm=b_norm=b_hat[0];
    while (res_norm>tol && max_iter)
    {
        for(int i=0;i<m;i++)
        {
            kernel(vecs[i],vecs[i+1]);

            
            last_h=calc(i);
            
            for(int j=0;j<i;j++)
            {
                tmp_0=cos[j]*H[i][j]-sin[j]*H[i][j+1];
                H[i][j+1]=sin[j]*H[i][j]+cos[j]*H[i][j+1];
                H[i][j]=tmp_0;
            }
            
            tmp_0=sqrt(H[i][i]*H[i][i]+last_h*last_h);
            
            cos[i]=H[i][i]/tmp_0;
            sin[i]=-last_h/tmp_0;
            
            H[i][i]=cos[i]*H[i][i]-sin[i]*last_h;
            b_hat[i+1]=sin[i]*b_hat[i];
            b_hat[i]*=cos[i];
            
            res_norm=fabs(b_hat[i+1]);
            
            if(res_norm<tol)
            {
                norm=solve_y(i+1,x);
                return true;
            }
        }
        
        
        
        norm=solve_y(m,x);
        
        max_iter--;
        
        if(max_iter==0)
            return fabs(b_hat[m]);
        restart();
    }
    
    
    return false;
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<class T,class C0>
T __GMRES<T,C0>::calc(int ivec)
{
    for(int j=0;j<ivec+1;j++)
        ans_lcl[j]=0.0;
    
    for(int i=0;i<n;i++)
        for(int j=0;j<ivec+1;j++)
            ans_lcl[j]+=Q[ivec+1][i]*Q[j][i];
    
    MPI_Allreduce(ans_lcl,H[ivec],ivec+1,Vec<type0>::MPI_T,MPI_SUM,world);
    
    T ans_lcl_=0.0,ans;
    for(int i=0;i<n;i++)
    {
        for(int j=0;j<ivec+1;j++)
            Q[ivec+1][i]-=H[ivec][j]*Q[j][i];
        
        ans_lcl_+=Q[ivec+1][i]*Q[ivec+1][i];
    }
    
    MPI_Allreduce(&ans_lcl_,&ans,1,Vec<type0>::MPI_T,MPI_SUM,world);
    ans=1.0/sqrt(ans);
    
    for(int i=0;i<n;i++)
        Q[ivec+1][i]*=ans;
    
    return 1.0/ans;
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<class T,class C0>
T __GMRES<T,C0>::solve_y(int nvecs,type0* x)
{
    /*
    printf("H={");
    for(int i=0;i<nvecs;i++)
    {
        printf("{");
        for(int j=0;j<nvecs;j++)
        {
            if(i>j)
                printf("0.0");
            else
                printf("%0.10lf",H[j][i]);
            if(j!=nvecs-1)
                printf(",");
        }
        printf("}");
        
        if(i!=nvecs-1)
            printf(",");
    }
    printf("};\n\n");
    printf("b={");
    for(int i=0;i<nvecs;i++)
    {
        printf("%0.10lf",b_hat[i]);
        if(i!=nvecs-1)
            printf(",");
    }
    printf("};\n\n");*/
        
    for(int i=nvecs-1;i>-1;i--)
    {
        y[i]=b_hat[i];
        for(int j=i+1;j<nvecs;j++)
            y[i]-=H[j][i]*y[j];
        y[i]/=H[i][i];
    }
    /*
    printf("y={");
    for(int i=0;i<nvecs;i++)
    {
        printf("%0.10lf",y[i]);
        if(i!=nvecs-1)
            printf(",");
    }
    printf("};\n\n");*/
    
    /*
    printf("r={");
    for(int i=0;i<nvecs;i++)
    {
        type0 r=b_hat[i];
        for(int j=i;j<nvecs;j++)
            r-=H[j][i]*y[j];
        printf("%0.10lf",r);
        if(i!=nvecs-1)
            printf(",");
    }
    printf("};\n\n");
    */
    
    
    
    for(int i=0;i<n;i++)
        for(int ivec=0;ivec<nvecs;ivec++)
            x[i]+=y[ivec]*Q[ivec][i];
    
    T norm=0.0;
    for(int ivec=0;ivec<nvecs;ivec++)
        norm+=y[ivec]*y[ivec];
    
    return sqrt(norm);
}

/*--------------------------------------------
 
 --------------------------------------------*/
template<class T,class C0>
void __GMRES<T,C0>::restart()
{
    T tmp_0=fabs(b_hat[m]);
    
    
    if(b_hat[m]>0.0)
        b_hat[m]=1.0;
    else
        b_hat[m]=-1.0;
    
    for(int i=m-1;i>-1;i--)
    {
        b_hat[i]=sin[i]*b_hat[i+1];
        b_hat[i+1]*=cos[i];
    }
    
    for(int i=0;i<n;i++)
    {
        Q[0][i]*=b_hat[0];
        for(int ivec=1;ivec<m+1;ivec++)
            Q[0][i]+=b_hat[ivec]*Q[ivec][i];
    }
    
    b_hat[0]=tmp_0;
    for(int i=1;i<m+1;i++)
        b_hat[i]=0.0;
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<class T,class C0>
void __GMRES<T,C0>::start(T* b)
{
    T ans_lcl=0.0,norm,inv_norm;
    for(int i=0;i<n;i++)
        ans_lcl+=b[i]*b[i];
    
    MPI_Allreduce(&ans_lcl,&norm,1,Vec<type0>::MPI_T,MPI_SUM,world);
    norm=sqrt(norm);
    inv_norm=1.0/norm;
    
    for(int i=0;i<n;i++)
        Q[0][i]=inv_norm*b[i];
    
    b_hat[0]=norm;
}

/*--------------------------------------------
 
 --------------------------------------------*/
namespace MAPP_NS
{
    class GMRES
    {
    private:
        const int m;
        const int n;
        type0* Q;
        type0** A_hat;
        type0* Ax_hat;
        type0(* cos_sin)[2];
        
        
        type0* x_hat;
        

        type0 calc(type0*,type0*);
        type0 calc(int,type0*,type0*);
        
        type0 solve_y(int,type0*);
        
        
        MPI_Comm& world;
    protected:
    public:
        GMRES(Atoms*,int,int);
        ~GMRES();
        template<class KERNEL>
        bool solve(KERNEL&,Vec<type0>*,type0,type0&,Vec<type0>*);
        
    };
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<class KERNEL>
bool GMRES::solve(KERNEL& A,Vec<type0>* Ax,type0 tol,type0& norm,Vec<type0>* x)
{
    type0* __Ax=Ax->begin();
    type0* __x=x->begin();
    calc(__Ax,__x);
    
    for(int i=0;i<m;i++)
    {
        A(x,Ax);
        if(calc(i,__Ax,__x)<tol)
        {
            norm=solve_y(i+1,__x);
            return true;
        }
    }
    
    norm=solve_y(m,__x);
    return false;
}














#endif
