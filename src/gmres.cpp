#include "gmres.h"
#include "memory.h"
/*--------------------------------------------
 
 --------------------------------------------*/
__GMRES::__GMRES(Atoms* __atoms,int __m,int __dim):
world(__atoms->world),
m(__m),
dim(__dim),
n(__dim*__atoms->natms)
{
    
    
    H=new type0*[m+1];
    *H=new type0[m*(m+1)/2+1];
    for(int i=1;i<m+1;i++)
        H[i]=H[i-1]+i;
    
    Memory::alloc(b_hat,m+1);
    Memory::alloc(cos,m+1);
    Memory::alloc(sin,m+1);
    Memory::alloc(y,m+1);
    Memory::alloc(ans_lcl,m+1);
    Memory::alloc(Q,(m+1)*n);
}
/*--------------------------------------------
 
 --------------------------------------------*/
__GMRES::~__GMRES()
{
    Memory::dealloc(Q);
    Memory::dealloc(ans_lcl);
    Memory::dealloc(y);
    Memory::dealloc(sin);
    Memory::dealloc(cos);
    Memory::dealloc(b_hat);
    delete [] *H;
    delete [] H;
}
/*--------------------------------------------
 
 --------------------------------------------*/
type0 __GMRES::calc(type0* Ax,type0* x)
{
    type0 norm_sq_lcl=0.0,norm;
    for(int i=0;i<n;i++)
        norm_sq_lcl+=Ax[i]*Ax[i];
    
    MPI_Allreduce(&norm_sq_lcl,&norm,1,MPI_TYPE0,MPI_SUM,world);
    norm=sqrt(norm);
    type0 norm_inv=1.0/norm;
    
    type0* q=Q;
    for(int i=0;i<n;i++,q+=m+1)
        x[i]=*q=norm_inv*Ax[i];
    
    b_hat[0]=norm;
    return norm;
}
/*--------------------------------------------
 
 --------------------------------------------*/
type0 __GMRES::calc(int iter,type0* Ax,type0* x)
{
    type0* h=H[iter];
    int ivec=iter+1;
    for(int j=0;j<ivec;j++)
        ans_lcl[j]=0.0;
    
    type0* q=Q;
    for(int i=0;i<n;i++)
    {
        type0 __ax=Ax[i];
        for(int j=0;j<ivec;j++)
            ans_lcl[j]+=__ax*q[j];
        q+=m+1;
    }
    
    
    MPI_Allreduce(ans_lcl,h,ivec,MPI_TYPE0,MPI_SUM,world);
    type0 norm_sq_lcl=0.0,norm;
    
    q=Q;
    for(int i=0;i<n;i++,q+=m+1)
    {
        q[ivec]=Ax[i];
        for(int j=0;j<ivec;j++)
            q[ivec]-=h[j]*q[j];
        norm_sq_lcl+=q[ivec]*q[ivec];
    }
    
    MPI_Allreduce(&norm_sq_lcl,&norm,1,MPI_TYPE0,MPI_SUM,world);
    norm=sqrt(norm);
    type0 norm_inv=1.0/norm;
    
    
    q=Q+ivec;
    for(int i=0;i<n;i++,q+=m+1)
        x[i]=*q*=norm_inv;
    
    
    type0 tmp;
    for(int i=0;i<iter;i++)
    {
        tmp=cos[i]*h[i]-sin[i]*h[i+1];
        h[i+1]=sin[i]*h[i]+cos[i]*h[i+1];
        h[i]=tmp;
    }
    
    tmp=sqrt(h[iter]*h[iter]+norm*norm);
    cos[iter]=h[iter]/tmp;
    sin[iter]=-norm/tmp;
   
    h[iter]=cos[iter]*h[iter]-sin[iter]*norm;
    b_hat[iter+1]=sin[iter]*b_hat[iter];
    b_hat[iter]*=cos[iter];

    
    
    return fabs(b_hat[iter+1]);
    
}
/*--------------------------------------------
 
 --------------------------------------------*/
type0 __GMRES::solve_y(int nvecs,type0* x)
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
     printf("};\n\n");
    */
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
    
    
    type0* q=Q;
    for(int i=0;i<n;i++,q+=m+1)
    {
        x[i]=0.0;
        for(int ivec=0;ivec<nvecs;ivec++)
            x[i]+=y[ivec]*q[ivec];
    }
    
    type0 norm=0.0;
    for(int ivec=0;ivec<nvecs;ivec++)
        norm+=y[ivec]*y[ivec];
    
    return sqrt(norm);
}














