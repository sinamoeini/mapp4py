#include "gmres.h"
#include "memory.h"
/*--------------------------------------------
 
 --------------------------------------------*/
GMRES::GMRES(Atoms* __atoms,int __m,int __dim):
m(__m),
n(__dim*__atoms->natms_lcl),
world(__atoms->world)
{
    
    
    A_hat=new type0*[m];
    *A_hat=new type0[m*(m+1)/2];
    for(int i=1;i<m;i++)
        A_hat[i]=A_hat[i-1]+i;
    
    Memory::alloc(Ax_hat,m+1);
    Memory::alloc(cos_sin,m+1);
    Memory::alloc(x_hat,m+1);
    Memory::alloc(Q,(m+1)*n);
}
/*--------------------------------------------
 
 --------------------------------------------*/
GMRES::~GMRES()
{
    Memory::dealloc(Q);
    Memory::dealloc(x_hat);
    Memory::dealloc(cos_sin);
    Memory::dealloc(Ax_hat);
    delete [] *A_hat;
    delete [] A_hat;
}
/*--------------------------------------------
 
 --------------------------------------------*/
type0 GMRES::calc(type0* Ax,type0* x)
{
    type0 norm_sq_lcl=0.0,norm;
    for(int i=0;i<n;i++)
        norm_sq_lcl+=Ax[i]*Ax[i];
    
    MPI_Allreduce(&norm_sq_lcl,&norm,1,Vec<type0>::MPI_T,MPI_SUM,world);
    norm=sqrt(norm);
    type0 norm_inv=1.0/norm;
    
    type0* q=Q;
    for(int i=0;i<n;i++,q+=m+1)
        x[i]=*q=norm_inv*Ax[i];
    
    Ax_hat[0]=norm;
    return norm;
}
/*--------------------------------------------
 
 --------------------------------------------*/
type0 GMRES::calc(int iter,type0* RESTRICT Ax,type0* RESTRICT x)
{
    type0* h=A_hat[iter];
    int ivec=iter+1;
    for(int j=0;j<ivec;j++)
        x_hat[j]=0.0;
    
    type0* q=Q;
    for(int i=0;i<n;i++)
    {
        type0 __ax=Ax[i];
        for(int j=0;j<ivec;j++)
            x_hat[j]+=__ax*q[j];
        q+=m+1;
    }
    
    
    MPI_Allreduce(x_hat,h,ivec,Vec<type0>::MPI_T,MPI_SUM,world);
    type0 norm_sq_lcl=0.0,norm;
    
    q=Q;
    for(int i=0;i<n;i++,q+=m+1)
    {
        q[ivec]=Ax[i];
        for(int j=0;j<ivec;j++)
            q[ivec]-=h[j]*q[j];
        norm_sq_lcl+=q[ivec]*q[ivec];
    }
    
    MPI_Allreduce(&norm_sq_lcl,&norm,1,Vec<type0>::MPI_T,MPI_SUM,world);
    norm=sqrt(norm);
    type0 norm_inv=1.0/norm;
    
    
    q=Q+ivec;
    for(int i=0;i<n;i++,q+=m+1)
        x[i]=*q*=norm_inv;
    
    
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

    
    
    return fabs(Ax_hat[iter+1]);
    
}
/*--------------------------------------------
 
 --------------------------------------------*/
type0 GMRES::solve_y(int nvecs,type0* x)
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
        x_hat[i]=Ax_hat[i];
        for(int j=i+1;j<nvecs;j++)
            x_hat[i]-=A_hat[j][i]*x_hat[j];
        x_hat[i]/=A_hat[i][i];
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
            x[i]+=x_hat[ivec]*q[ivec];
    }
    
    type0 norm=0.0;
    for(int ivec=0;ivec<nvecs;ivec++)
        norm+=x_hat[ivec]*x_hat[ivec];
    
    return sqrt(norm);
}














