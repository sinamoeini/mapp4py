#include "xmath.h"
#define TOLERANCE 1.0e-10
#include <cstring>
#include "memory.h"
using namespace MAPP_NS;
/*--------------------------------------------
 return the list of all possible groups of 
 integers that their products are equal to a 
 specific number
 --------------------------------------------*/
void XMath::fac_list(int no,int dim,int*& list,int& list_sz)
{
    list_sz=0;
    list=NULL;
    int* tmp_list=new int[dim];
    fac_list_rec(no,dim,0,tmp_list,list,list_sz);
    delete [] tmp_list;
}
/*--------------------------------------------
 
 --------------------------------------------*/
void XMath::fac_list_rec(int no,int idim,int ipos,
int*& tmp_list,int*& list,int& list_sz)
{
    
    if(idim>1)
    {
        for(int i=1;i<=no;i++)
            if(no%i==0)
            {
                tmp_list[ipos]=i;
                fac_list_rec(no/i,idim-1,ipos+1,tmp_list,list,list_sz);
            }
    }
    else
    {
        tmp_list[ipos]=no;
        int* list_=new int[(list_sz+1)*(ipos+1)];
        memcpy(list_,list,list_sz*(ipos+1)*sizeof(int));
        memcpy(list_+(ipos+1)*list_sz,tmp_list,(ipos+1)*sizeof(int));
        delete [] list;
        list=list_;
        list_sz++;
    }
}
/*--------------------------------------------
 autogrid the domain
 --------------------------------------------*/
void XMath::square2lo_tri(type0(&H_old)[3][3],type0(&H_new)[3][3])
{
    type0** Q=new type0*[3];
    for(int idim=0;idim<3;idim++)
        Q[idim]=new type0[3];

    type0 H0H0;
    type0 H0H1;
    H0H0=H0H1=0.0;
    for(int i=0;i<3;i++)
    {
        H0H0+=H_old[0][i]*H_old[0][i];
        H0H1+=H_old[0][i]*H_old[1][i];
    }
    Q[0][0]=H_old[0][0];
    Q[0][1]=H_old[0][1];
    Q[0][2]=H_old[0][2];
    Q[1][0]=H0H0*H_old[1][0]-H0H1*H_old[0][0];
    Q[1][1]=H0H0*H_old[1][1]-H0H1*H_old[0][1];
    Q[1][2]=H0H0*H_old[1][2]-H0H1*H_old[0][2];
    Q[2][0]=H_old[0][1]*H_old[1][2]-H_old[0][2]*H_old[1][1];
    Q[2][1]=H_old[0][2]*H_old[1][0]-H_old[0][0]*H_old[1][2];
    Q[2][2]=H_old[0][0]*H_old[1][1]-H_old[0][1]*H_old[1][0];
    for(int i=0;i<3;i++)
    {
        H0H0=0.0;
        for(int j=0;j<3;j++)
            H0H0+=Q[i][j]*Q[i][j];
        
        H0H1=sqrt(H0H0);
        for(int j=0;j<3;j++)
            Q[i][j]/=H0H1;
    }
    
    for(int i=0;i<3;i++)
        for(int j=0;j<3;j++)
            H_new[i][j]=0.0;
    
    for(int i=0;i<3;i++)
        for(int j=0;j<i+1;j++)
        {
            for(int k=0;k<3;k++)
                H_new[i][j]+=H_old[i][k]*Q[j][k];
        }
    
    for(int i=0;i<3;i++)
        delete [] Q[i];
    delete [] Q;
}
/*--------------------------------------------
 inversion funtion to calculate B whenever H
 is changed
 --------------------------------------------*/
void XMath::invert(type0** A,type0** Ainv,int dim)
{
    if(dim==0)
        return;
    
    
    type0** ATA=new type0*[dim];
    for(int idim=0;idim<dim;idim++)
        ATA[idim]=new type0[dim];
    type0* c=new type0[dim];
    type0* x=new type0[dim];
    type0* g=new type0[dim];
    type0* g0=new type0[dim];
    type0* h=new type0[dim];
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
    
    for(int i=0;i<dim;i++)
        delete [] ATA[i];
    delete [] ATA;
    delete [] c;
    delete [] x;
    delete [] g;
    delete [] g0;
    delete [] h;
}
/*--------------------------------------------
 inversion funtion to calculate B whenever H
 is changed
 test: 
     XMath* xmath=new XMath();
    
    int dim=4;
    type0** A;
    type0** A_inv;
    type0** I;
    CREATE2D(A,dim,dim);
    CREATE2D(A_inv,dim,dim);
    CREATE2D(I,dim,dim);
    
    for(int i=0;i<dim;i++)
        for(int j=0;j<dim;j++)
            A[i][j]=0.0;

    A[0][0]=1.25;
    A[1][0]=23.45; A[1][1]=0.87;
    A[2][0]=33.9;  A[2][1]=52.08; A[2][2]=7.85;
    A[3][0]=1.025; A[3][1]=75.9;  A[3][2]=9.06; A[3][3]=12.35;
    
    xmath->invert_lower_triangle(A,A_inv,dim);
    
    for(int i=0;i<dim;i++)
    {
        for(int j=0;j<dim;j++)
        {
            I[i][j]=0.0;
            for(int k=0;k<dim;k++)
                I[i][j]+=A[i][k]*A_inv[k][j];
            
            printf("%lf ",I[i][j]);
        }
        printf("\n");
    }
    
    DEL_2D(A);
    DEL_2D(A_inv);
    DEL_2D(I);
 --------------------------------------------*/
/*
void XMath::invert_lower_triangle(type0**& A,type0**& A_inv,int& dim)
{
    for(int i=0;i<dim;i++)
    {
        A_inv[i][i]=1.0/A[i][i];
        for(int j=i+1;j<dim;j++)
            A_inv[i][j]=0.0;
        for(int j=i-1;j>-1;j--)
        {
            A_inv[i][j]=0.0;
            for(int k=j;k<i;k++)
                A_inv[i][j]-=A[i][k]*A_inv[k][j];
            A_inv[i][j]*=A_inv[i][i];
        }
    }
}*/
/*--------------------------------------------
 construct a legendre polynomial of degree n
 --------------------------------------------*/
void XMath::quadrature_lg(int n,type0* x,type0* w)
{
    int m=n/2+1;
    int iter,ord,icurs,jcurs;
    int max_iter=50;
    type0 a,u0,inv_u0,f,up,df,tmp0,tol;
    type0 ii,jj,del_u0;
    type0* p_coef=new type0[m];
    type0* dp_coef=new type0[m];
    
    up=0.0;
    for(int i=0;i<m;i++)
        p_coef[i]=dp_coef[i]=0.0;
    
    p_coef[0]=1.0;
    
    for(int i=1;i<n+1;i++)
    {
        ii=static_cast<type0>(i);
        if(i%2==0)
        {
            m=i/2+1;
            p_coef[m-1]=(2.0-1.0/ii)*p_coef[m-2];
            for(int j=m-2;j>0;j--)
            {
                jj=static_cast<type0>(j);
                p_coef[j]*=-(2.0*jj+1)/ii;
                p_coef[j]+=(1.0+(2*jj-1.0)/ii)*p_coef[j-1];
            }
            p_coef[0]*=-1.0/ii;
        }
        else
        {
            m=(i+1)/2;
            for(int j=0;j<m-1;j++)
            {
                jj=static_cast<type0>(j);
                p_coef[j]*=1.0+2.0*jj/ii;
                p_coef[j]-=2.0*(jj+1.0)*p_coef[j+1]/ii;
            }
            //2.0-1.0/ii=(1+2*(m-1)/i)
            p_coef[m-1]*=2.0-1.0/ii;
        }
    }
    
    m=n/2+1;
    
    for(int i=1;i<m;i++)
    {
        dp_coef[i-1]=p_coef[i]*static_cast<type0>(i);
    }
    
    tol=std::numeric_limits<type0>::epsilon();
    ord=m;
    icurs=n-1;
    a=p_coef[m-1];
    
    for(int i=0;i<m;i++)
        p_coef[i]/=a;
    for(int i=1;i<m;i++)
        dp_coef[i-1]/=a;
    
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
    delete [] p_coef;
    delete [] dp_coef;
    
    
    if(n%2==0)
    {
        icurs++;
        jcurs=icurs-1;
        for(int i=icurs;i<n;i++)
        {
            tmp0=a;
            for(int j=icurs;j<n;j++)
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
        for(int i=icurs;i<n;i++)
            tmp0*=-x[i]*x[i];
        w[icurs-1]=2.0/(tmp0*tmp0);
        
        jcurs=icurs-2;
        for(int i=icurs;i<n;i++)
        {
            tmp0=a*x[i];
            for(int j=icurs;j<n;j++)
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
/*--------------------------------------------
 construct a legendre polynomial of degree n
 --------------------------------------------*/
void XMath::quadrature_hg(int n,type0* x,type0* w)
{
    int m=n/2+1;
    int iter,ord,icurs;
    int max_iter=50;
    type0 a,u0,inv_u0,f,up,df,tmp0,tol,tmp1;
    type0 ii,del_u0;

    
    type0* p_1=new type0[m];
    type0* p_2=new type0[m];
    type0* p_coef=new type0[m];
    type0* dp_coef=new type0[m];
    
    up=0.0;
    for(int i=0;i<m;i++)
        p_1[i]=p_2[i]=p_coef[i]=dp_coef[i]=0.0;
    
    p_1[0]=1.0;
    p_coef[0]=1.0;

    
    
    for(int i=1;i<n+1;i++)
    {
        ii=static_cast<type0>(i);
        m=i/2+1;
        if(i%2==0)
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
    
    m=n/2+1;

    
    for(int i=1;i<m;i++)
    {
        dp_coef[i-1]=p_coef[i]*static_cast<type0>(i);
    }
    
    tol=std::numeric_limits<type0>::epsilon();
    ord=m;
    icurs=n/2;
    if(n%2==1)
    {
        x[icurs]=0.0;
        icurs++;
    }
    
    a=p_coef[m-1];
    
    for(int i=0;i<m;i++)
        p_coef[i]/=a;
    for(int i=1;i<m;i++)
        dp_coef[i-1]/=a;
    
    

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
        for(int i=0;i<m;i++)
        {
            tmp1+=p_2[i]*tmp0;
            tmp0*=u0;
        }
        
        if(n%2==0)
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
    
    tmp0=sqrt(M_PI)/static_cast<type0>(n);
    for(int i=0;i<n-1;i++)
    {
        tmp0*=static_cast<type0>(2*(i+1));
    }
    if(n%2==1)
    {
        w[n/2]=tmp0/(p_2[0]*p_2[0]);
    }


    for(int i=0;i<n/2;i++)
    {
        icurs--;
        w[icurs]*=tmp0;
        w[i]=w[icurs];
        x[i]=-x[icurs];
    }
    
    delete [] p_1;
    delete [] p_2;
    delete [] p_coef;
    delete [] dp_coef;

}
/*--------------------------------------------
 calculates square root of 3x3 matrix
 ref: L. P. Franca
 An Algorithm to Compute The Square Root of
 a 3x3 Positive Definite Matrix
 Computers Math. Applic. Vol. 18, No. 5,
 pp. 459-466, 1989
 --------------------------------------------*/
bool XMath::Msqrt(type0(&A)[3][3],type0(&Asq)[3][3])
{
    type0 IA=0;
    for(int i=0;i<3;i++)
        IA+=A[i][i];
    type0 IIA=0;
    for(int i=0;i<3;i++)
        for(int j=0;j<3;j++)
            IIA-=A[i][j]*A[j][i];
    
    IIA+=IA*IA;
    IIA*=0.5;
    type0 IIIA=M3DET(A);
    type0 k=IA*IA-3*IIA;
    if(k<0.0) return false;
    
    if(k<TOLERANCE)
    {
        if(IA<=0.0)
            return false;
        M3ZERO(Asq);
        for(int i=0;i<3;i++)
            Asq[i][i]=sqrt(IA/3.0);
        return true;
    }
    
    type0 l=IA*(IA*IA -4.5*IIA)+13.5*IIIA;
    
    type0 temp=l/(k*sqrt(k));
    if(temp>1.0||temp<-1.0)
        return false;
    
    type0 phi=acos(temp);
    type0 lambda=sqrt((1.0/3.0)*(IA+2*sqrt(k)*cos(phi/3.0)));
    
    type0 IIIAsq=sqrt(IIIA);
    type0 y=-lambda*lambda+IA+2*IIIAsq/lambda;
    if(y<0.0) return false;
    type0 IAsq=lambda+sqrt(y);
    type0 IIAsq=0.5*(IAsq*IAsq-IA);
    
    type0 coef0=IAsq*IIAsq-IIIAsq;
    if(coef0==0)
        return false;
    coef0=1.0/coef0;
    
    M3ZERO(Asq);
    
    for(int i=0;i<3;i++)
        for(int j=0;j<3;j++)
            for(int k=0;k<3;k++)
                Asq[i][j]-=coef0*A[i][k]*A[k][j];
    
    type0 coef1=coef0*(IAsq*IAsq-IIAsq);
    for(int i=0;i<3;i++)
        for(int j=0;j<3;j++)
            Asq[i][j]+=coef1*A[i][j];
    
    type0 coef2=coef0*IAsq*IIIAsq;
    for(int i=0;i<3;i++)
        Asq[i][i]+=coef2;
    return true;
}
/*--------------------------------------------
 calculates square root of 3x3 matrix
 ref: L. P. Franca
 An Algorithm to Compute The Square Root of
 a 3x3 Positive Definite Matrix
 Computers Math. Applic. Vol. 18, No. 5,
 pp. 459-466, 1989
 --------------------------------------------*/
bool XMath::Msqrt(type0(&A)[2][2],type0(&Asq)[2][2])
{
    type0 tau=A[0][0]+A[1][1];
    type0 delta=A[0][0]*A[1][1]-A[0][1]*A[1][0];
    if(delta<0.0) return false;
    type0 s=sqrt(delta);
    if(tau+2.0*s<=0.0) return false;
    type0 t=sqrt(tau+2.0*s);
    Asq[0][0]=(A[0][0]+s)/t;
    Asq[0][1]=A[0][1]/t;
    Asq[1][0]=A[1][0]/t;
    Asq[1][1]=(A[1][1]+s)/t;
    
    
    return true;
}
