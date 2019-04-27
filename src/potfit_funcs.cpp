#include "potfit_funcs.h"
#include "memory.h"
#include "random.h"
#define NMAX_ALPHA_SHRNK_ATTMPS 10
#include <cmath>
using namespace MAPP_NS;
/*--------------------------------------------
 
 --------------------------------------------*/
type0 PotFitAux::find_max_alpha(const type0 xlo,const type0 xhi,bool dof,type0 max_dx,type0 x,type0 h)
{
    
    if(dof==false || h==0.0) return std::numeric_limits<type0>::infinity();
    type0 max_alpha=fabs(max_dx/h);
    
    if(h>0.0 && std::isnan(xhi)==false)
    {
        max_alpha=MIN(max_alpha,(xhi-x)/h);
        for(int it=0;it<NMAX_ALPHA_SHRNK_ATTMPS && x+max_alpha*h>xhi;it++)
            max_alpha=nextafter(max_alpha,0.0);
        if(x+max_alpha*h>xhi)  return 0.0;
    }
    else if(h<0.0 && std::isnan(xlo)==false)
    {
        max_alpha=MIN(max_alpha,(xlo-x)/h);
        for(int it=0;it<NMAX_ALPHA_SHRNK_ATTMPS && x+max_alpha*h<xlo;it++)
            max_alpha=nextafter(max_alpha,0.0);
        if(x+max_alpha*h<xlo) return 0.0;
        
    }
    
    return max_alpha;
}
/*--------------------------------------------
 
 --------------------------------------------*/
type0 PotFitAux::find_max_alpha(const type0 xlo,const type0 xhi,type0 max_dx,type0 x,type0 h)
{
    if(h==0.0) return std::numeric_limits<type0>::infinity();
    type0 max_alpha=fabs(max_dx/h);
    
    if(h>0.0 && std::isnan(xhi)==false)
    {
        max_alpha=MIN(max_alpha,(xhi-x)/h);
        for(int it=0;it<NMAX_ALPHA_SHRNK_ATTMPS && x+max_alpha*h>xhi;it++)
            max_alpha=nextafter(max_alpha,0.0);
        if(x+max_alpha*h>xhi)  return 0.0;
    }
    else if(h<0.0 && std::isnan(xlo)==false)
    {
        max_alpha=MIN(max_alpha,(xlo-x)/h);
        for(int it=0;it<NMAX_ALPHA_SHRNK_ATTMPS && x+max_alpha*h<xlo;it++)
            max_alpha=nextafter(max_alpha,0.0);
        if(x+max_alpha*h<xlo) return 0.0;
    }
    
    return max_alpha;
}
/*----------------------------------------------------------------------------------------
     _____    _   _   _____
    |  _  \  | | | | /  _  \
    | |_| |  | |_| | | | | |
    |  _  /  |  _  | | | | |
    | | \ \  | | | | | |_| |
    |_|  \_\ |_| |_| \_____/
 
 ----------------------------------------------------------------------------------------*/
PotFitPairFunc::PotFitPairFunc(const char* __name):
name(__name),nvars(0),rc(0.0),vars(NULL),dvars_lcl(NULL),dvars(NULL),
A_name(std::string("A_")+std::string(name)),
dA_name_max(std::string("dA_")+std::string(name)+std::string("_max")),
A_name_dof(std::string("A_")+std::string(name)+std::string("_dof"))
{
}
/*--------------------------------------------
 
 --------------------------------------------*/
PotFitPairFunc::~PotFitPairFunc()
{
}
/*--------------------------------------------
 
 --------------------------------------------*/
void PotFitPairFunc::random_neigh(Random* random)
{
    for(size_t i=0;i<nvars;i++)
        hvars[i]=(2.0*random->uniform()-1.0)*dvars_max[i];
    
}
/*--------------------------------------------
 
 --------------------------------------------*/
int PotFitPairFunc::set_A(PyObject* val)
{
    VarAPI<type0*> var(A_name.c_str());
    if(var.set(val)!=0) return -1;
    if(var.__var__.size!=nvars)
    {
        PyErr_Format(PyExc_TypeError,"size mismatch in %s",A_name.c_str());
        return -1;
    }
    memcpy(vars,var.val,nvars*sizeof(type0));
    return 0;
}
/*--------------------------------------------
 
 --------------------------------------------*/
PyObject* PotFitPairFunc::get_A()
{
    size_t* szp=&nvars;
    return var<type0*>::build(vars,&szp);
}
/*--------------------------------------------
 
 --------------------------------------------*/
PyObject* PotFitPairFunc::get_dA()
{
    size_t* szp=&nvars;
    return var<type0*>::build(dvars,&szp);
}
/*--------------------------------------------
 
 --------------------------------------------*/
int PotFitPairFunc::set_dA_max(PyObject* val)
{
    VarAPI<type0*> var(dA_name_max.c_str());
    if(var.set(val)!=0) return -1;
    if(var.__var__.size!=nvars)
    {
        PyErr_Format(PyExc_TypeError,"size mismatch in %s",dA_name_max.c_str());
        return -1;
    }
    memcpy(dvars_max,var.val,nvars*sizeof(type0));
    return 0;
}
/*--------------------------------------------
 
 --------------------------------------------*/
PyObject* PotFitPairFunc::get_dA_max()
{
    size_t* szp=&nvars;
    return var<type0*>::build(dvars_max,&szp);
}
/*--------------------------------------------
 
 --------------------------------------------*/
int PotFitPairFunc::set_A_dof(PyObject* val)
{
    VarAPI<bool*> var(A_name_dof.c_str());
    if(var.set(val)!=0) return -1;
    if(var.__var__.size!=nvars)
    {
        PyErr_Format(PyExc_TypeError,"size mismatch in %s",A_name_dof.c_str());
        return -1;
    }
    memcpy(dofs,var.val,nvars*sizeof(bool));
    return 0;
}
/*--------------------------------------------
 
 --------------------------------------------*/
PyObject* PotFitPairFunc::get_A_dof()
{
    size_t* szp=&nvars;
    return var<bool*>::build(dofs,&szp);
}
/*--------------------------------------------
 
 --------------------------------------------*/
PotFitRhoSpl::PotFitRhoSpl(const char* __name):
PotFitPairFunc(__name),
R(NULL)
{
}
/*--------------------------------------------
 
 --------------------------------------------*/
PotFitRhoSpl::~PotFitRhoSpl()
{
    Memory::dealloc(R);
}
/*--------------------------------------------
 
 --------------------------------------------*/
type0 PotFitRhoSpl::sort(type0* A,type0* R,size_t sz)
{
    if(sz==0) return 0.0;
    size_t* key=NULL;
    Memory::alloc(key,sz);
    for(size_t i=0;i<sz;i++) key[i]=i;
    
    XMath::quicksort(key,key+sz,
    [&R](size_t* rank_i,size_t* rank_j){return (R[*rank_i]>R[*rank_j]);},
    [&A,&R](size_t* rank_i,size_t* rank_j)
    {
        std::swap(R[*rank_i],R[*rank_j]);
        std::swap(A[*rank_i],A[*rank_j]);
    });
    
    Memory::dealloc(key);
    return R[0];
    
}
/*--------------------------------------------
 
 --------------------------------------------*/
type0 PotFitRhoSpl::F(type0 r)
{
    if(r>=rc) return 0.0;
    type0 ans=0.0;
    for(size_t i=0;i<nvars && r<R[i];i++)
        ans+=vars[i]*Algebra::pow<3>(R[i]-r);
    return ans;
}
/*--------------------------------------------
 
 --------------------------------------------*/
type0 PotFitRhoSpl::dF(type0 r)
{
    if(r>=rc) return 0.0;
    type0 ans=0.0;
    for(size_t i=0;i<nvars && r<R[i];i++)
        ans+=-3.0*vars[i]*Algebra::pow<2>(R[i]-r);
    return ans;
}
/*--------------------------------------------
 
 --------------------------------------------*/
type0 PotFitRhoSpl::ddF(type0 r)
{
    if(r>=rc) return 0.0;
    type0 ans=0.0;
    for(size_t i=0;i<nvars && r<R[i];i++)
    ans+=6.0*vars[i]*(R[i]-r);
    return ans;
}
/*--------------------------------------------
 
 --------------------------------------------*/
void PotFitRhoSpl::DF(type0 coef,type0 r,type0* dv)
{
    if(r>=rc) return;
    for(size_t i=0;i<nvars && r<R[i];i++)
        dv[i]+=coef*Algebra::pow<3>(R[i]-r);
}
/*--------------------------------------------
 
 --------------------------------------------*/
void PotFitRhoSpl::DdF(type0 coef,type0 r,type0* dv)
{
    if(r>=rc) return;
    for(size_t i=0;i<nvars && r<R[i];i++)
        dv[i]+=-coef*3.0*Algebra::pow<2>(R[i]-r);
}
/*--------------------------------------------
 
 --------------------------------------------*/
int PotFitRhoSpl::quadratic(type0* c,type0* ans)
{
    type0 delta_sq=c[0]*c[0]-4.0*c[1];
    if(delta_sq<0.0) return 0;
    if(delta_sq==0.0)
    {
        ans[0]=-0.5*c[0];
        return 1;
    }
    
    type0 delta=sqrt(delta_sq);
    ans[0]=0.5*(-c[0]+delta);
    ans[1]=0.5*(-c[0]-delta);
    return 2;
}
/*--------------------------------------------
 
 --------------------------------------------*/
int PotFitRhoSpl::cubic(type0* c,type0* ans)
{
    type0 p=c[1]-c[0]*c[0]/3.0;
    type0 q=2.0*c[0]*c[0]*c[0]/27.0-c[0]*c[1]/3.0+c[2];
    
    if(p==0.0)
    {
        if(q>=0.0) ans[0]=-pow(q,1.0/3.0)-c[0]/3.0;
        else ans[0]=pow(-q,1.0/3.0)-c[0]/3.0;
        return 1;
    }
    if(q==0.0)
    {
        int nroots=0;
        ans[nroots++]=-c[0]/3.0;
        if(p<0.0)
        {
            ans[nroots++]=sqrt(-p)-c[0]/3.0;
            ans[nroots++]=-sqrt(-p)-c[0]/3.0;
        }
        return nroots;
    }
    
    type0 t=sqrt(fabs(p)/3.0);
    type0 g=1.5*q/(p*t);
    if(p>0.0)
    {
        ans[0]=-2.0*t*sinh(asinh(g)/3.0)-c[0]/3.0;
        return 1;
    }
    
    
    if(4.0*p*p*p+27.0*q*q<0.0)
    {
        type0 theta=acos(g)/3.0;
        ans[0]=2.0*t*cos(theta)-c[0]/3.0;
        ans[1]=2.0*t*cos(theta+2.0*M_PI/3.0)-c[0]/3.0;
        ans[2]=2.0*t*cos(theta-2.0*M_PI/3.0)-c[0]/3.0;
        return 3;
    }
    
    if(q>0.0)
    {
        ans[0]=-2.0*t*cosh(acosh(-g)/3.0)-c[0]/3.0;
        return 1;
    }
    
    ans[0]=2.0*t*cosh(acosh(g)/3.0)-c[0]/3.0;
    return 1;
}
/*--------------------------------------------
 
 --------------------------------------------*/
int PotFitRhoSpl::quartic(type0* c,type0* ans)
{
    type0 p=c[1]-0.375*c[0]*c[0];
    type0 q=0.125*c[0]*(c[0]*c[0]-4.0*c[1])+c[2];
    type0 r=0.25*p*p+0.01171875*c[0]*c[0]*c[0]*c[0]-0.0625*c[0]*c[0]*c[1]+0.25*c[0]*c[2]-c[3];
    if(q==0)
    {
        
        if(r<0.0) return 0;
        int nroots=0;
        
        if(r-0.5*p>=0.0)
        {
            type0 absy=sqrt(r-0.5*p);
            ans[nroots++]=absy-0.25*c[0];
            ans[nroots++]=-absy-0.25*c[0];
        }
        if(r+0.5*p>=0.0)
        {
            type0 absy=sqrt(r+0.5*p);
            ans[nroots++]=absy-0.25*c[0];
            ans[nroots++]=-absy-0.25*c[0];
        }
        
        return nroots;
    }
    type0 __c[3]={p,r,-0.125*q*q};
    type0 __ans[3];
    cubic(__c,__ans);
    type0 m=__ans[0];
    
    if(m<0.0) return 0;
    type0 sqrt_2m=sqrt(2.0*m);
    int nroots=0;
    if(-m-p+q/sqrt_2m>=0.0)
    {
        type0 delta=sqrt(2.0*(-m-p+q/sqrt_2m));
        ans[nroots++]=0.5*(-sqrt_2m+delta)-0.25*c[0];
        ans[nroots++]=0.5*(-sqrt_2m-delta)-0.25*c[0];
    }
    
    if(-m-p-q/sqrt_2m>=0.0)
    {
        type0 delta=sqrt(2.0*(-m-p-q/sqrt_2m));
        ans[nroots++]=0.5*(sqrt_2m+delta)-0.25*c[0];
        ans[nroots++]=0.5*(sqrt_2m-delta)-0.25*c[0];
    }
    
    return nroots;
}
/*--------------------------------------------
 
 --------------------------------------------*/
void PotFitRhoSpl::slpine_max_alpha(
type0 (&c0)[4],type0 (&dc)[4],
const type0& x_lo,const type0& x_hi,type0& max_alpha)
{

    type0 __c[5]={
        c0[0]*dc[1]-c0[1]*dc[0],
        2.0*(c0[0]*dc[2]-c0[2]*dc[0]),
        3.0*(c0[0]*dc[3]-c0[3]*dc[0])+(c0[1]*dc[2]-c0[2]*dc[1]),
        2.0*(c0[1]*dc[3]-c0[3]*dc[1]),
        c0[2]*dc[3]-c0[3]*dc[2]};
    
    type0 ans[4];
    int n;
    if(__c[0]!=0.0)
    {
        __c[1]/=__c[0];
        __c[2]/=__c[0];
        __c[3]/=__c[0];
        __c[4]/=__c[0];
        n=quartic(__c+1,ans);
    }
    else if(__c[1]!=0.0)
    {
        __c[2]/=__c[1];
        __c[3]/=__c[1];
        __c[4]/=__c[1];
        n=cubic(__c+2,ans);
    }
    else if(__c[2]!=0.0)
    {
        __c[3]/=__c[2];
        __c[4]/=__c[2];
        n=quadratic(__c+3,ans);
    }
    else if(__c[3]!=0.0)
    {
        n=1;
        ans[0]=-__c[4]/__c[3];
    }
    else
    {
        if(c0[0]*dc[0]<0.0)
            max_alpha=MIN(max_alpha,-c0[0]/dc[0]);
        else
            max_alpha=0.0;
        return;
    }
    

    type0 f0;
    type0 df;
    for(int i=0;i<n;i++)
    {
        if(x_lo<ans[i] && ans[i]<x_hi &&(df=((dc[0]*ans[i]+dc[1])*ans[i]+dc[2])*ans[i]+dc[3])<0.0)
        {
            f0=((c0[0]*ans[i]+c0[1])*ans[i]+c0[2])*ans[i]+c0[3];
            max_alpha=MIN(max_alpha,-f0/df);
        }
    }
    return;
    
}
/*--------------------------------------------
 
 --------------------------------------------*/
void PotFitRhoSpl::slpine_max_alpha(type0* h,
type0 (&c0)[4],type0 (&dc)[4],
const type0& x_lo,const type0& x_hi,type0& max_alpha)
{
    
    // when we come here alpha is not infinity
    int n=0;
    type0 __max_alpha;
    type0 c[4];
    type0 __c[2];
    type0 ans[2];
    type0 df;
    bool trial=true;
    for(int it=0;it<NMAX_ALPHA_SHRNK_ATTMPS && trial;it++)
    {
        Algebra::V_eq<4>(c0,c);
        Algebra::V_add_x_mul_V<4>(max_alpha,dc,c);
        
        if(c[0]!=0.0)
        {
            __c[0]=2.0*c[1]/(3.0*c[0]);
            __c[1]=c[2]/(3.0*c[0]);
            n=quadratic(__c,ans);
        }
        else if(c[1]!=0.0)
        {
            ans[0]=-c[2]/c[1];
            n=1;
        }
        else
        {
            trial=false;
            continue;
        }
        trial=false;
        for(int i=0;i<n;i++)
        {
            df=F(max_alpha,h,ans[i]);
            if((x_lo<=ans[i] && ans[i]<x_hi) && df<0.0)
            {
                trial=true;
                __max_alpha=-F(ans[i])/df;
                if(__max_alpha>=max_alpha)
                    max_alpha=nextafter(max_alpha,0.0);
                else
                    max_alpha=__max_alpha;
            }
        }
    }
    
    if(trial) max_alpha=0.0;
    
}
/*--------------------------------------------
 
 --------------------------------------------*/
type0 PotFitRhoSpl::F(type0* h,type0 r)
{
    if(r>=rc) return 0.0;
    type0 ans=0.0;
    for(size_t i=0;i<nvars && r<R[i];i++)
        ans+=h[i]*Algebra::pow<3>(R[i]-r);
    return ans;
}
/*--------------------------------------------
 
 --------------------------------------------*/
type0 PotFitRhoSpl::F(type0 alpha,type0* h,type0 r)
{
    if(r>=rc) return 0.0;
    type0 ans=0.0;
    for(size_t i=0;i<nvars && r<R[i];i++)
        ans+=(vars[i]+alpha*h[i])*Algebra::pow<3>(R[i]-r);
    return ans;
}
/*--------------------------------------------
 
 --------------------------------------------*/
void PotFitRhoSpl::find_max_alpha(type0& max_alpha)
{
    
    type0 c0[4]={0.0,0.0,0.0,0.0};
    type0 dc[4]={0.0,0.0,0.0,0.0};
    type0 x_hi,x_lo,f_lo,df_lo;
    bool dof=false;
    for(size_t i=0;i<nvars;i++)
    {
        if(max_alpha==0.0) continue;
        
        c0[0]+=-vars[i];
        c0[1]+=vars[i]*3.0*R[i];
        c0[2]+=-vars[i]*3.0*R[i]*R[i];
        c0[3]+=vars[i]*R[i]*R[i]*R[i];
        if(dofs[i]) dof=true;
        if(!dof) continue;
        dc[0]+=-hvars[i];
        dc[1]+=hvars[i]*3.0*R[i];
        dc[2]+=-hvars[i]*3.0*R[i]*R[i];
        dc[3]+=hvars[i]*R[i]*R[i]*R[i];
        x_hi=R[i];
        if(i!=nvars-1) x_lo=R[i+1];
        else x_lo=0.0;
        
        df_lo=F(hvars,x_lo);
        if(df_lo!=0.0)
            max_alpha=MIN(max_alpha,dvars_max[i]/fabs(df_lo));
        
        f_lo=F(x_lo);
        if((df_lo<0.0 && std::isinf(max_alpha)) || f_lo+max_alpha*df_lo<0.0)
        {
            max_alpha=MIN(max_alpha,-f_lo/df_lo);
            
            for(int it=0;it<NMAX_ALPHA_SHRNK_ATTMPS&&F(max_alpha,hvars,x_lo)<0.0;it++)
                max_alpha=nextafter(max_alpha,0.0);
            if(F(max_alpha,hvars,x_lo)<0.0) max_alpha=0.0;
        }
        
        if(std::isinf(max_alpha))
            slpine_max_alpha(c0,dc,x_lo,x_hi,max_alpha);
        
        slpine_max_alpha(hvars,c0,dc,x_lo,x_hi,max_alpha);
    }
}
/*--------------------------------------------
 
 --------------------------------------------*/
void PotFitRhoSpl::random_neigh(Random* random)
{
    for(size_t i=0;i<nvars;i++)
        hvars[i]=(2.0*random->uniform()-1.0)*dvars_max[i];
    /*
    type0 x_hi,x_lo;
    for(size_t i=0;i<nvars;i++)
    {
        x_hi=R[i];
        if(i!=nvars-1) x_lo=R[i+1];
        else x_lo=0.0;
        hvars[i]/=Algebra::pow<3>(x_hi-x_lo);
        
        for(size_t j=i+1;j<nvars;j++)
        {
            if(j!=nvars-1) x_lo=R[j+1];
            else x_lo=0.0;
            hvars[j]-=hvars[i]*Algebra::pow<3>(x_hi-x_lo);
        }
    }*/
}
/*--------------------------------------------
 
 --------------------------------------------*/
bool PotFitRhoSpl::validate()
{
    type0 c0[4]={0.0,0.0,0.0,0.0};
    type0 c1[3]={0.0,0.0,0.0};
    type0 ans[3];
    bool dof=false;
    type0 x_hi,x_lo;
    int nroots;
    
    for(size_t i=0;i<nvars;i++)
    {
        c0[0]+=-vars[i];
        c0[1]+=vars[i]*3.0*R[i];
        c0[2]+=-vars[i]*3.0*R[i]*R[i];
        c0[3]+=vars[i]*R[i]*R[i]*R[i];
        if(dofs[i]) dof=true;
        if(!dof) continue;
        x_hi=R[i];
        if(i!=nvars-1) x_lo=R[i+1];
        else x_lo=0.0;
        nroots=0;
        if(c0[0]!=0.0)
        {
            c1[0]=c0[1]/c0[0];
            c1[1]=c0[2]/c0[0];
            c1[2]=c0[3]/c0[0];
            nroots=cubic(c1,ans);
        }
        else if(c0[1]!=0.0)
        {
            c1[0]=c0[2]/c0[1];
            c1[1]=c0[3]/c0[1];
            nroots=quadratic(c1,ans);
        }
        else if(c0[2]!=0.0)
        {
            nroots=1;
            ans[0]=-c0[3]/c0[2];
        }
        else if(c0[3]<0.0) return false;
        
        for(int j=0;j<nroots;j++)
            if(x_lo<=ans[j] && ans[j]<x_hi) return false;
        
    }
    return true;
}
/*--------------------------------------------
 
 --------------------------------------------*/
int PotFitRhoSpl::set_init(PyObject* val,type0*& data,size_t& data_sz)
{
    VarAPI<type0* [2]> var(name);
    var.val[0]=NULL;
    var.val[1]=NULL;
    if(var.set(val)!=0) return -1;
    if(var.__var__[0].size!=var.__var__[1].size)
    {
        PyErr_Format(PyExc_TypeError,"size mismatch in %s",name);
        return -1;
    }
    nvars=var.__var__[0].size;
    Memory::dealloc(R);
    R=NULL;
    Memory::alloc(R,nvars);
    Memory::grow(data,data_sz,data_sz+nvars);
    memcpy(data+data_sz,var.val[0],nvars*sizeof(type0));
    memcpy(R,var.val[1],nvars*sizeof(type0));
    rc=PotFitRhoSpl::sort(data+data_sz,R,nvars);
    data_sz+=nvars;
    return 0;
}
/*--------------------------------------------
 
 --------------------------------------------*/
int PotFitRhoSpl::set(PyObject* val)
{
    VarAPI<type0*> var(name);
    if(var.set(val)!=0) return -1;
    if(var.__var__.size!=nvars)
    {
        PyErr_Format(PyExc_TypeError,"size mismatch in %s",name);
        return -1;
    }
    memcpy(vars,var.val,nvars*sizeof(type0));
    return 0;
}
/*--------------------------------------------
 
 --------------------------------------------*/
PotFitRhoAng::PotFitRhoAng(const char* __name):
PotFitPairFunc(__name)
{
    rc=2.8;
}
/*--------------------------------------------
 
 --------------------------------------------*/
PotFitRhoAng::~PotFitRhoAng()
{
}
/*--------------------------------------------
 
 --------------------------------------------*/
type0 PotFitRhoAng::F(type0 r)
{
    if(r>=rc) return 0.0;
    type0 y=exp(-vars[1]*r);
    return vars[0]*exp(1.0/(r-rc))*Algebra::pow<6>(r)*y*(1.0+512.0*y);
}
/*--------------------------------------------
 
 --------------------------------------------*/
type0 PotFitRhoAng::dF(type0 r)
{
    if(r>=rc) return 0.0;
    type0 y=exp(-vars[1]*r);
    type0 rho=vars[0]*exp(1.0/(r-rc))*Algebra::pow<6>(r)*y*(1.0+512.0*y);
    return rho*(6.0/r-1.0/((r-rc)*(r-rc))-vars[1]*(1.0+y/(0.001953125+y)));
}
/*--------------------------------------------
 
 --------------------------------------------*/
type0 PotFitRhoAng::ddF(type0 r)
{
    if(r>=rc) return 0.0;
    type0 y=exp(-vars[1]*r);
    type0 z=(r-rc)*(r-rc);
    type0 x0=1.0+512.0*y;
    type0 x1=1.0+1024.0*y;
    type0 x2=1.0+2048.0*y;
    
    
    return
    vars[0]*r*r*r*r*y*exp(1.0/(r-rc))*
    ((30.0+r*(r+r*r-11.0*z-rc*rc)/(z*z))*x0+
     2.0*vars[1]*r*(r/z-6.0)*x1+
    vars[1]*vars[1]*r*r*x2);
}
/*--------------------------------------------
 
 --------------------------------------------*/
void PotFitRhoAng::DF(type0 coef,type0 r,type0* dv)
{
    if(r>=rc) return;
    type0 y=exp(-vars[1]*r);
    type0 c0=exp(1.0/(r-rc))*Algebra::pow<6>(r);
    dv[0]+=coef*c0*y*(1.0+512.0*y);
    dv[1]+=-coef*vars[0]*r*c0*y*(1.0+1024.0*y);
}
/*--------------------------------------------
 
 --------------------------------------------*/
void PotFitRhoAng::DdF(type0 coef,type0 r,type0* dv)
{
    if(r>=rc) return;
    type0 y=exp(-vars[1]*r);
    dv[0]+=coef*(6.0/r-1.0/((r-rc)*(r-rc))-vars[1]*(1.0+y/(0.001953125+y)))*exp(1.0/(r-rc))*Algebra::pow<6>(r)*y*(1.0+512.0*y);
    type0 c0=exp(1.0/(r-rc))*Algebra::pow<7>(r)*y*(1.0+1024.0*y);
    dv[1]+=coef*vars[0]*c0*(vars[1]*(1.0+y/(0.0009765625+y))-7.0/r+1.0/Algebra::pow<2>(r-rc));
}
/*--------------------------------------------
 
 --------------------------------------------*/
void PotFitRhoAng::find_max_alpha(type0& max_alpha)
{
    max_alpha=MIN(max_alpha,PotFitAux::find_max_alpha(0.0,std::numeric_limits<type0>::quiet_NaN(),dvars_max[0],vars[0],hvars[0]));
    max_alpha=MIN(max_alpha,PotFitAux::find_max_alpha(0.0,std::numeric_limits<type0>::quiet_NaN(),dvars_max[1],vars[1],hvars[1]));

}
/*--------------------------------------------
 
 --------------------------------------------*/
bool PotFitRhoAng::validate()
{
    if(vars[0]<0.0) return false;
    if(vars[1]<0.0) return false;
    return true;
}
/*--------------------------------------------
 
 --------------------------------------------*/
int PotFitRhoAng::set_init(PyObject* val,type0*& data,size_t& data_sz)
{
    
    VarAPI<type0 [2]> var(name);
    if(var.set(val)!=0) return -1;
    nvars=2;
    Memory::grow(data,data_sz,data_sz+nvars);
    memcpy(data+data_sz,var.val,nvars*sizeof(type0));
    data_sz+=nvars;
    return 0;
}
/*--------------------------------------------
 
 --------------------------------------------*/
int PotFitRhoAng::set(PyObject* val)
{
    VarAPI<type0[2]> var(name);
    if(var.set(val)!=0) return -1;
    memcpy(vars,var.val,nvars*sizeof(type0));
    return 0;
}
/*----------------------------------------------------------------------------------------
     _____   _   _   _
    |  _  \ | | | | | |
    | |_| | | |_| | | |
    |  ___/ |  _  | | |
    | |     | | | | | |
    |_|     |_| |_| |_|
 
 ----------------------------------------------------------------------------------------*/
PotFitPhiSpl::PotFitPhiSpl(const char* __name):
PotFitPairFunc(__name),
R(NULL)
{
}
/*--------------------------------------------
 
 --------------------------------------------*/
PotFitPhiSpl::~PotFitPhiSpl()
{
    Memory::dealloc(R);
}
/*--------------------------------------------
 
 --------------------------------------------*/
type0 PotFitPhiSpl::F(type0 r)
{
    if(r>=rc) return 0.0;
    type0 ans=0.0;
    for(size_t i=0;i<nvars && r<R[i];i++)
        ans+=vars[i]*Algebra::pow<3>(R[i]-r);
    return ans;
}
/*--------------------------------------------
 
 --------------------------------------------*/
type0 PotFitPhiSpl::dF(type0 r)
{
    if(r>=rc) return 0.0;
    type0 ans=0.0;
    for(size_t i=0;i<nvars && r<R[i];i++)
        ans+=-3.0*vars[i]*Algebra::pow<2>(R[i]-r);
    return ans;
}
/*--------------------------------------------
 
 --------------------------------------------*/
type0 PotFitPhiSpl::ddF(type0 r)
{
    if(r>=rc) return 0.0;
    type0 ans=0.0;
    for(size_t i=0;i<nvars && r<R[i];i++)
    ans+=6.0*vars[i]*(R[i]-r);
    return ans;
}
/*--------------------------------------------
 
 --------------------------------------------*/
void PotFitPhiSpl::DF(type0 coef,type0 r,type0* dv)
{
    if(r>=rc) return;
    for(size_t i=0;i<nvars && r<R[i];i++)
        dv[i]+=coef*Algebra::pow<3>(R[i]-r);
}
/*--------------------------------------------
 
 --------------------------------------------*/
void PotFitPhiSpl::DdF(type0 coef,type0 r,type0* dv)
{
    if(r>=rc) return;
    for(size_t i=0;i<nvars && r<R[i];i++)
        dv[i]+=-coef*3.0*Algebra::pow<2>(R[i]-r);
}
/*--------------------------------------------
  
  --------------------------------------------*/
type0 PotFitPhiSpl::F(type0* h,type0 r)
{
    if(r>=rc) return 0.0;
    type0 ans=0.0;
    for(size_t i=0;i<nvars && r<R[i];i++)
        ans+=h[i]*Algebra::pow<3>(R[i]-r);
    return ans;
}
/*--------------------------------------------
 
 --------------------------------------------*/
void PotFitPhiSpl::find_max_alpha(type0& max_alpha)
{
    type0 x_lo,df_lo;
    bool dof=false;
    for(size_t i=0;i<nvars;i++)
    {
        if(max_alpha==0.0) continue;
        if(hvars[i]!=0.0) dof=true;
        if(!dof) continue;

        if(i!=nvars-1) x_lo=R[i+1];
        else x_lo=0.0;
        
        df_lo=F(hvars,x_lo);
        if(df_lo!=0.0)
            max_alpha=MIN(max_alpha,dvars_max[i]/fabs(df_lo));
    }
}
/*--------------------------------------------
 
 --------------------------------------------*/
void PotFitPhiSpl::random_neigh(Random* random)
{
    for(size_t i=0;i<nvars;i++)
        hvars[i]=(2.0*random->uniform()-1.0)*dvars_max[i];
    /*
    type0 x_hi,x_lo;
    for(size_t i=0;i<nvars;i++)
    {
        x_hi=R[i];
        if(i!=nvars-1) x_lo=R[i+1];
        else x_lo=0.0;
        hvars[i]/=Algebra::pow<3>(x_hi-x_lo);
        
        for(size_t j=i+1;j<nvars;j++)
        {
            if(j!=nvars-1) x_lo=R[j+1];
            else x_lo=0.0;
            hvars[j]-=hvars[i]*Algebra::pow<3>(x_hi-x_lo);
        }
    }*/
}
/*--------------------------------------------
 
 --------------------------------------------*/
bool PotFitPhiSpl::validate()
{
    return true;
}
/*--------------------------------------------
 
 --------------------------------------------*/
int PotFitPhiSpl::set_init(PyObject* val,type0*& data,size_t& data_sz)
{
    VarAPI<type0* [2]> var(name);
    var.val[0]=NULL;
    var.val[1]=NULL;
    
    if(var.set(val)!=0) return -1;
    if(var.__var__[0].size!=var.__var__[1].size)
    {
        PyErr_Format(PyExc_TypeError,"size mismatch in %s",name);
        return -1;
    }
    nvars=var.__var__[0].size;
    Memory::dealloc(R);
    R=NULL;
    Memory::alloc(R,nvars);
    Memory::grow(data,data_sz,data_sz+nvars);
    memcpy(data+data_sz,var.val[0],nvars*sizeof(type0));
    memcpy(R,var.val[1],nvars*sizeof(type0));
    rc=PotFitRhoSpl::sort(data+data_sz,R,nvars);
    data_sz+=nvars;
    return 0;
}
/*--------------------------------------------
 
 --------------------------------------------*/
int PotFitPhiSpl::set(PyObject* val)
{
    VarAPI<type0*> var(name);
    if(var.set(val)!=0) return -1;
    if(var.__var__.size!=nvars)
    {
        PyErr_Format(PyExc_TypeError,"size mismatch in %s",name);
        return -1;
    }
    memcpy(vars,var.val,nvars*sizeof(type0));
    return 0;
}
/*----------------------------------------------------------------------------------------
 _____       ___  ___   _____   _____   _____   _____   _____   _____
| ____|     /   |/   | |  _  \ | ____| |  _  \ |  _  \ | ____| |  _  \
| |__      / /|   /| | | |_| | | |__   | | | | | | | | | |__   | | | |
|  __|    / / |__/ | | |  _  { |  __|  | | | | | | | | |  __|  | | | |
| |___   / /       | | | |_| | | |___  | |_| | | |_| | | |___  | |_| |
|_____| /_/        |_| |_____/ |_____| |_____/ |_____/ |_____| |_____/
 
 ----------------------------------------------------------------------------------------*/
PotFitEmbFunc::PotFitEmbFunc(const char* __name):
name(__name),nvars(0),vars(NULL),dvars_lcl(NULL),dvars(NULL),
A_name(std::string("A_")+std::string(name)),
dA_name_max(std::string("dA_")+std::string(name)+std::string("_max")),
A_name_dof(std::string("A_")+std::string(name)+std::string("_dof"))
{
}
/*--------------------------------------------
 
 --------------------------------------------*/
PotFitEmbFunc::~PotFitEmbFunc()
{
}
/*--------------------------------------------
 
 --------------------------------------------*/
void PotFitEmbFunc::random_neigh(Random* random)
{
    for(size_t i=0;i<nvars;i++)
        hvars[i]=(2.0*random->uniform()-1.0)*dvars_max[i];
}
/*--------------------------------------------
 
 --------------------------------------------*/
int PotFitEmbFunc::set_A(PyObject* val)
{
    VarAPI<type0*> var(A_name.c_str());
    if(var.set(val)!=0) return -1;
    if(var.__var__.size!=nvars)
    {
        PyErr_Format(PyExc_TypeError,"size mismatch in %s",A_name.c_str());
        return -1;
    }
    memcpy(vars,var.val,nvars*sizeof(type0));
    return 0;
}
/*--------------------------------------------
 
 --------------------------------------------*/
PyObject* PotFitEmbFunc::get_A()
{
    size_t* szp=&nvars;
    return var<type0*>::build(vars,&szp);
}
/*--------------------------------------------
 
 --------------------------------------------*/
PyObject* PotFitEmbFunc::get_dA()
{
    size_t* szp=&nvars;
    return var<type0*>::build(dvars,&szp);
}
/*--------------------------------------------
 
 --------------------------------------------*/
int PotFitEmbFunc::set_dA_max(PyObject* val)
{
    VarAPI<type0*> var(dA_name_max.c_str());
    if(var.set(val)!=0) return -1;
    if(var.__var__.size!=nvars)
    {
        PyErr_Format(PyExc_TypeError,"size mismatch in %s",dA_name_max.c_str());
        return -1;
    }
    memcpy(dvars_max,var.val,nvars*sizeof(type0));
    return 0;
}
/*--------------------------------------------
 
 --------------------------------------------*/
PyObject* PotFitEmbFunc::get_dA_max()
{
    size_t* szp=&nvars;
    return var<type0*>::build(dvars_max,&szp);
}
/*--------------------------------------------
 
 --------------------------------------------*/
int PotFitEmbFunc::set_A_dof(PyObject* val)
{
    VarAPI<bool*> var(A_name_dof.c_str());
    if(var.set(val)!=0) return -1;
    if(var.__var__.size!=nvars)
    {
        PyErr_Format(PyExc_TypeError,"size mismatch in %s",A_name_dof.c_str());
        return -1;
    }
    memcpy(dofs,var.val,nvars*sizeof(bool));
    return 0;
}
/*--------------------------------------------
 
 --------------------------------------------*/
PyObject* PotFitEmbFunc::get_A_dof()
{
    size_t* szp=&nvars;
    return var<bool*>::build(dofs,&szp);
}
/*--------------------------------------------
 
 --------------------------------------------*/
PotFitEmbAck::PotFitEmbAck(const char* __name):
PotFitEmbFunc(__name)
{
}
/*--------------------------------------------
 
 --------------------------------------------*/
PotFitEmbAck::~PotFitEmbAck()
{
}
/*--------------------------------------------
 
 --------------------------------------------*/
type0 PotFitEmbAck::F(type0 rho)
{
    return vars[0]*sqrt(rho)+(vars[1]+vars[2]*rho*rho)*rho*rho+vars[3]*rho;
}
/*--------------------------------------------
 
 --------------------------------------------*/
type0 PotFitEmbAck::dF(type0 rho)
{
    if(rho==0.0) return 0.0;
    return 0.5*vars[0]/sqrt(rho)+2.0*rho*(vars[1]+2.0*vars[2]*rho*rho)+vars[3];
}
/*--------------------------------------------
 
 --------------------------------------------*/
type0 PotFitEmbAck::ddF(type0 rho)
{
    if(rho==0.0) return 0.0;
    
    return -0.25*vars[0]/(rho*sqrt(rho))+2.0*vars[1]+12.0*vars[2]*rho*rho;
}
/*--------------------------------------------
 
 --------------------------------------------*/
void PotFitEmbAck::DF(type0 coef,type0 rho,type0* dv)
{
    if(rho==0.0) return;
    dv[0]+=coef*sqrt(rho);
    dv[1]+=coef*rho*rho;
    dv[2]+=coef*rho*rho*rho*rho;
    dv[3]+=coef*rho;
}
/*--------------------------------------------
 
 --------------------------------------------*/
void PotFitEmbAck::DdF(type0 coef,type0 rho,type0* dv)
{
    if(rho==0.0) return;
    dv[0]+=coef*0.5/sqrt(rho);
    dv[1]+=coef*2.0*rho;
    dv[2]+=coef*4.0*rho*rho*rho;
    dv[3]+=coef;
}
/*--------------------------------------------
 
 --------------------------------------------*/
void PotFitEmbAck::find_max_alpha(type0& max_alpha)
{
    max_alpha=MIN(max_alpha,PotFitAux::find_max_alpha(std::numeric_limits<type0>::quiet_NaN(),0.0,dvars_max[0],vars[0],hvars[0]));
    max_alpha=MIN(max_alpha,PotFitAux::find_max_alpha(std::numeric_limits<type0>::quiet_NaN(),std::numeric_limits<type0>::quiet_NaN(),dvars_max[1],vars[1],hvars[1]));
    max_alpha=MIN(max_alpha,PotFitAux::find_max_alpha(std::numeric_limits<type0>::quiet_NaN(),std::numeric_limits<type0>::quiet_NaN(),dvars_max[2],vars[2],hvars[2]));
    max_alpha=MIN(max_alpha,PotFitAux::find_max_alpha(std::numeric_limits<type0>::quiet_NaN(),std::numeric_limits<type0>::quiet_NaN(),dvars_max[3],vars[3],hvars[3]));
}
/*--------------------------------------------
 
 --------------------------------------------*/
bool PotFitEmbAck::validate()
{
    if(vars[0]>0.0) return false;
    if(vars[2]>0.0 && 0.875*pow(-2.0*vars[0],4.0/7.0)*pow(vars[2],3.0/7.0)+vars[1]<0.0) return false;
    return true;
}
/*--------------------------------------------
 
 --------------------------------------------*/
int PotFitEmbAck::set_init(PyObject* val,type0*& data,size_t& data_sz)
{
    
    VarAPI<type0 [4]> var("F");
    if(var.set(val)!=0) return -1;
    nvars=4;
    Memory::grow(data,data_sz,data_sz+nvars);
    memcpy(data+data_sz,var.val,nvars*sizeof(type0));
    data_sz+=nvars;
    return 0;
}
/*--------------------------------------------
 
 --------------------------------------------*/
int PotFitEmbAck::set(PyObject* val)
{
    VarAPI<type0[4]> var(name);
    if(var.set(val)!=0) return -1;
    memcpy(vars,var.val,nvars*sizeof(type0));
    return 0;
}
/*--------------------------------------------
 
 --------------------------------------------*/
PotFitEmbAng::PotFitEmbAng(const char* __name):
PotFitEmbFunc(__name)
{
    type0 t=3.0*sqrt(5.0)*cos((M_PI-acos(0.4*sqrt(5.0)))/3.0)-2.5;
    lim=pow(t,2.0/3.0)*((t+5.0)*t-5.0)/(9.0*t*Algebra::pow<3>(1.0+t));
}
/*--------------------------------------------
 
 --------------------------------------------*/
PotFitEmbAng::~PotFitEmbAng()
{
}
/*--------------------------------------------
 
 --------------------------------------------*/
type0 PotFitEmbAng::F(type0 rho)
{
    type0 rho_b=rho/vars[3];
    type0 rho_b_2_3=pow(rho_b,2.0/3.0);
    //return (vars[0]*rho_b+vars[1]+vars[2]*rho_b_2_3/(1.0+rho_b))*rho_b;
    return (vars[0]*rho_b+vars[1]+vars[2]*lim*rho_b+vars[2]*rho_b_2_3/(1.0+rho_b))*rho_b;
}
/*--------------------------------------------
 
 --------------------------------------------*/
type0 PotFitEmbAng::dF(type0 rho)
{
    if(rho==0.0) return 0.0;
    type0 rho_b=rho/vars[3];
    type0 rho_b_2_3=pow(rho_b,2.0/3.0);
    //return (2.0*vars[0]*rho_b+vars[1]+vars[2]*rho_b_2_3*(5.0+2.0*rho_b)/(3.0*Algebra::pow<2>(1.0+rho_b)))/vars[3];
    //printf("<<<%lf\n",-rho_b*(2.0*vars[0]*rho_b+vars[1]+2.0*vars[2]*lim*rho_b+vars[2]*rho_b_2_3*(5.0+2.0*rho_b)/(3.0*Algebra::pow<2>(1.0+rho_b)))/vars[3]);
    return (2.0*vars[0]*rho_b+vars[1]+2.0*vars[2]*lim*rho_b+vars[2]*rho_b_2_3*(5.0+2.0*rho_b)/(3.0*Algebra::pow<2>(1.0+rho_b)))/vars[3];
}
/*--------------------------------------------
 
 --------------------------------------------*/
type0 PotFitEmbAng::ddF(type0 rho)
{
    if(rho==0.0) return 0.0;
    type0 rho_b=rho/vars[3];
    type0 rho_b_2_3=pow(rho_b,2.0/3.0);
    //return 2.0*(vars[0]-vars[2]*rho_b_2_3*((rho_b+5.0)*rho_b-5.0)/(9.0*rho_b*Algebra::pow<3>(1.0+rho_b)))/(vars[3]*vars[3]);
    return 2.0*(vars[0]-vars[2]*lim-vars[2]*rho_b_2_3*((rho_b+5.0)*rho_b-5.0)/(9.0*rho_b*Algebra::pow<3>(1.0+rho_b)))/(vars[3]*vars[3]);
}
/*--------------------------------------------
 
 --------------------------------------------*/
void PotFitEmbAng::DF(type0 coef,type0 rho,type0* dv)
{
    if(rho==0.0) return;
    type0 rho_b=rho/vars[3];
    type0 rho_b_2_3=pow(rho_b,2.0/3.0);
    dv[0]+=coef*rho_b*rho_b;
    dv[1]+=coef*rho_b;
    //dv[2]+=rho_b*rho_b_2_3/(1.0+rho_b);
    dv[2]+=coef*(lim*rho_b*rho_b+rho_b*rho_b_2_3/(1.0+rho_b));
    //dv[3]+=coef*(-rho_b*(2.0*vars[0]*rho_b+vars[1]+vars[2]*rho_b_2_3*(5.0+2.0*rho_b)/(3.0*Algebra::pow<2>(1.0+rho_b)))/vars[3]);
    dv[3]+=-rho_b*coef*(2.0*vars[0]*rho_b+vars[1]+2.0*vars[2]*lim*rho_b+vars[2]*rho_b_2_3*(5.0+2.0*rho_b)/(3.0*Algebra::pow<2>(1.0+rho_b)))/vars[3];
    //printf("=====%lf\n",(-rho_b*(2.0*vars[0]*rho_b+vars[1]+vars[2]*rho_b_2_3*(5.0+2.0*rho_b)/(3.0*Algebra::pow<2>(1.0+rho_b)))/vars[3]));
}
/*--------------------------------------------
 
 --------------------------------------------*/
void PotFitEmbAng::DdF(type0 coef,type0 rho,type0* dv)
{
    if(rho==0.0) return;
    type0 rho_b=rho/vars[3];
    type0 rho_b_2_3=pow(rho_b,2.0/3.0);
    dv[0]+=coef*2.0*rho_b/vars[3];
    dv[1]+=coef*1.0/vars[3];
    //dv[2]+=(rho_b_2_3*(5.0+2.0*rho_b)/(3.0*(1.0+rho_b)*(1.0+rho_b)))/vars[3];
    dv[2]+=coef*(2.0*lim*rho_b+rho_b_2_3*(5.0+2.0*rho_b)/(3.0*(1.0+rho_b)*(1.0+rho_b)))/vars[3];
    dv[3]+=coef*
    (vars[0]+
    vars[1]*4.0*rho_b+
    vars[2]*rho_b_2_3*((4.0*rho_b+11.0)*rho_b+25.0)/(9.0*Algebra::pow<3>((1.0+rho_b))))/(vars[3]*vars[3]);
}
/*--------------------------------------------
 
 --------------------------------------------*/
void PotFitEmbAng::find_max_alpha(type0& max_alpha)
{
    max_alpha=MIN(max_alpha,PotFitAux::find_max_alpha(0.0,std::numeric_limits<type0>::quiet_NaN(),dvars_max[2],vars[2],hvars[2]));
    max_alpha=MIN(max_alpha,PotFitAux::find_max_alpha(0.0,std::numeric_limits<type0>::quiet_NaN(),dvars_max[3],vars[3],hvars[3]));
    max_alpha=MIN(max_alpha,PotFitAux::find_max_alpha(std::numeric_limits<type0>::quiet_NaN(),std::numeric_limits<type0>::quiet_NaN(),dvars_max[1],vars[1],hvars[1]));
    max_alpha=MIN(max_alpha,PotFitAux::find_max_alpha(0.0,std::numeric_limits<type0>::quiet_NaN(),dvars_max[0],vars[0],hvars[0]));
}
/*--------------------------------------------
 
 --------------------------------------------*/
bool PotFitEmbAng::validate()
{
    if(vars[3]<0.0) return false;
    if(vars[0]<0.0) return false;
    if(vars[2]<0.0) return false;
    return true;
}
/*--------------------------------------------
 
 --------------------------------------------*/
int PotFitEmbAng::set_init(PyObject* val,type0*& data,size_t& data_sz)
{
    
    VarAPI<type0 [4]> var(name);
    if(var.set(val)!=0) return -1;
    nvars=4;
    Memory::grow(data,data_sz,data_sz+nvars);
    memcpy(data+data_sz,var.val,nvars*sizeof(type0));
    data_sz+=nvars;
    return 0;
}
/*--------------------------------------------
 
 --------------------------------------------*/
int PotFitEmbAng::set(PyObject* val)
{
    VarAPI<type0[4]> var(name);
    if(var.set(val)!=0) return -1;
    memcpy(vars,var.val,nvars*sizeof(type0));
    return 0;
}


