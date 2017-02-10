#include <stdlib.h>
#include "ff_eam_dmd.h"
#include "neighbor_dmd.h"
#include "elements.h"
#include "memory.h"
#include "xmath.h"
#include "atoms_dmd.h"
#include "MAPP.h"
#include "dynamic_dmd.h"
#include "read_eam.h"
#include <limits>
#define PI_IN_SQ 0.564189583547756286948079451561
using namespace MAPP_NS;
/*--------------------------------------------
 constructor
 --------------------------------------------*/
ForceFieldEAMDMD::
ForceFieldEAMDMD(AtomsDMD*& atoms,
type0 __dr,type0 __drho,size_t __nr,size_t __nrho,
type0(***&& __r_phi_arr)[4],type0(***&& __rho_arr)[4],type0(**&& __F_arr)[5],
type0**&& __cut):
ForceFieldDMD(atoms),
dr(__dr),
drho(__drho),
nr(__nr),
nrho(__nrho),
r_phi_arr(__r_phi_arr),
r_rho_arr(__rho_arr),
F_arr(__F_arr),
rho_phi(NULL),
drho_phi_dr(NULL),
drho_phi_dalpha(NULL),
max_pairs(0),
psi_IJ(NULL),
psi_JI(NULL),
phi_IJ(NULL),
psi_r_IJ(NULL),
psi_r_JI(NULL),
phi_r_IJ(NULL),
psi_alpha_IJ(NULL),
psi_alpha_JI(NULL),
phi_alpha_IJ(NULL),
n_phi_psi(0),
phi_psi_cmp(NULL),
phi_psi_sz(NULL),
M_IJ(NULL),
M_N_sz_sz(0),
t_ptr(NULL),
s_ptr(NULL),
crd_ptr(NULL),
mu_ptr(NULL),
dE_ptr(NULL),
E_ptr(NULL),
cv_ptr(NULL),
x_tmp_ptr(NULL),
N(atoms->N)
{
    
    __r_phi_arr=NULL;
    __rho_arr=NULL;
    __F_arr=NULL;
    
    dr_inv=1.0/dr;
    drho_inv=1.0/drho;
    rho_max=static_cast<type0>(nrho)*drho;

    stride=nelems*(nelems+1)/2;
    phi_psi_sz_sz=0;
    n_phi_psi=0;
    M_N_sz_sz=0;
    kbT=beta=-1.0;
    
    Memory::alloc(type2rho_pair_ij,nelems,nelems);
    Memory::alloc(type2rho_pair_ji,nelems,nelems);
    Memory::alloc(type2phi_pair_ij,nelems,nelems);
    
    Memory::alloc(xi,N);
    Memory::alloc(wi_0,N);
    Memory::alloc(wi_1,N);
    Memory::alloc(g_fac,nelems);
    Memory::alloc(c_0,nelems);
    Memory::alloc(c_1,nelems);
    
    memcpy(xi,atoms->xi,N*sizeof(type0));
    memcpy(wi_0,atoms->wi,N*sizeof(type0));
    for(int i=0;i<N;i++)
        wi_1[i]=wi_0[i]*xi[i];

    stride=0;
    for(elem_type ielem=0;ielem<nelems;ielem++)
        for(elem_type jelem=0;jelem<nelems;jelem++)
        {
            /*
            for(size_t i=0;i<nr;i++)
                Algebra::Do<4>::func([this,&i,&ielem,&jelem](int j){r_rho_arr[ielem][jelem][i][j]*=static_cast<type0>(i)*dr;});
             */
            if(ielem==jelem)
            {
                type2rho_pair_ij[ielem][jelem]=type2rho_pair_ji[jelem][ielem]=stride++;
                type2phi_pair_ij[ielem][jelem]=stride++;
            }
            else
            {
                type2rho_pair_ij[ielem][jelem]=stride++;
                type2rho_pair_ji[jelem][ielem]=stride++;
                type2phi_pair_ij[ielem][jelem]=stride++;
            }
        }
    
    
    for(elem_type ielem=0;ielem<nelems;ielem++)
    {
        for(size_t i=0;i<nr;i++)
            r_rho_arr[ielem][0][i][0]*=static_cast<type0>(i)*dr;
        ReadEAM::interpolate(nr,dr,r_rho_arr[ielem][0]);
        for(elem_type jelem=1;jelem<nelems;jelem++)
        {
            if(r_rho_arr[ielem][jelem]!=r_rho_arr[ielem][jelem-1])
            {
                for(size_t i=0;i<nr;i++)
                    r_rho_arr[ielem][jelem][i][0]*=static_cast<type0>(i)*dr;
                ReadEAM::interpolate(nr,dr,r_rho_arr[ielem][jelem]);
            }
        }
    }
    

    type0 tmp;
    for(size_t i=0;i<nelems;i++)
        for(size_t j=0;j<i+1;j++)
        {
            cut[i][j]=cut[j][i]=__cut[i][j];
            tmp=__cut[i][j];
            cut_sq[i][j]=cut_sq[j][i]=tmp*tmp;
        }
    Memory::dealloc(__cut);
}
/*--------------------------------------------
 destructor
 --------------------------------------------*/
ForceFieldEAMDMD::~ForceFieldEAMDMD()
{
    Memory::dealloc(F_arr);
    Memory::dealloc(r_rho_arr);
    Memory::dealloc(r_phi_arr);

    Memory::dealloc(c_1);
    Memory::dealloc(c_0);
    Memory::dealloc(g_fac);
    Memory::dealloc(wi_1);
    Memory::dealloc(wi_0);
    Memory::dealloc(xi);
    
    Memory::dealloc(type2phi_pair_ij);
    Memory::dealloc(type2rho_pair_ji);
    Memory::dealloc(type2rho_pair_ij);
}
/*--------------------------------------------
 python constructor
 --------------------------------------------*/
void ForceFieldEAMDMD::ml_new(PyMethodDef& method_0,PyMethodDef& method_1,PyMethodDef& method_2)
{
    method_0.ml_flags=METH_VARARGS | METH_KEYWORDS;
    method_0.ml_name="ff_eam_funcfl";
    method_0.ml_doc="I will add doc here";
    method_0.ml_meth=(PyCFunction)(PyCFunctionWithKeywords)
    [](PyObject* self,PyObject* args,PyObject* kwds)->PyObject*
    {
        AtomsDMD::Object* __self=reinterpret_cast<AtomsDMD::Object*>(self);
        size_t& nelems=__self->atoms->elements->nelems;
        
        FuncAPI<std::string*> f("ff_eam_funcfl",{"funcfl_files"});
        f.noptionals=1;
        f.var<0>().dynamic_size(nelems);
        if(f(args,kwds)) return NULL;
        
        size_t nr,nrho;
        type0 dr,drho;
        type0** r_c;
        type0(** F)[5]=NULL;
        type0(*** r_phi)[4]=NULL;
        type0(*** rho)[4]=NULL;
        try
        {
            ReadEAM::funcfl(nelems,f.val<0>(),dr,drho,nr,nrho,r_phi,rho,F,r_c);
        }
        catch(char* err_msg)
        {
            PyErr_SetString(PyExc_TypeError,err_msg);
            delete [] err_msg;
            return NULL;
        }
        
        delete __self->ff;
        __self->ff=new ForceFieldEAMDMD(__self->atoms,dr,drho,nr,nrho,std::move(r_phi),std::move(rho),std::move(F),std::move(r_c));
        Py_RETURN_NONE;
    };
    
    method_1.ml_flags=METH_VARARGS | METH_KEYWORDS;
    method_1.ml_name="ff_eam_setfl";
    method_1.ml_doc="I will add doc here";
    method_1.ml_meth=(PyCFunction)(PyCFunctionWithKeywords)
    [](PyObject* self,PyObject* args,PyObject* kwds)->PyObject*
    {
        AtomsDMD::Object* __self=reinterpret_cast<AtomsDMD::Object*>(self);
        size_t& nelems=__self->atoms->elements->nelems;
        FuncAPI<std::string> f("ff_eam_setfl",{"setfl_file"});
        f.noptionals=1;
        if(f(args,kwds)) return NULL;
        
        
        size_t nr,nrho;
        type0 dr,drho;
        type0** r_c;
        type0(** F)[5]=NULL;
        type0(*** r_phi)[4]=NULL;
        type0(*** rho)[4]=NULL;
        try
        {
            ReadEAM::setfl(nelems,__self->atoms->elements->names,f.val<0>(),dr,drho,nr,nrho,r_phi,rho,F,r_c);
        }
        catch(char* err_msg)
        {
            PyErr_SetString(PyExc_TypeError,err_msg);
            delete [] err_msg;
            return NULL;
        }
        
        delete __self->ff;
        __self->ff=new ForceFieldEAMDMD(__self->atoms,dr,drho,nr,nrho,std::move(r_phi),std::move(rho),std::move(F),std::move(r_c));
        Py_RETURN_NONE;
    };
    
    method_2.ml_flags=METH_VARARGS | METH_KEYWORDS;
    method_2.ml_name="ff_eam_fs";
    method_2.ml_doc="I will add doc here";
    method_2.ml_meth=(PyCFunction)(PyCFunctionWithKeywords)
    [](PyObject* self,PyObject* args,PyObject* kwds)->PyObject*
    {
        AtomsDMD::Object* __self=reinterpret_cast<AtomsDMD::Object*>(self);
        size_t& nelems=__self->atoms->elements->nelems;
        FuncAPI<std::string> f("ff_eam_fs",{"fs_file"});
        f.noptionals=1;
        if(f(args,kwds)) return NULL;
        
        size_t nr,nrho;
        type0 dr,drho;
        type0** r_c;
        type0(** F)[5]=NULL;
        type0(*** r_phi)[4]=NULL;
        type0(*** rho)[4]=NULL;
        try
        {
            ReadEAM::fs(nelems,__self->atoms->elements->names,f.val<0>(),dr,drho,nr,nrho,r_phi,rho,F,r_c);
        }
        catch(char* err_msg)
        {
            PyErr_SetString(PyExc_TypeError,err_msg);
            delete [] err_msg;
            return NULL;
        }
        
        delete __self->ff;
        __self->ff=new ForceFieldEAMDMD(__self->atoms,dr,drho,nr,nrho,std::move(r_phi),std::move(rho),std::move(F),std::move(r_c));
        Py_RETURN_NONE;
    };
}
/*--------------------------------------------
 force calculation
 --------------------------------------------*/
void ForceFieldEAMDMD::force_calc()
{
    
    
    if(max_pairs<neighbor->no_pairs)
    {
        delete [] rho_phi;
        delete [] drho_phi_dr;
        delete [] drho_phi_dalpha;
        
        max_pairs=neighbor->no_pairs;
        size_t no_0=max_pairs*3;
        Memory::alloc(rho_phi,no_0);
        Memory::alloc(drho_phi_dr,no_0);
        Memory::alloc(drho_phi_dalpha,no_0);
    }
    
    for(int i=0;i<max_pairs*3;i++) rho_phi[i]=drho_phi_dr[i]=drho_phi_dalpha[i]=0.0;
    
    
    
    
    
    type0 r,r_inv;
    size_t m;
    type0* coef;
    type0 tmp0,tmp1;
    type0 fpair,apair;
    
    type0 p,c_iv;
    
    type0 dx_ij[__dim__];
    
    type0 const* c=atoms->c->begin();

    
    
    type0* dE=dE_ptr->begin();
    type0* mu=mu_ptr->begin();
    type0* E=E_ptr->begin();
    const int n=atoms->natms*c_dim;
    for(int i=0;i<n;i++) E[i]=0.0;
    
    elem_type const* elem_vec=atoms->elem->begin();
    
    type0 alpha_ij,rsq;
    elem_type ielem,jelem;
    
    
    
    const type0* x=atoms->x->begin();
    const type0* alpha=atoms->alpha->begin();
    int** neighbor_list=neighbor->neighbor_list;
    int* neighbor_list_size=neighbor->neighbor_list_size;
    size_t istart=0;
    for(int i=0;i<n;i++)
    {
        if(c[i]<0.0) continue;
        ielem=elem_vec[i];
        const int neigh_sz=neighbor_list_size[i];
        for(int j,__j=0;__j<neigh_sz;__j++,istart+=3)
        {
            j=neighbor_list[i][__j];
            jelem=elem_vec[j];
            rsq=Algebra::RSQ<__dim__>(x+(i/c_dim)*__dim__,x+(j/c_dim)*__dim__);
            r=sqrt(rsq);
            alpha_ij=sqrt(alpha[i]*alpha[i]+alpha[j]*alpha[j]);
            if(r-alpha_ij*xi[N-1]>=cut[ielem][jelem]) continue;
            type0 upper=(r+cut[ielem][jelem])/alpha_ij;
            type0 lower=(r-cut[ielem][jelem])/alpha_ij;
            type0 __r,p,tmp0,tmp1;
            type0* coef;

            
            r_inv=1.0/r;
            
            type0 __rho_phi[3]{[0 ... 2]=0.0};
            type0 __drho_phi_dr[3]{[0 ... 2]=0.0};
            type0 __drho_phi_dalpha[3]{[0 ... 2]=0.0};
            for(int l=0;l<N;l++)
            {
                if(xi[l]<=lower && xi[l]>=upper) continue;
                
                __r=r-xi[l]*alpha_ij;
                p=fabs(__r)*dr_inv;
                m=static_cast<size_t>(p);
                m=MIN(m,nr-2);
                p-=m;
                p=MIN(p,1.0);
                
                coef=r_phi_arr[ielem][jelem][m];
                tmp0=((coef[3]*p+coef[2])*p+coef[1])*p+coef[0];
                tmp1=((3.0*coef[3]*p+2.0*coef[2])*p+coef[1])*dr_inv;
                if(__r<0.0) tmp0*=-1.0;
                
                __rho_phi[0]+=wi_0[l]*tmp0;
                __drho_phi_dr[0]+=wi_0[l]*tmp1;
                __drho_phi_dalpha[0]+=wi_1[l]*tmp1;
                
                coef=r_rho_arr[ielem][jelem][m];
                tmp0=((coef[3]*p+coef[2])*p+coef[1])*p+coef[0];
                tmp1=((3.0*coef[3]*p+2.0*coef[2])*p+coef[1])*dr_inv;
                if(__r<0.0) tmp0*=-1.0;
                
                __rho_phi[1]+=wi_0[l]*tmp0;
                __drho_phi_dr[1]+=wi_0[l]*tmp1;
                __drho_phi_dalpha[1]+=wi_1[l]*tmp1;
                
                coef=r_rho_arr[jelem][ielem][m];
                tmp0=((coef[3]*p+coef[2])*p+coef[1])*p+coef[0];
                tmp1=((3.0*coef[3]*p+2.0*coef[2])*p+coef[1])*dr_inv;
                if(__r<0.0) tmp0*=-1.0;
                
                __rho_phi[2]+=wi_0[l]*tmp0;
                __drho_phi_dr[2]+=wi_0[l]*tmp1;
                __drho_phi_dalpha[2]+=wi_1[l]*tmp1;
            }
            
            tmp0=PI_IN_SQ*r_inv;
            
            Algebra::Do<3>::func([&__rho_phi,&__drho_phi_dr,&__drho_phi_dalpha,&tmp0,&r_inv,&alpha_ij,&istart,this]
            (int i)
            {
                __rho_phi[i]*=tmp0;
                __drho_phi_dr[i]*=tmp0;
                __drho_phi_dr[i]-=__rho_phi[i]*r_inv;
                __drho_phi_dr[i]*=r_inv;
                __drho_phi_dalpha[i]*=-tmp0/alpha_ij;
                
                rho_phi[istart+i]=__rho_phi[i];
                drho_phi_dr[istart+i]=__drho_phi_dr[i];
                drho_phi_dalpha[istart+i]=__drho_phi_dalpha[i];
            });
      
            
            E[i]+=c[j]*__rho_phi[2];
            
            if(j<n)
            {
                E[j]+=c[i]*__rho_phi[1];
                nrgy_strss_lcl[0]+=c[i]*c[j]*__rho_phi[0];
            }
            else
                nrgy_strss_lcl[0]+=0.5*c[i]*c[j]*__rho_phi[0];
        }
        
        
        
        p=E[i]*drho_inv;
        m=static_cast<size_t>(p);
        m=MIN(m,nrho-2);
        p-=m;
        p=MIN(p,1.0);
        coef=F_arr[ielem][m];
        
        tmp0=(((coef[4]*p+coef[3])*p+coef[2])*p+coef[1])*p+coef[0];
        tmp1=(((4.0*coef[4]*p+3.0*coef[3])*p+2.0*coef[2])*p+coef[1])*drho_inv;
        
        if(E[i]>rho_max) tmp0+=tmp1*(E[i]-rho_max);
        
        E[i]=tmp0;
        dE[i]=tmp1;
        mu[i]=tmp0;
        if(c[i]!=0.0)
            nrgy_strss_lcl[0]+=c[i]*(tmp0+c_0[ielem]-3.0*kbT*log(alpha[i]));
        
        nrgy_strss_lcl[0]+=kbT*calc_ent(c[i]);
    }
    

    
    const int __natms=atoms->natms;
    
    for(int i=0;i<__natms;i++)
    {
        c_iv=1.0;
        for(int ic=0;ic<c_dim;ic++)
            if(c[i*c_dim+ic]>0.0)
                c_iv-=c[i*c_dim+ic];
        nrgy_strss_lcl[0]+=kbT*calc_ent(c_iv);
    }
    
    dynamic->update(dE_ptr);

    
    type0* fvec=f->begin();
    type0* f_alphavec=f_alpha->begin();
    type0 f_i[__dim__]={[0 ... __dim__-1]=0.0};
    type0 x_i[__dim__];
    istart=0;
    for(int i=0;i<n;i++)
    {
        
        
        if(i%c_dim==0)
        {
            Algebra::zero<__dim__>(f_i);
            Algebra::V_eq<__dim__>(x+(i/c_dim)*__dim__,x_i);
        }
        
        if(c[i]<0.0) continue;
        
        type0 c_i=c[i];
        type0 dE_i=dE[i];
        type0 f_alpha_i=0.0;
        type0 mu_i=0.0;
        const int neigh_sz=neighbor_list_size[i];
        for(int j,__j=0;__j<neigh_sz;__j++,istart+=3)
        {
            j=neighbor_list[i][__j];
            Algebra::DX<__dim__>(x_i,x+(j/c_dim)*__dim__,dx_ij);
            
            fpair=-(drho_phi_dr[istart+2]*dE_i+drho_phi_dr[istart+1]*dE[j]+drho_phi_dr[istart])*c_i*c[j];
            apair=-(drho_phi_dalpha[istart+2]*dE_i+drho_phi_dalpha[istart+1]*dE[j]+drho_phi_dalpha[istart])*c_i*c[j];
            mu_i+=c[j]*(rho_phi[istart]+rho_phi[istart+1]*dE[j]);
            if(j<n) mu[j]+=c_i*(rho_phi[istart]+rho_phi[istart+2]*dE_i);
            
            Algebra::V_add_x_mul_V<__dim__>(fpair,dx_ij,f_i);
            f_alpha_i+=alpha[i]*apair;
            
            if(j<n)
            {
                Algebra::V_add_x_mul_V<__dim__>(-fpair,dx_ij,fvec+(j/c_dim)*x_dim);
                f_alphavec[j]+=alpha[j]*apair;
            }
            else
                fpair*=0.5;
            Algebra::DyadicV(-fpair,dx_ij,&nrgy_strss_lcl[1]);
        }
        
        f_alpha_i+=3.0*kbT*c[i]/alpha[i];
        f_alphavec[i]+=f_alpha_i;
        mu[i]+=mu_i;
        
        if(i%c_dim==0)
            Algebra::V_add<__dim__>(f_i,fvec+__dim__*(i/c_dim));
    }
    
}
/*--------------------------------------------
 energy calculation
 --------------------------------------------*/
void ForceFieldEAMDMD::energy_calc()
{
    type0* E=E_ptr->begin();
    const int n=atoms->natms*c_dim;
    for(int i=0;i<n;i++) E[i]=0.0;
    
    elem_type const* elem_vec=atoms->elem->begin();
    
    type0 alpha_ij,rsq,r,r_inv,p,tmp0,tmp1;
    elem_type ielem,jelem;
    size_t m;
    
    const type0* x=atoms->x->begin();
    const type0* alpha=atoms->alpha->begin();
    int** neighbor_list=neighbor->neighbor_list;
    int* neighbor_list_size=neighbor->neighbor_list_size;
    type0 const* c=atoms->c->begin();
    type0 const* coef;

    size_t istart=0;
    for(int i=0;i<n;i++)
    {
        ielem=elem_vec[i];
        const int neigh_sz=neighbor_list_size[i];
        for(int j,__j=0;__j<neigh_sz;__j++,istart+=3)
        {
            j=neighbor_list[i][__j];
            jelem=elem_vec[j];
            rsq=Algebra::RSQ<__dim__>(x+(i/c_dim)*__dim__,x+(j/c_dim)*__dim__);
            r=sqrt(rsq);
            alpha_ij=sqrt(alpha[i]*alpha[i]+alpha[j]*alpha[j]);
            if(r-alpha_ij*xi[N-1]>=cut[ielem][jelem]) continue;
            
            type0 upper=(r+cut[ielem][jelem])/alpha_ij;
            type0 lower=(r-cut[ielem][jelem])/alpha_ij;
            type0 __r,p,tmp0,tmp1;
            type0* coef;
            
            r_inv=1.0/r;
            
            type0 __arr[3]{[0 ... 2]=0.0};
            for(int l=0;l<N;l++)
            {
                if(xi[l]<=lower && xi[l]>=upper) continue;
                
                __r=r-xi[l]*alpha_ij;
                p=fabs(__r)*dr_inv;
                m=static_cast<size_t>(p);
                m=MIN(m,nr-2);
                p-=m;
                p=MIN(p,1.0);
                
                coef=r_phi_arr[ielem][jelem][m];
                tmp0=((coef[3]*p+coef[2])*p+coef[1])*p+coef[0];
                tmp1=((3.0*coef[3]*p+2.0*coef[2])*p+coef[1])*dr_inv;
                if(__r<0.0) tmp0*=-1.0;
                __arr[0]+=wi_0[l]*tmp0;

                
                coef=r_rho_arr[ielem][jelem][m];
                tmp0=((coef[3]*p+coef[2])*p+coef[1])*p+coef[0];
                tmp1=((3.0*coef[3]*p+2.0*coef[2])*p+coef[1])*dr_inv;
                if(__r<0.0) tmp0*=-1.0;
                __arr[1]+=wi_0[l]*tmp0;
                
                coef=r_rho_arr[jelem][ielem][m];
                tmp0=((coef[3]*p+coef[2])*p+coef[1])*p+coef[0];
                tmp1=((3.0*coef[3]*p+2.0*coef[2])*p+coef[1])*dr_inv;
                if(__r<0.0) tmp0*=-1.0;
                __arr[2]+=wi_0[l]*tmp0;
            }
            
            tmp0=PI_IN_SQ*r_inv;
            
            Algebra::Do<3>::func([&__arr,&tmp0,this](int i)
            {
                __arr[i]*=tmp0;
            });
            
            //printf("%e %e %.12e %.12e %.12e\n",r,alpha_ij,__drho_phi_dr[0],__drho_phi_dr[2],__drho_phi_dr[1]);
            E[i]+=c[j]*__arr[2];
            
            if(j<n)
            {
                E[j]+=c[i]*__arr[1];
                nrgy_strss_lcl[0]+=c[i]*c[j]*__arr[0];
            }
            else
                nrgy_strss_lcl[0]+=0.5*c[i]*c[j]*__arr[0];
        }
        
        
        if(c[i]<0.0) continue;
        tmp0=E[i];
        p=tmp0*drho_inv;
        m=static_cast<size_t>(p);
        m=MIN(m,nrho-2);
        p-=m;
        p=MIN(p,1.0);
        coef=F_arr[ielem][m];
        tmp1=(((coef[4]*p+coef[3])*p+coef[2])*p+coef[1])*p+coef[0];
        if(E[i]>rho_max)
            tmp1+=(((4.0*coef[4]*p+3.0*coef[3])*p+2.0*coef[2])*p+coef[1])*drho_inv*(tmp0-rho_max);
        
        if(c[i]!=0.0)
            nrgy_strss_lcl[0]+=c[i]*(tmp1+c_0[ielem]-3.0*kbT*log(alpha[i]));
        nrgy_strss_lcl[0]+=kbT*calc_ent(c[i]);
    }
    const int __natms=atoms->natms;
    type0 c_iv;
    for(int i=0;i<__natms;i++)
    {
        c_iv=1.0;
        for(int ic=0;ic<c_dim;ic++)
            if(c[i*c_dim+ic]>0.0)
                c_iv-=c[i*c_dim+ic];
        nrgy_strss_lcl[0]+=kbT*calc_ent(c_iv);
    }
}
/*--------------------------------------------
 init
 --------------------------------------------*/
void ForceFieldEAMDMD::init()
{
    setup();
    set_temp(300.0);
    cv_ptr=new Vec<type0>(atoms,1);
    E_ptr=new Vec<type0>(atoms,c_dim);
    dE_ptr=new Vec<type0>(atoms,c_dim);
    mu_ptr=new Vec<type0>(atoms,c_dim,"mu");

    crd_ptr=new Vec<type0>(atoms,c_dim);
    s_ptr=new Vec<type0>(atoms,c_dim);
    x_tmp_ptr=new Vec<type0>(atoms,1);
    t_ptr=new Vec<type0>(atoms,c_dim);
    
    type0* mu=mu_ptr->begin();
    const int n=atoms->natms*c_dim;
    for(int i=0;i<n;i++) mu[i]=0.0;
}
/*--------------------------------------------
 fin
 --------------------------------------------*/
void ForceFieldEAMDMD::fin()
{
    
    Memory::dealloc(rho_phi);
    Memory::dealloc(drho_phi_dr);
    Memory::dealloc(drho_phi_dalpha);
    max_pairs=0;
    
    Memory::dealloc(phi_psi_sz);
    phi_psi_sz_sz=0;

    
    Memory::dealloc(psi_IJ);
    Memory::dealloc(psi_JI);
    Memory::dealloc(phi_IJ);
    Memory::dealloc(psi_r_IJ);
    Memory::dealloc(psi_r_JI);
    Memory::dealloc(phi_r_IJ);
    Memory::dealloc(psi_alpha_IJ);
    Memory::dealloc(psi_alpha_JI);
    Memory::dealloc(phi_alpha_IJ);
    n_phi_psi=0;
    
    
    Memory::dealloc(M_IJ);
    M_N_sz_sz=0;
    
    delete t_ptr;
    delete s_ptr;
    delete crd_ptr;
    delete mu_ptr;
    delete dE_ptr;
    delete E_ptr;
    delete cv_ptr;
    delete x_tmp_ptr;
    t_ptr=s_ptr=crd_ptr=mu_ptr=dE_ptr=E_ptr=cv_ptr=x_tmp_ptr=NULL;
}
/*--------------------------------------------
 set the temperature in the simulation
 --------------------------------------------*/
void ForceFieldEAMDMD::set_temp(type0 T)
{
    type0 kb=8.617332478e-5;
    type0 hbar=6.5821192815e-16;
    type0 mass;
    type0 deb_l;
    
    for(size_t i=0;i<nelems;i++)
    {
        mass=atoms->elements->masses[i];
        c_1[i]=sqrt(0.5*kb*T/mass)/M_PI;
        mass*=1.0364269184093291236e-28;
        deb_l=hbar*hbar*2.0/(mass*kb*T);
        c_0[i]=1.5*kb*T*(log(deb_l)-1.0);
    }
    
    kbT=kb*T;
    beta=1.0/kbT;
}
/*--------------------------------------------
 return M_{ij}^{\alpha}
 --------------------------------------------*/
inline type0 ForceFieldEAMDMD::calc_ent(type0 x)
{
    type0 ans=x*log(x);
    if(isnan(ans))
        return 0.0;
    return ans;
}
/*--------------------------------------------
 claculate F and dF and dFF
 --------------------------------------------*/
void ForceFieldEAMDMD::dc()
{
    calc_mu();

    type0* xvec=atoms->x->begin();
    type0* c=atoms->c->begin();
    type0* c_dvec=atoms->c_d->begin();
    elem_type* elem_vec=atoms->elem->begin();
    type0* cv=cv_ptr->begin();
    type0* mu=mu_ptr->begin();
    const int n=atoms->natms*c_dim;
    for(int i=0;i<n;i++)
        c_dvec[i]=0.0;
    
    type0 gamma_i,gamma_j,x_ij_sq,alpha_Q_sq;
    type0* xi;
    type0* xj;
    type0 mu_ji,alpha_i,alpha_j,Qi;
    
    int iatm,jatm,icmp,jcmp;
    type0 d_i,d_j,dc_ij;
    
    int** neighbor_list_2nd=neighbor_dmd->neighbor_list_2nd;
    int* neighbor_list_size_2nd=neighbor_dmd->neighbor_list_size_2nd;
    
    for(int ic_dim=0,ielem;ic_dim<n;ic_dim++)
    {
        ielem=elem_vec[ic_dim];
        icmp=ic_dim%c_dim;
        iatm=ic_dim/c_dim;
        xi=xvec+x_dim*iatm;
        alpha_i=xi[3+icmp];
        
        for(int j=0,jc_dim;j<neighbor_list_size_2nd[ic_dim];j++)
        {
            
            jc_dim=neighbor_list_2nd[ic_dim][j];
            jcmp=jc_dim%c_dim;
            jatm=jc_dim/c_dim;
            xj=xvec+x_dim*jatm;
            alpha_j=xj[3+icmp];
            mu_ji=beta*(mu[jc_dim]-mu[ic_dim]);
            
            x_ij_sq=(xi[0]-xj[0])*(xi[0]-xj[0])
            +(xi[1]-xj[1])*(xi[1]-xj[1])
            +(xi[2]-xj[2])*(xi[2]-xj[2]);
            gamma_i=x_ij_sq/(alpha_i*alpha_i);
            gamma_j=x_ij_sq/(alpha_j*alpha_j);
            
            calc_Q(ielem,gamma_i,gamma_j,mu_ji,Qi,alpha_Q_sq);
            
            d_i=d_j=-c_1[ielem]*x_ij_sq*alpha_Q_sq;
            d_i/=alpha_i*alpha_i*alpha_i;
            d_j/=alpha_j*alpha_j*alpha_j;
            
            dc_ij=d_i*exp(-Qi)*c[ic_dim]*cv[jatm]-d_j*exp(-Qi+mu_ji)*c[jc_dim]*cv[iatm];

            c_dvec[ic_dim]+=dc_ij;
            if(jc_dim<n)
            {
                c_dvec[jc_dim]-=dc_ij;
            }
        }
    }
}
/*--------------------------------------------
 calculate norm of d^2c/dt^2
 --------------------------------------------*/
type0 ForceFieldEAMDMD::ddc_norm()
{
    const int n=atoms->natms*c_dim;
    type0* tmp=NULL;
    Memory::alloc(tmp,n);
    for(int i=0;i<n;i++) tmp[i]=0;
    update_J(1.0,tmp,tmp);
    delete [] tmp;
    
    type0* c_dvec=atoms->c_d->begin();
    type0* s=s_ptr->begin();
    memcpy(s,c_dvec,sizeof(type0)*n);
    operator()(s_ptr,t_ptr);
    type0* t=t_ptr->begin();
    type0* c=atoms->c->begin();
    type0 ans_lcl=0.0,ans,tmp0;
    
    for(int ic_dim=0;ic_dim<phi_psi_sz_sz;ic_dim++)
        if(c[ic_dim]>=0.0)
        {
            tmp0=t[ic_dim]+s[ic_dim];
            ans_lcl+=tmp0*tmp0;
        }
    
    MPI_Allreduce(&ans_lcl,&ans,1,Vec<type0>::MPI_T,MPI_SUM,world);
    return sqrt(ans);
}
/*--------------------------------------------
 calculate norm of d^2c/dt^2
 --------------------------------------------*/
void ForceFieldEAMDMD::ddc(type0* ddc_)
{
}
/*--------------------------------------------
 calculate mu crd dE ddE and local energy
 --------------------------------------------*/
void ForceFieldEAMDMD::calc_mu()
{
    type0* c=atoms->c->begin();
    elem_type* elem_vec=atoms->elem->begin();
    type0* xvec=atoms->x->begin();
    type0* E=E_ptr->begin();
    type0* mu=mu_ptr->begin();
    type0* dE=dE_ptr->begin();
    type0* cv=cv_ptr->begin();
    type0* ddE=E;
    const int n=atoms->natms*c_dim;
    for(int i=0;i<n;i++)
        mu[i]=E[i]=0.0;
    int ielem;
    size_t m;
    type0 p,tmp0,tmp1,tmp2;
    type0* coef;
    
    nrgy_strss_lcl[0]=0.0;
    const int nall=atoms->natms+atoms->natms_ph;
    for(int i=0;i<nall;i++)
    {
        cv[i]=1.0;
        for(int j=0;j<c_dim;j++)
            if(c[i*c_dim+j]>=0.0)
                cv[i]-=c[i*c_dim+j];
    }
    
    int istart=0;
    for(int ic_dim=0;ic_dim<n;ic_dim++)
    {
        for(int j=0,jc_dim;j<phi_psi_sz[ic_dim];j++)
        {
            jc_dim=phi_psi_cmp[istart];
            E[ic_dim]+=c[jc_dim]*psi_JI[istart];
            mu[ic_dim]+=c[jc_dim]*phi_IJ[istart];
            
            
            if(jc_dim<phi_psi_sz_sz)
            {
                E[jc_dim]+=c[ic_dim]*psi_IJ[istart];
                mu[jc_dim]+=c[ic_dim]*phi_IJ[istart];
                nrgy_strss_lcl[0]+=c[ic_dim]*c[jc_dim]*phi_IJ[istart];
            }
            else
                nrgy_strss_lcl[0]+=0.5*c[ic_dim]*c[jc_dim]*phi_IJ[istart];
                
            istart++;
        }
        if(c[ic_dim]>=0.0)
        {
            ielem=elem_vec[ic_dim];
            p=E[ic_dim]*drho_inv;
            m=static_cast<size_t>(p);
            m=MIN(m,nrho-2);
            p-=m;
            p=MIN(p,1.0);
            coef=F_arr[ielem][m];
            tmp0=(((coef[4]*p+coef[3])*p+coef[2])*p+coef[1])*p+coef[0];
            tmp1=(((4.0*coef[4]*p+3.0*coef[3])*p+2.0*coef[2])*p+coef[1])*drho_inv;
            tmp2=(((12.0*coef[4]*p+6.0*coef[3])*p+2.0*coef[2]))*drho_inv*drho_inv;
            if(E[ic_dim]>rho_max)
            {
                tmp0+=tmp1*(E[ic_dim]-rho_max);
                tmp2=0.0;
            }
            
            dE[ic_dim]=tmp1;
            ddE[ic_dim]=tmp2;
            mu[ic_dim]+=tmp0;
            tmp2=tmp0+c_0[ielem];
            if(c[ic_dim]!=0.0)
                tmp2-=3.0*kbT*log(xvec[3*(ic_dim/c_dim +1)+ic_dim]);
            nrgy_strss_lcl[0]+=c[ic_dim]*tmp2+kbT*calc_ent(c[ic_dim]);
        }
        if(ic_dim%c_dim==c_dim-1)
            nrgy_strss_lcl[0]+=kbT*calc_ent(cv[ic_dim/c_dim]);
    }
    
    dynamic->update(dE_ptr);
    
    istart=0;
    for(int ic_dim=0;ic_dim<phi_psi_sz_sz;ic_dim++)
    {
        for(int j=0,jc_dim;j<phi_psi_sz[ic_dim];j++)
        {
            jc_dim=phi_psi_cmp[istart];
            mu[ic_dim]+=psi_IJ[istart]*c[jc_dim]*dE[jc_dim];
            
            if(jc_dim<phi_psi_sz_sz)
            {
                mu[jc_dim]+=psi_JI[istart]*c[ic_dim]*dE[ic_dim];
            }
            istart++;
        }
    }
    
    dynamic->update(mu_ptr);
}
/*--------------------------------------------
 create the sparse matrices
 --------------------------------------------*/
void ForceFieldEAMDMD::init_static()
{
    delete [] phi_psi_sz;
    
    delete [] psi_IJ;
    delete [] psi_JI;
    delete [] phi_IJ;
    delete [] psi_r_IJ;
    delete [] psi_r_JI;
    delete [] phi_r_IJ;
    delete [] psi_alpha_IJ;
    delete [] psi_alpha_JI;
    delete [] phi_alpha_IJ;

    delete [] M_IJ;
    
    size_t no_pairs=neighbor->no_pairs;
    phi_psi_sz_sz=atoms->natms*c_dim;
    
    Memory::alloc(phi_psi_sz,phi_psi_sz_sz);
    Memory::alloc(phi_psi_cmp,c_dim*c_dim*no_pairs);
    Memory::alloc(psi_IJ,c_dim*c_dim*no_pairs);
    Memory::alloc(psi_JI,c_dim*c_dim*no_pairs);
    Memory::alloc(phi_IJ,c_dim*c_dim*no_pairs);
    Memory::alloc(psi_r_IJ,c_dim*c_dim*no_pairs);
    Memory::alloc(psi_r_JI,c_dim*c_dim*no_pairs);
    Memory::alloc(phi_r_IJ,c_dim*c_dim*no_pairs);
    Memory::alloc(psi_alpha_IJ,c_dim*c_dim*no_pairs);
    Memory::alloc(psi_alpha_JI,c_dim*c_dim*no_pairs);
    Memory::alloc(phi_alpha_IJ,c_dim*c_dim*no_pairs);
    
    for(int i=0;i<phi_psi_sz_sz;i++)
        phi_psi_sz[i]=0;
    
    
    type0* c=atoms->c->begin();
    elem_type* elem_vec=atoms->elem->begin();
    int** neighbor_list=neighbor->neighbor_list;
    int* neighbor_list_size=neighbor->neighbor_list_size;
    
    int istart=0;
    int iistart=0;
    size_t isps=0;
    type0 psi_ij,psi_ji,phi_ij;
    type0 psi_r_ij,psi_r_ji,phi_r_ij;
    type0 psi_alpha_ij,psi_alpha_ji,phi_alpha_ij;
    const int natms=atoms->natms;
    for(int iatm=0;iatm<natms;iatm++)
    {
        
        for(int ic_dim=iatm*c_dim,ielem;ic_dim<(iatm+1)*c_dim && c[ic_dim]>=0.0;ic_dim++)
        {
            ielem=elem_vec[ic_dim];
            
            iistart=istart;
            for(int j=0,jatm;j<neighbor_list_size[iatm];j++)
            {
                jatm=neighbor_list[iatm][j];

                for(int jc_dim=jatm*c_dim,jelem;jc_dim<(jatm+1)*c_dim && c[jc_dim]>=0.0;jc_dim++)
                {
                    jelem=elem_vec[jc_dim];
                    
                    psi_ij=rho_phi[iistart+type2rho_pair_ij[ielem][jelem]];
                    psi_ji=rho_phi[iistart+type2rho_pair_ji[jelem][ielem]];
                    phi_ij=rho_phi[iistart+type2phi_pair_ij[ielem][jelem]];
                    
                    psi_r_ij=drho_phi_dr[iistart+type2rho_pair_ij[ielem][jelem]];
                    psi_r_ji=drho_phi_dr[iistart+type2rho_pair_ji[jelem][ielem]];
                    phi_r_ij=drho_phi_dr[iistart+type2phi_pair_ij[ielem][jelem]];
                    
                    psi_alpha_ij=drho_phi_dalpha[iistart+type2rho_pair_ij[ielem][jelem]];
                    psi_alpha_ji=drho_phi_dalpha[iistart+type2rho_pair_ji[jelem][ielem]];
                    phi_alpha_ij=drho_phi_dalpha[iistart+type2phi_pair_ij[ielem][jelem]];
                    
                    if(psi_ij!=0.0 || psi_ji!=0.0 || phi_ij!=0.0 ||
                       psi_r_ij!=0.0 || psi_r_ji!=0.0 || phi_r_ij!=0.0 ||
                       psi_alpha_ij!=0.0 || psi_alpha_ji!=0.0 || phi_alpha_ij!=0.0)
                    {
                        phi_psi_sz[ic_dim]++;
                        phi_psi_cmp[isps]=jc_dim;
                        
                        psi_IJ[isps]=psi_ij;
                        psi_JI[isps]=psi_ji;
                        phi_IJ[isps]=phi_ij;
                        
                        psi_r_IJ[isps]=psi_r_ij;
                        psi_r_JI[isps]=psi_r_ji;
                        phi_r_IJ[isps]=phi_r_ij;
                        
                        psi_alpha_IJ[isps]=psi_alpha_ij;
                        psi_alpha_JI[isps]=psi_alpha_ji;
                        phi_alpha_IJ[isps]=phi_alpha_ij;
                        
                        isps++;
                    }
                }
                iistart+=stride;
            }
        }
        istart=iistart;
    }
    if(max_pairs)
    {
        delete [] rho_phi;
        delete [] drho_phi_dr;
        delete [] drho_phi_dalpha;
        max_pairs=0;
    }
    
    GROW(phi_psi_cmp,c_dim*c_dim*no_pairs,isps);
    GROW(psi_IJ,c_dim*c_dim*no_pairs,isps);
    GROW(psi_JI,c_dim*c_dim*no_pairs,isps);
    GROW(phi_IJ,c_dim*c_dim*no_pairs,isps);
    GROW(psi_r_IJ,c_dim*c_dim*no_pairs,isps);
    GROW(psi_r_JI,c_dim*c_dim*no_pairs,isps);
    GROW(phi_r_IJ,c_dim*c_dim*no_pairs,isps);
    GROW(psi_alpha_IJ,c_dim*c_dim*no_pairs,isps);
    GROW(psi_alpha_JI,c_dim*c_dim*no_pairs,isps);
    GROW(phi_alpha_IJ,c_dim*c_dim*no_pairs,isps);
    n_phi_psi=isps;

    
    
    
    int* neighbor_list_size_2nd=neighbor_dmd->neighbor_list_size_2nd;
    M_N_sz_sz=0;
    const int n=atoms->natms*c_dim;
    for(int i=0;i<n;i++)
        M_N_sz_sz+=neighbor_list_size_2nd[i];
    
    Memory::alloc(M_IJ,3*M_N_sz_sz);
    
}
/*--------------------------------------------
 create the sparse matrices
 --------------------------------------------*/
type0 ForceFieldEAMDMD::update_J(type0 alpha,type0* a,type0* g)
{
    calc_mu();
    
    type0* c=atoms->c->begin();
    type0* xvec=atoms->x->begin();
    elem_type* elem_vec=atoms->elem->begin();
    type0* mu=mu_ptr->begin();
    type0* cv=cv_ptr->begin();
    
    int istart;
    
    
    type0* c_dvec=atoms->c_d->begin();
    for(int i=0;i<phi_psi_sz_sz;i++)
    {
        c_dvec[i]=0.0;
        if(c[i]>=0.0)
            g[i]=c[i]+a[i];
        else
            g[i]=0.0;
    }
    

    
    type0 iota=log(alpha);
    alpha_tmp=alpha;
    
    type0 gamma_i,gamma_j,x_ij_sq,alpha_Q_sq;
    type0 theta_i,theta_j;

    type0* xi;
    type0* xj;
    type0 ans_lcl=0.0,dc_ij,dg_ij,mu_ji,alpha_i,alpha_j;
    type0 d_i,d_j,Qi;
    type0 exp_mod_Qi,exp_mod_Qj;
    int iatm,jatm,icmp,jcmp;
    int** neighbor_list_2nd=neighbor_dmd->neighbor_list_2nd;
    int* neighbor_list_size_2nd=neighbor_dmd->neighbor_list_size_2nd;
    
    istart=0;
    for(int ic_dim=0,ielem;ic_dim<phi_psi_sz_sz;ic_dim++)
    {
        ielem=elem_vec[ic_dim];
        icmp=ic_dim%c_dim;
        iatm=ic_dim/c_dim;
        xi=xvec+x_dim*iatm;
        alpha_i=xi[3+icmp];
        
        for(int j=0,jc_dim;j<neighbor_list_size_2nd[ic_dim];j++)
        {
            jc_dim=neighbor_list_2nd[ic_dim][j];
            jcmp=jc_dim%c_dim;
            jatm=jc_dim/c_dim;
            xj=xvec+x_dim*jatm;
            alpha_j=xj[3+icmp];
            mu_ji=beta*(mu[jc_dim]-mu[ic_dim]);
            
            x_ij_sq=(xi[0]-xj[0])*(xi[0]-xj[0])
            +(xi[1]-xj[1])*(xi[1]-xj[1])
            +(xi[2]-xj[2])*(xi[2]-xj[2]);
            gamma_i=x_ij_sq/(alpha_i*alpha_i);
            gamma_j=x_ij_sq/(alpha_j*alpha_j);
            
            calc_Q(ielem,gamma_i,gamma_j,mu_ji,Qi,alpha_Q_sq,theta_i);
            
            exp_mod_Qi=exp(-Qi+iota);
            exp_mod_Qj=exp(-Qi+mu_ji+iota);
            
            d_i=d_j=-c_1[ielem]*x_ij_sq*alpha_Q_sq;
            d_i/=alpha_i*alpha_i*alpha_i;
            d_j/=alpha_j*alpha_j*alpha_j;
            theta_j=theta_i+1.0;
            
            dc_ij=d_i*exp(-Qi)*c[ic_dim]*cv[jatm]-d_j*exp(-Qi+mu_ji)*c[jc_dim]*cv[iatm];
            dg_ij=d_i*exp_mod_Qi*c[ic_dim]*cv[jatm]-d_j*exp_mod_Qj*c[jc_dim]*cv[iatm];

            M_IJ[istart++]=kbT*(d_i*theta_i*exp_mod_Qi*c[ic_dim]*cv[jatm]-d_j*theta_j*exp_mod_Qj*c[jc_dim]*cv[iatm]);
            M_IJ[istart++]=d_i*exp_mod_Qi;
            M_IJ[istart++]=d_j*exp_mod_Qj;
            g[ic_dim]-=dg_ij;
            c_dvec[ic_dim]+=dc_ij;
            
            if(jc_dim<phi_psi_sz_sz)
            {
                g[jc_dim]+=dg_ij;
                c_dvec[jc_dim]-=dc_ij;
            }
        }
        if(c[ic_dim]>=0.0)
        {
            
            ans_lcl+=g[ic_dim]*g[ic_dim];
        }
    }
    type0 ans;
    MPI_Allreduce(&ans_lcl,&ans,1,Vec<type0>::MPI_T,MPI_SUM,world);
    return sqrt(ans);
    
}
/*--------------------------------------------
 create the sparse matrices
 --------------------------------------------*/
void ForceFieldEAMDMD::operator()(Vec<type0>* x_ptr,Vec<type0>* Ax_ptr)
{
    dynamic->update(x_ptr);
    
    type0 tmp0;
    type0* c=atoms->c->begin();
    type0* ddE=E_ptr->begin();
    type0* dE=dE_ptr->begin();
    type0* x=x_ptr->begin();
    type0* Ax=Ax_ptr->begin();
    type0* b0=crd_ptr->begin();
    
    for(int ic_dim=0;ic_dim<phi_psi_sz_sz;ic_dim++)
        b0[ic_dim]=Ax[ic_dim]=0.0;
    
    int istart=0;
    
    for(int ic_dim=0;ic_dim<phi_psi_sz_sz;ic_dim++)
        for(int j=0,jc_dim;j<phi_psi_sz[ic_dim];j++)
        {
            jc_dim=phi_psi_cmp[istart];
            Ax[ic_dim]+=c[ic_dim]*ddE[ic_dim]*psi_JI[istart]*x[jc_dim];
            if(jc_dim<phi_psi_sz_sz)
                Ax[jc_dim]+=c[jc_dim]*ddE[jc_dim]*psi_IJ[istart]*x[ic_dim];
            istart++;
        }
    
    dynamic->update(Ax_ptr);
    
    istart=0;
    for(int ic_dim=0;ic_dim<phi_psi_sz_sz;ic_dim++)
        for(int j=0,jc_dim;j<phi_psi_sz[ic_dim];j++)
        {
            jc_dim=phi_psi_cmp[istart];
            tmp0=psi_IJ[istart]*dE[ic_dim]+psi_JI[istart]*dE[jc_dim]+phi_IJ[istart];
            
            b0[ic_dim]+=psi_IJ[istart]*Ax[jc_dim]+tmp0*x[jc_dim];
            if(jc_dim<phi_psi_sz_sz)
                b0[jc_dim]+=psi_JI[istart]*Ax[ic_dim]+tmp0*x[ic_dim];
            
            istart++;
        }
    
    dynamic->update(crd_ptr);
    
    type0* x_tmp=x_tmp_ptr->begin();
    type0* cv=cv_ptr->begin();
    const int nall=atoms->natms+atoms->natms_ph;
    for(int i=0;i<nall;i++)
    {
        x_tmp[i]=0.0;
        for(int j=0;j<c_dim;j++)
            if(c[i*c_dim+j]>=0.0)
                x_tmp[i]+=x[i*c_dim+j];
    }

    for(int ic_dim=0;ic_dim<phi_psi_sz_sz;ic_dim++)
        Ax[ic_dim]=-x[ic_dim];
    
    
    int** neighbor_list_2nd=neighbor_dmd->neighbor_list_2nd;
    int* neighbor_list_size_2nd=neighbor_dmd->neighbor_list_size_2nd;
    
    istart=0;
    for(int ic_dim=0;ic_dim<phi_psi_sz_sz;ic_dim++)
        for(int j=0,jc_dim;j<neighbor_list_size_2nd[ic_dim];j++)
        {
            jc_dim=neighbor_list_2nd[ic_dim][j];
            
            tmp0=M_IJ[istart++]*(b0[jc_dim]-b0[ic_dim]);
            tmp0+=M_IJ[istart++]*(cv[jc_dim/c_dim]*x[ic_dim]-c[ic_dim]*x_tmp[jc_dim/c_dim]);
            tmp0-=M_IJ[istart++]*(cv[ic_dim/c_dim]*x[jc_dim]-c[jc_dim]*x_tmp[ic_dim/c_dim]);
            Ax[ic_dim]+=tmp0;
            
            if(jc_dim<phi_psi_sz_sz)
                Ax[jc_dim]-=tmp0;
        }
}
/*--------------------------------------------
 force calculation
 --------------------------------------------*/
void ForceFieldEAMDMD::
force_calc_static(bool st_clc)
{
    type0* xvec=atoms->x->begin();
    type0* fvec=f->begin();
    type0* dE=dE_ptr->begin();
    type0* c=atoms->c->begin();
    
    int i_comp,I_comp;
    int j_comp,J_comp;
    type0 dx0,dx1,dx2;
    type0 fpair,apair;
    
    if(st_clc)
        for(int i=1;i<7;i++)
            nrgy_strss_lcl[i]=0.0;
   
    int istart=0;
    
    for(int ic_dim=0;ic_dim<phi_psi_sz_sz;ic_dim++)
    {
        i_comp=x_dim*(ic_dim/c_dim);
        I_comp=3*(ic_dim/c_dim)+ic_dim;
    
        for(int j=0,jc_dim;j<phi_psi_sz[ic_dim];j++)
        {
            jc_dim=phi_psi_cmp[istart];
            fpair=-(psi_r_JI[istart]*dE[ic_dim]+psi_r_IJ[istart]*dE[jc_dim]+phi_r_IJ[istart])*c[ic_dim]*c[jc_dim];
            apair=-(psi_alpha_JI[istart]*dE[ic_dim]+psi_alpha_IJ[istart]*dE[jc_dim]+phi_alpha_IJ[istart])*c[ic_dim]*c[jc_dim];
            
            j_comp=x_dim*(jc_dim/c_dim);
            J_comp=3*(jc_dim/c_dim)+jc_dim;
            
            dx0=xvec[i_comp]-xvec[j_comp];
            dx1=xvec[i_comp+1]-xvec[j_comp+1];
            dx2=xvec[i_comp+2]-xvec[j_comp+2];
            
            fvec[i_comp]+=dx0*fpair;
            fvec[i_comp+1]+=dx1*fpair;
            fvec[i_comp+2]+=dx2*fpair;
            fvec[I_comp]+=apair*xvec[I_comp];
            
            if(jc_dim<phi_psi_sz_sz)
            {
                fvec[j_comp]-=dx0*fpair;
                fvec[j_comp+1]-=dx1*fpair;
                fvec[j_comp+2]-=dx2*fpair;
                fvec[J_comp]+=apair*xvec[J_comp];
            }
            
            if(jc_dim>=phi_psi_sz_sz)
                fpair*=0.5;
            
            if(st_clc)
            {
                nrgy_strss_lcl[1]-=fpair*dx0*dx0;
                nrgy_strss_lcl[2]-=fpair*dx1*dx1;
                nrgy_strss_lcl[3]-=fpair*dx2*dx2;
                nrgy_strss_lcl[4]-=fpair*dx1*dx2;
                nrgy_strss_lcl[5]-=fpair*dx2*dx0;
                nrgy_strss_lcl[6]-=fpair*dx0*dx1;
            }
            
            istart++;
        }
        
        if(c[ic_dim]>0.0)
            fvec[I_comp]+=3.0*kbT*c[ic_dim]/xvec[I_comp];
    }
    
    if(st_clc)
    {
        for(int i=0;i<7;i++)
            nrgy_strss[i]=0.0;
        
        MPI_Allreduce(nrgy_strss_lcl,nrgy_strss,7,Vec<type0>::MPI_T,MPI_SUM,world);
    }
    else
        MPI_Allreduce(nrgy_strss_lcl,nrgy_strss,1,Vec<type0>::MPI_T,MPI_SUM,world);
}

/*--------------------------------------------
 gamma_i = x_ij*x_ij/(alpha_i*alpha_i)
 gamma_j = x_ij*x_ij/(alpha_j*alpha_j)
 mu_ji = beta*(mu_j-mu_i)+3*log(alpha_j/alpha_i)
 
 Q=beta*U
 
 alpha_Q_sq=alpha_U_sq/(x_ij*x_ij)
 --------------------------------------------*/
inline void ForceFieldEAMDMD::
calc_Q(int& ielem,type0& gamma_i,type0& gamma_j,type0& mu_ji
,type0& Q,type0& alpha_Q_sq)
{
    type0 z0;
    type0 a=5.0*(gamma_i-gamma_j-6.0*mu_ji);
    type0 xi=-2.0*(gamma_i+gamma_j);
    type0 d=sqrt((xi-a)*(xi-a)-8.0*a*gamma_i);
    z0=-4.0*gamma_i/((xi-a)-d);
    
    /*
    if(z0<0.0 || z0>1.0)
    {
        Error::abort("could not find z0");
    }
    */
    
    Q=z0*z0*(1.0-z0)*(1.0-z0)*((1.0-z0)*gamma_i+z0*gamma_j);
    Q+=z0*z0*z0*(6.0*z0*z0-15.0*z0+10.0)*mu_ji;
    Q*=g_fac[ielem];
    Q+=0.5*(1.0-g_fac[ielem])*mu_ji;
    alpha_Q_sq=1.0/((1.0-z0)*gamma_i+z0*gamma_j);
    
}
/*--------------------------------------------
 gamma_i = x_ij*x_ij/(alpha_i*alpha_i)
 gamma_j = x_ij*x_ij/(alpha_j*alpha_j)
 mu_ji = beta*(mu_j-mu_i)+3*log(alpha_j/alpha_i)
 
 Q=beta*U
 dQ=dU/d((mu_j-mu_i)+3*log(alpha_j/alpha_i)/beta)
 
 alpha_Q_sq=alpha_U_sq/(x_ij*x_ij)
 dalpha_Q_sq=1/(beta*x_ij*x_ij)
 dalpha_Q_sq/d((mu_j-mu_i)+3*log(alpha_j/alpha_i)/beta)
 --------------------------------------------*/
inline void ForceFieldEAMDMD::
calc_Q(int& ielem,type0& gamma_i,type0& gamma_j,type0& mu_ji
,type0& Q,type0& alpha_Q_sq,type0& theta)
{
    type0 z0,dz0,dQ;
    type0 a=5.0*(gamma_i-gamma_j-6.0*mu_ji);
    type0 xi=-2.0*(gamma_i+gamma_j);
    type0 d=sqrt((xi-a)*(xi-a)-8.0*a*gamma_i);
    z0=-4.0*gamma_i/((xi-a)-d);
    dz0=30.0*z0*(1.0-z0)/d;
    
    /*
    if(z0<0.0 || z0>1.0)
    {
        Error::abort("could not find z0");
    }
    */
    dQ=z0*z0*z0*(6.0*z0*z0-15.0*z0+10.0);
    Q=z0*z0*(1.0-z0)*(1.0-z0)*((1.0-z0)*gamma_i+z0*gamma_j);
    Q+=dQ*mu_ji;
    Q*=g_fac[ielem];
    Q+=0.5*(1.0-g_fac[ielem])*mu_ji;
    dQ*=g_fac[ielem];
    dQ+=0.5*(1.0-g_fac[ielem]);
    
    alpha_Q_sq=1.0/((1.0-z0)*gamma_i+z0*gamma_j);
    theta=alpha_Q_sq*(gamma_i-gamma_j)*dz0-dQ;
}









