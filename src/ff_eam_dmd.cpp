#include "ff_eam_dmd.h"
#include <stdlib.h>
#include "neighbor_dmd.h"
#include "elements.h"
#include "memory.h"
#include "xmath.h"
#include "atoms_dmd.h"
#include "MAPP.h"
#include "dynamic_dmd.h"
#include "import_eam.h"
#include <limits>
#define PI_IN_SQ 0.564189583547756286948079451561
using namespace MAPP_NS;
/*--------------------------------------------
 constructor
 --------------------------------------------*/
ForceFieldEAMDMD::
ForceFieldEAMDMD(AtomsDMD* atoms,
type0 __dr,type0 __drho,size_t __nr,size_t __nrho,
type0(***&& __r_phi_arr)[4],type0(***&& __rho_arr)[4],type0(**&& __F_arr)[5],
type0**&& __cut,type0*&& __r_crd):
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
M_IJ(NULL),
vec0(NULL),
vec1(NULL),
vec2(NULL),
vec3(NULL),
mu_ptr(NULL),
dE_ptr(NULL),
ddE_ptr(NULL),
cv_ptr(NULL),
N(atoms->N)
{
    
    __r_phi_arr=NULL;
    __rho_arr=NULL;
    __F_arr=NULL;
    
    dr_inv=1.0/dr;
    drho_inv=1.0/drho;
    rho_max=static_cast<type0>(nrho)*drho;

    kbT=beta=-1.0;
    

    
    Memory::alloc(xi,N);
    Memory::alloc(wi_0,N);
    Memory::alloc(wi_1,N);
    Memory::alloc(c_0,nelems);
    Memory::alloc(c_1,nelems);
    
    memcpy(xi,atoms->xi,N*sizeof(type0));
    memcpy(wi_0,atoms->wi,N*sizeof(type0));
    for(int i=0;i<N;i++)
        wi_1[i]=wi_0[i]*xi[i];


    for(elem_type elem_i=0;elem_i<nelems;elem_i++)
    {
        r_crd[elem_i]=__r_crd[elem_i];
    }
    Memory::dealloc(__r_crd);
    
    for(elem_type elem_i=0;elem_i<nelems;elem_i++)
    {
        for(size_t i=0;i<nr;i++)
            r_rho_arr[elem_i][0][i][0]*=static_cast<type0>(i)*dr;
        ImportEAM::interpolate(nr,dr,r_rho_arr[elem_i][0]);
        for(elem_type elem_j=1;elem_j<nelems;elem_j++)
        {
            if(r_rho_arr[elem_i][elem_j]!=r_rho_arr[elem_i][elem_j-1])
            {
                for(size_t i=0;i<nr;i++)
                    r_rho_arr[elem_i][elem_j][i][0]*=static_cast<type0>(i)*dr;
                ImportEAM::interpolate(nr,dr,r_rho_arr[elem_i][elem_j]);
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
    Memory::dealloc(wi_1);
    Memory::dealloc(wi_0);
    Memory::dealloc(xi);

}
/*--------------------------------------------
 python constructor
 --------------------------------------------*/
void ForceFieldEAMDMD::ml_new(PyMethodDef& method_0,PyMethodDef& method_1,PyMethodDef& method_2)
{
    method_0.ml_flags=METH_VARARGS | METH_KEYWORDS;
    method_0.ml_name="ff_eam_funcfl";
    method_0.ml_meth=(PyCFunction)(PyCFunctionWithKeywords)
    [](PyObject* self,PyObject* args,PyObject* kwds)->PyObject*
    {
        AtomsDMD::Object* __self=reinterpret_cast<AtomsDMD::Object*>(self);
        size_t& nelems=__self->atoms->elements.nelems;
        
        FuncAPI<std::string*,type0*> f("ff_eam_funcfl",{"funcfl_files","r_crd"});
        f.v<0>().dynamic_size(nelems);
        f.v<1>().dynamic_size(nelems);
        f.logics<1>()[0]=VLogics("gt",0.0);
        if(f(args,kwds)) return NULL;
        
        size_t nr,nrho;
        type0 dr,drho;
        type0** r_c;
        type0(** F)[5]=NULL;
        type0(*** r_phi)[4]=NULL;
        type0(*** rho)[4]=NULL;
        try
        {
            ImportEAM::funcfl(nelems,f.val<0>(),dr,drho,nr,nrho,r_phi,rho,F,r_c);
        }
        catch(char* err_msg)
        {
            PyErr_SetString(PyExc_TypeError,err_msg);
            delete [] err_msg;
            return NULL;
        }
        
        for(size_t i=0;i<nelems;i++)
            if(f.v<1>()[i]>r_c[i][i])
            {
                Memory::dealloc(r_c);
                Memory::dealloc(F);
                Memory::dealloc(r_phi);
                Memory::dealloc(rho);
                PyErr_Format(PyExc_TypeError,"r_crd[%zu] should be less than r_c[%zu][%zu]",i,i,i);
                return NULL;
            }
        
        delete __self->ff;
        __self->ff=new ForceFieldEAMDMD(__self->atoms,dr,drho,nr,nrho,std::move(r_phi),std::move(rho),std::move(F),std::move(r_c),std::move(f.val<1>()));
        Py_RETURN_NONE;
    };
    method_0.ml_doc=(char*)R"---(
    ff_eam_funcfl(funcfl_files)
   
    Tabulated EAM force field given by FuncFL file/s
    
    Assigns EAM force field to system
    
    Parameters
    ----------
    funcfl_files : string[nelems]
        list of relative paths to DYNAMO files with FuncFL format
    
    Returns
    -------
    None
   
    Notes
    -----
    This is tabulated form of Embedded Atom Method (EAM) potential
    
    
    Examples
    --------
    Ni
    
    ::
     
        >>> from mapp import dmd
        >>> sim=dmd.cfg("configs/Ni-DMD.cfg")
        >>> sim.ff_eam_funcfl("potentials/niu3.eam")
    
    

    )---";
    
    method_1.ml_flags=METH_VARARGS | METH_KEYWORDS;
    method_1.ml_name="ff_eam_setfl";
    method_1.ml_meth=(PyCFunction)(PyCFunctionWithKeywords)
    [](PyObject* self,PyObject* args,PyObject* kwds)->PyObject*
    {
        AtomsDMD::Object* __self=reinterpret_cast<AtomsDMD::Object*>(self);
        size_t& nelems=__self->atoms->elements.nelems;
        FuncAPI<std::string,type0*> f("ff_eam_setfl",{"setfl_file","r_crd"});
        f.v<1>().dynamic_size(nelems);
        f.logics<1>()[0]=VLogics("gt",0.0);
        if(f(args,kwds)) return NULL;
        
        
        size_t nr,nrho;
        type0 dr,drho;
        type0** r_c;
        type0(** F)[5]=NULL;
        type0(*** r_phi)[4]=NULL;
        type0(*** rho)[4]=NULL;
        try
        {
            ImportEAM::setfl(nelems,__self->atoms->elements.names,f.val<0>(),dr,drho,nr,nrho,r_phi,rho,F,r_c);
        }
        catch(char* err_msg)
        {
            PyErr_SetString(PyExc_TypeError,err_msg);
            delete [] err_msg;
            return NULL;
        }
        
        for(size_t i=0;i<nelems;i++)
            if(f.v<1>()[i]>r_c[i][i])
            {
                Memory::dealloc(r_c);
                Memory::dealloc(F);
                Memory::dealloc(r_phi);
                Memory::dealloc(rho);
                PyErr_Format(PyExc_TypeError,"r_crd[%zu] should be less than r_c[%zu][%zu]",i,i,i);
                return NULL;
            }
        
        delete __self->ff;
        __self->ff=new ForceFieldEAMDMD(__self->atoms,dr,drho,nr,nrho,std::move(r_phi),std::move(rho),std::move(F),std::move(r_c),std::move(f.val<1>()));
        Py_RETURN_NONE;
    };
    method_1.ml_doc=(char*)R"---(
    ff_eam_setfl(setfl_file)
   
    Tabulated EAM force field given by a single SetFL file
    
    Assigns EAM force field to system
    
    Parameters
    ----------
    setfl_file : string
        relative path to DYNAMO file with SetFL format
    
    Returns
    -------
    None
   
    Notes
    -----
    This is tabulated form of Embedded Atom Method (EAM) potential
    
    
    Examples
    --------
    Cu
    
    ::
     
        >>> from mapp import dmd
        >>> sim=dmd.cfg("configs/Cu-DMD.cfg")
        >>> sim.ff_eam_setfl("potentials/Cu_mishin.eam.alloy")
    
    

    )---";
    
    method_2.ml_flags=METH_VARARGS | METH_KEYWORDS;
    method_2.ml_name="ff_eam_fs";
    method_2.ml_meth=(PyCFunction)(PyCFunctionWithKeywords)
    [](PyObject* self,PyObject* args,PyObject* kwds)->PyObject*
    {
        AtomsDMD::Object* __self=reinterpret_cast<AtomsDMD::Object*>(self);
        size_t& nelems=__self->atoms->elements.nelems;
        FuncAPI<std::string,type0*> f("ff_eam_fs",{"fs_file","r_crd"});
        f.v<1>().dynamic_size(nelems);
        f.logics<1>()[0]=VLogics("gt",0.0);
        if(f(args,kwds)) return NULL;
        
        size_t nr,nrho;
        type0 dr,drho;
        type0** r_c;
        type0(** F)[5]=NULL;
        type0(*** r_phi)[4]=NULL;
        type0(*** rho)[4]=NULL;
        try
        {
            ImportEAM::fs(nelems,__self->atoms->elements.names,f.val<0>(),dr,drho,nr,nrho,r_phi,rho,F,r_c);
        }
        catch(char* err_msg)
        {
            PyErr_SetString(PyExc_TypeError,err_msg);
            delete [] err_msg;
            return NULL;
        }
        
        for(size_t i=0;i<nelems;i++)
            if(f.v<1>()[i]>r_c[i][i])
            {
                Memory::dealloc(r_c);
                Memory::dealloc(F);
                Memory::dealloc(r_phi);
                Memory::dealloc(rho);
                PyErr_Format(PyExc_TypeError,"r_crd[%zu] should be less than r_c[%zu][%zu]",i,i,i);
                return NULL;
            }
        
        delete __self->ff;
        __self->ff=new ForceFieldEAMDMD(__self->atoms,dr,drho,nr,nrho,std::move(r_phi),std::move(rho),std::move(F),std::move(r_c),std::move(f.val<1>()));
        Py_RETURN_NONE;
    };
    method_2.ml_doc=(char*)R"---(
    ff_eam_fs(fs_file)
   
    Tabulated Finnis-Sinclair EAM
    
    Assigns Finnis-Sinclair EAM force field to system. For explanation of the parameter see the Notes section.
    
    Parameters
    ----------
    fs_file : string
        relative path to DYNAMO file with fs format
    
    Returns
    -------
    None
   
    Notes
    -----
    This is tabulated form of Finnis-Sinclair Embedded Atom Method (EAM) potential
    
    
    Examples
    --------
    Iron Hydrogrn mixture
    ::
     
        >>> from mapp import dmd
        >>> sim=dmd.cfg("configs/FeH-DMD.cfg")
        >>> sim.ff_eam_fs("potentials/FeH.eam.fs")
    
    

    )---";
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
    
    type0 p,cv_i;
    
    type0 dx_ij[__dim__];
    
    type0 const* c=atoms->c->begin();

    
    
    type0* dE=dE_ptr->begin();
    type0* mu=mu_ptr->begin();
    type0* rho=E_ptr->begin();
    const int n=atoms->natms_lcl*c_dim;
    for(int i=0;i<n;i++) rho[i]=0.0;
    
    elem_type const* elem_vec=atoms->elem->begin();
    
    type0 alpha_ij,rsq;
    elem_type elem_i,elem_j;
    
    
    
    const type0* x=atoms->x->begin();
    const type0* alpha=atoms->alpha->begin();
    int** neighbor_list=neighbor->neighbor_list;
    int* neighbor_list_size=neighbor->neighbor_list_size;
    size_t istart=0;
    for(int i=0;i<n;i++)
    {
        type0 c_i=c[i];
        if(c_i<0.0) continue;
        elem_i=elem_vec[i];
        const int neigh_sz=neighbor_list_size[i];
        for(int j,__j=0;__j<neigh_sz;__j++,istart+=3)
        {
            j=neighbor_list[i][__j];
            elem_j=elem_vec[j];
            rsq=Algebra::RSQ<__dim__>(x+(i/c_dim)*__dim__,x+(j/c_dim)*__dim__);
            r=sqrt(rsq);
            alpha_ij=sqrt(alpha[i]*alpha[i]+alpha[j]*alpha[j]);
            if(r-alpha_ij*xi[N-1]>=cut[elem_i][elem_j]) continue;
            type0 upper=(r+cut[elem_i][elem_j])/alpha_ij;
            type0 lower=(r-cut[elem_i][elem_j])/alpha_ij;
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
                
                coef=r_phi_arr[elem_i][elem_j][m];
                tmp0=((coef[3]*p+coef[2])*p+coef[1])*p+coef[0];
                tmp1=((3.0*coef[3]*p+2.0*coef[2])*p+coef[1])*dr_inv;
                if(__r<0.0) tmp0*=-1.0;
                
                __rho_phi[0]+=wi_0[l]*tmp0;
                __drho_phi_dr[0]+=wi_0[l]*tmp1;
                __drho_phi_dalpha[0]+=wi_1[l]*tmp1;
                
                coef=r_rho_arr[elem_i][elem_j][m];
                tmp0=((coef[3]*p+coef[2])*p+coef[1])*p+coef[0];
                tmp1=((3.0*coef[3]*p+2.0*coef[2])*p+coef[1])*dr_inv;
                if(__r<0.0) tmp0*=-1.0;
                
                __rho_phi[1]+=wi_0[l]*tmp0;
                __drho_phi_dr[1]+=wi_0[l]*tmp1;
                __drho_phi_dalpha[1]+=wi_1[l]*tmp1;
                
                coef=r_rho_arr[elem_j][elem_i][m];
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
      
            
            rho[i]+=c[j]*__rho_phi[2];
            
            if(j<n)
            {
                rho[j]+=c_i*__rho_phi[1];
                nrgy_strss_lcl[0]+=c_i*c[j]*__rho_phi[0];
            }
            else
                nrgy_strss_lcl[0]+=0.5*c_i*c[j]*__rho_phi[0];
        }
        
        
        
        p=rho[i]*drho_inv;
        m=static_cast<size_t>(p);
        m=MIN(m,nrho-2);
        p-=m;
        p=MIN(p,1.0);
        coef=F_arr[elem_i][m];
        
        tmp0=(((coef[4]*p+coef[3])*p+coef[2])*p+coef[1])*p+coef[0];
        tmp1=(((4.0*coef[4]*p+3.0*coef[3])*p+2.0*coef[2])*p+coef[1])*drho_inv;
        
        if(rho[i]>rho_max) tmp0+=tmp1*(rho[i]-rho_max);
        
        rho[i]=tmp0;
        dE[i]=tmp1;
        mu[i]=tmp0;
        if(c_i!=0.0)
            nrgy_strss_lcl[0]+=c_i*(tmp0+c_0[elem_i]-3.0*kbT*log(alpha[i]));
        
        nrgy_strss_lcl[0]+=kbT*calc_ent(c_i);
    }
    
    const int __natms=atoms->natms_lcl;
    
    for(int i=0;i<__natms;i++)
    {
        cv_i=1.0;
        for(int ic=0;ic<c_dim;ic++)
            if(c[i*c_dim+ic]>0.0)
                cv_i-=c[i*c_dim+ic];
        nrgy_strss_lcl[0]+=kbT*calc_ent(cv_i);
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
        
        type0 c_i=c[i];
        type0 alpha_i=alpha[i];
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
            f_alpha_i+=alpha_i*apair;
            
            if(j<n)
            {
                Algebra::V_add_x_mul_V<__dim__>(-fpair,dx_ij,fvec+(j/c_dim)*__dim__);
                f_alphavec[j]+=alpha[j]*apair;
            }
            else
                fpair*=0.5;
            Algebra::DyadicV(-fpair,dx_ij,&nrgy_strss_lcl[1]);
        }
        
        f_alpha_i+=3.0*kbT*c_i/alpha_i;
        f_alphavec[i]+=f_alpha_i;
        mu[i]+=mu_i;
        
        if((i+1)%c_dim==0)
            Algebra::V_add<__dim__>(f_i,fvec+__dim__*(i/c_dim));
    }    
}
/*--------------------------------------------
 energy calculation
 --------------------------------------------*/
void ForceFieldEAMDMD::energy_calc()
{
    type0* rho=ddE_ptr->begin();
    const int n=atoms->natms_lcl*c_dim;
    for(int i=0;i<n;i++) rho[i]=0.0;
    
    elem_type const* elem_vec=atoms->elem->begin();
    
    type0 alpha_ij,rsq,r,r_inv,p,tmp0,tmp1;
    elem_type elem_i,elem_j;
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
        elem_i=elem_vec[i];
        const int neigh_sz=neighbor_list_size[i];
        for(int j,__j=0;__j<neigh_sz;__j++,istart+=3)
        {
            j=neighbor_list[i][__j];
            elem_j=elem_vec[j];
            rsq=Algebra::RSQ<__dim__>(x+(i/c_dim)*__dim__,x+(j/c_dim)*__dim__);
            r=sqrt(rsq);
            alpha_ij=sqrt(alpha[i]*alpha[i]+alpha[j]*alpha[j]);
            if(r-alpha_ij*xi[N-1]>=cut[elem_i][elem_j]) continue;
            
            type0 upper=(r+cut[elem_i][elem_j])/alpha_ij;
            type0 lower=(r-cut[elem_i][elem_j])/alpha_ij;
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
                
                coef=r_phi_arr[elem_i][elem_j][m];
                tmp0=((coef[3]*p+coef[2])*p+coef[1])*p+coef[0];
                tmp1=((3.0*coef[3]*p+2.0*coef[2])*p+coef[1])*dr_inv;
                if(__r<0.0) tmp0*=-1.0;
                __arr[0]+=wi_0[l]*tmp0;

                
                coef=r_rho_arr[elem_i][elem_j][m];
                tmp0=((coef[3]*p+coef[2])*p+coef[1])*p+coef[0];
                tmp1=((3.0*coef[3]*p+2.0*coef[2])*p+coef[1])*dr_inv;
                if(__r<0.0) tmp0*=-1.0;
                __arr[1]+=wi_0[l]*tmp0;
                
                coef=r_rho_arr[elem_j][elem_i][m];
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
            
            rho[i]+=c[j]*__arr[2];
            
            if(j<n)
            {
                rho[j]+=c[i]*__arr[1];
                nrgy_strss_lcl[0]+=c[i]*c[j]*__arr[0];
            }
            else
                nrgy_strss_lcl[0]+=0.5*c[i]*c[j]*__arr[0];
        }
        
        
        if(c[i]<0.0) continue;
        tmp0=rho[i];
        p=tmp0*drho_inv;
        m=static_cast<size_t>(p);
        m=MIN(m,nrho-2);
        p-=m;
        p=MIN(p,1.0);
        coef=F_arr[elem_i][m];
        tmp1=(((coef[4]*p+coef[3])*p+coef[2])*p+coef[1])*p+coef[0];
        if(rho[i]>rho_max)
            tmp1+=(((4.0*coef[4]*p+3.0*coef[3])*p+2.0*coef[2])*p+coef[1])*drho_inv*(tmp0-rho_max);
        
        if(c[i]!=0.0)
            nrgy_strss_lcl[0]+=c[i]*(tmp1+c_0[elem_i]-3.0*kbT*log(alpha[i]));
        nrgy_strss_lcl[0]+=kbT*calc_ent(c[i]);
    }
    const int __natms=atoms->natms_lcl;
    type0 cv_i;
    for(int i=0;i<__natms;i++)
    {
        cv_i=1.0;
        for(int ic=0;ic<c_dim;ic++)
            if(c[i*c_dim+ic]>0.0)
                cv_i-=c[i*c_dim+ic];
        nrgy_strss_lcl[0]+=kbT*calc_ent(cv_i);
    }
}
/*--------------------------------------------
 init
 --------------------------------------------*/
void ForceFieldEAMDMD::init()
{
    pre_init();
    set_temp(300.0);
    
    
    mu_ptr=new DMDVec<type0>(atoms,0.0,"mu");
    cv_ptr=new Vec<type0>(atoms,1);
    E_ptr=new Vec<type0>(atoms,c_dim);
    dE_ptr=new Vec<type0>(atoms,c_dim);
    ddE_ptr=new Vec<type0>(atoms,c_dim);
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
    
    delete ddE_ptr;
    delete dE_ptr;
    delete E_ptr;
    delete cv_ptr;
    delete mu_ptr;
    
    dE_ptr=ddE_ptr=cv_ptr=vec0=vec1=vec2=vec3=NULL;
    mu_ptr=NULL;
    post_fin();
}
/*--------------------------------------------
 create the sparse matrices
 --------------------------------------------*/
void ForceFieldEAMDMD::init_static()
{
    Memory::alloc(M_IJ,3*neighbor->no_pairs_2nd);
    vec0=new Vec<type0>(atoms,c_dim);
    vec1=new Vec<type0>(atoms,1);
    vec2=new Vec<type0>(atoms,c_dim);
    vec3=new Vec<type0>(atoms,c_dim);
}
/*--------------------------------------------
 create the sparse matrices
 --------------------------------------------*/
void ForceFieldEAMDMD::fin_static()
{
    
    delete vec3;
    delete vec2;
    delete vec1;
    delete vec0;
    Memory::dealloc(M_IJ);
}
/*--------------------------------------------
 set the temperature in the simulation
 --------------------------------------------*/
void ForceFieldEAMDMD::set_temp(type0 T)
{
    type0 kb=atoms->kB;
    type0 h=atoms->h;
    type0 mass;
    type0 deb_l;
    
    for(size_t i=0;i<nelems;i++)
    {
        mass=atoms->elements.masses[i];
        c_1[i]=sqrt(0.5*kb*T/mass)/M_PI;
        deb_l=h*h/(2.0*M_PI*M_PI*mass*kb*T);
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
    if(x==0.0) return 0.0;
    return x*log(x);
}
/*--------------------------------------------
 external vecs
 atoms->x;
 atoms->alpha;
 atoms->c;
 atoms->elem;
 
 internal vecs
 mu_ptr;
 cv_ptr;
 --------------------------------------------*/
void ForceFieldEAMDMD::dc()
{
    /*
    calc_mu();
    
    //external vecs
    const type0* x=atoms->x->begin();
    const type0* alpha=atoms->alpha->begin();
    const type0* c=atoms->c->begin();
    const elem_type* elem_vec=atoms->elem->begin();
    type0* c_d=atoms->c_d->begin();
    
    
    //internal vecs
    type0* cv=cv_ptr->begin();
    type0* mu=mu_ptr->begin();
    
    const int n=atoms->natms_lcl*c_dim;
    for(int i=0;i<n;i++) c_d[i]=0.0;
    
    int** neighbor_list=neighbor->neighbor_list_2nd;
    int* neighbor_list_size=neighbor->neighbor_list_size_2nd;
    elem_type elem_i;
    type0 rsq,alpha_i,alpha_j,alpha_i_sq,alpha_j_sq,gamma_i,gamma_j,mu_i,mu_ji;
    type0 Qi,alpha_Q_sq,d_i,d_j;
    type0 dc_ij;
    size_t istart=0;
    for(int i=0;i<n;i++)
    {
        if(c[i]<0.0) continue;
        elem_i=elem_vec[i];
        alpha_i=alpha[i];
        alpha_i_sq=alpha_i*alpha_i;
        mu_i=mu[i];
        const int neigh_sz=neighbor_list_size[i];
        for(int j,__j=0;__j<neigh_sz;__j++,istart+=3)
        {
            j=neighbor_list[i][__j];
            alpha_j=alpha[j];
            alpha_j_sq=alpha_j*alpha_j;
            rsq=Algebra::RSQ<__dim__>(x+(i/c_dim)*__dim__,x+(j/c_dim)*__dim__);
            gamma_i=rsq/alpha_i_sq;
            gamma_j=rsq/alpha_j_sq;
            mu_ji=beta*(mu[j]-mu_i);
            calc_Q(gamma_i,gamma_j,mu_ji,Qi,alpha_Q_sq);
            
            d_i=d_j=-c_1[elem_i]*rsq*alpha_Q_sq;
            d_i/=alpha_i_sq*alpha_i;
            d_j/=alpha_j_sq*alpha_j;
            
            dc_ij=d_i*exp(-Qi)*c[i]*cv[j/c_dim]-d_j*exp(-Qi+mu_ji)*c[j]*cv[i/c_dim];
            
            
            c_d[i]+=dc_ij;
            if(j<n)
                c_d[j]-=dc_ij;
        }
    }
     */
}
/*--------------------------------------------
 calculate norm of d^2c/dt^2
 
 external vecs
 atoms->c_d;
 atoms->c;
 
 internal vecs
 s_ptr;
 t_ptr;
 --------------------------------------------*/
type0 ForceFieldEAMDMD::ddc_norm()
{
    const int n=atoms->natms_lcl*c_dim;
    type0* tmp=NULL;
    Memory::alloc(tmp,n);
    for(int i=0;i<n;i++) tmp[i]=0;
    update_J(1.0,tmp,tmp);
    delete [] tmp;
    
    //external vec
    type0* c_dvec=atoms->c_d->begin();
    
    //internal vec
    type0* s=vec2->begin();
    
    memcpy(s,c_dvec,sizeof(type0)*n);
    
    operator()(vec2,vec3);
    
    //internal vec
    type0* t=vec3->begin();
    
    //external vec
    type0* c=atoms->c->begin();
    type0 ans_lcl=0.0,ans;
    for(int i=0;i<n;i++)
        if(c[i]>=0.0)
            ans_lcl+=(t[i]-s[i])*(t[i]-s[i]);
    
    MPI_Allreduce(&ans_lcl,&ans,1,Vec<type0>::MPI_T,MPI_SUM,world);
    return sqrt(ans);
}
/*--------------------------------------------
 calculate cv
 calculate mu
 
 external vecs
 atoms->alpha;
 atoms->c;
 atoms->elem;

 internal vecs
 cv_ptr;
 E_ptr;
 dE_ptr;
 mu_ptr;
 --------------------------------------------*/
void ForceFieldEAMDMD::calc_mu()
{
    /*--------------------------------
     2 things are hapenning here:
     0. calculate c_v for local atoms
     1. calculate the non interacting 
     part of energy
     --------------------------------*/
    //external vecs
    const type0* alpha=atoms->alpha->begin();
    type0 const* c=atoms->c->begin();
    elem_type const* elem_vec=atoms->elem->begin();
    
    //internal vecs
    type0* cv=cv_ptr->begin();
    
    int natms_lcl=atoms->natms_lcl;
    type0 en=0.0;
    for(int i=0;i<natms_lcl;i++)
    {
        type0 cv_i=1.0;
        for(int j=0;j<c_dim;j++)
            if(c[i*c_dim+j]>0.0)
            {
                en+=kbT*calc_ent(c[i*c_dim+j])
                +c[i*c_dim+j]*
                (c_0[elem_vec[i*c_dim+j]]-3.0*kbT*log(alpha[i*c_dim+j]));
                cv_i-=c[i*c_dim+j];
            }
        cv[i]=cv_i;
        en+=kbT*calc_ent(cv_i);
    }
    
    nrgy_strss_lcl[0]+=en;
    
    
    /*--------------------------------
     we will calculate the rest of c_v
     here
     --------------------------------*/
    int nall=natms_lcl+atoms->natms_ph;
    for(int i=natms_lcl;i<nall;i++)
    {
        type0 cv_i=1.0;
        for(int j=0;j<c_dim;j++)
            if(c[i*c_dim+j]>0.0)
                cv_i-=c[i*c_dim+j];
        cv[i]=cv_i;
    }
    
    
    
    /*--------------------------------
     calculate the electron densities
     --------------------------------*/
    //internal vecs
    type0* rho=E_ptr->begin();
    type0* mu=mu_ptr->begin();
    
    int n=atoms->natms_lcl*c_dim;
    for(int i=0;i<n;i++) mu[i]=rho[i]=0.0;
    
    size_t m;
    type0 p,__E,__dE,__rho;
    type0* coef;
    
    //internal vecs
    type0* dE=dE_ptr->begin();
    type0* ddE=ddE_ptr->begin();
    int** neighbor_list=neighbor->neighbor_list;
    int* neighbor_list_size=neighbor->neighbor_list_size;
    size_t istart=0;
    for(int i=0;i<n;i++)
    {
        if(c[i]<0.0) continue;
        const int neigh_sz=neighbor_list_size[i];
        for(int j,__j=0;__j<neigh_sz;__j++,istart+=3)
        {
            j=neighbor_list[i][__j];
            rho[i]+=c[j]*rho_phi[istart+2];

            if(j<n)
                rho[j]+=c[i]*rho_phi[istart+1];
        }
        
        __rho=rho[i];
        p=__rho*drho_inv;
        m=static_cast<size_t>(p);
        m=MIN(m,nrho-2);
        p-=m;
        p=MIN(p,1.0);
        coef=F_arr[elem_vec[i]][m];
        
        __E=(((coef[4]*p+coef[3])*p+coef[2])*p+coef[1])*p+coef[0];
        __dE=(((4.0*coef[4]*p+3.0*coef[3])*p+2.0*coef[2])*p+coef[1])*drho_inv;
        ddE[i]=(((12.0*coef[4]*p+6.0*coef[3])*p+2.0*coef[2]))*drho_inv*drho_inv;
        if(__rho>rho_max)
            __E+=__dE*(__rho-rho_max);
        
        nrgy_strss_lcl[0]+=c[i]*__E;
        rho[i]=__E;
        mu[i]+=__E;
        dE[i]=__dE;
    }
    
    dynamic->update(dE_ptr);
    
    
    istart=0;
    for(int i=0;i<n;i++)
    {
        const int neigh_sz=neighbor_list_size[i];
        for(int j,__j=0;__j<neigh_sz;__j++,istart+=3)
        {
            j=neighbor_list[i][__j];
            mu[i]+=c[j]*(rho_phi[istart]+rho_phi[istart+1]*dE[j]);
            if(j<n)
            {
                mu[j]+=c[i]*(rho_phi[istart]+rho_phi[istart+2]*dE[i]);
                nrgy_strss_lcl[0]+=c[i]*c[j]*rho_phi[istart];
            }
            else
                nrgy_strss_lcl[0]+=0.5*c[i]*c[j]*rho_phi[istart];
        }
    }
    dynamic->update(mu_ptr);
}
/*--------------------------------------------
 force calculation
 --------------------------------------------*/
void ForceFieldEAMDMD::force_calc_static()
{
    int n=atoms->natms_lcl*c_dim;
    type0 const* E=E_ptr->begin();
    type0 const* c=atoms->c->begin();
    type0 const* cv=cv_ptr->begin();
    type0 const* alpha=atoms->alpha->begin();
    elem_type const* elem=atoms->elem->begin();
    for(int i=0;i<n;i++)
    {
        if(i%c_dim==0)
            nrgy_strss_lcl[0]+=kbT*calc_ent(cv[i]);
        if(c[i]>=0.0)
            nrgy_strss_lcl[0]+=c[i]*(E[i]+c_0[elem[i]]-3.0*kbT*log(alpha[i]))+kbT*calc_ent(c[i]);
    }
    
    
    
    type0 const* dE=dE_ptr->begin();
    type0 const* x=atoms->x->begin();
    

    
    type0 x_i[__dim__];
    type0 dx_ij[__dim__];
    
    type0 fpair;
    int** neighbor_list=neighbor->neighbor_list;
    int* neighbor_list_size=neighbor->neighbor_list_size;
    size_t istart=0;
    for(int i=0;i<n;i++)
    {
        if(i%c_dim==0) Algebra::V_eq<__dim__>(x+(i/c_dim)*__dim__,x_i);
        type0 c_i=c[i];
        type0 dE_i=dE[i];
        const int neigh_sz=neighbor_list_size[i];
        
        for(int j,__j=0;__j<neigh_sz;__j++,istart+=3)
        {
            j=neighbor_list[i][__j];
            Algebra::DX<__dim__>(x_i,x+(j/c_dim)*__dim__,dx_ij);
            fpair=-(drho_phi_dr[istart+2]*dE_i+drho_phi_dr[istart+1]*dE[j]+drho_phi_dr[istart])*c_i*c[j];
            if(j<n)
                nrgy_strss_lcl[0]+=rho_phi[istart]*c_i*c[j];
            else
            {
                fpair*=0.5;
                nrgy_strss_lcl[0]+=0.5*rho_phi[istart]*c_i*c[j];
            }
            
            Algebra::DyadicV(-fpair,dx_ij,&nrgy_strss_lcl[1]);
        }
    }
}
/*--------------------------------------------
 create the sparse matrices
 
 external vecs
 atoms->x;
 atoms->alpha;
 atoms->c;
 atoms->elem;
 atoms->c_d;
 
 
 internal vecs
 type0* mu=mu_ptr->begin();
 type0* cv=cv_ptr->begin();
 
 F=c-scalar*c_d(c)+a
 --------------------------------------------*/
type0 ForceFieldEAMDMD::update_J(type0 scalar,type0* a,type0* F)
{
    type0 iota=log(scalar);
    
    calc_mu();
    
    //external vecs
    const type0* x=atoms->x->begin();
    const type0* alpha=atoms->alpha->begin();
    const type0* c=atoms->c->begin();
    const elem_type* elem_vec=atoms->elem->begin();
    type0* c_d=atoms->c_d->begin();
    
    
    //internal vecs
    type0* cv=cv_ptr->begin();
    type0* mu=mu_ptr->begin();
    
    const int n=atoms->natms_lcl*c_dim;
    for(int i=0;i<n;i++) F[i]=c_d[i]=0.0;
    
    
    type0 rsq,alpha_i,alpha_j,alpha_sq_i,alpha_sq_j,gamma_i,gamma_j,mu_i,mu_ji;
    type0 Q_i,gamma_ij_inv,theta_i,theta_j,exp_mod_Q_i,exp_mod_Q_j,d_i,d_j;
    type0 dc_ij,dg_ij,c_1_ielem;
    type0 ans_lcl=0.0;
    size_t istart=0;
    int** neighbor_list=neighbor->neighbor_list_2nd;
    int* neighbor_list_size=neighbor->neighbor_list_size_2nd;
    for(int i=0;i<n;i++)
    {
        if(c[i]<0.0) continue;
        alpha_i=alpha[i];
        alpha_sq_i=alpha_i*alpha_i;
        mu_i=mu[i];
        c_1_ielem=c_1[elem_vec[i]];
        const int neigh_sz=neighbor_list_size[i];
        for(int j,__j=0;__j<neigh_sz;__j++,istart+=3)
        {
            j=neighbor_list[i][__j];
            alpha_j=alpha[j];
            alpha_sq_j=alpha_j*alpha_j;
            rsq=Algebra::RSQ<__dim__>(x+(i/c_dim)*__dim__,x+(j/c_dim)*__dim__);
            gamma_i=rsq/alpha_sq_i;
            gamma_j=rsq/alpha_sq_j;
            mu_ji=beta*(mu[j]-mu_i);
            calc_Q(gamma_i,gamma_j,mu_ji,Q_i,gamma_ij_inv,theta_i);
            
            exp_mod_Q_i=exp(-Q_i+iota);
            exp_mod_Q_j=exp(-Q_i+mu_ji+iota);

            d_i=d_j=-c_1_ielem*rsq*gamma_ij_inv;
            d_i/=alpha_sq_i*alpha_i;
            d_j/=alpha_sq_j*alpha_j;
            theta_j=theta_i+1.0;
            
            dc_ij=exp(-Q_i)*(d_i*c[i]*cv[j/c_dim]-d_j*exp(mu_ji)*c[j]*cv[i/c_dim]);
            dg_ij=exp(-Q_i+iota)*(d_j*exp(mu_ji)*c[j]*cv[i/c_dim]-d_i*c[i]*cv[j/c_dim]);
            
            M_IJ[istart]=beta*exp(-Q_i+iota)*(d_i*theta_i*c[i]*cv[j/c_dim]-d_j*theta_j*exp(mu_ji)*c[j]*cv[i/c_dim]);
            M_IJ[istart+1]=d_i*exp_mod_Q_i;
            M_IJ[istart+2]=d_j*exp_mod_Q_j;
            
            F[i]+=dg_ij;
            c_d[i]+=dc_ij;

            if(j<n)
            {
                F[j]-=dg_ij;
                c_d[j]-=dc_ij;
            }
        }
        F[i]+=c[i]+a[i];
        ans_lcl+=F[i]*F[i];
    }
    
    type0 ans;
    MPI_Allreduce(&ans_lcl,&ans,1,Vec<type0>::MPI_T,MPI_SUM,world);
    return sqrt(ans);
}
/*--------------------------------------------
 create the sparse matrices
 
 external vecs
 atoms->c;
 
 internal vecs
 E_ptr;
 dE_ptr;
 vec0;
 vec1;
 cv_ptr;
 
 input vecs
 x_ptr;
 Ax_ptr;
 
 F_i=c_i-beta*(d c_i)/dt
 J_ij=d F_i/(d c_j)
 Ax_i=-J_ij*x_j;
 --------------------------------------------*/
void ForceFieldEAMDMD::operator()(Vec<type0>* x_ptr,Vec<type0>* Ax_ptr)
{
    dynamic->update(x_ptr);
    
    //external vecs
    type0* c=atoms->c->begin();
    
    //internal vecs
    type0* ddE=ddE_ptr->begin();
    type0* dE=dE_ptr->begin();
    type0* b0=vec0->begin();
    
    //input vecs
    type0* x=x_ptr->begin();
    type0* Ax=Ax_ptr->begin();
    
    int** neighbor_list=neighbor->neighbor_list;
    int* neighbor_list_size=neighbor->neighbor_list_size;
    size_t istart=0;
    const int n=atoms->natms_lcl*c_dim;
    for(int i=0;i<n;i++) b0[i]=Ax[i]=0.0;
    for(int i=0;i<n;i++)
    {
        const int neigh_sz=neighbor_list_size[i];
        type0 ci_dEEi=c[i]*ddE[i];
        for(int j,__j=0;__j<neigh_sz;__j++,istart+=3)
        {
            j=neighbor_list[i][__j];
            Ax[i]+=ci_dEEi*rho_phi[istart+2]*x[j];
            if(j<n)
                Ax[j]+=c[j]*ddE[j]*rho_phi[istart+1]*x[i];
        }
    }
    dynamic->update(Ax_ptr);
    
    type0 tmp0;
    istart=0;
    for(int i=0;i<n;i++)
    {
        const int neigh_sz=neighbor_list_size[i];
        for(int j,__j=0;__j<neigh_sz;__j++,istart+=3)
        {
            j=neighbor_list[i][__j];
            tmp0=rho_phi[istart]+dE[i]*rho_phi[istart+2]+dE[j]*rho_phi[istart+1];
            b0[i]+=rho_phi[istart+1]*Ax[j]+tmp0*x[j];
            if(j<n)
                b0[j]+=rho_phi[istart+2]*Ax[i]+tmp0*x[i];
        }
    }
    
    dynamic->update(vec0);
    
    //internal vecs
    type0* xv=vec1->begin();
    type0* cv=cv_ptr->begin();
    
    const int nall=atoms->natms_lcl+atoms->natms_ph;
    for(int i=0;i<nall;i++)
    {
        xv[i]=0.0;
        for(int j=0;j<c_dim;j++)
            if(c[i*c_dim+j]>=0.0)
                xv[i]+=x[i*c_dim+j];
    }

    for(int i=0;i<n;i++)
        Ax[i]=x[i];
    
    
    neighbor_list=neighbor->neighbor_list_2nd;
    neighbor_list_size=neighbor->neighbor_list_size_2nd;
    istart=0;
    for(int i=0;i<n;i++)
    {
        const int neigh_sz=neighbor_list_size[i];
        for(int j,__j=0;__j<neigh_sz;__j++,istart+=3)
        {
            j=neighbor_list[i][__j];
            
            tmp0=M_IJ[istart]*(b0[j]-b0[i]);
            tmp0+=M_IJ[istart+1]*(cv[j/c_dim]*x[i]-c[i]*xv[j/c_dim]);
            tmp0-=M_IJ[istart+2]*(cv[i/c_dim]*x[j]-c[j]*xv[i/c_dim]);
            Ax[i]-=tmp0;
            
            if(j<n) Ax[j]+=tmp0;
        }
    }
    
}
/*--------------------------------------------
 gamma_i = x_ij*x_ij/(alpha_i*alpha_i)
 gamma_j = x_ij*x_ij/(alpha_j*alpha_j)
 mu_ji = beta*(mu_j-mu_i)
 
 Q=beta*U
 
 alpha_Q_sq=alpha_U_sq/(x_ij*x_ij)
 --------------------------------------------*/
inline void ForceFieldEAMDMD::
calc_Q(type0& gamma_i,type0& gamma_j,type0& mu_ji
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
    alpha_Q_sq=1.0/((1.0-z0)*gamma_i+z0*gamma_j);
}
/*--------------------------------------------
 gamma_i = x_ij*x_ij/(alpha_i*alpha_i)
 gamma_j = x_ij*x_ij/(alpha_j*alpha_j)
 mu_ji = beta*(mu_j-mu_i)
 
 Q=beta*U
 
 dz0: derivative of z0 w.r.t. mu_ji
 dQ:  derivative of Q  w.r.t. mu_ji
 
 alpha_Q_sq=alpha_U_sq/(x_ij*x_ij)
 dalpha_Q_sq/dmu_ji
 --------------------------------------------*/
inline void ForceFieldEAMDMD::
calc_Q(type0& gamma_i,type0& gamma_j,type0& mu_ji
,type0& Q,type0& alpha_Q_sq,type0& theta)
{
    type0 z0,dz0,dQ;
    type0 a=5.0*(gamma_i-gamma_j-6.0*mu_ji);
    type0 xi=-2.0*(gamma_i+gamma_j);
    type0 delta_sqrt=sqrt((xi-a)*(xi-a)-8.0*a*gamma_i);
    z0=-4.0*gamma_i/((xi-a)-delta_sqrt);
    dz0=30.0*z0*(1.0-z0)/delta_sqrt;
    /*
     here 
    
     
     */
    dQ=z0*z0*z0*(6.0*z0*z0-15.0*z0+10.0);
    Q=z0*z0*(1.0-z0)*(1.0-z0)*((1.0-z0)*gamma_i+z0*gamma_j);
    Q+=dQ*mu_ji;
    
    alpha_Q_sq=1.0/((1.0-z0)*gamma_i+z0*gamma_j);
    theta=alpha_Q_sq*(gamma_i-gamma_j)*dz0-dQ;

}
