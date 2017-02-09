#include "ff_fs.h"
#include "neighbor_md.h"
#include "atoms_md.h"
#include "elements.h"
#include "memory.h"
#include "dynamic_md.h"
using namespace MAPP_NS;
/*--------------------------------------------
 Finnis-Sinclair (FS) potential
 ref:
 T. T. Lau, C. J. Forst, X. Lin, J. D. Gale,
 S. Yip, & K. J. Van Vliet
 Many-Body Potential for Point Defect Clusters
 in Fe-C Alloys
 Phys. Rev. Lett. Vol. 98, pp. 215501, 2007
 --------------------------------------------*/
/*--------------------------------------------
 constructor
 --------------------------------------------*/
ForceFieldFS::ForceFieldFS(AtomsMD*& __atoms,type0*&& __A,
type0**&& __t1,type0**&& __t2,type0**&& __k1,type0**&& __k2,
type0**&& __k3,type0**&& r_c_phi,type0**&& r_c_rho):
ForceFieldMD(__atoms),
A(__A),
t1(__t1),
t2(__t2),
k1(__k1),
k2(__k2),
k3(__k3),
cut_sq_phi(r_c_phi),
cut_sq_rho(r_c_rho)
{
    __A=NULL;
    __t1=NULL;
    __t2=NULL;
    __k1=NULL;
    __k2=NULL;
    __k3=NULL;
    r_c_phi=NULL;
    r_c_rho=NULL;
    
    type0 tmp;
    for(size_t i=0;i<nelems;i++)
        for(size_t j=0;j<i+1;j++)
        {
            cut[i][j]=cut[j][i]=MAX(cut_sq_phi[i][j],cut_sq_rho[i][j]);
            cut_sq[i][j]=cut_sq[j][i]=cut[i][j]*cut[i][j];
            tmp=cut_sq_phi[i][j];
            cut_sq_phi[i][j]=cut_sq_phi[j][i]=tmp*tmp;
            tmp=cut_sq_rho[i][j];
            cut_sq_rho[i][j]=cut_sq_rho[j][i]=tmp*tmp;
        }
}
/*--------------------------------------------
 destructor
 --------------------------------------------*/
ForceFieldFS::~ForceFieldFS()
{
    Memory::dealloc(k3);
    Memory::dealloc(k2);
    Memory::dealloc(k1);
    Memory::dealloc(t2);
    Memory::dealloc(t1);
    Memory::dealloc(A);
    Memory::dealloc(cut_sq_phi);
    Memory::dealloc(cut_sq_rho);
}
/*--------------------------------------------
 python constructor
 --------------------------------------------*/
void ForceFieldFS::ml_new(PyMethodDef& tp_methods)
{
    tp_methods.ml_flags=METH_VARARGS | METH_KEYWORDS;
    tp_methods.ml_name="ff_fs";
    tp_methods.ml_doc="I will add doc here";
    
    tp_methods.ml_meth=(PyCFunction)(PyCFunctionWithKeywords)
    [](PyObject* self,PyObject* args,PyObject* kwds)->PyObject*
    {
        AtomsMD::Object* __self=reinterpret_cast<AtomsMD::Object*>(self);
        size_t& nelems=__self->atoms->elements->nelems;
        
        FuncAPI<type0*,type0**,type0**,symm<type0**>,symm<type0**>,symm<type0**>,symm<type0**>,symm<type0**>>
        f("ff_lj",{"A","t1","t2","k1","k2","k3","r_c_phi","r_c_rho"});

        f.var<0>().dynamic_size(nelems);
        f.logics<0>()[0]=VLogics("gt",0.0);
        f.var<1>().dynamic_size(nelems,nelems);
        f.var<2>().dynamic_size(nelems,nelems);
        f.var<3>().dynamic_size(nelems,nelems);
        f.var<4>().dynamic_size(nelems,nelems);
        f.var<5>().dynamic_size(nelems,nelems);
        f.var<6>().dynamic_size(nelems,nelems);
        f.logics<6>()[0]=VLogics("gt",0.0);
        f.var<7>().dynamic_size(nelems,nelems);
        f.logics<7>()[0]=VLogics("gt",0.0);
        
        if(f(args,kwds)) return NULL;
        
        delete __self->ff;
        __self->ff=new ForceFieldFS(__self->atoms,f.mov<0>(),f.mov<1>(),f.mov<2>(),
        f.mov<3>(),f.mov<4>(),f.mov<5>(),f.mov<6>(),f.mov<7>());
        Py_RETURN_NONE;
    };
}
/*--------------------------------------------
 force and energy calculation
 --------------------------------------------*/
void ForceFieldFS::force_calc()
{
    type0* xvec=atoms->x->begin();
    type0* fvec=f->begin();
    type0* rho=rho_ptr->begin();
    elem_type* evec=elem->begin();
    
    int iatm,jatm;
    int icomp,jcomp;
    elem_type ielem,jelem;
    type0 dx0,dx1,dx2,rsq,en;
    type0 dr_rho,dr_phi,r,rho_coef,phi_coef;
    
    int** neighbor_list=neighbor->neighbor_list;
    int* neighbor_list_size=neighbor->neighbor_list_size;
    const int natms=atoms->natms;
    for(int i=0;i<natms;i++)
        rho[i]=0.0;
    
    for(iatm=0;iatm<natms;iatm++)
    {
        ielem=evec[iatm];
        icomp=3*iatm;
        for(int j=0;j<neighbor_list_size[iatm];j++)
        {
            jatm=neighbor_list[iatm][j];
            jelem=evec[jatm];
            
            jcomp=3*jatm;
            dx0=xvec[icomp]-xvec[jcomp];
            dx1=xvec[icomp+1]-xvec[jcomp+1];
            dx2=xvec[icomp+2]-xvec[jcomp+2];
            rsq=dx0*dx0+dx1*dx1+dx2*dx2;
            
            if(rsq<cut_sq_rho[ielem][jelem])
            {
                r=sqrt(rsq);
                dr_rho=r-sqrt(cut_sq_rho[ielem][jelem]);
                rho[iatm]+=dr_rho*dr_rho*(t1[jelem][ielem]
                +t2[jelem][ielem]*dr_rho);
                
                if(jatm<natms)
                    rho[jatm]+=dr_rho*dr_rho*(t1[ielem][jelem]
                    +t2[ielem][jelem]*dr_rho);
                
            }
        }
    }
    
    dynamic->update(rho_ptr);
    
    for(iatm=0;iatm<natms;iatm++)
    {
        ielem=evec[iatm];
        icomp=3*iatm;
        
        if(rho[iatm]>0.0)
            nrgy_strss_lcl[0]+=-A[ielem]*sqrt(rho[iatm]);
        
        for(int j=0;j<neighbor_list_size[iatm];j++)
        {
            jatm=neighbor_list[iatm][j];
            jelem=evec[jatm];
            
            jcomp=3*jatm;
            dx0=xvec[icomp]-xvec[jcomp];
            dx1=xvec[icomp+1]-xvec[jcomp+1];
            dx2=xvec[icomp+2]-xvec[jcomp+2];
            rsq=dx0*dx0+dx1*dx1+dx2*dx2;
            
            if(rsq<cut_sq[ielem][jelem])
            {
                r=sqrt(rsq);
                
                if(rsq<cut_sq_rho[ielem][jelem])
                {
                    dr_rho=r-sqrt(cut_sq_rho[ielem][jelem]);
                    
                    rho_coef=0.0;
                    
                    if(rho[iatm]>0.0)
                        rho_coef+=A[ielem]*(dr_rho*(1.0*t1[jelem][ielem]
                        +1.5*t2[jelem][ielem]*dr_rho))/sqrt(rho[iatm]);
                    
                    if(rho[jatm]>0.0)
                        rho_coef+=A[jelem]*(dr_rho*(1.0*t1[ielem][jelem]
                        +1.5*t2[ielem][jelem]*dr_rho))/sqrt(rho[jatm]);
                    
                    rho_coef*=1.0/r;
                    
                    fvec[icomp]+=dx0*rho_coef;
                    fvec[icomp+1]+=dx1*rho_coef;
                    fvec[icomp+2]+=dx2*rho_coef;
                    
                    if(jatm<natms)
                    {
                        fvec[jcomp]-=dx0*rho_coef;
                        fvec[jcomp+1]-=dx1*rho_coef;
                        fvec[jcomp+2]-=dx2*rho_coef;
                    }
                    
                    if(jatm>=natms)
                        rho_coef*=0.5;
                    
                    nrgy_strss_lcl[1]-=rho_coef*dx0*dx0;
                    nrgy_strss_lcl[2]-=rho_coef*dx0*dx1;
                    nrgy_strss_lcl[3]-=rho_coef*dx0*dx2;
                    nrgy_strss_lcl[4]-=rho_coef*dx1*dx1;
                    nrgy_strss_lcl[5]-=rho_coef*dx1*dx2;
                    nrgy_strss_lcl[6]-=rho_coef*dx2*dx2;
                }
                
                if(rsq<cut_sq_phi[ielem][jelem])
                {
                    dr_phi=r-sqrt(cut_sq_phi[ielem][jelem]);
                    phi_coef=2.0*dr_phi*(k1[ielem][jelem]+k2[ielem][jelem]*r+k3[ielem][jelem]*rsq)
                    +dr_phi*dr_phi*(k2[ielem][jelem]+2.0*k3[ielem][jelem]*r);
                    phi_coef*=-1.0/r;
                    en=dr_phi*dr_phi*(k1[ielem][jelem]+k2[ielem][jelem]*r+k3[ielem][jelem]*rsq);
                    
                    fvec[icomp]+=dx0*phi_coef;
                    fvec[icomp+1]+=dx1*phi_coef;
                    fvec[icomp+2]+=dx2*phi_coef;
                    
                    if(jatm<natms)
                    {
                        fvec[jcomp]-=dx0*phi_coef;
                        fvec[jcomp+1]-=dx1*phi_coef;
                        fvec[jcomp+2]-=dx2*phi_coef;
                    }
                    
                    if(jatm>=natms)
                    {
                        phi_coef*=0.5;
                        en*=0.5;
                    }
                    
                    nrgy_strss_lcl[0]+=en;
                    nrgy_strss_lcl[1]-=phi_coef*dx0*dx0;
                    nrgy_strss_lcl[2]-=phi_coef*dx0*dx1;
                    nrgy_strss_lcl[3]-=phi_coef*dx0*dx2;
                    nrgy_strss_lcl[4]-=phi_coef*dx1*dx1;
                    nrgy_strss_lcl[5]-=phi_coef*dx1*dx2;
                    nrgy_strss_lcl[6]-=phi_coef*dx2*dx2;
                }
            }
        }
    }
}
/*--------------------------------------------
 only energy calculation this is useful for
 minimization/linesearch methods that do not
 use derivatives of energy
 --------------------------------------------*/
void ForceFieldFS::energy_calc()
{
    type0* xvec=atoms->x->begin();
    type0* rho=rho_ptr->begin();
    elem_type* evec=elem->begin();
    
    int iatm,jatm;
    int icomp,jcomp;
    elem_type ielem,jelem;
    type0 dx0,dx1,dx2,rsq;
    type0 dr_rho,dr_phi,r;
    
    int** neighbor_list=neighbor->neighbor_list;
    int* neighbor_list_size=neighbor->neighbor_list_size;
    const int natms=atoms->natms;
    for(int i=0;i<natms;i++)
        rho[i]=0.0;
    
    for(iatm=0;iatm<natms;iatm++)
    {
        ielem=evec[iatm];
        icomp=__dim__*iatm;
        for(int j=0;j<neighbor_list_size[iatm];j++)
        {
            jatm=neighbor_list[iatm][j];
            jelem=evec[jatm];
            
            jcomp=3*jatm;
            dx0=xvec[icomp]-xvec[jcomp];
            dx1=xvec[icomp+1]-xvec[jcomp+1];
            dx2=xvec[icomp+2]-xvec[jcomp+2];
            rsq=dx0*dx0+dx1*dx1+dx2*dx2;

            
            if(rsq<cut_sq[ielem][jelem])
            {
                r=sqrt(rsq);
                if(rsq<cut_sq_rho[ielem][jelem])
                {
                    dr_rho=r-sqrt(cut_sq_rho[ielem][jelem]);
                    rho[iatm]+=dr_rho*dr_rho*(t1[jelem][ielem]
                    +t2[jelem][ielem]*dr_rho);
                    
                    if(jatm<natms)
                        rho[jatm]+=dr_rho*dr_rho*(t1[ielem][jelem]
                        +t2[ielem][jelem]*dr_rho);
                }
                
                if(rsq<cut_sq_phi[ielem][jelem])
                {
                    dr_phi=r-sqrt(cut_sq_phi[ielem][jelem]);
                    if(jatm<natms)
                        nrgy_strss_lcl[0]+=dr_phi*dr_phi*(k1[ielem][jelem]+k2[ielem][jelem]*r+k3[ielem][jelem]*rsq);
                    else
                        nrgy_strss_lcl[0]+=0.5*dr_phi*dr_phi*(k1[ielem][jelem]+k2[ielem][jelem]*r+k3[ielem][jelem]*rsq);
                }
            }
        }
        
        if(rho[iatm]>0.0)
            nrgy_strss_lcl[0]+=-A[ielem]*sqrt(rho[iatm]);
    }
}
/*--------------------------------------------
 initiate before a run
 --------------------------------------------*/
void ForceFieldFS::init()
{
    setup();
    rho_ptr=new Vec<type0>(atoms,1);
}
/*--------------------------------------------
 after a run
 --------------------------------------------*/
void ForceFieldFS::fin()
{
    delete rho_ptr;
}
/*--------------------------------------------
 init xchng
 --------------------------------------------*/
void ForceFieldFS::init_xchng()
{
}
/*--------------------------------------------
 fin xchng
 --------------------------------------------*/
void ForceFieldFS::fin_xchng()
{
}
/*--------------------------------------------
 pre xchng energy
 --------------------------------------------*/
void ForceFieldFS::pre_xchng_energy(GCMC*)
{
}
/*--------------------------------------------
 xchng energy
 --------------------------------------------*/
type0 ForceFieldFS::xchng_energy(GCMC*)
{
    return 0.0;
}
/*--------------------------------------------
 post xchng energy
 --------------------------------------------*/
void ForceFieldFS::post_xchng_energy(GCMC*)
{
}
