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
ForceFieldFS::ForceFieldFS(AtomsMD* __atoms,type0*&& __A,
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
    
    tp_methods.ml_meth=(PyCFunction)(PyCFunctionWithKeywords)
    [](PyObject* self,PyObject* args,PyObject* kwds)->PyObject*
    {
        AtomsMD::Object* __self=reinterpret_cast<AtomsMD::Object*>(self);
        size_t& nelems=__self->atoms->elements.nelems;
        
        FuncAPI<type0*,type0**,type0**,symm<type0**>,symm<type0**>,symm<type0**>,symm<type0**>,symm<type0**>>
        f("ff_fs",{"A","t1","t2","k1","k2","k3","r_c_phi","r_c_rho"});

        f.v<0>().dynamic_size(nelems);
        f.logics<0>()[0]=VLogics("gt",0.0);
        f.v<1>().dynamic_size(nelems,nelems);
        f.v<2>().dynamic_size(nelems,nelems);
        f.v<3>().dynamic_size(nelems,nelems);
        f.v<4>().dynamic_size(nelems,nelems);
        f.v<5>().dynamic_size(nelems,nelems);
        f.v<6>().dynamic_size(nelems,nelems);
        f.logics<6>()[0]=VLogics("gt",0.0);
        f.v<7>().dynamic_size(nelems,nelems);
        f.logics<7>()[0]=VLogics("gt",0.0);
        
        if(f(args,kwds)) return NULL;
        
        delete __self->ff;
        __self->ff=new ForceFieldFS(__self->atoms,f.mov<0>(),f.mov<1>(),f.mov<2>(),
        f.mov<3>(),f.mov<4>(),f.mov<5>(),f.mov<6>(),f.mov<7>());
        Py_RETURN_NONE;
    };
    
    tp_methods.ml_doc=(char*)R"---(
    ff_fs(A,t1,t2,k1,k2,k3,r_c_phi,r_c_rho)
   
    Finnis-Sinclair EAM
    
    Assigns Finnis-Sinclair EAM force field to system. For explanation of the parameter see the Notes section.
    
    Parameters
    ----------
    A : double[nelems]
        :math:`A`
    t1 : double[nelems][nelems]
        :math:`t_1`
    t2 : double[nelems][nelems]
        :math:`t_2`
    k1 : symmetric double[nelems][nelems]
        :math:`k_1`
    k2 : symmetric double[nelems][nelems]
        :math:`k_2`
    k3 : symmetric double[nelems][nelems]
        :math:`k_3`
    r_c_phi : symmetric double[nelems][nelems]
        :math:`r_{c,\phi}`
    r_c_rho : symmetric double[nelems][nelems]
        :math:`r_{c,\rho}`

    
    Returns
    -------
    None
   
    Notes
    -----
    This is the analytical form of Finnis-Sinclair Embedded Atom Method (EAM) potential
    
    .. math::
        U=\sum_{i}\left( -A_{\alpha}\sqrt{\sum_{j\neq i} \rho_{\beta\alpha}(r_{ij})}  + \frac{1}{2}\sum_{j\neq i} \phi_{\beta\alpha}(r_{ij}) \right),
    
    where
    
    .. math::
        \rho_{\beta\alpha}(r)=
        \left\{\begin{array}{ll}
        t^{\alpha\beta}_1(r-r^{\alpha\beta}_{c,\rho})^2+t^{\alpha\beta}_2(r-r^{\alpha\beta}_{c,\rho})^3, \quad & r<r^{\alpha\beta}_{c,\rho}\\
        0 & r>r^{\alpha\beta}_{c,\rho}\
        \end{array}\right.

            
    and
        
    .. math::
        \phi_{\beta\alpha}(r)=
        \left\{\begin{array}{ll}
        (r-r^{\alpha\beta}_{c,\phi})^2(k^{\alpha\beta}_1+k^{\alpha\beta}_2 r+k^{\alpha\beta}_3 r^2), \quad & r<r^{\alpha\beta}_{c,\phi}\\
        0 & r>r^{\alpha\beta}_{c,\phi}\
        \end{array}\right.

    
    Examples
    --------
    Iron Carbon mixture
    ::
     
        >>> from mapp import md
        >>> sim=md.cfg("configs/Cementite.cfg")
        >>> sim.ff_fs(A=[1.8289905,2.9588787],t1=[[1.0,10.024001],[10.482408,0.0]],
                      t2=[[0.504238,1.638980],[3.782595,-7.329211]],
                      k1=[[1.237115],[8.972488,22.061824]],
                      k2=[[-0.35921],[-4.086410,-17.468518]],
                      k3=[[-0.038560],[1.483233,4.812639]],
                      r_c_phi=[[3.40],[2.468801,2.875598]],
                      r_c_rho=[[3.569745],[2.545937,2.892070]])

    )---";
}
/*--------------------------------------------
 force and energy calculation
 --------------------------------------------*/
void ForceFieldFS::force_calc()
{
    type0* xvec=atoms->x->begin();
    type0* fvec=f->begin();
    type0* rho=rho_ptr->begin();
    elem_type* evec=atoms->elem->begin();
    
    int iatm,jatm;
    int icomp,jcomp;
    elem_type ielem,jelem;
    type0 dx0,dx1,dx2,rsq,en;
    type0 dr_rho,dr_phi,r,rho_coef,phi_coef;
    
    int** neighbor_list=neighbor->neighbor_list;
    int* neighbor_list_size=neighbor->neighbor_list_size;
    const int natms_lcl=atoms->natms_lcl;
    for(int i=0;i<natms_lcl;i++)
        rho[i]=0.0;
    
    for(iatm=0;iatm<natms_lcl;iatm++)
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
                
                if(jatm<natms_lcl)
                    rho[jatm]+=dr_rho*dr_rho*(t1[ielem][jelem]
                    +t2[ielem][jelem]*dr_rho);
                
            }
        }
    }
    
    dynamic->update(rho_ptr);
    
    for(iatm=0;iatm<natms_lcl;iatm++)
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
                    
                    if(jatm<natms_lcl)
                    {
                        fvec[jcomp]-=dx0*rho_coef;
                        fvec[jcomp+1]-=dx1*rho_coef;
                        fvec[jcomp+2]-=dx2*rho_coef;
                    }
                    
                    if(jatm>=natms_lcl)
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
                    
                    if(jatm<natms_lcl)
                    {
                        fvec[jcomp]-=dx0*phi_coef;
                        fvec[jcomp+1]-=dx1*phi_coef;
                        fvec[jcomp+2]-=dx2*phi_coef;
                    }
                    
                    if(jatm>=natms_lcl)
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
    elem_type* evec=atoms->elem->begin();
    
    int iatm,jatm;
    int icomp,jcomp;
    elem_type ielem,jelem;
    type0 dx0,dx1,dx2,rsq;
    type0 dr_rho,dr_phi,r;
    
    int** neighbor_list=neighbor->neighbor_list;
    int* neighbor_list_size=neighbor->neighbor_list_size;
    const int natms_lcl=atoms->natms_lcl;
    for(int i=0;i<natms_lcl;i++)
        rho[i]=0.0;
    
    for(iatm=0;iatm<natms_lcl;iatm++)
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
                    
                    if(jatm<natms_lcl)
                        rho[jatm]+=dr_rho*dr_rho*(t1[ielem][jelem]
                        +t2[ielem][jelem]*dr_rho);
                }
                
                if(rsq<cut_sq_phi[ielem][jelem])
                {
                    dr_phi=r-sqrt(cut_sq_phi[ielem][jelem]);
                    if(jatm<natms_lcl)
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
    pre_init();
    rho_ptr=new Vec<type0>(atoms,1);
}
/*--------------------------------------------
 after a run
 --------------------------------------------*/
void ForceFieldFS::fin()
{
    delete rho_ptr;
    post_fin();
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
