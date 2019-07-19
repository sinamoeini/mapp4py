#include "ff_fs.h"
#include "neighbor_md.h"
#include "atoms_md.h"
#include "elements.h"
#include "memory.h"
#include "dynamic_md.h"
#include "gcmc.h"
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
cut_phi(r_c_phi),
cut_rho(r_c_rho),
t1(__t1),
t2(__t2),
k1(__k1),
k2(__k2),
k3(__k3),
A(__A)
{
    gcmc_n_cutoff=2;
    gcmc_n_vars=2;
    gcmc_tag_enabled=true;
    
    __A=NULL;
    __t1=NULL;
    __t2=NULL;
    __k1=NULL;
    __k2=NULL;
    __k3=NULL;
    r_c_phi=NULL;
    r_c_rho=NULL;
    
    Memory::alloc(cut_sq_phi,nelems,nelems);
    Memory::alloc(cut_sq_rho,nelems,nelems);
    type0 tmp;
    for(size_t i=0;i<nelems;i++)
        for(size_t j=0;j<i+1;j++)
        {
            cut[i][j]=cut[j][i]=MAX(cut_phi[i][j],cut_rho[i][j]);
            cut_sq[i][j]=cut_sq[j][i]=cut[i][j]*cut[i][j];
            tmp=cut_phi[i][j];
            cut_phi[j][i]=tmp;
            cut_sq_phi[i][j]=cut_sq_phi[j][i]=tmp*tmp;
            tmp=cut_rho[i][j];
            cut_rho[j][i]=tmp;
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
    Memory::dealloc(cut_phi);
    Memory::dealloc(cut_rho);
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
    
    tp_methods.ml_meth=(PyCFunction)(PyCFunctionWithKeywords)(
    [](PyObject* self,PyObject* args,PyObject* kwds)->PyObject*
    {
        AtomsMD::Object* __self=reinterpret_cast<AtomsMD::Object*>(self);

        
        FuncAPI<type0*,sq<type0**>,sq<type0**>,symm<type0**>,symm<type0**>,
        symm<type0**>,symm<type0**>,symm<type0**>,std::string*>
        f("ff_fs",{"A","t1","t2","k1","k2","k3","r_c_phi","r_c_rho","elems"});
        f.noptionals=1;
        f.logics<0>()[0]=VLogics("gt",0.0);
        f.logics<6>()[0]=VLogics("ge",0.0);
        f.logics<7>()[0]=VLogics("ge",0.0);
        
        const std::string* names=__self->atoms->elements.names;
        const size_t nelems=__self->atoms->elements.nelems;
        if(f(args,kwds)) return NULL;
        if(f.remap<8,0,1,2,3,4,5,6,7>("elements present in system",names,nelems))
            return NULL;
        
        
        delete __self->ff;
        __self->ff=new ForceFieldFS(__self->atoms,f.mov<0>(),f.mov<1>(),f.mov<2>(),
        f.mov<3>(),f.mov<4>(),f.mov<5>(),f.mov<6>(),f.mov<7>());
        Py_RETURN_NONE;
    });
    
    tp_methods.ml_doc=(char*)R"---(
    ff_fs(A,t1,t2,k1,k2,k3,r_c_phi,r_c_rho,elems=None)
   
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
    elems : string[nelems]
        mapping elements
    
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
        >>> sim.ff_fs(A=[1.8289905,2.9588787],
                      t1=[[1.0,10.024001],[10.482408,0.0]],
                      t2=[[0.504238,1.638980],[3.782595,-7.329211]],
                      k1=[[1.237115],[8.972488,22.061824]],
                      k2=[[-0.35921],[-4.086410,-17.468518]],
                      k3=[[-0.038560],[1.483233,4.812639]],
                      r_c_phi=[[3.40],[2.468801,2.875598]],
                      r_c_rho=[[3.569745],[2.545937,2.892070]],
                      elems=['Fe','C'])

    )---";
}
/*--------------------------------------------
 force and energy calculation
 --------------------------------------------*/
void ForceFieldFS::__force_calc()
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
                dr_rho=r-cut_rho[ielem][jelem];
                rho[iatm]+=dr_rho*dr_rho*(t1[jelem][ielem]
                +t2[jelem][ielem]*dr_rho);
                
                if(jatm<natms_lcl)
                    rho[jatm]+=dr_rho*dr_rho*(t1[ielem][jelem]
                    +t2[ielem][jelem]*dr_rho);
                
            }
        }
    }
    
    update(rho_ptr);

    for(iatm=0;iatm<natms_lcl;iatm++)
    {
        ielem=evec[iatm];
        icomp=3*iatm;
        
        if(rho[iatm]>0.0)
            __vec_lcl[0]+=-A[ielem]*sqrt(rho[iatm]);
        
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
                    dr_rho=r-cut_rho[ielem][jelem];
                    
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
                    
                    __vec_lcl[1]-=rho_coef*dx0*dx0;
                    __vec_lcl[2]-=rho_coef*dx0*dx1;
                    __vec_lcl[3]-=rho_coef*dx0*dx2;
                    __vec_lcl[4]-=rho_coef*dx1*dx1;
                    __vec_lcl[5]-=rho_coef*dx1*dx2;
                    __vec_lcl[6]-=rho_coef*dx2*dx2;
                }
                
                if(rsq<cut_sq_phi[ielem][jelem])
                {
                    dr_phi=r-cut_phi[ielem][jelem];
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
                    
                    __vec_lcl[0]+=en;
                    __vec_lcl[1]-=phi_coef*dx0*dx0;
                    __vec_lcl[2]-=phi_coef*dx0*dx1;
                    __vec_lcl[3]-=phi_coef*dx0*dx2;
                    __vec_lcl[4]-=phi_coef*dx1*dx1;
                    __vec_lcl[5]-=phi_coef*dx1*dx2;
                    __vec_lcl[6]-=phi_coef*dx2*dx2;
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
void ForceFieldFS::__energy_calc()
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
                    dr_rho=r-cut_rho[ielem][jelem];
                    rho[iatm]+=dr_rho*dr_rho*(t1[jelem][ielem]
                    +t2[jelem][ielem]*dr_rho);
                    
                    if(jatm<natms_lcl)
                        rho[jatm]+=dr_rho*dr_rho*(t1[ielem][jelem]
                        +t2[ielem][jelem]*dr_rho);
                }
                
                if(rsq<cut_sq_phi[ielem][jelem])
                {
                    dr_phi=r-cut_phi[ielem][jelem];
                    if(jatm<natms_lcl)
                        __vec_lcl[0]+=dr_phi*dr_phi*(k1[ielem][jelem]+k2[ielem][jelem]*r+k3[ielem][jelem]*rsq);
                    else
                        __vec_lcl[0]+=0.5*dr_phi*dr_phi*(k1[ielem][jelem]+k2[ielem][jelem]*r+k3[ielem][jelem]*rsq);
                }
            }
        }
        
        if(rho[iatm]>0.0)
            __vec_lcl[0]+=-A[ielem]*sqrt(rho[iatm]);
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
    F_ptr=new Vec<type0>(atoms,1);
    rho_xchng_ptr=new Vec<type0>(atoms,1);
    F_xchng_ptr=new Vec<type0>(atoms,1);
    
    const type0* x=atoms->x->begin();
    type0* rho=rho_ptr->begin();
    type0* F=F_ptr->begin();
    type0* rho_xchng=rho_xchng_ptr->begin();
    type0* F_xchng=F_xchng_ptr->begin();
    elem_type* evec=atoms->elem->begin();
    type0 rsq,r;
    
    elem_type ielem,jelem;
    
    
    int** neighbor_list=neighbor->neighbor_list;
    int* neighbor_list_size=neighbor->neighbor_list_size;
    const int natms_lcl=atoms->natms_lcl;
    for(int iatm=0;iatm<natms_lcl;iatm++) rho[iatm]=rho_xchng[iatm]=F_xchng[iatm]=0.0;
    
    type0 dr_rho;
    for(int iatm=0;iatm<natms_lcl;iatm++)
    {
        ielem=evec[iatm];
        const int list_size=neighbor_list_size[iatm];
        for(int j=0,jatm;j<list_size;j++)
        {
            jatm=neighbor_list[iatm][j];
            jelem=evec[jatm];
            rsq=Algebra::RSQ<__dim__>(x+iatm*__dim__,x+jatm*__dim__);
            if(rsq>=cut_sq[ielem][jelem]) continue;
            r=sqrt(rsq);
            
            
            dr_rho=r-cut_rho[ielem][jelem];
            if(dr_rho <0.0)
            {
                rho[iatm]+=dr_rho*dr_rho*(t1[jelem][ielem]
                        +t2[jelem][ielem]*dr_rho);
                
                
                if(jatm<natms_lcl)
                    rho[jatm]+=dr_rho*dr_rho*(t1[ielem][jelem]
                            +t2[ielem][jelem]*dr_rho);
            }
            
        }
        F[iatm]=rho[iatm]<0.0 ? 0.0:(-A[ielem]*sqrt(rho[iatm]));
    }

}
/*--------------------------------------------
 fin xchng
 --------------------------------------------*/
void ForceFieldFS::fin_xchng()
{
    delete F_xchng_ptr;
    delete rho_xchng_ptr;
    delete F_ptr;
}
/*--------------------------------------------
 pre xchng energy
 --------------------------------------------*/
void ForceFieldFS::pre_xchng_energy(GCMC* gcmc)
{
    int& icomm=gcmc->icomm;
    type0 en;
    type0 rho_iatm_lcl;
    
    int& iatm=gcmc->iatm;
    elem_type& ielem=gcmc->ielem;
    
    int& jatm=gcmc->jatm;
    elem_type& jelem=gcmc->jelem;
    
    
    type0&rsq=gcmc->rsq;
    
    type0 r;
    type0 c0=1.0,en0=0.0;
    
    type0* rho=rho_ptr->begin();
    type0* rho_xchng=rho_xchng_ptr->begin();
    int* tag=gcmc->tag_vec_p->begin();
    const int natms_lcl=atoms->natms_lcl;
    for(int i=0;i<natms_lcl;i++)
        rho_xchng[i]=rho[i];
    type0 dr_rho,dr_phi;
    for(gcmc->reset_icomm();icomm!=-1;gcmc->next_icomm())
    {
        c0=1.0;
        if(gcmc->xchng_mode==DEL_MODE) c0=-1.0;
        else if(gcmc->xchng_mode==NOEX_MODE) continue;
        
        en=rho_iatm_lcl=0.0;
        for(gcmc->reset_iatm();iatm!=-1;gcmc->next_iatm())
            for(gcmc->reset_jatm();jatm!=-1;gcmc->next_jatm())
            {
                r=sqrt(rsq);
                
                dr_rho=r-cut_rho[ielem][jelem];
                dr_phi=r-cut_phi[ielem][jelem];
                
                
                if(jatm<natms_lcl)
                {
                    if(dr_phi<0.0)
                        en+=dr_phi*dr_phi*(k1[ielem][jelem]+k2[ielem][jelem]*r+k3[ielem][jelem]*rsq);
                    if(dr_rho<0.0)
                        rho_xchng[jatm]+=c0*dr_rho*dr_rho*(t1[ielem][jelem]
                        +t2[ielem][jelem]*dr_rho);
                }
                else if(dr_phi<0.0)
                    en+=0.5*dr_phi*dr_phi*(k1[ielem][jelem]+k2[ielem][jelem]*r+k3[ielem][jelem]*rsq);
                
                if(dr_rho<0.0)
                    rho_iatm_lcl+=dr_rho*dr_rho*(t1[jelem][ielem]+t2[jelem][ielem]*dr_rho);
            }
        
        
        
        type0* F=F_ptr->begin();
        type0* F_xchng=F_xchng_ptr->begin();
        elem_type* evec=atoms->elem->begin();
        
        type0 tmp0;
        en0=0.0;
        for(int i=0;i<natms_lcl;i++)
            if(tag[i]==icomm)
            {
                tmp0=rho_xchng[i];
                
                F_xchng[i]=rho_xchng[i]<0.0 ? 0.0:(-A[evec[i]]*sqrt(rho_xchng[i]));
                en0+=F_xchng[i]-F[i];
            }
        
        en+=en0*c0;
        
        gcmc->lcl_vars[0]=en;
        gcmc->lcl_vars[1]=rho_iatm_lcl;
        
    }
}
/*--------------------------------------------
 xchng energy
 --------------------------------------------*/
type0 ForceFieldFS::xchng_energy(GCMC* gcmc)
{
    int& icomm=gcmc->icomm;
    for(gcmc->reset_icomm();icomm!=-1;gcmc->next_icomm())
        MPI_Reduce(gcmc->lcl_vars,gcmc->vars,2,Vec<type0>::MPI_T,MPI_SUM,gcmc->curr_root,*gcmc->curr_comm);
    
    
    if(gcmc->im_root)
    {
        //restart the comms
        gcmc->reset_icomm();
        if(gcmc->xchng_mode==NOEX_MODE) return 0.0;
        return (gcmc->vars[1]<0.0 ? 0.0:(-A[gcmc->ielem]*sqrt(gcmc->vars[1])))+gcmc->vars[0];
        
    }
    return 0.0;
}
/*--------------------------------------------
 post xchng energy
 --------------------------------------------*/
void ForceFieldFS::post_xchng_energy(GCMC* gcmc)
{
    int* tag=gcmc->tag_vec_p->begin();
    type0* rho=rho_ptr->begin();
    type0* F=F_ptr->begin();
    
    type0* rho_xchng=rho_xchng_ptr->begin();
    type0* F_xchng=F_xchng_ptr->begin();
    const int natms_lcl=atoms->natms_lcl;
    for(int i=0;i<natms_lcl;i++)
        if(tag[i]==0)
        {
            rho[i]=rho_xchng[i];
            F[i]=F_xchng[i];
        }
    
    if(gcmc->im_root && gcmc->xchng_mode==INS_MODE && gcmc->root_succ)
    {
        rho[natms_lcl-1]=gcmc->vars[1];
        F[natms_lcl-1]=rho[natms_lcl-1]<0.0 ? 0.0:(-A[gcmc->ielem]*sqrt(rho[natms_lcl-1]));
    }
}
