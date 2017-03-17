/*--------------------------------------------
 Created by Sina on 07/15/13.
 Copyright (c) 2013 MIT. All rights reserved.
 --------------------------------------------*/
#include "ff_lj.h"
#include "neighbor_md.h"
#include "atoms_md.h"
#include "elements.h"
#include "memory.h"
#include "pgcmc.h"
#include "api.h"
using namespace MAPP_NS;
/*--------------------------------------------
 constructor
 --------------------------------------------*/
ForceFieldLJ::
ForceFieldLJ(AtomsMD* atoms,type0**&& __epsilon,
type0**&& __sigma,type0**&& __cut,bool __shift):
ForceFieldMD(atoms),
epsilon(__epsilon),
sigma(__sigma),
shift(__shift),
offset(NULL)
{
    gcmc_n_cutoff=1;
    gcmc_n_vars=1;
    gcmc_tag_enabled=false;
    
    __epsilon=NULL;
    __sigma=NULL;

    for(size_t i=0;i<nelems;i++)
        for(size_t j=0;j<i+1;j++)
        {
            cut[i][j]=cut[j][i]=__cut[i][j];
            cut_sq[i][j]=cut_sq[j][i]=__cut[i][j]*__cut[i][j];
        }
    Memory::dealloc(__cut);
    
    Memory::alloc(offset,nelems,nelems);
    for(size_t ielem=0;ielem<nelems;ielem++)
        for(size_t jelem=0;jelem<nelems;jelem++)
            offset[ielem][jelem]=0.0;
    if(shift)
    {
        type0 sig2,sig6,sig12;
        for(size_t i=0;i<nelems;i++)
            for(size_t j=0;j<i+1;j++)
            {
                sig2=sigma[i][j]*sigma[i][j]/cut_sq[i][j];
                sig6=sig2*sig2*sig2;
                sig12=sig6*sig6;
                offset[i][j]=-4.0*epsilon[i][j]*(sig12-sig6);
                offset[j][i]=offset[i][j];
            }
    }
    
}
/*--------------------------------------------
 destructor
 --------------------------------------------*/
ForceFieldLJ::~ForceFieldLJ()
{
    Memory::dealloc(sigma);
    Memory::dealloc(epsilon);
    Memory::dealloc(offset);
}
/*--------------------------------------------
 python constructor
 --------------------------------------------*/
void ForceFieldLJ::ml_new(PyMethodDef& tp_methods)
{
    tp_methods.ml_flags=METH_VARARGS | METH_KEYWORDS;
    tp_methods.ml_name="ff_lj";
    tp_methods.ml_meth=(PyCFunction)(PyCFunctionWithKeywords)
    [](PyObject* self,PyObject* args,PyObject* kwds)->PyObject*
    {
        AtomsMD::Object* __self=reinterpret_cast<AtomsMD::Object*>(self);
        size_t& nelems=__self->atoms->elements.nelems;
        
        FuncAPI<symm<type0**>,symm<type0**>,symm<type0**>,bool> f("ff_lj",{"eps","sigma","r_c","shift"});
        f.noptionals=1;
        f.v<0>().dynamic_size(nelems,nelems);
        f.logics<0>()[0]=VLogics("gt",0.0);
        f.v<1>().dynamic_size(nelems,nelems);
        f.logics<1>()[0]=VLogics("gt",0.0);
        f.v<2>().dynamic_size(nelems,nelems);
        f.logics<2>()[0]=VLogics("gt",0.0);
        f.val<3>()=false;
        
        if(f(args,kwds)) return NULL;
        
        delete __self->ff;
        __self->ff=new ForceFieldLJ(__self->atoms,std::move(f.val<0>()),std::move(f.val<1>()),std::move(f.val<2>()),f.val<3>());
        Py_RETURN_NONE;
    };

    tp_methods.ml_doc=R"---(
    ff_lj(eps,sigma,r_c,shift=False)
   
    Lennard-Jones potential
    
    
    see Notes section below
    
    Parameters
    ----------
    eps : symmetric double[nelems][nelems]
        :math:`\epsilon`
    sigma : symmetric double[nelems][nelems]
        :math:`\sigma`
    r_c : symmetric double[nelems][nelems]
        :math:`r_c`
    shift : bool
        shift the tail if set to True
    
    Returns
    -------
    None
   
    Notes
    -----
    This is the famous Lennard Jones potential
    
    .. math::
        U=\frac{1}{2}\sum_{i}\sum_{j\neq i}
        \left\{\begin{array}{ll}
        4\epsilon_{\alpha\beta}\biggl[\left( \frac{\sigma_{\alpha\beta}}{r_{ij}}\right)^{12}-\left( \frac{\sigma_{\alpha\beta}}{r_{ij}}\right)^6\biggr] &r_{ij}<r^{\alpha\beta}_c\\
        0 &r_{ij}>r^{\alpha\beta}_c
        \end{array}\right.
    
    Examples
    --------
    Kob-Anderson potential
    ::
     
        >>> from mapp import md
        >>> sim=md.cfg("configs/KA.cfg")
        >>> sim.ff_lj(sigma=[[1.0],[0.8,0.88]],
                      eps=[[1.0],[1.5,0.5]],
                      r_c=[[2.5],[2.0,2.2]],
                      shift=False)
        
    )---";
}
/*--------------------------------------------
 initiate before a run
 --------------------------------------------*/
void ForceFieldLJ::init()
{
    pre_init();
}
/*--------------------------------------------
 after a run
 --------------------------------------------*/
void ForceFieldLJ::fin()
{
    post_fin();
}
/*--------------------------------------------
 init xchng
 --------------------------------------------*/
void ForceFieldLJ::init_xchng()
{

}
/*--------------------------------------------
 fin xchng
 --------------------------------------------*/
void ForceFieldLJ::fin_xchng()
{

}
/*--------------------------------------------
 pre xchng energy
 --------------------------------------------*/
void ForceFieldLJ::pre_xchng_energy(GCMC* gcmc)
{
    
    int& icomm=gcmc->icomm;
    type0 en;
    
    int& iatm=gcmc->iatm;
    elem_type& ielem=gcmc->ielem;
    
    int& jatm=gcmc->jatm;
    elem_type& jelem=gcmc->jelem;
    
    type0&rsq=gcmc->rsq;
    
    type0 sig2,sig6;
    
    const int natms_lcl=atoms->natms_lcl;
    
    for(gcmc->reset_icomm();icomm!=-1;gcmc->next_icomm())
    {
        en=0.0;
        for(gcmc->reset_iatm();iatm!=-1;gcmc->next_iatm())
            for(gcmc->reset_jatm();jatm!=-1;gcmc->next_jatm())
            {
                sig2=sigma[ielem][jelem]*sigma[ielem][jelem]/rsq;
                sig6=sig2*sig2*sig2;
                
                if(jatm<natms_lcl)
                    en+=4.0*epsilon[ielem][jelem]*sig6*(sig6-1.0)+offset[ielem][jelem];
                else
                    en+=2.0*epsilon[ielem][jelem]*sig6*(sig6-1.0)+offset[ielem][jelem];
            }
        
        gcmc->lcl_vars[0]=en;
        
    }
}
/*--------------------------------------------
 xchng energy
 --------------------------------------------*/
type0 ForceFieldLJ::xchng_energy(GCMC* gcmc)
{
    int& icomm=gcmc->icomm;
    for(gcmc->reset_icomm();icomm!=-1;gcmc->next_icomm())
        MPI_Reduce(gcmc->lcl_vars,gcmc->vars,1,Vec<type0>::MPI_T,MPI_SUM,gcmc->curr_root,*gcmc->curr_comm);
    if(gcmc->im_root)
        return gcmc->vars[0];
    
    return 0.0;
}
/*--------------------------------------------
 post xchng energy
 --------------------------------------------*/
void ForceFieldLJ::post_xchng_energy(GCMC*)
{
}
/*--------------------------------------------
 force and energy calculation
 --------------------------------------------*/
#include "xmath.h"
void ForceFieldLJ::force_calc()
{
    const type0* x=atoms->x->begin();
    type0* fvec=f->begin();
    elem_type* evec=atoms->elem->begin();
    
    elem_type ielem,jelem;
    type0 rsq;
    
    
    
    int** neighbor_list=neighbor->neighbor_list;
    int* neighbor_list_size=neighbor->neighbor_list_size;
    
    type0 x_i[__dim__];
    type0 dx_ij[__dim__];
    const int natms_lcl=atoms->natms_lcl;
    for(int iatm=0;iatm<natms_lcl;iatm++)
    {
        ielem=evec[iatm];
        Algebra::V_eq<__dim__>(x+iatm*__dim__,x_i);
        type0 f_i[__dim__]{[0 ... __dim__-1]=0.0};
        const int list_size=neighbor_list_size[iatm];
        for(int j=0,jatm;j<list_size;j++)
        {
            jatm=neighbor_list[iatm][j];
            jelem=evec[jatm];
            rsq=Algebra::DX_RSQ(x_i,x+jatm*__dim__,dx_ij);
            if(rsq>=cut_sq[ielem][jelem]) continue;
            
            type0 sig2=sigma[ielem][jelem]*sigma[ielem][jelem]/rsq;
            type0 sig6=sig2*sig2*sig2;
            
            type0 eps=epsilon[ielem][jelem];
            type0 fpair=24.0*eps*sig6*(2.0*sig6-1.0)/rsq;
            type0 en=4.0*eps*sig6*(sig6-1.0)+offset[ielem][jelem];
            
            Algebra::V_add_x_mul_V<__dim__>(fpair,dx_ij,f_i);
            
            if(jatm<natms_lcl)
                Algebra::V_add_x_mul_V<__dim__>(-fpair,dx_ij,fvec+__dim__*jatm);
            else
            {
                fpair*=0.5;
                en*=0.5;
            }
            nrgy_strss_lcl[0]+=en;
            Algebra::DyadicV(-fpair,dx_ij,&nrgy_strss_lcl[1]);
        }
        
        Algebra::V_add<__dim__>(f_i,fvec+iatm*__dim__);
    }
}
/*--------------------------------------------
 only energy calculation this is useful for
 minimization/linesearch methods that do not
 use derivatives of energy
 --------------------------------------------*/
void ForceFieldLJ::energy_calc()
{
    const type0* x=atoms->x->begin();
    elem_type* evec=atoms->elem->begin();
    
    elem_type ielem,jelem;
    
    int** neighbor_list=neighbor->neighbor_list;
    int* neighbor_list_size=neighbor->neighbor_list_size;
    
    type0 xi[__dim__];
    const int natms_lcl=atoms->natms_lcl;
    for(int iatm=0;iatm<natms_lcl;iatm++)
    {
        ielem=evec[iatm];
        Algebra::V_eq<__dim__>(x+iatm*__dim__,xi);
        const int list_size=neighbor_list_size[iatm];
        for(int j=0,jatm;j<list_size;j++)
        {
            jatm=neighbor_list[iatm][j];
            jelem=evec[jatm];
            type0 rsq=Algebra::RSQ<__dim__>(xi,x+jatm*__dim__);
            
            if(rsq>=cut_sq[ielem][jelem]) continue;
            type0 sig2=sigma[ielem][jelem]*sigma[ielem][jelem]/rsq;
            type0 sig6=sig2*sig2*sig2;
            
            if(jatm<natms_lcl)
                nrgy_strss_lcl[0]+=4.0*epsilon[ielem][jelem]*sig6*(sig6-1.0)+offset[ielem][jelem];
            else
                nrgy_strss_lcl[0]+=2.0*epsilon[ielem][jelem]*sig6*(sig6-1.0)+0.5*offset[ielem][jelem];
        }
    }

}
