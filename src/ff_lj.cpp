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
ForceFieldLJ(AtomsMD*& atoms,type0**&& __epsilon,
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
    tp_methods.ml_doc="I will add doc here";
    
    tp_methods.ml_meth=(PyCFunction)(PyCFunctionWithKeywords)
    [](PyObject* self,PyObject* args,PyObject* kwds)->PyObject*
    {
        AtomsMD::Object* __self=reinterpret_cast<AtomsMD::Object*>(self);
        size_t& nelems=__self->atoms->elements->nelems;
        
        FuncAPI<symm<type0**>,symm<type0**>,symm<type0**>,bool> f("ff_lj",{"eps","sigma","r_c","shift"});
        f.noptionals=1;
        f.var<0>().dynamic_size(nelems,nelems);
        f.logics<0>()[0]=VLogics("gt",0.0);
        f.var<1>().dynamic_size(nelems,nelems);
        f.logics<1>()[0]=VLogics("gt",0.0);
        f.var<2>().dynamic_size(nelems,nelems);
        f.logics<2>()[0]=VLogics("gt",0.0);
        f.val<3>()=false;
        
        if(f(args,kwds)) return NULL;
        
        delete __self->ff;
        __self->ff=new ForceFieldLJ(__self->atoms,std::move(f.val<0>()),std::move(f.val<1>()),std::move(f.val<2>()),f.val<3>());
        Py_RETURN_NONE;
    };
}
/*--------------------------------------------
 initiate before a run
 --------------------------------------------*/
void ForceFieldLJ::init()
{
    setup();
}
/*--------------------------------------------
 after a run
 --------------------------------------------*/
void ForceFieldLJ::fin()
{
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
    
    const int natms=atoms->natms;
    
    for(gcmc->reset_icomm();icomm!=-1;gcmc->next_icomm())
    {
        en=0.0;
        for(gcmc->reset_iatm();iatm!=-1;gcmc->next_iatm())
            for(gcmc->reset_jatm();jatm!=-1;gcmc->next_jatm())
            {
                sig2=sigma[ielem][jelem]*sigma[ielem][jelem]/rsq;
                sig6=sig2*sig2*sig2;
                
                if(jatm<natms)
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
        MPI_Reduce(gcmc->lcl_vars,gcmc->vars,1,MPI_TYPE0,MPI_SUM,gcmc->curr_root,*gcmc->curr_comm);
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
    elem_type* evec=elem->begin();
    
    elem_type ielem,jelem;
    type0 rsq;
    
    
    
    int** neighbor_list=neighbor->neighbor_list;
    int* neighbor_list_size=neighbor->neighbor_list_size;
    
    type0 xi[__dim__];
    type0 dxij[__dim__];
    const int natms=atoms->natms;
    for(int iatm=0;iatm<natms;iatm++)
    {
        ielem=evec[iatm];
        Algebra::V_eq<__dim__>(x+iatm*__dim__,xi);
        type0 fi[__dim__]{[0 ... __dim__-1]=0.0};
        const int list_size=neighbor_list_size[iatm];
        for(int j=0,jatm;j<list_size;j++)
        {
            jatm=neighbor_list[iatm][j];
            jelem=evec[jatm];
            rsq=Algebra::DX_RSQ(xi,x+jatm*__dim__,dxij);
            if(rsq>=cut_sq[ielem][jelem]) continue;
            
            type0 sig2=sigma[ielem][jelem]*sigma[ielem][jelem]/rsq;
            type0 sig6=sig2*sig2*sig2;
            
            type0 eps=epsilon[ielem][jelem];
            type0 fpair=24.0*eps*sig6*(2.0*sig6-1.0)/rsq;
            type0 en=4.0*eps*sig6*(sig6-1.0)+offset[ielem][jelem];
            
            Algebra::V_add_x_mul_V<__dim__>(fpair,dxij,fi);
            
            if(jatm<natms)
                Algebra::V_add_x_mul_V<__dim__>(-fpair,dxij,fvec+__dim__*jatm);
            else
            {
                fpair*=0.5;
                en*=0.5;
            }
            nrgy_strss_lcl[0]+=en;
            Algebra::DyadicV(-fpair,dxij,&nrgy_strss_lcl[1]);
        }
        
        Algebra::V_add<__dim__>(fi,fvec+iatm*__dim__);
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
    elem_type* evec=elem->begin();
    
    elem_type ielem,jelem;
    
    int** neighbor_list=neighbor->neighbor_list;
    int* neighbor_list_size=neighbor->neighbor_list_size;
    
    type0 xi[__dim__];
    const int natms=atoms->natms;
    for(int iatm=0;iatm<natms;iatm++)
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
            
            if(jatm<natms)
                nrgy_strss_lcl[0]+=4.0*epsilon[ielem][jelem]*sig6*(sig6-1.0)+offset[ielem][jelem];
            else
                nrgy_strss_lcl[0]+=2.0*epsilon[ielem][jelem]*sig6*(sig6-1.0)+0.5*offset[ielem][jelem];
        }
    }

}
