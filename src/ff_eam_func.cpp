#include "ff_eam_func.h"
#include "mapp_func.h"
#include "neighbor_md.h"
#include "atoms_md.h"
#include "elements.h"
#include "memory.h"
#include "dynamic_md.h"
#include "gcmc.h"
using namespace MAPP_NS;
/*--------------------------------------------
 This is for my personal and fitting of iron
 and hydrogen
 --------------------------------------------*/
/*--------------------------------------------
 constructor
 --------------------------------------------*/
ForceFieldEAMFunc::ForceFieldEAMFunc(AtomsMD* __atoms,
EAMFunc* __eam_func,
elem_type*&& __elem_map):
ForceFieldMD(__atoms),
eam_func(__eam_func),
elem_map(__elem_map)
{
    gcmc_n_cutoff=2;
    gcmc_n_vars=2;
    gcmc_tag_enabled=true;
    
    __elem_map=NULL;
    for(size_t i=0;i<nelems;i++)
        for(size_t j=0;j<i+1;j++)
        {
            cut[i][j]=cut[j][i]=eam_func->rc[elem_map[i]][elem_map[j]];
            cut_sq[i][j]=cut_sq[j][i]=cut[i][j]*cut[i][j];
        }
}
/*--------------------------------------------
 destructor
 --------------------------------------------*/
ForceFieldEAMFunc::~ForceFieldEAMFunc()
{
    delete [] elem_map;
    delete eam_func;
}
/*--------------------------------------------
 force and energy calculation
 --------------------------------------------*/
void ForceFieldEAMFunc::__force_calc()
{
    const type0* x=atoms->x->begin();
    elem_type* evec=atoms->elem->begin();
    type0* rho=rho_ptr->begin();
    elem_type ielem,jelem;
    type0 rsq,r,phi,fpair;
    
    int** neighbor_list=neighbor->neighbor_list;
    int* neighbor_list_size=neighbor->neighbor_list_size;
    
    type0 dx_ij[__dim__];
    const int natms_lcl=atoms->natms_lcl;
    for(int iatm=0;iatm<natms_lcl;iatm++) rho[iatm]=0.0;
    for(int iatm=0;iatm<natms_lcl;iatm++)
    {
        ielem=elem_map[evec[iatm]];
        const int list_size=neighbor_list_size[iatm];
        for(int j=0,jatm;j<list_size;j++)
        {
            jatm=neighbor_list[iatm][j];
            jelem=elem_map[evec[jatm]];
            rsq=Algebra::RSQ<__dim__>(x+iatm*__dim__,x+jatm*__dim__);
            if(rsq>=cut_sq[ielem][jelem]) continue;
            r=sqrt(rsq);
            phi=eam_func->phi(ielem,jelem,r);
            rho[iatm]+=eam_func->rho(jelem,ielem,r);
            
            if(jatm<natms_lcl)
                rho[jatm]+=eam_func->rho(ielem,jelem,r);
            else
                phi*=0.5;
            __vec_lcl[0]+=phi;
        }
        __vec_lcl[0]+=eam_func->F(ielem,rho[iatm]);
    }
    
#ifdef NEW_UPDATE
    update(rho_ptr);
#else
    dynamic->update(rho_ptr);
#endif
    
    type0* fvec=f->begin();
    for(int iatm=0;iatm<natms_lcl;iatm++)
    {
        ielem=elem_map[evec[iatm]];
        const int list_size=neighbor_list_size[iatm];
        for(int j=0,jatm;j<list_size;j++)
        {
            jatm=neighbor_list[iatm][j];
            jelem=elem_map[evec[jatm]];
            rsq=Algebra::DX_RSQ(x+iatm*__dim__,x+jatm*__dim__,dx_ij);
            if(rsq>=cut_sq[ielem][jelem]) continue;
            r=sqrt(rsq);
            fpair=eam_func->fpair(ielem,jelem,rho[iatm],rho[jatm],r);
            
            
            Algebra::V_add_x_mul_V<__dim__>(fpair,dx_ij,fvec+iatm*__dim__);
            if(jatm<natms_lcl)
                Algebra::V_add_x_mul_V<__dim__>(-fpair,dx_ij,fvec+jatm*__dim__);
            else
                fpair*=0.5;
            
            Algebra::DyadicV<__dim__>(-fpair,dx_ij,&__vec_lcl[1]);
            
        }
    }
    
}
/*--------------------------------------------
 only energy calculation this is useful for
 minimization/linesearch methods that do not
 use derivatives of energy
 --------------------------------------------*/
void ForceFieldEAMFunc::__energy_calc()
{
    const type0* x=atoms->x->begin();
    elem_type* evec=atoms->elem->begin();
    type0* rho=rho_ptr->begin();
    elem_type ielem,jelem;
    type0 rsq,r,phi;
    
    int** neighbor_list=neighbor->neighbor_list;
    int* neighbor_list_size=neighbor->neighbor_list_size;
    

    const int natms_lcl=atoms->natms_lcl;
    for(int iatm=0;iatm<natms_lcl;iatm++) rho[iatm]=0.0;
    for(int iatm=0;iatm<natms_lcl;iatm++)
    {
        ielem=elem_map[evec[iatm]];
        const int list_size=neighbor_list_size[iatm];
        for(int j=0,jatm;j<list_size;j++)
        {
            jatm=neighbor_list[iatm][j];
            jelem=elem_map[evec[jatm]];
            rsq=Algebra::RSQ<__dim__>(x+iatm*__dim__,x+jatm*__dim__);
            if(rsq>=cut_sq[ielem][jelem]) continue;
            r=sqrt(rsq);
            phi=eam_func->phi(ielem,jelem,r);
            rho[iatm]+=eam_func->rho(jelem,ielem,r);
            
            if(jatm<natms_lcl)
                rho[jatm]+=eam_func->rho(ielem,jelem,r);
            else
                phi*=0.5;
            __vec_lcl[0]+=phi;
        }
        __vec_lcl[0]+=eam_func->F(ielem,rho[iatm]);
    }
}
/*--------------------------------------------
 initiate before a run
 --------------------------------------------*/
void ForceFieldEAMFunc::init()
{
    pre_init();
    rho_ptr=new Vec<type0>(atoms,1,"rho");
    
}
/*--------------------------------------------
 after a run
 --------------------------------------------*/
void ForceFieldEAMFunc::fin()
{
    delete rho_ptr;
    post_fin();
}
/*--------------------------------------------
 init xchng
 --------------------------------------------*/
void ForceFieldEAMFunc::init_xchng()
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
    
    for(int iatm=0;iatm<natms_lcl;iatm++)
    {
        ielem=elem_map[evec[iatm]];
        const int list_size=neighbor_list_size[iatm];
        for(int j=0,jatm;j<list_size;j++)
        {
            jatm=neighbor_list[iatm][j];
            jelem=elem_map[evec[jatm]];
            rsq=Algebra::RSQ<__dim__>(x+iatm*__dim__,x+jatm*__dim__);
            if(rsq>=cut_sq[ielem][jelem]) continue;
            r=sqrt(rsq);
            
            rho[iatm]+=eam_func->rho(jelem,ielem,r);
            
            if(jatm<natms_lcl)
                rho[jatm]+=eam_func->rho(ielem,jelem,r);
        }
        F[iatm]=eam_func->F(ielem,rho[iatm]);
    }
}
/*--------------------------------------------
 fin xchng
 --------------------------------------------*/
void ForceFieldEAMFunc::fin_xchng()
{
    delete F_xchng_ptr;
    delete rho_xchng_ptr;
    delete F_ptr;
}
/*--------------------------------------------
 pre xchng energy
 --------------------------------------------*/
void ForceFieldEAMFunc::pre_xchng_energy(GCMC* gcmc)
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

                if(jatm<natms_lcl)
                {
                    en+=eam_func->phi(elem_map[ielem],elem_map[jelem],r);
                    rho_xchng[jatm]+=c0*eam_func->rho(elem_map[ielem],elem_map[jelem],r);
                }
                else
                    en+=0.5*eam_func->phi(elem_map[ielem],elem_map[jelem],r);
                
                rho_iatm_lcl+=eam_func->rho(elem_map[jelem],elem_map[ielem],r);
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
                F_xchng[i]=eam_func->F(elem_map[evec[i]],rho_xchng[i]);
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
type0 ForceFieldEAMFunc::xchng_energy(GCMC* gcmc)
{
    int& icomm=gcmc->icomm;
    for(gcmc->reset_icomm();icomm!=-1;gcmc->next_icomm())
        MPI_Reduce(gcmc->lcl_vars,gcmc->vars,2,Vec<type0>::MPI_T,MPI_SUM,gcmc->curr_root,*gcmc->curr_comm);
 
    
    if(gcmc->im_root)
    {
        //restart the comms
        gcmc->reset_icomm();
        if(gcmc->xchng_mode==NOEX_MODE) return 0.0;
        return eam_func->F(elem_map[gcmc->ielem],gcmc->vars[1])+gcmc->vars[0];
    }
    return 0.0;
}
/*--------------------------------------------
 post xchng energy
 --------------------------------------------*/
void ForceFieldEAMFunc::post_xchng_energy(GCMC* gcmc)
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
        F[natms_lcl-1]=eam_func->F(elem_map[gcmc->ielem],rho[natms_lcl-1]);
    }
}
/*--------------------------------------------
 python constructor
 --------------------------------------------*/
#include <dlfcn.h>
void ForceFieldEAMFunc::ml_new(PyMethodDef& tp_methods)
{
    tp_methods.ml_flags=METH_VARARGS | METH_KEYWORDS;
    tp_methods.ml_name="ff_eam_func";
    
    tp_methods.ml_meth=(PyCFunction)(PyCFunctionWithKeywords)(
    [](PyObject* self,PyObject* args,PyObject* kwds)->PyObject*
    {
        
        FuncAPI<std::string> f("ff_eam_func",{"so_file"});
        if(f(args,kwds)) return NULL;

        void* so=dlopen(f.val<0>().c_str(),RTLD_LAZY);
        create_eam_ff_t* ff_eam=reinterpret_cast<create_eam_ff_t*>( dlsym(so,"create_eam_ff"));
        EAMFunc* ff_eam_func=ff_eam();
        AtomsMD::Object* __self=reinterpret_cast<AtomsMD::Object*>(self);
        elem_type* __elem_map=NULL;
        try
        {
            __elem_map=ff_eam_func->remap(__self->atoms->elements.names,__self->atoms->elements.__nelems);
        }
        catch(std::string& err_msg)
        {
            delete ff_eam_func;
            PyErr_SetString(PyExc_TypeError,err_msg.c_str());
            return NULL;
        }
                
        delete __self->ff;
        __self->ff=new ForceFieldEAMFunc(__self->atoms,ff_eam_func,std::move(__elem_map));
        Py_RETURN_NONE;
    });
    
    tp_methods.ml_doc=(char*)"";
}
