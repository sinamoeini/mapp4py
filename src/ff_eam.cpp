#include "ff_eam.h"
#include "neighbor_md.h"
#include "atoms_md.h"
#include "elements.h"
#include "memory.h"
#include "api.h"
#include "import_eam.h"
#include "dynamic_md.h"
#include "gcmc.h"
using namespace MAPP_NS;
/*--------------------------------------------
 constructor
 --------------------------------------------*/
ForceFieldEAM::ForceFieldEAM(AtomsMD* atoms,
type0 __dr,type0 __drho,size_t __nr,size_t __nrho,
type0(***&& __r_phi_arr)[7],type0(***&& __rho_arr)[7],type0(**&& __F_arr)[7],
type0**&& __cut):
ForceFieldMD(atoms),
max_pairs(0),
drhoi_dr(NULL),
drhoj_dr(NULL),
dr(__dr),
drho(__drho),
nr(__nr),
nrho(__nrho),
r_phi_arr(__r_phi_arr),
rho_arr(__rho_arr),
F_arr(__F_arr)
{
    gcmc_n_cutoff=2;
    gcmc_n_vars=2;
    gcmc_tag_enabled=true;
    
    __r_phi_arr=NULL;
    __rho_arr=NULL;
    __F_arr=NULL;
    
    for(size_t i=0;i<nelems;i++)
        for(size_t j=0;j<i+1;j++)
        {
            cut[i][j]=cut[j][i]=__cut[i][j];
            cut_sq[i][j]=cut_sq[j][i]=__cut[i][j]*__cut[i][j];
        }
    Memory::dealloc(__cut);
    
    dr_inv=1.0/dr;
    drho_inv=1.0/drho;
    rho_max=static_cast<type0>(nrho)*drho;
}
/*--------------------------------------------
 destructor
 --------------------------------------------*/
ForceFieldEAM::~ForceFieldEAM()
{
    Memory::dealloc(drhoj_dr);
    Memory::dealloc(drhoi_dr);
    Memory::dealloc(F_arr);
    Memory::dealloc(rho_arr);
    Memory::dealloc(r_phi_arr);
}
/*--------------------------------------------
 python constructor
 --------------------------------------------*/
void ForceFieldEAM::ml_new(PyMethodDef& method_0,PyMethodDef& method_1,PyMethodDef& method_2)
{
    method_0.ml_flags=METH_VARARGS | METH_KEYWORDS;
    method_0.ml_name="ff_eam_funcfl";
    method_0.ml_meth=(PyCFunction)(PyCFunctionWithKeywords)
    [](PyObject* self,PyObject* args,PyObject* kwds)->PyObject*
    {
        AtomsMD::Object* __self=reinterpret_cast<AtomsMD::Object*>(self);
        size_t& nelems=__self->atoms->elements.nelems;
        
        FuncAPI<std::string*> f("ff_eam_funcfl",{"funcfl_files"});
        f.v<0>().dynamic_size(nelems);
        if(f(args,kwds)) return NULL;
        
        size_t nr,nrho;
        type0 dr,drho;
        type0** r_c;
        type0(** F)[7]=NULL;
        type0(*** r_phi)[7]=NULL;
        type0(*** rho)[7]=NULL;
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
        
        delete __self->ff;
        __self->ff=new ForceFieldEAM(__self->atoms,dr,drho,nr,nrho,std::move(r_phi),std::move(rho),std::move(F),std::move(r_c));
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
     
        >>> from mapp import md
        >>> sim=md.cfg("configs/Ni.cfg")
        >>> sim.ff_eam_funcfl("potentials/niu3.eam")
    
    

    )---";
    
    method_1.ml_flags=METH_VARARGS | METH_KEYWORDS;
    method_1.ml_name="ff_eam_setfl";
    method_1.ml_meth=(PyCFunction)(PyCFunctionWithKeywords)
    [](PyObject* self,PyObject* args,PyObject* kwds)->PyObject*
    {
        AtomsMD::Object* __self=reinterpret_cast<AtomsMD::Object*>(self);
        size_t& nelems=__self->atoms->elements.nelems;
        FuncAPI<std::string> f("ff_eam_setfl",{"setfl_file"});
        if(f(args,kwds)) return NULL;
        
        
        size_t nr,nrho;
        type0 dr,drho;
        type0** r_c;
        type0(** F)[7]=NULL;
        type0(*** r_phi)[7]=NULL;
        type0(*** rho)[7]=NULL;
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
        
        delete __self->ff;
        __self->ff=new ForceFieldEAM(__self->atoms,dr,drho,nr,nrho,std::move(r_phi),std::move(rho),std::move(F),std::move(r_c));
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
     
        >>> from mapp import md
        >>> sim=md.cfg("configs/Cu.cfg")
        >>> sim.ff_eam_setfl("potentials/Cu_mishin.eam.alloy")

    
    
    )---";
    
    
    method_2.ml_flags=METH_VARARGS | METH_KEYWORDS;
    method_2.ml_name="ff_eam_fs";
    method_2.ml_meth=(PyCFunction)(PyCFunctionWithKeywords)
    [](PyObject* self,PyObject* args,PyObject* kwds)->PyObject*
    {
        AtomsMD::Object* __self=reinterpret_cast<AtomsMD::Object*>(self);
        size_t& nelems=__self->atoms->elements.nelems;
        FuncAPI<std::string> f("ff_eam_fs",{"fs_file"});
        if(f(args,kwds)) return NULL;
        
        size_t nr,nrho;
        type0 dr,drho;
        type0** r_c;
        type0(** F)[7]=NULL;
        type0(*** r_phi)[7]=NULL;
        type0(*** rho)[7]=NULL;
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
        
        delete __self->ff;
        __self->ff=new ForceFieldEAM(__self->atoms,dr,drho,nr,nrho,std::move(r_phi),std::move(rho),std::move(F),std::move(r_c));
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
     
        >>> from mapp import md
        >>> sim=md.cfg("configs/FeH.cfg")
        >>> sim.ff_eam_fs("potentials/FeH.eam.fs")
    
    

    )---";
    
}
/*--------------------------------------------
 force calculation
 --------------------------------------------*/
void ForceFieldEAM::force_calc()
{
    if(max_pairs<neighbor->no_pairs)
    {
        Memory::dealloc(drhoj_dr);
        Memory::dealloc(drhoi_dr);
        max_pairs=neighbor->no_pairs;
        Memory::alloc(drhoi_dr,max_pairs);
        Memory::alloc(drhoj_dr,max_pairs);
    }
    
    type0* xvec=atoms->x->begin();
    type0* fvec=f->begin();
    elem_type* evec=atoms->elem->begin();
    type0* rho=rho_ptr->begin();
    
    int iatm,jatm;
    int icomp,jcomp;
    elem_type ielem,jelem;
    type0 dx0,dx1,dx2,rsq,z2p,z2;
    type0 r,p,r_inv,fpair,tmp0,tmp1;
    type0 drho_i_dr,drho_j_dr,dphi_dr;
    type0 rho_i,rho_j,phi;
    size_t m,istart;
    type0* coef;
    
    int** neighbor_list=neighbor->neighbor_list;
    int* neighbor_list_size=neighbor->neighbor_list_size;

    const int natms_lcl=atoms->natms_lcl;
    for(iatm=0;iatm<natms_lcl;iatm++) rho[iatm]=0.0;
    
    istart=0;
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
            drhoi_dr[istart]=drhoj_dr[istart]=0.0;
            
            if(rsq>=cut_sq[ielem][jelem])
            {
                istart++;
                continue;
            }
            
            r=sqrt(rsq);
            r_inv=1.0/r;
            p=r*dr_inv;
            m=static_cast<size_t>(p);
            m=MIN(m,nr-2);
            
            p-=m;
            p=MIN(p,1.0);
            
            coef=rho_arr[jelem][ielem][m];
            rho_i=((coef[3]*p+coef[2])*p+coef[1])*p+coef[0];
            drho_i_dr=(coef[6]*p+coef[5])*p+coef[4];
            coef=rho_arr[ielem][jelem][m];
            rho_j=((coef[3]*p+coef[2])*p+coef[1])*p+coef[0];
            drho_j_dr=(coef[6]*p+coef[5])*p+coef[4];
            
            coef=r_phi_arr[ielem][jelem][m];
            z2p=(coef[6]*p + coef[5])*p+coef[4];
            z2=((coef[3]*p+coef[2])*p+coef[1])*p+coef[0];
            
            phi=z2*r_inv;
            dphi_dr=z2p*r_inv-phi*r_inv;
            
            rho[iatm]+=rho_i;
            if(jatm<natms_lcl)
                rho[jatm]+=rho_j;
            
            fpair=-dphi_dr*r_inv;
            
            fvec[icomp]+=fpair*dx0;
            fvec[icomp+1]+=fpair*dx1;
            fvec[icomp+2]+=fpair*dx2;
            
            if(jatm<natms_lcl)
            {
                fvec[jcomp]-=fpair*dx0;
                fvec[jcomp+1]-=fpair*dx1;
                fvec[jcomp+2]-=fpair*dx2;
            }
            
            if(jatm>=natms_lcl)
            {
                fpair*=0.5;
                phi*=0.5;
            }
            
            nrgy_strss_lcl[0]+=phi;
            nrgy_strss_lcl[1]-=fpair*dx0*dx0;
            nrgy_strss_lcl[2]-=fpair*dx0*dx1;
            nrgy_strss_lcl[3]-=fpair*dx0*dx2;
            nrgy_strss_lcl[4]-=fpair*dx1*dx1;
            nrgy_strss_lcl[5]-=fpair*dx1*dx2;
            nrgy_strss_lcl[6]-=fpair*dx2*dx2;
            
            drhoi_dr[istart]=-drho_i_dr*r_inv;
            drhoj_dr[istart]=-drho_j_dr*r_inv;
            
            istart++;
        }
        p=rho[iatm]*drho_inv;
        m=static_cast<size_t>(p);
        m=MIN(m,nrho-2);
        p-=m;
        p=MIN(p,1.0);
        coef=F_arr[ielem][m];
        tmp1=(coef[6]*p+coef[5])*p+coef[4];
        tmp0=((coef[3]*p+coef[2])*p+coef[1])*p+coef[0];
        if(rho[iatm]>rho_max)
            tmp0+=tmp1*(rho[iatm]-rho_max);
        nrgy_strss_lcl[0]+=tmp0;
        rho[iatm]=tmp1;
    }
    
    dynamic->update(rho_ptr);

    
    istart=0;
    for(iatm=0;iatm<natms_lcl;iatm++)
    {
        ielem=evec[iatm];
        icomp=3*iatm;
        for(int j=0;j<neighbor_list_size[iatm];j++)
        {
            if(drhoi_dr[istart]!=0.0 || drhoj_dr[istart]!=0.0)
            {
                jatm=neighbor_list[iatm][j];
                jelem=evec[jatm];
                jcomp=3*jatm;
                
                fpair=rho[iatm]*drhoi_dr[istart]+rho[jatm]*drhoj_dr[istart];
                
                dx0=xvec[icomp]-xvec[jcomp];
                dx1=xvec[icomp+1]-xvec[jcomp+1];
                dx2=xvec[icomp+2]-xvec[jcomp+2];
                
                fvec[icomp]+=dx0*fpair;
                fvec[icomp+1]+=dx1*fpair;
                fvec[icomp+2]+=dx2*fpair;
                
                if(jatm<natms_lcl)
                {
                    fvec[jcomp]-=dx0*fpair;
                    fvec[jcomp+1]-=dx1*fpair;
                    fvec[jcomp+2]-=dx2*fpair;
                }
                
                if(jatm>=natms_lcl)
                    fpair*=0.5;
                
                nrgy_strss_lcl[1]-=fpair*dx0*dx0;
                nrgy_strss_lcl[2]-=fpair*dx0*dx1;
                nrgy_strss_lcl[3]-=fpair*dx0*dx2;
                nrgy_strss_lcl[4]-=fpair*dx1*dx1;
                nrgy_strss_lcl[5]-=fpair*dx1*dx2;
                nrgy_strss_lcl[6]-=fpair*dx2*dx2;
            }
            istart++;
        }
    }
}
/*--------------------------------------------
 energy calculation
 --------------------------------------------*/
void ForceFieldEAM::energy_calc()
{
    type0* xvec=atoms->x->begin();
    type0* rho=rho_ptr->begin();
    elem_type* evec=atoms->elem->begin();
    
    int iatm,jatm;
    int icomp,jcomp;
    elem_type ielem,jelem;
    type0 dx0,dx1,dx2,rsq;
    type0 r,p,phi,tmp0;
    size_t m;
    type0* coef;
    
    int** neighbor_list=neighbor->neighbor_list;
    int* neighbor_list_size=neighbor->neighbor_list_size;
    const int natms_lcl=atoms->natms_lcl;
    for(iatm=0;iatm<natms_lcl;iatm++) rho[iatm]=0.0;
    
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
            
            if(rsq<cut_sq[ielem][jelem])
            {
                r=sqrt(rsq);
                
                p=r*dr_inv;
                m=static_cast<size_t>(p);
                m=MIN(m,nr-2);
                p-=m;
                p=MIN(p,1.0);

                coef=r_phi_arr[ielem][jelem][m];
                phi=(((coef[3]*p+coef[2])*p+coef[1])*p+coef[0])/r;
                
                coef=rho_arr[jelem][ielem][m];
                rho[iatm]+=((coef[3]*p+coef[2])*p+coef[1])*p+coef[0];
                if(jatm<natms_lcl)
                {
                    coef=rho_arr[ielem][jelem][m];
                    rho[jatm]+=((coef[3]*p+coef[2])*p+coef[1])*p+coef[0];
                    nrgy_strss_lcl[0]+=phi;
                }
                else
                    nrgy_strss_lcl[0]+=0.5*phi;
                

            }
        }
        
        p=rho[iatm]*drho_inv;
        m=static_cast<size_t>(p);
        m=MIN(m,nrho-2);
        p-=m;
        p=MIN(p,1.0);
        coef=F_arr[ielem][m];
        tmp0=((coef[3]*p+coef[2])*p+coef[1])*p+coef[0];
        if(rho[iatm]>rho_max)
            tmp0+=((coef[6]*p+coef[5])*p+coef[4])*(rho[iatm]-rho_max);
        nrgy_strss_lcl[0]+=tmp0;

    }
}
/*--------------------------------------------
 pre gcmc energy claculate the increase or
 decrease in electron density
 --------------------------------------------*/
void ForceFieldEAM::pre_xchng_energy(GCMC* gcmc)
{
    int& icomm=gcmc->icomm;
    type0 en;
    type0 rho_iatm_lcl;
    
    int& iatm=gcmc->iatm;
    elem_type& ielem=gcmc->ielem;
    
    int& jatm=gcmc->jatm;
    elem_type& jelem=gcmc->jelem;
    
    type0&rsq=gcmc->rsq;
    
    type0 r,p;
    size_t m;
    type0* coef;
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
        if(gcmc->xchng_mode==DEL_MODE)
            c0=-1.0;
        
        en=rho_iatm_lcl=0.0;
        for(gcmc->reset_iatm();iatm!=-1;gcmc->next_iatm())
            for(gcmc->reset_jatm();jatm!=-1;gcmc->next_jatm())
            {
                r=sqrt(rsq);
                
                p=r*dr_inv;
                m=static_cast<size_t>(p);
                m=MIN(m,nr-2);
                p-=m;
                p=MIN(p,1.0);
                
                coef=r_phi_arr[ielem][jelem][m];
                
                if(jatm<natms_lcl)
                {
                    en+=(((coef[3]*p+coef[2])*p+coef[1])*p+coef[0])/r;
                    coef=rho_arr[ielem][jelem][m];
                    rho_xchng[jatm]+=c0*(((coef[3]*p+coef[2])*p+coef[1])*p+coef[0]);
                }
                else
                    en+=0.5*(((coef[3]*p+coef[2])*p+coef[1])*p+coef[0])/r;
                
                
                coef=rho_arr[jelem][ielem][m];
                rho_iatm_lcl+=((coef[3]*p+coef[2])*p+coef[1])*p+coef[0];
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
                p=tmp0*drho_inv;
                m=static_cast<size_t>(p);
                m=MIN(m,nrho-2);
                p-=m;
                p=MIN(p,1.0);
                coef=F_arr[evec[i]][m];
                F_xchng[i]=((coef[3]*p+coef[2])*p+coef[1])*p+coef[0];
                if(tmp0>rho_max)
                    F_xchng[i]+=((coef[6]*p+coef[5])*p+coef[4])*(tmp0-rho_max);
                en0+=F_xchng[i]-F[i];
            }
        
        en+=en0*c0;
        
        gcmc->lcl_vars[0]=en;
        gcmc->lcl_vars[1]=rho_iatm_lcl;
        
    }
}
/*--------------------------------------------
 calculate the energy if I am root if not
 pass the lcl variables
 --------------------------------------------*/
type0 ForceFieldEAM::xchng_energy(GCMC* gcmc)
{
    int& icomm=gcmc->icomm;
    for(gcmc->reset_icomm();icomm!=-1;gcmc->next_icomm())
        MPI_Reduce(gcmc->lcl_vars,gcmc->vars,2,Vec<type0>::MPI_T,MPI_SUM,gcmc->curr_root,*gcmc->curr_comm);
 
    
    if(gcmc->im_root)
    {
        //restart the comms
        gcmc->reset_icomm();
        type0 rho_iatm,en;
        size_t m;
        type0* coef;
        
        en=gcmc->vars[0];
        rho_iatm=gcmc->vars[1];
        type0 p=rho_iatm*drho_inv;
        m=static_cast<size_t>(p);
        m=MIN(m,nrho-2);
        p-=m;
        p=MIN(p,1.0);
        coef=F_arr[gcmc->ielem][m];
        type0 F_iatm=((coef[3]*p+coef[2])*p+coef[1])*p+coef[0];
        if(rho_iatm>rho_max)
            F_iatm+=((coef[6]*p+coef[5])*p+coef[4])*(rho_iatm-rho_max);
        
        en+=F_iatm;
        return en;
    }
    return 0.0;
}
/*--------------------------------------------
 calculate the energy if I am root if not
 pass the lcl variables
 --------------------------------------------*/
void ForceFieldEAM::post_xchng_energy(GCMC* gcmc)
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
        type0 p=rho[natms_lcl-1]*drho_inv;
        size_t m=static_cast<size_t>(p);
        m=MIN(m,nrho-2);
        p-=m;
        p=MIN(p,1.0);
        type0* coef=F_arr[gcmc->ielem][m];
        F[natms_lcl-1]=((coef[3]*p+coef[2])*p+coef[1])*p+coef[0];
        if(rho[natms_lcl-1]>rho_max)
            F[natms_lcl-1]+=((coef[6]*p+coef[5])*p+coef[4])*(rho[natms_lcl-1]-rho_max);
    }
}
/*--------------------------------------------
 init before running
 --------------------------------------------*/
void ForceFieldEAM::init()
{
    pre_init();
    rho_ptr=new Vec<type0>(atoms,1);
    
}
/*--------------------------------------------
 fin after running
 --------------------------------------------*/
void ForceFieldEAM::fin()
{
    Memory::dealloc(drhoj_dr);
    Memory::dealloc(drhoi_dr);
    max_pairs=0;
    delete rho_ptr;
    post_fin();
}
/*--------------------------------------------
 init xchng
 --------------------------------------------*/
void ForceFieldEAM::init_xchng()
{
    F_ptr=new Vec<type0>(atoms,1);
    rho_xchng_ptr=new Vec<type0>(atoms,1);
    F_xchng_ptr=new Vec<type0>(atoms,1);
    
    type0* xvec=atoms->x->begin();
    type0* rho=rho_ptr->begin();
    type0* F=F_ptr->begin();
    type0* rho_xchng=rho_xchng_ptr->begin();
    type0* F_xchng=F_xchng_ptr->begin();
    elem_type* evec=atoms->elem->begin();
    
    int iatm,jatm;
    int icomp,jcomp;
    elem_type ielem,jelem;
    type0 dx0,dx1,dx2,rsq;
    type0 r,p;
    size_t m;
    type0* coef;
    
    int** neighbor_list=neighbor->neighbor_list;
    int* neighbor_list_size=neighbor->neighbor_list_size;
    const int natms_lcl=atoms->natms_lcl;
    for(iatm=0;iatm<natms_lcl;iatm++) rho[iatm]=rho_xchng[iatm]=F_xchng[iatm]=0.0;
    
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
            
            if(rsq<cut_sq[ielem][jelem])
            {
                r=sqrt(rsq);
                
                p=r*dr_inv;
                m=static_cast<size_t>(p);
                m=MIN(m,nr-2);
                p-=m;
                p=MIN(p,1.0);
                
                coef=r_phi_arr[ielem][jelem][m];
                
                coef=rho_arr[jelem][ielem][m];
                rho[iatm]+=((coef[3]*p+coef[2])*p+coef[1])*p+coef[0];
                if(jatm<natms_lcl)
                {
                    coef=rho_arr[ielem][jelem][m];
                    rho[jatm]+=((coef[3]*p+coef[2])*p+coef[1])*p+coef[0];
                }
            }
        }
        
        p=rho[iatm]*drho_inv;
        m=static_cast<size_t>(p);
        m=MIN(m,nrho-2);
        p-=m;
        p=MIN(p,1.0);
        ielem=evec[iatm];
        coef=F_arr[ielem][m];
        F[iatm]=((coef[3]*p+coef[2])*p+coef[1])*p+coef[0];
    }
}
/*--------------------------------------------
 fin xchng
 --------------------------------------------*/
void ForceFieldEAM::fin_xchng()
{
    delete F_xchng_ptr;
    delete rho_xchng_ptr;
    delete F_ptr;
}

 
