#ifndef __MAPP__ff_eam_potfit__
#define __MAPP__ff_eam_potfit__
#include "ff_md.h"
#include "xmath.h"
#include "potfit_funcs.h"
#include "atoms_md.h"
#include "neighbor_md.h"
#include "dynamic_md.h"
#include "memory.h"
namespace MAPP_NS
{

    
    template<size_t NELEMS>
    class ForceFieldEAMPotFit: public ForceFieldMD
    {
    private:
        
        type0 phi_calc(elem_type,elem_type,type0);
        type0 rho_calc(elem_type,elem_type,type0);
        type0 F_calc(elem_type,type0);
        type0 fpair_calc(elem_type,elem_type,type0,type0,type0);
        void reorder(std::string*,size_t);
        
        Vec<type0>* rho_vec_ptr;
        template<class T,size_t N>
        size_t find_uniq(T** orig,T**& uniq)
        {
            T* __uniq[N];
            memcpy(__uniq,orig,N*sizeof(T*));
            XMath::quicksort(__uniq,__uniq+N,
            [](T** rank_i,T** rank_j){return (*rank_i<*rank_j);},
            [](T** rank_i,T** rank_j){std::swap(*rank_i,*rank_j);});
            size_t sz=1;
            T* prev=*__uniq;
            for(size_t i=1;i<N;i++)
                if(__uniq[i]!=prev)
                {
                    prev=__uniq[i];
                    sz++;
                }

            Memory::alloc(uniq,sz);
            *uniq=prev=*__uniq;
            sz=1;
            for(size_t i=1;i<N;i++)
                if(__uniq[i]!=prev)
                {
                    uniq[sz]=__uniq[i];
                    prev=__uniq[i];
                    sz++;
                }
            
            return sz;
        }
        
    protected:
        void force_calc();
        void energy_calc();
        void pre_xchng_energy(GCMC*){};
        type0 xchng_energy(GCMC*){return 0.0;};
        void post_xchng_energy(GCMC*){};
    public:
        ForceFieldEAMPotFit(AtomsMD*,
        std::string(&)[NELEMS],
        type0*&,size_t,
        PotFitPairFunc*(&)[NELEMS][NELEMS],PotFitPairFunc*(&)[NELEMS][NELEMS],PotFitEmbFunc*(&)[NELEMS]);
        virtual ~ForceFieldEAMPotFit();
        void set_cutoff();
        void init();
        void fin();
        void init_xchng(){};
        void fin_xchng(){};
        void energy_gradient();
        void force_gradient();
        type0 mean_rho(elem_type);
        
        size_t nvs;
        type0* vs;
        type0* v0s;
        type0* dvs;
        type0* dvs_lcl;
        bool* dofs;
        type0* fvs;
        type0* f0vs;
        type0* hvs;
        type0* dvs_max;
        
        
        std::string names[NELEMS];
        PotFitPairFunc* rho_ptr[NELEMS][NELEMS];
        PotFitPairFunc* phi_ptr[NELEMS][NELEMS];
        PotFitEmbFunc* F_ptr[NELEMS];
        
        PotFitPairFunc** uniq_rho_ptr;
        PotFitPairFunc** uniq_phi_ptr;
        PotFitEmbFunc** uniq_F_ptr;
        size_t nuniq_phi_ptrs;
        size_t nuniq_rho_ptrs;
        size_t nuniq_F_ptrs;
        
    };
    
    class ForceFieldEAMPotFitAckOgata
    {
    private:
    protected:
    public:
        static constexpr size_t nelems=2;
        static ForceFieldEAMPotFit<nelems>* get_new_ff(class AtomsMD*,PyObject*,PyObject*);
        static void ml_new(PyMethodDef&);
    };
}
/*--------------------------------------------
 constructor
 --------------------------------------------*/
template<size_t NELEMS>
ForceFieldEAMPotFit<NELEMS>::ForceFieldEAMPotFit(AtomsMD* __atoms,
std::string (&__names)[NELEMS],
type0*& __vs,size_t __nvs,
PotFitPairFunc*(& __phi_ptr)[NELEMS][NELEMS],
PotFitPairFunc*(&__rho_ptr)[NELEMS][NELEMS],
PotFitEmbFunc*(&__F_ptr)[NELEMS]):
ForceFieldMD(__atoms),vs(__vs),nvs(__nvs),
uniq_rho_ptr(NULL),
uniq_phi_ptr(NULL),
uniq_F_ptr(NULL)
{
    memcpy(&phi_ptr[0][0],&__phi_ptr[0][0],NELEMS*NELEMS*sizeof(PotFitPairFunc*));
    memcpy(&rho_ptr[0][0],&__rho_ptr[0][0],NELEMS*NELEMS*sizeof(PotFitPairFunc*));
    memcpy(F_ptr,__F_ptr,NELEMS*sizeof(PotFitEmbFunc*));
    for(size_t i=0;i<NELEMS;i++) names[i]=__names[i];
    reorder(__atoms->elements.names, __atoms->elements.nelems);
    set_cutoff();
    Memory::alloc(v0s,nvs);
    Memory::alloc(dvs,nvs);
    Memory::alloc(dvs_lcl,nvs);
    Memory::alloc(dofs,nvs);
    Memory::alloc(fvs,nvs);
    Memory::alloc(f0vs,nvs);
    Memory::alloc(hvs,nvs);
    Memory::alloc(dvs_max,nvs);

#define POTFIT_OFFSET(A) \
offset=A->vars-vs; \
A->dvars_lcl=dvs_lcl+offset; \
A->dofs=dofs+offset; \
A->dvars_max=dvs_max+offset; \
A->hvars=hvs+offset
    
    ptrdiff_t offset;
    for(size_t i=0;i<NELEMS;i++)
    {
        POTFIT_OFFSET(F_ptr[i]);
        
        for(size_t j=0;j<NELEMS;j++)
        {
            POTFIT_OFFSET(phi_ptr[i][j]);
            POTFIT_OFFSET(rho_ptr[i][j]);
        }
    }
    
    for(size_t i=0;i<nvs;i++)
    {
        dofs[i]=true;
        dvs_max[i]=1.0;
    }
    
    
    
    
    nuniq_rho_ptrs=nuniq_phi_ptrs=nuniq_F_ptrs=0;
    nuniq_phi_ptrs=NELEMS*(NELEMS+1)/2;
    nuniq_rho_ptrs=nuniq_F_ptrs=NELEMS;
    Memory::alloc(uniq_phi_ptr,nuniq_phi_ptrs);
    Memory::alloc(uniq_rho_ptr,nuniq_rho_ptrs);
    Memory::alloc(uniq_F_ptr,nuniq_F_ptrs);
    
    nuniq_rho_ptrs=nuniq_phi_ptrs=nuniq_F_ptrs=0;
    for(size_t i=0;i<NELEMS;i++)
    {
        uniq_rho_ptr[nuniq_rho_ptrs++]=rho_ptr[i][0];
        uniq_F_ptr[nuniq_F_ptrs++]=F_ptr[i];
        
        for(size_t j=0;j<i+1;j++)
            uniq_phi_ptr[nuniq_phi_ptrs++]=phi_ptr[i][j];
    }
}
/*--------------------------------------------
 destructor
 --------------------------------------------*/
template<size_t NELEMS>
ForceFieldEAMPotFit<NELEMS>::~ForceFieldEAMPotFit()
{
    Memory::dealloc(dvs_max);
    Memory::dealloc(hvs);
    Memory::dealloc(f0vs);
    Memory::dealloc(fvs);
    Memory::dealloc(dofs);
    Memory::dealloc(dvs_lcl);
    Memory::dealloc(dvs);
    Memory::dealloc(v0s);
    Memory::dealloc(vs);
    for(size_t i=0;i<nuniq_phi_ptrs;i++) delete uniq_phi_ptr[i];
    for(size_t i=0;i<nuniq_rho_ptrs;i++) delete uniq_rho_ptr[i];
    for(size_t i=0;i<nuniq_F_ptrs;i++) delete uniq_F_ptr[i];
    *uniq_phi_ptr=NULL;
    *uniq_rho_ptr=NULL;
    *uniq_F_ptr=NULL;
    Memory::dealloc(uniq_phi_ptr);
    Memory::dealloc(uniq_rho_ptr);
    Memory::dealloc(uniq_F_ptr);
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<size_t NELEMS>
void ForceFieldEAMPotFit<NELEMS>::reorder(std::string* __names,size_t __nelems)
{
    if(__nelems>NELEMS) throw 1;
    
    size_t old_2_new[NELEMS];
    for(size_t i=0;i<NELEMS;i++) old_2_new[i]=i;
    for(size_t i=0;i<__nelems;i++)
    {
        bool fnd=false;
        for(size_t j=i;j<NELEMS&&!fnd;j++)
            if(strcmp(__names[i].c_str(),names[old_2_new[j]].c_str())==0)
            {
                std::swap(old_2_new[i],old_2_new[j]);
                fnd=true;
            }
        if(!fnd)
            throw 2;
    }
    
    
    PotFitPairFunc* __rho_ptr[NELEMS][NELEMS];
    PotFitPairFunc* __phi_ptr[NELEMS][NELEMS];
    PotFitEmbFunc* __F_ptr[NELEMS];
    std::string ___names[NELEMS];
    for(size_t i=0;i<NELEMS;i++)
    {
        ___names[i]=names[i];
        __F_ptr[i]=F_ptr[i];
        for(size_t j=0;j<NELEMS;j++)
        {
            __phi_ptr[i][j]=phi_ptr[i][j];
            __rho_ptr[i][j]=rho_ptr[i][j];
        }
    }
    
    
    for(size_t i=0;i<NELEMS;i++)
    {
        names[i]=___names[old_2_new[i]];
        F_ptr[i]=__F_ptr[old_2_new[i]];
        for(size_t j=0;j<NELEMS;j++)
        {
            phi_ptr[i][j]=__phi_ptr[old_2_new[i]][old_2_new[j]];
            rho_ptr[i][j]=__rho_ptr[old_2_new[i]][old_2_new[j]];
        }
    }
    
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<size_t NELEMS>
void ForceFieldEAMPotFit<NELEMS>::set_cutoff()
{
    for(size_t i=0;i<nelems;i++)
        for(size_t j=0;j<i+1;j++)
        {
            cut[i][j]=cut[j][i]=MAX(phi_ptr[i][j]->rc,MAX(rho_ptr[i][j]->rc,rho_ptr[j][i]->rc));
            cut_sq[i][j]=cut_sq[j][i]=cut[i][j]*cut[i][j];
        }
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<size_t NELEMS>
void ForceFieldEAMPotFit<NELEMS>::init()
{
    pre_init();
    rho_vec_ptr=new Vec<type0>(atoms,1,"rho");
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<size_t NELEMS>
void ForceFieldEAMPotFit<NELEMS>::fin()
{
    delete rho_vec_ptr;
    post_fin();
}
/*--------------------------------------------
 force and energy calculation
 --------------------------------------------*/
template<size_t NELEMS>
void ForceFieldEAMPotFit<NELEMS>::force_calc()
{
    type0* xvec=atoms->x->begin();
    type0* fvec=f->begin();
    type0* rho=rho_vec_ptr->begin();
    elem_type* evec=atoms->elem->begin();
    
    int iatm,jatm;
    elem_type ielem,jelem;
    type0 r,rsq;
    type0 fpair;
    type0 dx_ij[__dim__];
    int** neighbor_list=neighbor->neighbor_list;
    int* neighbor_list_size=neighbor->neighbor_list_size;
    const int natms_lcl=atoms->natms_lcl;
    for(int i=0;i<natms_lcl;i++)
        rho[i]=0.0;
    
    for(iatm=0;iatm<natms_lcl;iatm++)
    {
        ielem=evec[iatm];
        for(int j=0;j<neighbor_list_size[iatm];j++)
        {
            jatm=neighbor_list[iatm][j];
            jelem=evec[jatm];
            rsq=Algebra::RSQ<__dim__>(xvec+iatm*__dim__,xvec+jatm*__dim__);
            if(rsq>=cut_sq[ielem][jelem]) continue;
            r=sqrt(rsq);
            

            rho[iatm]+=rho_calc(jelem,ielem,r);
            if(jatm<natms_lcl)
            {
                rho[jatm]+=rho_calc(ielem,jelem,r);
                __vec_lcl[0]+=phi_calc(ielem,jelem,r);
            }
            else
                __vec_lcl[0]+=0.5*phi_calc(ielem,jelem,r);
        }
        __vec_lcl[0]+=F_calc(ielem,rho[iatm]);
    }
    
    
#ifdef NEW_UPDATE
    update(rho_vec_ptr);
#else
    dynamic->update(rho_vec_ptr);
#endif
    

    for(iatm=0;iatm<natms_lcl;iatm++)
    {
        ielem=evec[iatm];

        for(int j=0;j<neighbor_list_size[iatm];j++)
        {
            jatm=neighbor_list[iatm][j];
            jelem=evec[jatm];

            rsq=Algebra::DX_RSQ<__dim__>(xvec+iatm*__dim__,xvec+jatm*__dim__,dx_ij);
            if(rsq>=cut_sq[ielem][jelem]) continue;
            
            r=sqrt(rsq);
            
            
            fpair=fpair_calc(ielem,jelem,rho[iatm],rho[jatm],r);
            
            if(fpair==0.0) continue;
            
            Algebra::V_add_x_mul_V<__dim__>(fpair,dx_ij,fvec+iatm*__dim__);
            if(jatm<natms_lcl)
                Algebra::V_add_x_mul_V<__dim__>(-fpair,dx_ij,fvec+jatm*__dim__);
            else
                fpair*=0.5;
            
            Algebra::DyadicV<__dim__>(-fpair,dx_ij,&__vec_lcl[1]);
        }
    }
    type0 f_sum_lcl[__dim__];
    Algebra::zero<__dim__>(f_sum_lcl);
    for(int i=0;i<natms_lcl;i++)
        Algebra::V_add<__dim__>(fvec+__dim__*i,f_sum_lcl);

    type0 f_corr[__dim__];
    MPI_Allreduce(f_sum_lcl,f_corr,__dim__,Vec<type0>::MPI_T,MPI_SUM,world);
    type0 a=-1.0/static_cast<type0>(atoms->natms);
    Algebra::Do<__dim__>::func([&a,&f_corr](int i){f_corr[i]*=a;});
    for(int i=0;i<natms_lcl;i++)
        Algebra::V_add<__dim__>(f_corr,fvec+__dim__*i);
}
/*--------------------------------------------
 energy calculation
 --------------------------------------------*/
template<size_t NELEMS>
void ForceFieldEAMPotFit<NELEMS>::energy_calc()
{
    type0* xvec=atoms->x->begin();
    type0* rho=rho_vec_ptr->begin();
    elem_type* evec=atoms->elem->begin();
    
    int iatm,jatm;
    elem_type ielem,jelem;
    type0 r,rsq;
    
    int** neighbor_list=neighbor->neighbor_list;
    int* neighbor_list_size=neighbor->neighbor_list_size;
    const int natms_lcl=atoms->natms_lcl;
    for(int i=0;i<natms_lcl;i++)
        rho[i]=0.0;
    
    for(iatm=0;iatm<natms_lcl;iatm++)
    {
        ielem=evec[iatm];
        for(int j=0;j<neighbor_list_size[iatm];j++)
        {
            jatm=neighbor_list[iatm][j];
            jelem=evec[jatm];
            rsq=Algebra::RSQ<__dim__>(xvec+iatm*__dim__,xvec+jatm*__dim__);
            if(rsq>=cut_sq[ielem][jelem]) continue;
            r=sqrt(rsq);
            
            
            rho[iatm]+=rho_calc(jelem,ielem,r);
            if(jatm<natms_lcl)
            {
                rho[jatm]+=rho_calc(ielem,jelem,r);
                __vec_lcl[0]+=phi_calc(ielem,jelem,r);
            }
            else
                __vec_lcl[0]+=0.5*phi_calc(ielem,jelem,r);
        }
        // add the embedded energy here
        __vec_lcl[0]+=F_calc(ielem,rho[iatm]);
    }
}
/*--------------------------------------------
 energy calculation
 --------------------------------------------*/
template<size_t NELEMS>
void ForceFieldEAMPotFit<NELEMS>::energy_gradient()
{
    for(size_t i=0;i<nvs;i++) dvs_lcl[i]=0.0;
    type0* xvec=atoms->x->begin();
    type0* rho=rho_vec_ptr->begin();
    elem_type* evec=atoms->elem->begin();
    
    int iatm,jatm;
    elem_type ielem,jelem;
    type0 r,rsq,coef;
    
    int** neighbor_list=neighbor->neighbor_list;
    int* neighbor_list_size=neighbor->neighbor_list_size;
    const int natms_lcl=atoms->natms_lcl;
    for(int i=0;i<natms_lcl;i++)
        rho[i]=0.0;
    
    for(iatm=0;iatm<natms_lcl;iatm++)
    {
        ielem=evec[iatm];
        for(int j=0;j<neighbor_list_size[iatm];j++)
        {
            jatm=neighbor_list[iatm][j];
            jelem=evec[jatm];
            rsq=Algebra::RSQ<__dim__>(xvec+iatm*__dim__,xvec+jatm*__dim__);
            if(rsq>=cut_sq[ielem][jelem]) continue;
            r=sqrt(rsq);
            coef=jatm<natms_lcl ? 1.0:0.5;
            phi_ptr[ielem][jelem]->DF(coef,r);
            
            rho[iatm]+=rho_calc(jelem,ielem,r);
            if(jatm<natms_lcl)
            {
                rho[jatm]+=rho_calc(ielem,jelem,r);
                __vec_lcl[0]+=phi_calc(ielem,jelem,r);
            }
            else
                __vec_lcl[0]+=0.5*phi_calc(ielem,jelem,r);
        }

        __vec_lcl[0]+=F_calc(ielem,rho[iatm]);
        F_ptr[ielem]->DF(rho[iatm]);
    }
    
    type0 dFj,dFi;
    
    for(iatm=0;iatm<natms_lcl;iatm++)
    {
        ielem=evec[iatm];
        dFi=F_ptr[ielem]->dF(rho[iatm]);
        for(int j=0;j<neighbor_list_size[iatm];j++)
        {
            jatm=neighbor_list[iatm][j];
            jelem=evec[jatm];
            rsq=Algebra::RSQ<__dim__>(xvec+iatm*__dim__,xvec+jatm*__dim__);
            if(rsq>=cut_sq[ielem][jelem]) continue;
            r=sqrt(rsq);
            
            rho_ptr[jelem][ielem]->DF(dFi,r);
            if(jatm<natms_lcl)
            {
                dFj=F_ptr[jelem]->dF(rho[jatm]);
                rho_ptr[ielem][jelem]->DF(dFj,r);
            }
        }
    }
    MPI_Allreduce(dvs_lcl,dvs,static_cast<int>(nvs),Vec<type0>::MPI_T,MPI_SUM,world);
    MPI_Allreduce(__vec_lcl,__vec,1,Vec<type0>::MPI_T,MPI_SUM,world);
    for(size_t i=0;i<nvs;i++)
        if(!dofs[i]) dvs[i]=0.0;
}
/*--------------------------------------------
 force_norm gradient 
 --------------------------------------------*/
template<size_t NELEMS>
void ForceFieldEAMPotFit<NELEMS>::force_gradient()
{
    for(size_t i=0;i<nvs;i++) dvs_lcl[i]=0.0;
    type0* xvec=atoms->x->begin();
    type0* rho=rho_vec_ptr->begin();
    elem_type* evec=atoms->elem->begin();
    
    int iatm,jatm;
    elem_type ielem,jelem;
    type0 r,rsq,coef;
    
    int** neighbor_list=neighbor->neighbor_list;
    int* neighbor_list_size=neighbor->neighbor_list_size;
    const int natms_lcl=atoms->natms_lcl;
    for(int i=0;i<natms_lcl;i++)
        rho[i]=0.0;
    
    for(iatm=0;iatm<natms_lcl;iatm++)
    {
        ielem=evec[iatm];
        for(int j=0;j<neighbor_list_size[iatm];j++)
        {
            jatm=neighbor_list[iatm][j];
            jelem=evec[jatm];
            rsq=Algebra::RSQ<__dim__>(xvec+iatm*__dim__,xvec+jatm*__dim__);
            if(rsq>=cut_sq[ielem][jelem]) continue;
            r=sqrt(rsq);
            coef=jatm<natms_lcl ? 1.0:0.5;
            phi_ptr[ielem][jelem]->DF(coef,r);
            
            rho[iatm]+=rho_calc(jelem,ielem,r);
            if(jatm<natms_lcl)
            {
                rho[jatm]+=rho_calc(ielem,jelem,r);
                __vec_lcl[0]+=phi_calc(ielem,jelem,r);
            }
            else
                __vec_lcl[0]+=0.5*phi_calc(ielem,jelem,r);
        }

        __vec_lcl[0]+=F_calc(ielem,rho[iatm]);
        F_ptr[ielem]->DF(rho[iatm]);
    }
    
    type0 dFj,dFi;
    
    for(iatm=0;iatm<natms_lcl;iatm++)
    {
        ielem=evec[iatm];
        dFi=F_ptr[ielem]->dF(rho[iatm]);
        for(int j=0;j<neighbor_list_size[iatm];j++)
        {
            jatm=neighbor_list[iatm][j];
            jelem=evec[jatm];
            rsq=Algebra::RSQ<__dim__>(xvec+iatm*__dim__,xvec+jatm*__dim__);
            if(rsq>=cut_sq[ielem][jelem]) continue;
            r=sqrt(rsq);
            
            rho_ptr[jelem][ielem]->DF(dFi,r);
            if(jatm<natms_lcl)
            {
                dFj=F_ptr[jelem]->dF(rho[jatm]);
                rho_ptr[ielem][jelem]->DF(dFj,r);
            }
        }
    }
    MPI_Allreduce(dvs_lcl,dvs,static_cast<int>(nvs),Vec<type0>::MPI_T,MPI_SUM,world);
    MPI_Allreduce(__vec_lcl,__vec,1,Vec<type0>::MPI_T,MPI_SUM,world);
    for(size_t i=0;i<nvs;i++)
        if(!dofs[i]) dvs[i]=0.0;
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<size_t NELEMS>
type0 ForceFieldEAMPotFit<NELEMS>::phi_calc(elem_type ielem,elem_type jelem,type0 r)
{
    return phi_ptr[ielem][jelem]->F(r);
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<size_t NELEMS>
type0 ForceFieldEAMPotFit<NELEMS>::mean_rho(elem_type ielem)
{
    elem_type* evec=atoms->elem->begin();
    type0* rho=rho_vec_ptr->begin();
    const int natms_lcl=atoms->natms_lcl;
    type0 __tmp_lcl[2]={0.0,0.0};
    type0 __tmp[2]={0.0,0.0};
    for(int i=0;i<natms_lcl;i++)
        if(evec[i]==ielem)
        {
            __tmp_lcl[0]+=rho[i];
            __tmp_lcl[1]++;
        }
    MPI_Allreduce(__tmp_lcl,__tmp,2,Vec<type0>::MPI_T,MPI_SUM,world);
    if(__tmp[1]==0.0) return 0.0;
    return __tmp[0]/__tmp[1];
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<size_t NELEMS>
type0 ForceFieldEAMPotFit<NELEMS>::rho_calc(elem_type ielem,elem_type jelem,type0 r)
{
    return rho_ptr[ielem][jelem]->F(r);
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<size_t NELEMS>
type0 ForceFieldEAMPotFit<NELEMS>::F_calc(elem_type ielem,type0 rho)
{
    return F_ptr[ielem]->F(rho);
    
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<size_t NELEMS>
type0 ForceFieldEAMPotFit<NELEMS>::fpair_calc(elem_type ielem,elem_type jelem,type0 rho_i,type0 rho_j,type0 r)
{
    return -(phi_ptr[ielem][jelem]->dF(r)+
             F_ptr[ielem]->dF(rho_i)*rho_ptr[jelem][ielem]->dF(r)+
             F_ptr[jelem]->dF(rho_j)*rho_ptr[ielem][jelem]->dF(r))/r;
}


#endif


