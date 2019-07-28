/*--------------------------------------------
 Created by Sina on 06/05/13.
 Copyright (c) 2013 MIT. All rights reserved.
 --------------------------------------------*/
#include "elements.h"
#include "mpi_compat.h"
#include "gcmc.h"
#include "memory.h"
#include "random.h"
#include "neighbor.h"
#include "ff_md.h"
#include "MAPP.h"
#include "atoms_md.h"
#include "comm.h"
#include "dynamic_md.h"
using namespace MAPP_NS;
/*--------------------------------------------
 constructor
 --------------------------------------------*/
GCMC::GCMC(AtomsMD*& __atoms,ForceFieldMD*&__ff,DynamicMD*& __dynamic,elem_type __gas_type,type0 __mu,type0 __T,int seed):
gas_type(__gas_type),
T(__T),
mu(__mu),
natms_lcl(__atoms->natms_lcl),
natms_ph(__atoms->natms_ph),
cut_sq(__ff->cut_sq),
s_lo(__atoms->comm.s_lo),
s_hi(__atoms->comm.s_hi),
dynamic(__dynamic),
world(__atoms->comm.world),
atoms(__atoms),
ff(__ff),
ielem(gas_type)
{
    random=new Random(seed);
    s_trials=new type0*[__dim__];
    *s_trials=NULL;
    del_ids=NULL;
    del_ids_sz=del_ids_cpcty=0;
    vars=lcl_vars=NULL;
}
/*--------------------------------------------
 destructor
 --------------------------------------------*/
GCMC::~GCMC()
{
    delete [] del_ids;
    delete [] s_trials;
    delete random;
}
/*--------------------------------------------
 
 --------------------------------------------*/
void GCMC::add_del_id(int* new_ids,int no)
{
    if(del_ids_sz+no>del_ids_cpcty)
    {
        int* del_ids_=new int[del_ids_sz+no];
        memcpy(del_ids_,del_ids,del_ids_sz*sizeof(int));
        del_ids_cpcty=del_ids_sz+no;
        delete [] del_ids;
        del_ids=del_ids_;
    }
    memcpy(del_ids+del_ids_sz,new_ids,sizeof(int)*no);
    del_ids_sz+=no;
}
/*--------------------------------------------
 
 --------------------------------------------*/
int GCMC::get_new_id()
{
    if(del_ids_sz)
    {
        del_ids_sz--;
        return del_ids[del_ids_sz];
    }
    else
    {
        max_id++;
        return max_id;
    }
}
/*--------------------------------------------
 
 --------------------------------------------*/
void GCMC::init()
{
    cut=ff->cut[ielem][0];
    for(size_t i=1;i<atoms->elements.nelems;i++)
        cut=MAX(cut,ff->cut[ielem][i]);
    
    gas_mass=atoms->elements.masses[gas_type];
    kbT=atoms->kB*T;
    beta=1.0/kbT;
    lambda=atoms->hP/sqrt(2.0*M_PI*kbT*gas_mass);
    sigma=sqrt(kbT/gas_mass);
    z_fac=1.0;
    for(int i=0;i<__dim__;i++) z_fac/=lambda;
    z_fac*=exp(beta*mu);
    vol=1.0;
    for(int i=0;i<__dim__;i++)vol*=atoms->H[i][i];
    
    id_type max_id_=0;
    id_type* id=atoms->id->begin();
    for(int i=0;i<natms_lcl;i++)
        max_id_=MAX(id[i],max_id_);
    MPI_Allreduce(&max_id_,&max_id,1,Vec<id_type>::MPI_T,MPI_MAX,world);
    for(int i=0;i<del_ids_sz;i++)
        max_id=MAX(max_id,del_ids[i]);
        
    ngas_lcl=0;
    elem_type* elem=atoms->elem->begin();
    for(int i=0;i<natms_lcl;i++)
        if(elem[i]==gas_type) ngas_lcl++;
    MPI_Allreduce(&ngas_lcl,&ngas,1,MPI_INT,MPI_SUM,world);
}
/*--------------------------------------------
 
 --------------------------------------------*/
void GCMC::fin()
{
    
}
/*--------------------------------------------
 
 --------------------------------------------*/
void GCMC::box_setup()
{
    int sz=0;
    max_ntrial_atms=1;
    for(int i=0;i<__dim__;i++)
    {
        type0 tmp=0.0;
        for(int j=i;j<__dim__;j++)
            tmp+=atoms->B[j][i]*atoms->B[j][i];
        cut_s[i]=sqrt(tmp)*cut;
        

        
        s_lo_ph[i]=s_lo[i]-cut_s[i];
        s_hi_ph[i]=s_hi[i]+cut_s[i];
        nimages_per_dim[i][0]=static_cast<int>(floor(s_hi_ph[i]));
        nimages_per_dim[i][1]=-static_cast<int>(floor(s_lo_ph[i]));
        max_ntrial_atms*=1+nimages_per_dim[i][0]+nimages_per_dim[i][1];
        sz+=1+nimages_per_dim[i][0]+nimages_per_dim[i][1];

    }
    
    *s_trials=new type0[sz];
    for(int i=1;i<__dim__;i++)
        s_trials[i]=s_trials[i-1]+1+nimages_per_dim[i-1][0]+nimages_per_dim[i-1][1];
}
/*--------------------------------------------
 
 --------------------------------------------*/
void GCMC::box_dismantle()
{
    delete [] *s_trials;
    *s_trials=NULL;
}



