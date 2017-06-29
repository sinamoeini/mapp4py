/*--------------------------------------------
 Created by Sina on 06/05/13.
 Copyright (c) 2013 MIT. All rights reserved.
 --------------------------------------------*/
#include "neighbor_dmd_sc.h"
#include "ff_dmd.h"
#include "atoms_dmd.h"
#include "memory.h"
#include "MAPP.h"
#include "xmath.h"
#include "cell.h"
using namespace MAPP_NS;
/*--------------------------------------------
 constructor
 --------------------------------------------*/
NeighborDMDSC::NeighborDMDSC(AtomsDMD* __atoms,
type0**& __cut_sk,type0*& __rsq_crd):
NeighborDMD(__atoms,__cut_sk,__rsq_crd)
{
}
/*--------------------------------------------
 destructor
 --------------------------------------------*/
NeighborDMDSC::~NeighborDMDSC()
{
}
/*--------------------------------------------
 
 --------------------------------------------*/
void NeighborDMDSC::mark_redndnt_ph(byte* mark)
{
    int natms_lcl=atoms->natms_lcl;
    int natms_ph=atoms->natms_ph;
    memset(mark,'0',natms_ph);
    for(int iatm=0;iatm<neighbor_list_size_size;iatm++)
        for(int j=0;j<neighbor_list_size[iatm];j++)
            if(neighbor_list[iatm][j]>=natms_lcl)
                mark[neighbor_list[iatm][j]-natms_lcl]='1';
}
/*--------------------------------------------
 
 --------------------------------------------*/
void NeighborDMDSC::rename_atoms(int* old_2_new)
{
    for(int iatm=0;iatm<neighbor_list_size_size;iatm++)
        for(int j=0;j<neighbor_list_size[iatm];j++)
            neighbor_list[iatm][j]=old_2_new[neighbor_list[iatm][j]];
}
/*--------------------------------------------
 initiation before DMD
 --------------------------------------------*/
void NeighborDMDSC::init()
{
    no_neigh_lists=0;
    cell=new Cell<1>(atoms);
}
/*--------------------------------------------
 finalize DMD
 --------------------------------------------*/
void NeighborDMDSC::fin()
{
    if(neighbor_list_size_size)
    {
        for(int i=0;i<neighbor_list_size_size;i++)
            delete [] neighbor_list[i];
        delete [] neighbor_list;
        delete [] neighbor_list_size;
    }
    neighbor_list_size_size=0;
    neighbor_list=NULL;
    neighbor_list_size=NULL;

    delete cell;    
}
/*--------------------------------------------
 create the neighbopr list
 s_or_x: 0 for x, 1 for s
 --------------------------------------------*/
void NeighborDMDSC::create_list(bool box_change)
{
    for(int i=0;i<neighbor_list_size_size;i++)
        delete [] neighbor_list[i];
    delete [] neighbor_list;
    delete [] neighbor_list_size;
    
    const int c_dim=atoms->c_dim;
    neighbor_list_size_size=atoms->natms_lcl;
    
    Memory::alloc(neighbor_list,neighbor_list_size_size);
    Memory::alloc(neighbor_list_size,neighbor_list_size_size);
    for(int i=0;i<neighbor_list_size_size;i++)
        neighbor_list_size[i]=0;
    
    
    elem_type* elem_vec=elem->begin();
    const type0* c=atoms->c->begin();
    const type0* x=atoms->x->begin();
    const type0* alpha=atoms->alpha->begin();
    type0 max_cut=atoms->max_cut;
    
    int* tmp_neigh_list;
    int tmp_neigh_list_sz=0;
    Memory::alloc(tmp_neigh_list,neighbor_list_size_size+(atoms->natms_ph));
    
    auto store=[this,&tmp_neigh_list_sz,&tmp_neigh_list,&c,&c_dim](int i)->void
    {
        neighbor_list[i]=NULL;
        if(tmp_neigh_list_sz)
        {
            neighbor_list[i]=new int[tmp_neigh_list_sz];
            memcpy(neighbor_list[i],tmp_neigh_list,tmp_neigh_list_sz*sizeof(int));
            neighbor_list_size[i]=tmp_neigh_list_sz;
            no_pairs+=tmp_neigh_list_sz;
            tmp_neigh_list_sz=0;
        }
    };
    
    no_pairs=0;
    if(pair_wise)
    {
        type0 max_alpha_i=0.0;
        cell->DO(box_change,[&max_alpha_i,this,&c_dim,&c,&alpha](int i)
        {
            type0 const * c_i=c+i*c_dim;
            type0 const * alpha_i=alpha+i*c_dim;
            for(int ic=0;ic<c_dim&& c_i[ic]>=0.0;ic++)
                max_alpha_i=MAX(max_alpha_i,alpha_i[ic]);
                
        }
                 ,
        [&tmp_neigh_list,&tmp_neigh_list_sz,&elem_vec,&x,&alpha,&c,&c_dim,&max_alpha_i,&max_cut,this](int i,int j)
        {
            if(i>=j) return;
            type0 r=sqrt(Algebra::RSQ<__dim__>(x+i*__dim__,x+j*__dim__));
            type0 const * c_j=c+j*c_dim;
            type0 const * alpha_j=alpha+j*c_dim;
            
            type0 max_alpha_j=0.0;
            for(int jc=0;jc<c_dim&& c_j[jc]>=0.0;jc++)
                max_alpha_j=MAX(max_alpha_j,alpha_j[jc]);
            if(r-scl*sqrt(max_alpha_i*max_alpha_i+max_alpha_j*max_alpha_j)<max_cut)
                tmp_neigh_list[tmp_neigh_list_sz++]=j;



        },store);
    }
    else
    {
        type0 max_alpha_i=0.0;
        cell->DO(box_change,[&max_alpha_i,this,&c_dim,&c,&alpha](int i)
        {
            type0 const * c_i=c+i*c_dim;
            type0 const * alpha_i=alpha+i*c_dim;
            for(int ic=0;ic<c_dim&& c_i[ic]>=0.0;ic++)
                max_alpha_i=MAX(max_alpha_i,alpha_i[ic]);
                
        }
                 ,
        [&tmp_neigh_list,&tmp_neigh_list_sz,&elem_vec,&x,&alpha,&c,&c_dim,&max_alpha_i,&max_cut,this](int i,int j)
        {
            if(i==j) return;
            type0 r=sqrt(Algebra::RSQ<__dim__>(x+i*__dim__,x+j*__dim__));
            type0 const * c_j=c+j*c_dim;
            type0 const * alpha_j=alpha+j*c_dim;
            
            type0 max_alpha_j=0.0;
            for(int jc=0;jc<c_dim&& c_j[jc]>=0.0;jc++)
                max_alpha_j=MAX(max_alpha_j,alpha_j[jc]);
            if(r-scl*sqrt(max_alpha_i*max_alpha_i+max_alpha_j*max_alpha_j)<max_cut)
                tmp_neigh_list[tmp_neigh_list_sz++]=j;



        },store);
        no_pairs/=2;
    }
   
    
    Memory::dealloc(tmp_neigh_list);
    
    no_neigh_lists++;
}
