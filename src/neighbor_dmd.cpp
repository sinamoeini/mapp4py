/*--------------------------------------------
 Created by Sina on 06/05/13.
 Copyright (c) 2013 MIT. All rights reserved.
 --------------------------------------------*/
#include "neighbor_dmd.h"
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
NeighborDMD::NeighborDMD(AtomsDMD* __atoms,
type0**& __cut_sk,type0*& __rsq_crd):
Neighbor(__atoms),
elem(__atoms->elem),
cut_sk(__cut_sk),
rsq_crd(__rsq_crd),
scl(__atoms->xi[__atoms->N-1]),
atoms(__atoms),
neighbor_list_2nd(NULL),
neighbor_list_size_2nd(NULL),
neighbor_list_size_size_2nd(0),
no_pairs_2nd(0)
{
}
/*--------------------------------------------
 destructor
 --------------------------------------------*/
NeighborDMD::~NeighborDMD()
{
}
/*--------------------------------------------
 
 --------------------------------------------*/
void NeighborDMD::mark_redndnt_ph(byte* mark)
{
    const int c_dim=atoms->c_dim;
    const int n=c_dim*atoms->natms_lcl;
    memset(mark,'0',atoms->natms_ph);
    for(int i=0;i<neighbor_list_size_size;i++)
        for(int j=0;j<neighbor_list_size[i];j++)
            if(neighbor_list[i][j]>=n)
                mark[(neighbor_list[i][j]-n)/c_dim]='1';
}
/*--------------------------------------------
 
 --------------------------------------------*/
void NeighborDMD::rename_atoms(int* old_2_new)
{
    const int c_dim=atoms->c_dim;
    for(int i=0;i<neighbor_list_size_size;i++)
        for(int j=0;j<neighbor_list_size[i];j++)
            neighbor_list[i][j]=old_2_new[neighbor_list[i][j]/c_dim]+(neighbor_list[i][j]%c_dim);
}
/*--------------------------------------------
 initiation before DMD
 --------------------------------------------*/
void NeighborDMD::init()
{
    no_neigh_lists=0;
    cell=new Cell<1>(atoms);
}
/*--------------------------------------------
 finalize DMD
 --------------------------------------------*/
void NeighborDMD::fin()
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
    
    if(neighbor_list_size_size_2nd)
    {
        for(int i=0;i<neighbor_list_size_size_2nd;i++)
            if(neighbor_list_size_2nd[i])
                delete [] neighbor_list_2nd[i];
        
        delete [] neighbor_list_2nd;
        delete [] neighbor_list_size_2nd;
    }
    neighbor_list_size_size_2nd=0;
    neighbor_list_2nd=NULL;
    neighbor_list_size_2nd=NULL;

    delete cell;    
}
/*--------------------------------------------
 create the neighbopr list
 s_or_x: 0 for x, 1 for s
 --------------------------------------------*/
void NeighborDMD::create_list(bool box_change)
{
    for(int i=0;i<neighbor_list_size_size;i++)
        delete [] neighbor_list[i];
    delete [] neighbor_list;
    delete [] neighbor_list_size;
    
    const int c_dim=atoms->c_dim;
    neighbor_list_size_size=atoms->natms_lcl*c_dim;
    
    Memory::alloc(neighbor_list,neighbor_list_size_size);
    Memory::alloc(neighbor_list_size,neighbor_list_size_size);
    for(int i=0;i<neighbor_list_size_size;i++)
        neighbor_list_size[i]=0;
    
    
    elem_type* elem_vec=elem->begin();
    const type0* c=atoms->c->begin();
    const type0* x=atoms->x->begin();
    const type0* alpha=atoms->alpha->begin();

    
    int** tmp_neigh_list;
    int* tmp_neigh_list_sz;
    Memory::alloc(tmp_neigh_list,c_dim,neighbor_list_size_size+(atoms->natms_ph)*c_dim);
    Memory::alloc(tmp_neigh_list_sz,c_dim);
    for(int ic=0;ic<c_dim;ic++) tmp_neigh_list_sz[ic]=0;
    
    auto store=[this,&tmp_neigh_list_sz,&tmp_neigh_list,&c,&c_dim](int i)->void
    {
        for(int ic=0;ic<c_dim;ic++)
        {
            neighbor_list[i*c_dim+ic]=NULL;
            if(tmp_neigh_list_sz[ic])
            {
                neighbor_list[i*c_dim+ic]=new int[tmp_neigh_list_sz[ic]];
                memcpy(neighbor_list[i*c_dim+ic],tmp_neigh_list[ic],tmp_neigh_list_sz[ic]*sizeof(int));
                neighbor_list_size[i*c_dim+ic]=tmp_neigh_list_sz[ic];
                no_pairs+=tmp_neigh_list_sz[ic];
                tmp_neigh_list_sz[ic]=0;
            }
        }
    };
    
    no_pairs=0;
    if(pair_wise)
    {
        cell->DO(box_change,[](int){},
        [&tmp_neigh_list,&tmp_neigh_list_sz,&elem_vec,&x,&alpha,&c,&c_dim,this](int i,int j)
        {
            if(i>=j) return;
            type0 r=sqrt(Algebra::RSQ<__dim__>(x+i*__dim__,x+j*__dim__));
            type0 const * c_i=c+i*c_dim;
            type0 const * c_j=c+j*c_dim;
            type0 const * alpha_i=alpha+i*c_dim;
            type0 const * alpha_j=alpha+j*c_dim;
            const elem_type* elem_i=elem_vec+i*c_dim;
            const elem_type* elem_j=elem_vec+j*c_dim;

            for(int ic=0;ic<c_dim&& c_i[ic]>=0.0;ic++)
                for(int jc=0;jc<c_dim&& c_j[jc]>=0.0;jc++)
                    if(r-scl*sqrt(alpha_i[ic]*alpha_i[ic]+alpha_j[jc]*alpha_j[jc])<cut_sk[elem_i[ic]][elem_j[jc]])
                        tmp_neigh_list[ic][tmp_neigh_list_sz[ic]++]=j*c_dim+jc;

        },store);
    }
    else
    {
        cell->DO(box_change,[](int){},
        [&tmp_neigh_list,&tmp_neigh_list_sz,&elem_vec,&x,&alpha,&c,&c_dim,this](int i,int j)
        {
            if(i==j) return;
            type0 r=sqrt(Algebra::RSQ<__dim__>(x+i*__dim__,x+j*__dim__));
            type0 const * c_i=c+i*c_dim;
            type0 const * c_j=c+j*c_dim;
            const elem_type* elem_i=elem_vec+i*c_dim;
            const elem_type* elem_j=elem_vec+j*c_dim;
            
            for(int ic=0;ic<c_dim&& c_i[ic]>=0.0;ic++)
                for(int jc=0;jc<c_dim&& c_j[jc]>=0.0;jc++)
                    if(r-scl*sqrt(alpha[ic]*alpha[ic]+alpha[jc]*alpha[jc])<cut_sk[elem_i[ic]][elem_j[jc]])
                        tmp_neigh_list[ic][tmp_neigh_list_sz[ic]++]=j*c_dim+jc;

        },store);
        no_pairs/=2;
    }
   
    
    Memory::dealloc(tmp_neigh_list);
    Memory::dealloc(tmp_neigh_list_sz);
    
    no_neigh_lists++;
}
/*--------------------------------------------
 initiation before DMD
 --------------------------------------------*/
void NeighborDMD::init_static()
{
    no_pairs_2nd=0;
}
/*--------------------------------------------
 finalize DMD
 --------------------------------------------*/
void NeighborDMD::fin_static()
{
    
    if(neighbor_list_size_size_2nd)
    {
        for(int i=0;i<neighbor_list_size_size_2nd;i++)
            delete [] neighbor_list_2nd[i];
        delete [] neighbor_list_2nd;
        delete [] neighbor_list_size_2nd;
    }
    neighbor_list_size_size_2nd=0;
    neighbor_list_2nd=NULL;
    neighbor_list_size_2nd=NULL;
    no_pairs_2nd=0;
}
/*--------------------------------------------
 create the neighbopr list
 s_or_x: 0 for x, 1 for s
 --------------------------------------------*/
void NeighborDMD::create_2nd_list()
{

    const int c_dim=atoms->c_dim;
    const int n=atoms->natms_lcl*c_dim;
    if(neighbor_list_size_size_2nd)
    {
        for(int i=0;i<neighbor_list_size_size_2nd;i++)
            delete [] neighbor_list_2nd[i];
        delete [] neighbor_list_2nd;
        delete [] neighbor_list_size_2nd;
    }
    neighbor_list_size_size_2nd=n;
    Memory::alloc(neighbor_list_2nd,neighbor_list_size_size_2nd);
    Memory::alloc(neighbor_list_size_2nd,neighbor_list_size_size_2nd);
    
    
    
    
    const type0* x=atoms->x->begin();
    elem_type* elem_vec=elem->begin();
    elem_type elem_i;
    type0 rsq_crd_ielem;
    
    int* tmp_neigh_list=NULL;
    Memory::alloc(tmp_neigh_list,atoms->natms_lcl+atoms->natms_ph);
    type0 x_i[__dim__];
    no_pairs_2nd=0;
    for(int i=0;i<n;i++)
    {
        elem_i=elem_vec[i];
        const int neigh_sz=neighbor_list_size[i];
        int __neigh_sz=0;
        Algebra::V_eq<__dim__>(x+(i/c_dim)*__dim__,x_i);
        rsq_crd_ielem=rsq_crd[elem_i];
        for(int j,__j=0;__j<neigh_sz;__j++)
        {
            j=neighbor_list[i][__j];
            if(elem_vec[j]!=elem_i) continue;
            if(Algebra::RSQ<__dim__>(x_i,x+(j/c_dim)*__dim__)>=rsq_crd_ielem) continue;
            tmp_neigh_list[__neigh_sz++]=j;
        }
        neighbor_list_2nd[i]=new int[__neigh_sz];
        memcpy(neighbor_list_2nd[i],tmp_neigh_list,__neigh_sz*sizeof(int));
        neighbor_list_size_2nd[i]=__neigh_sz;
        no_pairs_2nd+=__neigh_sz;
    }
    
    Memory::dealloc(tmp_neigh_list);
}

