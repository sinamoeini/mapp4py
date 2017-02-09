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
using namespace MAPP_NS;
/*--------------------------------------------
 constructor
 --------------------------------------------*/
NeighborDMD::NeighborDMD(AtomsDMD* __atoms,
type0**& __cut_sk,type0*& __rsq_crd):
elem(__atoms->elem),
c_vec(__atoms->c),
cut_sk(__cut_sk),
rsq_crd(__rsq_crd),
Neighbor(__atoms),
neighbor_list_2nd(NULL),
neighbor_list_size_2nd(NULL),
neighbor_list_size_size_2nd(0),
scl(__atoms->xi[__atoms->N-1]),
atoms_dmd(__atoms)
{
}
/*--------------------------------------------
 destructor
 --------------------------------------------*/
NeighborDMD::~NeighborDMD()
{
}
/*--------------------------------------------
 initiation before DMD
 --------------------------------------------*/
void NeighborDMD::init()
{
    Neighbor::init();
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
    
    Neighbor::fin();
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
    
    const int c_dim=c_vec->dim;
    neighbor_list_size_size=atoms->natms*c_dim;
    
    Memory::alloc(neighbor_list,neighbor_list_size_size);
    Memory::alloc(neighbor_list_size,neighbor_list_size_size);
    for(int i=0;i<neighbor_list_size_size;i++)
        neighbor_list_size[i]=0;
    
    
    elem_type* elem_vec=elem->begin();
    const type0* c=c_vec->begin();
    const type0* x=atoms->x->begin();
    const type0* alpha=atoms_dmd->alpha->begin();
    const int x_dim=atoms->x->dim;
    

    
    int** tmp_neigh_list;
    int* tmp_neigh_list_sz;
    Memory::alloc(tmp_neigh_list,c_dim,neighbor_list_size_size+atoms->natms_ph);
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
        [&tmp_neigh_list,&tmp_neigh_list_sz,&elem_vec,&x,&alpha,&x_dim,&c,&c_dim,this](int i,int j)
        {
            if(i>=j) return;
            type0 r=sqrt(Algebra::RSQ<__dim__>(x+i*x_dim,x+j*x_dim));
            type0 const * c_i=c+i*c_dim;
            type0 const * c_j=c+j*c_dim;
            type0 const * alpha_i=alpha+i*c_dim;
            type0 const * alpha_j=alpha+j*c_dim;
            const elem_type* elem_i=elem_vec+i*c_dim;
            const elem_type* elem_j=elem_vec+j*c_dim;

            for(int ic=0;ic<c_dim&& c_i[ic]>0.0;ic++)
                for(int jc=0;jc<c_dim&& c_j[jc]>0.0;jc++)
                    if(r-scl*sqrt(alpha_i[ic]*alpha_i[ic]-alpha_j[jc]*alpha_j[jc])<cut_sk[elem_i[ic]][elem_j[jc]])
                        tmp_neigh_list[ic][tmp_neigh_list_sz[ic]++]=j*c_dim+jc;

        },store);
    }
    else
    {
        cell->DO(box_change,[](int){},
        [&tmp_neigh_list,&tmp_neigh_list_sz,&elem_vec,&x,&alpha,&x_dim,&c,&c_dim,this](int i,int j)
        {
            if(i==j) return;
            type0 r=sqrt(Algebra::RSQ<__dim__>(x+i*x_dim,x+j*x_dim));
            type0 const * c_i=c+i*c_dim;
            type0 const * c_j=c+j*c_dim;
            const elem_type* elem_i=elem_vec+i*c_dim;
            const elem_type* elem_j=elem_vec+j*c_dim;
            
            for(int ic=0;ic<c_dim&& c_i[ic]>0.0;ic++)
                for(int jc=0;jc<c_dim&& c_j[jc]>0.0;jc++)
                    if(r-scl*sqrt(alpha[ic]*alpha[ic]-alpha[jc]*alpha[jc])<cut_sk[elem_i[ic]][elem_j[jc]])
                        tmp_neigh_list[ic][tmp_neigh_list_sz[ic]++]=j*c_dim+jc;

        },store);
        no_pairs/=2;
    }
   
    
    Memory::dealloc(tmp_neigh_list);
    Memory::dealloc(tmp_neigh_list_sz);
    
    no_neigh_lists++;
}
/*--------------------------------------------
 create the neighbopr list
 s_or_x: 0 for x, 1 for s
 --------------------------------------------*/
void NeighborDMD::create_2nd_list()
{
    type0 rsq;
    int iatm,jatm,icomp,iicomp;
    const int natms=atoms->natms;    
    const type0* x=atoms->x->begin();
    const type0* c=c_vec->begin();
    const int c_dim=c_vec->dim;
    const int x_dim=atoms->x->dim;
    elem_type* elem_vec=elem->begin();
    
    
    
    int** tmp_neigh_list=NULL;
    Memory::alloc(tmp_neigh_list,c_dim,atoms->natms+atoms->natms_ph);
    
    
    neighbor_list_size_size_2nd=natms*c_dim;
    Memory::alloc(neighbor_list_2nd,neighbor_list_size_size_2nd);
    Memory::alloc(neighbor_list_size_2nd,neighbor_list_size_size_2nd);
    for(int i=0;i<neighbor_list_size_size_2nd;i++)
        neighbor_list_size_2nd[i]=0;
    
    const type0* xi=x;
    const type0* xj;
    const type0* ci=c;
    const type0* cj;
    const elem_type* ielem=elem_vec;
    const elem_type* jelem;
    for(iatm=0;iatm<natms;iatm++,xi+=x_dim,ci+=c_dim,ielem+=c_dim)
    {
        icomp=(__dim__+c_dim)*iatm;
        iicomp=c_dim*iatm;
        for(int j=0;j<neighbor_list_size[iatm];j++)
        {
            jatm=neighbor_list[iatm][j];
            xj=x+x_dim*jatm;
            cj=c+c_dim*jatm;
            jelem=elem_vec+c_dim*jatm;
            
            rsq=Algebra::RSQ<__dim__>(xi,xj);

            for(int ic_dim=0;ic_dim<c_dim;ic_dim++)
                for(int jc_dim=0;jc_dim<c_dim;jc_dim++)
                    if(ielem[ic_dim]==jelem[jc_dim]
                       && ci[ic_dim]>=0.0
                       && cj[jc_dim]>=0.0
                       && rsq<rsq_crd[ielem[ic_dim]])
                    {
                        tmp_neigh_list[ic_dim][neighbor_list_size_2nd[c_dim*iatm+ic_dim]]=c_dim*jatm+jc_dim;
                        neighbor_list_size_2nd[c_dim*iatm+ic_dim]++;
                    }
            
        }
        
        for(int ic_dim=0;ic_dim<c_dim;ic_dim++)
        {
            neighbor_list_2nd[c_dim*iatm+ic_dim]=NULL;
            if(neighbor_list_size_2nd[c_dim*iatm+ic_dim])
            {
                Memory::alloc(neighbor_list_2nd[iicomp+ic_dim],neighbor_list_size_2nd[c_dim*iatm+ic_dim]);
                memcpy(neighbor_list_2nd[c_dim*iatm+ic_dim],tmp_neigh_list[ic_dim],neighbor_list_size_2nd[c_dim*iatm+ic_dim]*sizeof(int));
            }

        }
    }
    
    Memory::dealloc(tmp_neigh_list);
    
    no_neigh_lists++;
}

