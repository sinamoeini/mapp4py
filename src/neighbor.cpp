/*--------------------------------------------
 Created by Sina on 06/05/13.
 Copyright (c) 2013 MIT. All rights reserved.
 --------------------------------------------*/
#include "neighbor.h"
#include "ff.h"
#include "memory.h"
#include "comm.h"
#include "xmath.h"
using namespace MAPP_NS;
/*--------------------------------------------
 constructor
 --------------------------------------------*/
Neighbor::Neighbor(Atoms* atoms_):
atoms(atoms_),
pair_wise(true),
neighbor_list_size(NULL),
neighbor_list(NULL),
neighbor_list_size_size(0)
{
}
/*--------------------------------------------
 destructor
 --------------------------------------------*/
Neighbor::~Neighbor()
{
}
/*--------------------------------------------
 
 --------------------------------------------*/
void Neighbor::print_stats()
{
    /*-------temp_remove-------
    if(atoms->my_p==0)
    {
        fprintf(output,"neighbor stats:\n");
        fprintf(output,"no of neigh lists "
        "generated: %d\n",no_neigh_lists);
    }
     */
}
/*--------------------------------------------
 
 --------------------------------------------*/
void Neighbor::init()
{
    no_neigh_lists=0;
    cell=new Cell<1>(atoms);
}
/*--------------------------------------------
 
 --------------------------------------------*/
void Neighbor::fin()
{
    delete cell;
}
/*--------------------------------------------
 
 --------------------------------------------*/
void Neighbor::mark_redndnt_ph(byte* mark)
{
    int natms=atoms->natms;
    int natms_ph=atoms->natms_ph;
    memset(mark,'0',natms_ph);
    for(int iatm=0;iatm<neighbor_list_size_size;iatm++)
        for(int j=0;j<neighbor_list_size[iatm];j++)
            if(neighbor_list[iatm][j]>=natms)
                mark[neighbor_list[iatm][j]-natms]='1';
}
/*--------------------------------------------
 
 --------------------------------------------*/
void Neighbor::rename_atoms(int* old_2_new)
{
    for(int iatm=0;iatm<neighbor_list_size_size;iatm++)
        for(int j=0;j<neighbor_list_size[iatm];j++)
            neighbor_list[iatm][j]=old_2_new[neighbor_list[iatm][j]];
}
/*--------------------------------------------
 
 --------------------------------------------*/
/*
template<const int N,const int D,const int B>
class A
{
public:
    static constexpr int M=N/Algebra::__pow<B,D-1>::value;
    static constexpr int REM=N-M*Algebra::__pow<B,D-1>::value;
    static constexpr int L=M-(B-1)/2;
    static constexpr int V=(L==0 ? 0:(L>0 ? (L-1)*(L-1):(L+1)*(L+1)))+A<REM,D-1,B>::V;
};

template<const int N,const int B>
class A<N,1,B>
{
public:
    static constexpr int M=N;
    static constexpr int REM=N-M;
    static constexpr int L=M-(B-1)/2;
    static constexpr int V=(L==0 ? 0:(L>0 ? (L-1)*(L-1):(L+1)*(L+1)));
};

template<const int N,const int D,const int B>
class BB
{
public:
    static constexpr int value=((A<N,D,B>::V<((B-1)/2)*((B-1)/2))? 1:0)+BB<N-1,D,B>::value;
};
template<const int D,const int B>
class BB<0,D,B>
{
public:
    static constexpr int value=((A<0,D,B>::V<((B-1)/2)*((B-1)/2))? 1:0);
};
*/
