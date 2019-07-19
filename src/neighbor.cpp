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
//atoms(atoms_),
pair_wise(true),
neighbor_list(NULL),
neighbor_list_size(NULL),
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
