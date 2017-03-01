/*--------------------------------------------
 Created by Sina on 07/16/13.
 Copyright (c) 2013 MIT. All rights reserved.
 --------------------------------------------*/
#include "ff.h"
#include "MAPP.h"
#include "elements.h"
#include "timer.h"
#include "memory.h"
#include "atoms.h"
#include "comm.h"
#include "xmath.h"
using namespace MAPP_NS;
/*--------------------------------------------
 constructor
 --------------------------------------------*/
ForceField::ForceField(Atoms* __atoms):
nelems(__atoms->elements->nelems),
world(__atoms->comm.world)
{
    Memory::alloc(cut,nelems,nelems);
    Memory::alloc(cut_sq,nelems,nelems);
    Memory::alloc(cut_sk_sq,nelems,nelems);
    Memory::alloc(rsq_crd,nelems);
    f=new Vec<type0>(__atoms,__atoms->x->dim);
}
/*--------------------------------------------
 destructor
 --------------------------------------------*/
ForceField::~ForceField()
{
    delete f;
    Memory::dealloc(rsq_crd);
    Memory::dealloc(cut_sk_sq);
    Memory::dealloc(cut_sq);
    Memory::dealloc(cut);
}
/*--------------------------------------------
 
 --------------------------------------------*/
void ForceField::reset()
{
    type0* __f=f->begin();
    const int n=f->vec_sz*f->dim;
    for(int i=0;i<n;i++) __f[i]=0.0;
    for(int i=0;i<__nvoigt__+1;i++) nrgy_strss_lcl[i]=0.0;
}


