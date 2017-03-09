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
nelems(__atoms->elements.nelems),
world(__atoms->comm.world)
{
    Memory::alloc(cut,nelems,nelems);
    Memory::alloc(cut_sq,nelems,nelems);
    Memory::alloc(cut_sk_sq,nelems,nelems);
    Memory::alloc(rsq_crd,nelems);
}
/*--------------------------------------------
 destructor
 --------------------------------------------*/
ForceField::~ForceField()
{
    Memory::dealloc(rsq_crd);
    Memory::dealloc(cut_sk_sq);
    Memory::dealloc(cut_sq);
    Memory::dealloc(cut);
}



