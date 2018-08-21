/*--------------------------------------------
 Created by Sina on 07/15/13.
 Copyright (c) 2013 MIT. All rights reserved.
 --------------------------------------------*/
#include "ff_meam.h"
#include "neighbor_md.h"
#include "atoms_md.h"
#include "elements.h"
#include "memory.h"
#include "pgcmc.h"
#include "api.h"
using namespace MAPP_NS;
/*--------------------------------------------
 constructor
 --------------------------------------------*/
ForceFieldMEAM::
ForceFieldMEAM(AtomsMD* atoms):
ForceFieldMD(atoms)
{
    gcmc_n_cutoff=1;
    gcmc_n_vars=1;
    gcmc_tag_enabled=false;
    neighbor->pair_wise=false;
    
}
/*--------------------------------------------
 destructor
 --------------------------------------------*/
ForceFieldMEAM::~ForceFieldMEAM()
{
}
/*--------------------------------------------
 initiate before a run
 --------------------------------------------*/
void ForceFieldMEAM::init()
{
    pre_init();
}
/*--------------------------------------------
 after a run
 --------------------------------------------*/
void ForceFieldMEAM::fin()
{
    post_fin();
}
/*--------------------------------------------
 init xchng
 --------------------------------------------*/
void ForceFieldMEAM::init_xchng()
{

}
/*--------------------------------------------
 fin xchng
 --------------------------------------------*/
void ForceFieldMEAM::fin_xchng()
{

}
/*--------------------------------------------
 pre xchng energy
 --------------------------------------------*/
void ForceFieldMEAM::pre_xchng_energy(GCMC* gcmc)
{
    
}
/*--------------------------------------------
 xchng energy
 --------------------------------------------*/
type0 ForceFieldMEAM::xchng_energy(GCMC* gcmc)
{
    return 0.0;
}
/*--------------------------------------------
 post xchng energy
 --------------------------------------------*/
void ForceFieldMEAM::post_xchng_energy(GCMC*)
{
}
/*--------------------------------------------
 force and energy calculation
 --------------------------------------------*/
void ForceFieldMEAM::__force_calc()
{

    
    
    
    
    
}
/*--------------------------------------------
 only energy calculation this is useful for
 minimization/linesearch methods that do not
 use derivatives of energy
 --------------------------------------------*/
void ForceFieldMEAM::__energy_calc()
{
}
/*--------------------------------------------
 python constructor
 --------------------------------------------*/
void ForceFieldMEAM::ml_new(PyMethodDef& tp_methods)
{

}
