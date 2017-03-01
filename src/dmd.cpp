#include "dmd.h"
#include "atoms_dmd.h"
using namespace MAPP_NS;
/*--------------------------------------------
 
 --------------------------------------------*/
DMD::DMD():
max_nsteps(1000),
a_tol(sqrt(std::numeric_limits<type0>::epsilon())),
min_dt(std::numeric_limits<type0>::epsilon()),
c(NULL),
c_d(NULL)
{
}
/*--------------------------------------------
 
 --------------------------------------------*/
DMD::~DMD()
{
}
/*--------------------------------------------
 
 --------------------------------------------*/
void DMD::init_static()
{
    c_dim=atoms->c_dim;
    ncs=atoms->natms*c_dim;
    c=atoms->c->begin();
    c_d=atoms->c_d->begin();
}
/*--------------------------------------------
 
 --------------------------------------------*/
void DMD::fin_static()
{
    c=c_d=NULL;
}
