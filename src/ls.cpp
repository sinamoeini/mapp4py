#include "min_cg.h"
#include "min_cg_dmd.h"
#include "ls.h"
/*--------------------------------------------
 constructor
 --------------------------------------------*/
LineSearch::LineSearch():
golden(0.5+0.5*sqrt(5.0)),
epsilon_3_4(pow(std::numeric_limits<type0>::epsilon(),0.75)),
sqrt_epsilon(sqrt(std::numeric_limits<type0>::epsilon())),
prev_val(-1.0)
{    
}
/*--------------------------------------------
 destructor
 --------------------------------------------*/
LineSearch::~LineSearch()
{
    
}
/*--------------------------------------------
 constructor
 --------------------------------------------*/
LineSearchGoldenSection::LineSearchGoldenSection():
LineSearch()
{
    tol=sqrt(2.0*std::numeric_limits<type0>::epsilon());
    max_iter=5;
    brack=true;
}
/*--------------------------------------------
 destructor
 --------------------------------------------*/
LineSearchGoldenSection::~LineSearchGoldenSection()
{
}
/*--------------------------------------------
 
 --------------------------------------------*/
int LineSearchGoldenSection::line_min(MinCG* mincg,type0& nrgy,type0& alpha,int init_flag)
{
    return min(mincg,nrgy,alpha,init_flag);
}
/*--------------------------------------------
 
 --------------------------------------------*/
int LineSearchGoldenSection::line_min(MinCGDMD* mincg,type0& nrgy,type0& alpha,int init_flag)
{
    return min(mincg,nrgy,alpha,init_flag);
}
/*--------------------------------------------
 constructor
 --------------------------------------------*/
LineSearchBrent::LineSearchBrent():LineSearch()
{
    max_iter=5;
    tol=sqrt(2.0*std::numeric_limits<type0>::epsilon());
    zeps=std::numeric_limits<type0>::epsilon()*1.0e-3;
    brack=true;
}
/*--------------------------------------------
 destructor
 --------------------------------------------*/
LineSearchBrent::~LineSearchBrent()
{
}
/*--------------------------------------------
 
 --------------------------------------------*/
int LineSearchBrent::line_min(MinCG* mincg,type0& nrgy,type0& alpha,int init_flag)
{
    return min(mincg,nrgy,alpha,init_flag);
}
/*--------------------------------------------
 
 --------------------------------------------*/
int LineSearchBrent::line_min(MinCGDMD* mincg,type0& nrgy,type0& alpha,int init_flag)
{
    return min(mincg,nrgy,alpha,init_flag);
}
/*--------------------------------------------
 constructor
 --------------------------------------------*/
LineSearchBackTrack::LineSearchBackTrack():LineSearch()
{
    min_alpha=0.0;
    c=0.4;
    rho=0.5;
}
/*--------------------------------------------
 destructor
 --------------------------------------------*/
LineSearchBackTrack::~LineSearchBackTrack()
{
}
/*--------------------------------------------
 
 --------------------------------------------*/
int LineSearchBackTrack::line_min(MinCG* mincg,type0& nrgy,type0& alpha,int init_flag)
{
    return min(mincg,nrgy,alpha,init_flag);
}
/*--------------------------------------------
 
 --------------------------------------------*/
int LineSearchBackTrack::line_min(MinCGDMD* mincg,type0& nrgy,type0& alpha,int init_flag)
{
    return min(mincg,nrgy,alpha,init_flag);
}
