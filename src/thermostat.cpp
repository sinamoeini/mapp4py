#include "thermostat.h"
#include <math.h>
using namespace MAPP_NS;
/*--------------------------------------------
 
 --------------------------------------------*/
ThermostatNHC::ThermostatNHC():
t_relax(0),
freq_sq(0),
nlinks(0),
niters(0),
ddt(0),
ddt2(0),
ddt4(0),
eta_d(NULL)
{}
/*--------------------------------------------
 
 --------------------------------------------*/
ThermostatNHC::ThermostatNHC(const type0 __dt,
const type0 __t_relax,const int __nlinks,
const int __niters):
t_relax(__t_relax),
freq_sq(1.0/(__t_relax*__t_relax)),
nlinks(__nlinks),
niters(__niters),
eta_d(NULL),
ddt(__dt/static_cast<type0>(__niters)),
ddt2(0.5*__dt/static_cast<type0>(__niters)),
ddt4(0.25*__dt/static_cast<type0>(__niters))
{
    eta_d=new type0[nlinks];
    for(int ich=0;ich<nlinks;ich++) eta_d[ich]=0.0;
}
/*--------------------------------------------
 
 --------------------------------------------*/
ThermostatNHC::~ThermostatNHC()
{
    delete [] eta_d;
}
/*--------------------------------------------
 
 --------------------------------------------*/
type0 ThermostatNHC::operator()(type0 T_r,int ndof)
{
    type0 rescale=1.0;
    type0 prev_fac,tmp;
    if(nlinks==1)
    {
        for(int iter=0;iter<niters;iter++)
        {
            eta_d[0]+=ddt2*freq_sq*(T_r-1.0);
            prev_fac=exp(-ddt*eta_d[0]);
            rescale*=prev_fac;
            T_r*=prev_fac*prev_fac;
            eta_d[0]+=ddt2*freq_sq*(T_r-1.0);
        }
        return rescale;
    }
    
    
    auto f=[](type0 x)->type0
    {
        if(fabs(x)<0.015625)
            return 1.0+x/2.0+x*x/6.0+x*x*x/24.0;
        return (exp(x)-1.0)/x;
    };
    
    for(int iter=0;iter<niters;iter++)
    {
        eta_d[nlinks-1]+=ddt2*(eta_d[nlinks-2]*eta_d[nlinks-2]-freq_sq);
        prev_fac=-eta_d[nlinks-1];
        
        for(int ich=nlinks-2;ich>0;ich--)
        {
            eta_d[ich]+=ddt2*f(prev_fac*ddt2)*(prev_fac*eta_d[ich]+eta_d[ich-1]*eta_d[ich-1]-freq_sq);
            prev_fac=-eta_d[ich];
        }
        
        eta_d[0]+=ddt2*f(prev_fac*ddt2)*(prev_fac*eta_d[0]+freq_sq*(T_r-1.0));
        prev_fac=-eta_d[0];
        

        tmp=exp(ddt*prev_fac);
        rescale*=tmp;
        T_r*=tmp*tmp;
        prev_fac=(T_r-1.0)*freq_sq;

        
        eta_d[0]+=ddt2*f(-eta_d[1]*ddt2)*(-eta_d[1]*eta_d[0]+prev_fac);
        prev_fac=ndof*eta_d[0]*eta_d[0]-freq_sq;
        
        for(int ich=1;ich<nlinks-1;ich++)
        {
            eta_d[ich]+=ddt2*f(-eta_d[ich+1]*ddt2)*(-eta_d[ich+1]*eta_d[ich]+prev_fac);
            prev_fac=eta_d[ich]*eta_d[ich]-freq_sq;
        }
        
        eta_d[nlinks-1]+=ddt2*prev_fac;
    }
    
    
    /*
    for(int iter=0;iter<niters;iter++)
    {
        eta_d[nlinks-1]+=(eta_d[nlinks-2]*eta_d[nlinks-2]-freq_sq)*ddt2;
        prev_fac=exp(-ddt4*eta_d[nlinks-1]);
        
        for(int ich=nlinks-2;ich>0;ich--)
        {
            eta_d[ich]=(eta_d[ich]*prev_fac+(eta_d[ich-1]*eta_d[ich-1]-freq_sq)*ddt2)*prev_fac;
            prev_fac=exp(-ddt4*eta_d[ich]);
        }
        
        eta_d[0]=(eta_d[0]*prev_fac+(T_r-1.0)*freq_sq*ddt2)*prev_fac;
        prev_fac=exp(-ddt*eta_d[0]);
        
        
        rescale*=prev_fac;
        T_r*=prev_fac*prev_fac;
        prev_fac=(T_r-1.0)*freq_sq;
        
        
        tmp=exp(-ddt4*eta_d[1]);
        eta_d[0]=(eta_d[0]*tmp+prev_fac*ddt2)*tmp;
        prev_fac=ndof*eta_d[0]*eta_d[0]-freq_sq;
        
        for(int ich=1;ich<nlinks-1;ich++)
        {
            tmp=exp(-ddt4*eta_d[ich+1]);
            eta_d[ich]=(eta_d[ich]*tmp+prev_fac*ddt2)*tmp;
            prev_fac=eta_d[ich]*eta_d[ich]-freq_sq;
        }
        
        eta_d[nlinks-1]+=ddt2*prev_fac;
    }

    */
    return rescale;
}
