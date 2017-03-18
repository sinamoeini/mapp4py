#ifndef __MAPP__thermostat__
#define __MAPP__thermostat__
#include <stdio.h>
#include "global.h"
namespace MAPP_NS
{
    class ThermostatNHC
    {
    private:
        const type0 ddt;
        const type0 ddt2;
        const type0 ddt4;
        const type0 freq_sq;
        type0* eta_d;
    protected:
    public:
        const type0 t_relax;
        const int L;
        const int niters;
        ThermostatNHC(const type0,const type0,const int,const int);
        ThermostatNHC();
        ~ThermostatNHC();
        type0 operator()(type0,int);
    };
}
#endif
