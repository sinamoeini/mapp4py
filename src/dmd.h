#ifndef __MAPP__dmd__
#define __MAPP__dmd__
#include "api.h"
namespace MAPP_NS
{
    template<typename> class Vec;
    class DMD
    {
    private:
    protected:
    public:
        int c_dim;
        int ncs;
        
        
        int max_nsteps;
        
        type0 a_tol;
        type0 min_dt;
        
        type0 nc_dofs;
        
        type0 dt;
        type0 dt_min;
        type0 t_cur;
        type0 t_fin;
        
        
        type0* c;
        type0* c_d;
        
        class AtomsDMD* atoms;
        class ForceFieldDMD* ff;
        class DynamicDMD* dynamic;
        
        
        
        DMD();
        virtual ~DMD();
        void init_static();
        void fin_static();
        
    };
}
#endif 
