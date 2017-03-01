#ifndef __MAPP__dmd_implicit__
#define __MAPP__dmd_implicit__
#include "dmd.h"
namespace MAPP_NS
{
    template<typename,class> class GMRES;
    class DMDImplicit:public DMD
    {
    private:
    protected:
        bool nonlin();
        type0 update_c();
    public:
        int max_niters_nonlin;
        int max_niters_lin;
        
        int nnonlin_acc;
        int nnonlin_rej;
        int ninteg_acc;
        int ninteg_rej;
        int nintpol_acc;
        int nintpol_rej;
        
        type0 err_fac;
        type0 beta;
        
        type0* y_0;
        type0* c_0;
        type0* del_c;
        type0* F;
        type0* a;
        
        
        
        class GMRES<type0,ForceFieldDMD>* gmres;
        
        DMDImplicit();
        virtual ~DMDImplicit();
        void init_static();
        void fin_static();

        
    };
}
#endif 
