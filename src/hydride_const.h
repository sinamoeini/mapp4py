#ifndef __MAPP__hydride_const__
#define __MAPP__hydride_const__

#include "global.h"
#include "potfit.h"
namespace MAPP_NS
{
    class HydrideConst
    {
    private:
    protected:
    public:       
        HydrideConst(PotFit*,type0(*)[2],size_t,type0(*)[2],size_t,size_t*&&,size_t,size_t*&&,size_t
        ,type0*&&,type0**&&,type0*&&,type0**&&,type0**&&,type0**&&,type0***&&);
        ~HydrideConst();
        
        int natms;
        int nFe;
        int nH;
        size_t rho_sz;
        
        size_t FH_sz;
        
        
        type0* rho0s;
        type0* rhos;
        type0* Fs;
        type0* dFs;
        type0* ddFs;

        size_t phi_const_sz,phi_nconst_sz;
        size_t* phi_const_list;
        size_t* phi_nconst_list;
        
        PotFit*  potfit;
        //given and constant
        type0* phi_const0;//[phi_const_sz]
        //given and constant
        type0** A;//[phi_const_sz][phi_nconst_sz]
        //given and constant
        type0** B;//[natms][rho_sz]
        //given and constant
        type0*** C;//[phi_const_sz][natms][rho_sz]
        //given and constant
        type0** D0;//[phi_const_sz][natms]
        //given and constant
        type0** L;//[phi_const_sz][natms]
        
        type0** D;//[phi_const_sz][natms]
        type0** Dphi_Drho;//[phi_const_sz][rho_sz]
        type0** Dphi_DFH;//[phi_const_sz][3]
        void prep();
        void adj_deriv();
        void update();
        void readj_dof();
        
    };
}
#endif
