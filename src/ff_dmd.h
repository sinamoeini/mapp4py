#ifndef __MAPP__ff_dmd__
#define __MAPP__ff_dmd__
#include "ff.h"
namespace MAPP_NS
{
    template<class> class DMDVec;
    class ForceFieldDMD : public ForceField
    {
    friend class DynamicDMD;
    private:
    protected:
        class DynamicDMD* dynamic;
        class AtomsDMD* atoms;
        
        virtual void force_calc_static()=0;
        void pre_init();
        void post_fin();
    public:
        bool dof_empty,dof_alpha_empty,dof_c_empty;
        
        type0 max_cut;
        class NeighborDMD* neighbor;
        int c_dim;
        type0* rsq_crd;
        type0* r_crd;
        type0** cut_sk;
        Vec<type0>* f;
        DMDVec<type0>* f_alpha;
        
        ForceFieldDMD(class AtomsDMD*);
        virtual ~ForceFieldDMD();
        void reset();

        
        void force_calc_timer();
        type0 energy_calc_timer();
        
        type0 value_timer();
        void derivative_timer();
        void derivative_timer(type0(*&)[__dim__]);
        
        
        virtual void init_static()=0;
        virtual void fin_static()=0;
        void force_calc_static_timer();
        virtual void operator()(Vec<type0>*,Vec<type0>*)=0;
        virtual type0 update_J(type0,type0*,type0*)=0;
        virtual type0 ddc_norm()=0;
    };
}
#endif 
