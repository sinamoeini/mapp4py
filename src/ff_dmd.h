#ifndef __MAPP__ff_dmd__
#define __MAPP__ff_dmd__
#include "ff.h"
namespace MAPP_NS
{
    class ForceFieldDMD : public ForceField
    {
    friend class DynamicDMD;
    private:
    protected:
        class DynamicDMD* dynamic;
        class AtomsDMD* atoms;
        
        virtual void force_calc_static()=0;
    public:
        ForceFieldDMD(class AtomsDMD*);
        virtual ~ForceFieldDMD();
        type0* rsq_crd;
        type0* r_crd;
        void reset();
        int c_dim;
        
        void setup();
        
        class NeighborDMD* neighbor;
        
        type0 max_cut;
        bool dof_empty;
        type0** cut_sk;
        
        void force_calc_timer();
        type0 energy_calc_timer();
        
        type0 value_timer();
        void derivative_timer();
        void derivative_timer(type0(*&)[__dim__]);
        Vec<type0>* f_alpha;
        
        
        void force_calc_static_timer();
        virtual void init_static()=0;
        virtual void fin_static()=0;
        virtual void operator()(Vec<type0>*,Vec<type0>*)=0;
        virtual type0 update_J(type0,type0*,type0*)=0;
        virtual type0 ddc_norm()=0;
    };
}
#endif 
