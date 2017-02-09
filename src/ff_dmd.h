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
        class AtomsDMD*& atoms;
        class NeighborDMD* neighbor_dmd;
        
    public:
        ForceFieldDMD(class AtomsDMD*&);
        virtual ~ForceFieldDMD();
        type0* rsq_crd;
        type0* r_crd;
        void reset();
        int c_dim;
        int x_dim;
        
        void setup();
        type0 max_cut;
        type0** cut_sk;
        
        void force_calc_timer();
        type0 energy_calc_timer();
        
        type0 value_timer();
        void derivative_timer();
        void derivative_timer(type0(*&)[__dim__]);
        Vec<type0>* f_alpha;
    };
}
#endif 
