#ifndef __MAPP__ff_md__
#define __MAPP__ff_md__
#include "ff.h"
namespace MAPP_NS
{
    class ForceFieldMD : public ForceField
    {
        friend class DynamicMD;
    private:
    protected:
        class DynamicMD* dynamic;
        class AtomsMD* atoms;
        
        virtual void pre_xchng_energy(class GCMC*)=0;
        virtual type0 xchng_energy(class GCMC*)=0;
        virtual void post_xchng_energy(class GCMC*)=0;
        
        Vec<elem_type>*& elem;
    public:
        ForceFieldMD(class AtomsMD*);
        virtual ~ForceFieldMD();
        void setup();
        virtual void init_xchng()=0;
        virtual void fin_xchng()=0;
        class NeighborMD* neighbor;
        
        type0 max_cut;
        
        void pre_xchng_energy_timer(class GCMC*);
        type0 xchng_energy_timer(class GCMC*);
        void post_xchng_energy_timer(class GCMC*);
        
        void force_calc_timer();
        type0 energy_calc_timer();
        
        type0 value_timer();
        void derivative_timer();
        void derivative_timer(type0(*&)[__dim__]);
        
        int gcmc_n_vars;
        int gcmc_n_cutoff;
        bool gcmc_tag_enabled;
    };
}

#endif
