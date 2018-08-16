#ifndef __MAPP__ff_md__
#define __MAPP__ff_md__
#include "ff.h"
namespace MAPP_NS
{
    class ForceFieldMD : public ForceField
    {
    friend class DynamicMD;
    private:
        void reset();
    protected:
        class DynamicMD* dynamic;
        class AtomsMD* atoms;
        
        
        virtual void pre_xchng_energy(class GCMC*)=0;
        virtual type0 xchng_energy(class GCMC*)=0;
        virtual void post_xchng_energy(class GCMC*)=0;
    
        void pre_init();
        void post_fin();
        
    public:
        int gcmc_n_vars;
        int gcmc_n_cutoff;
        bool gcmc_tag_enabled;
        
        bool dof_empty;
        type0 max_cut;
        class NeighborMD* neighbor;
        Vec<type0>* f;
        type0 F_H[__dim__][__dim__];
        
        ForceFieldMD(class AtomsMD*);
        virtual ~ForceFieldMD();
        
        
        void force_calc_timer();

        type0 value_timer();
        type0* derivative_timer();
        //void derivative_timer(type0(*&)[__dim__]);
        //void derivative_timer(bool,type0(*&)[__dim__]);
        
        
        virtual void init_xchng()=0;
        virtual void fin_xchng()=0;
        void pre_xchng_energy_timer(class GCMC*);
        type0 xchng_energy_timer(class GCMC*);
        void post_xchng_energy_timer(class GCMC*);
        
        void calc_ndof();
        int nx_dof;
    };
}

#endif
