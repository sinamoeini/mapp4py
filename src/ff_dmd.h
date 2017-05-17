#ifndef __MAPP__ff_dmd__
#define __MAPP__ff_dmd__
#include "ff.h"
#include "min_vec.h"
namespace MAPP_NS
{
    template<class> class DMDVec;
    class ForceFieldDMD : public ForceField
    {
    friend class DynamicDMD;
    private:
        void reset();
        void reset_c_d();
    protected:
        class DynamicDMD* dynamic;
        class AtomsDMD* atoms;
        
        virtual void force_calc_static()=0;
        virtual void c_d_calc()=0;
        virtual void J(Vec<type0>*,Vec<type0>*)=0;
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
        DMDVec<type0>* c_d;
        
        ForceFieldDMD(class AtomsDMD*);
        virtual ~ForceFieldDMD();
                
        void force_calc_timer();
        type0 energy_calc_timer();
        
        type0 value_timer();
        void derivative_timer();
        void derivative_timer(type0(*&)[__dim__]);
        void derivative_timer(bool,type0(*&)[__dim__]);
        
        
        virtual void init_static()=0;
        virtual void fin_static()=0;
        void force_calc_static_timer();
        void c_d_calc_timer();
        type0 c_dd_norm_timer();
        void J_timer(Vec<type0>*,Vec<type0>*);
        
        
        virtual type0 prep(VecTens<type0,2>&)=0;
        virtual void J(VecTens<type0,2>&,VecTens<type0,2>&)=0;
        virtual void init_refine()=0;
        virtual void fin_refine()=0;
    };
}
#endif 
