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
        type0 __vec[__nvoigt__+2];
        type0 __vec_lcl[__nvoigt__+2];
        class DynamicDMD* dynamic;
#ifdef NEW_UPDATE
#else
        class DynamicDMD* dynamic;
#endif
        class AtomsDMD* atoms;
        
        void pre_init();
        void post_fin();
        void impose_dof(type0*,type0*);
        type0 norm_sq(type0*,type0*);
        
        virtual void __force_calc_static()=0;
        virtual void __c_d_calc()=0;
        virtual void __J(Vec<type0>*,Vec<type0>*)=0;
        virtual void __J(Vec<type0>*,Vec<type0>*,Vec<type0>*,Vec<type0>*)=0;
        virtual void __prepJ_n_res(Vec<type0>*,Vec<type0>*)=0;
    public:
        bool dof_empty,dof_alpha_empty,dof_c_empty;
        
        type0 max_cut;
        class NeighborDMD* neighbor;
        int c_dim;
        type0* rsq_crd;
        type0* r_crd;
        type0** cut_sk;
        type0* ave_mu;
        Vec<type0>* f;
        DMDVec<type0>* f_alpha;
        DMDVec<type0>* c_d;
        DMDVec<type0>* mu;
        type0 F_H[__dim__][__dim__];

        
        ForceFieldDMD(class AtomsDMD*);
        virtual ~ForceFieldDMD();
                
        void force_calc_timer();
        
        type0 value();
        type0* derivative();
        //void derivative_timer(type0(*&)[__dim__]);
        //void derivative_timer(bool,type0(*&)[__dim__]);
        
        
        virtual void init_static()=0;
        virtual void fin_static()=0;
        void force_calc_static();
        void c_d_calc();
        type0 c_dd_norm();
        void J(Vec<type0>*,Vec<type0>*);
        
        
        
        virtual type0 prepJ_n_res(Vec<type0>*,Vec<type0>*);
        type0* J(Vec<type0>*,Vec<type0>*,Vec<type0>*,Vec<type0>*);
        
        void calc_thermo();
        type0 err;
        type0 c_d_norm;
        
        void calc_ndof();
        int nx_dof,nalpha_dof,nc_dof;
    };
}
#endif 
