#ifndef __MAPP__ff_dmd__
#define __MAPP__ff_dmd__
#include "ff.h"
namespace MAPP_NS
{
    template<class> class DMDVec;
    class ForceFieldDMD : public ForceField
    {
    friend class DynamicDMD;
    template<bool,bool,bool>
    friend class NewDynamicDMD;
    friend class DAE;
    private:
        void reset();
        void reset_c_d();
        
#ifdef CONSERVATIVE
        type0* tot_mu_lcl;
        type0* tot_mu;
        type0* tot_nelms_lcl;
        type0* tot_nelms;
#endif
    protected:
        
        type0 kBT,beta;
        type0 __vec[__nvoigt__+3];
        type0 __vec_lcl[__nvoigt__+3];
        class AtomsDMD* atoms;
        
        void pre_init();
        void post_fin();
        void impose_dof(type0*,type0*);
        void impose_c_dof(type0*);
        void norm_sq(type0*,type0&,type0*,type0&);
        
        virtual void __force_calc_static()=0;
        virtual void __c_d_calc()=0;
        virtual void __J(Vec<type0>*,Vec<type0>*)=0;
        virtual void __J(Vec<type0>*,Vec<type0>*,Vec<type0>*,Vec<type0>*)=0;
        virtual void __prepJ_n_res(Vec<type0>*,Vec<type0>*)=0;
        virtual void __force_calc_gp(){};
        virtual void __energy_calc_gp(){};
    public:
        bool dof_empty,dof_alpha_empty,dof_c_empty;
        
        type0 max_cut;
        class NeighborDMD* neighbor;
        int c_dim;
        type0* mu_0;
        type0* lambda;
        type0* rsq_crd;
        type0* r_crd;
        type0* ext_mu;
        type0** cut_sk;
        type0* ave_mu;
        Vec<type0>* f;
        DMDVec<type0>* f_alpha;
        DMDVec<type0>* c_d;
        DMDVec<type0>* f_c;
        type0 F_H[__dim__][__dim__];

        
        ForceFieldDMD(class AtomsDMD*);
        virtual ~ForceFieldDMD();
                
        void force_calc_timer();
        
        type0 value();
        type0* derivative();
        
        template<bool=false>
        type0* new_derivative();
        template<bool=false>
        type0 new_value();
        
        
        
        virtual void init_static()=0;
        virtual void fin_static()=0;
        void force_calc_static();
        void c_d_calc();
        type0 c_dd_norm();
        void J(Vec<type0>*,Vec<type0>*);
        
        
        
        virtual void prepJ_n_res(Vec<type0>*,Vec<type0>*);
        type0* J(Vec<type0>*,Vec<type0>*,Vec<type0>*,Vec<type0>*);
        
        void Jnew(type0 (*)[__dim__],Vec<type0>*,Vec<type0>*,type0 (*)[__dim__],Vec<type0>*,Vec<type0>*);
        void Jnew(Vec<type0>*,Vec<type0>*,Vec<type0>*,Vec<type0>*);

        void calc_thermo();
        type0 err_sq_alpha;
        type0 err_sq_x;
        type0 c_d_norm;
        
        void calc_ndof();
        int nx_dof,nalpha_dof,nc_dof;
    };
    template<>
    type0* ForceFieldDMD::new_derivative<false>();
    template<>
    type0* ForceFieldDMD::new_derivative<true>();
    template<>
    type0 ForceFieldDMD::new_value<false>();
    template<>
    type0 ForceFieldDMD::new_value<true>();
}
#endif 
