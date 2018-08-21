#ifndef __MAPP__ff_eam_dmd_cluster__
#define __MAPP__ff_eam_dmd_cluster__
#include "ff_dmd.h"
namespace MAPP_NS
{
    template<class> class DMDVec;
    class ForceFieldEAMDMDCLUSTER : public ForceFieldDMD
    {
    private:
        size_t nr,nrho;
        type0 dr,drho,dr_inv,drho_inv,rho_max;
        type0(** F_arr)[5];
        type0(*** r_phi_arr)[4];
        type0(*** r_rho_arr)[4];

        
        
        Vec<type0>* E_ptr;
        Vec<type0>* dE_ptr;
        Vec<type0>* rho_ptr;
        Vec<type0>* b_ptr;
        Vec<type0>* d_ptr;
        Vec<type0>* theta_ptr;
        
        Vec<type0>* vec0;
        Vec<type0>* vec1;
        Vec<type0>* vec2;
        Vec<type0>* vec3;
        Vec<type0>* fcoef_ptr;
        
        
        /*--------------------------------------------*/
        type0 kbT,beta;
        type0* c_0;
        type0* c_1;
        type0* zeta;
        
        /*--------------------------------------------*/
        type0* rho_phi;
        type0* drho_phi_dr;
        type0* drho_phi_dalpha;
        type0* B_pair;
        bool* is_pair;
        size_t max_pairs;
        
        /*--------------------------------------------*/
        
        type0* xi;
        type0* wi_0;
        type0* wi_1;
        const int N;
        
        /*--------------------------------------------*/
        
        type0 calc_ent(type0);
        void set_temp();


        void calc_pair(type0,type0*,type0*,type0*&,type0*&,type0*&);
        
        /*--------------------------------------------*/
        void calc_Q(type0&,type0&,type0&,type0&,type0&);
        void calc_Q(type0&,type0&,type0&,type0&,type0&,type0&);
        /*--------------------------------------------*/
        
        void sc_loop();
        void _sc_loop();
        void __sc_loop();
        void ___sc_loop();

        
    protected:
        void __force_calc();
        void __energy_calc();
        void __force_calc_static(){};
        void __c_d_calc(){};
        void __J(Vec<type0>*,Vec<type0>*){};
        void __prepJ_n_res(Vec<type0>*,Vec<type0>*);
        void __J(Vec<type0>*,Vec<type0>*,Vec<type0>*,Vec<type0>*);
    public:
        
        void init_static(){};
        void fin_static(){};
        ForceFieldEAMDMDCLUSTER(class AtomsDMD*,
        type0,type0,size_t,size_t,
        type0(***&&)[4],type0(***&&)[4],type0(**&&)[5],
        type0**&&,type0*&&,type0*&&);
    
        ~ForceFieldEAMDMDCLUSTER();
        void init();
        void fin();
        
        
        
        
        type0* ave_mu(){return NULL;};
        static void ml_new(PyMethodDef&,PyMethodDef&,PyMethodDef&);
    };
}
#endif
 


