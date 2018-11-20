#ifndef __MAPP__ff_eam_dmd__
#define __MAPP__ff_eam_dmd__
#include "ff_dmd.h"
namespace MAPP_NS
{
    template<class> class DMDVec;
    class ForceFieldEAMDMD : public ForceFieldDMD
    {
    private:
        size_t nr,nrho;
        type0 dr,drho,dr_inv,drho_inv,rho_max;
        type0(** F_arr)[5];
        type0(*** r_phi_arr)[4];
        type0(*** r_rho_arr)[4];

        
        
        Vec<type0>* E_ptr;
        Vec<type0>* dE_ptr;
        Vec<type0>* ddE_ptr;
        Vec<type0>* cv_ptr;
        
        Vec<type0>* vec0;
        Vec<type0>* vec1;
        Vec<type0>* vec2;
        Vec<type0>* vec3;
        Vec<type0>* fcoef_ptr;
        
        
        /*--------------------------------------------*/
        type0* c_1;
        type0* zeta;
        
        /*--------------------------------------------*/
        type0* rho_phi;
        type0* drho_phi_dr;
        type0* drho_phi_dalpha;
        type0* ddrho_phi_drdr;
        type0* ddrho_phi_drdalpha;
        type0* ddrho_phi_dalphadalpha;
        size_t max_pairs0,max_pairs1;
        
        /*--------------------------------------------*/
        
        type0* xi;
        type0* wi_0;
        type0* wi_1;
        const int N;
        
        /*--------------------------------------------*/
        
        type0 calc_ent(type0 x)
        {
            return x==0.0 ? 0.0:x*log(x);
        }
        void calc_mu();
        void set_temp();


        type0* M_IJ;
        
        /*--------------------------------------------*/
        void calc_Q(type0&,type0&,type0&,type0&,type0&);
        void calc_Q(type0&,type0&,type0&,type0&,type0&,type0&);
        /*--------------------------------------------*/
        
        
    protected:
        void __force_calc_static();
        void __force_calc();
        void __energy_calc();
        
        void __force_calc_gp();
        void __energy_calc_gp();
        void __c_d_calc();
        void __J(Vec<type0>*,Vec<type0>*);
        void __prepJ_n_res(Vec<type0>*,Vec<type0>*);
        void __J(Vec<type0>*,Vec<type0>*,Vec<type0>*,Vec<type0>*);
    public:
        ForceFieldEAMDMD(class AtomsDMD*,
        type0,type0,size_t,size_t,
        type0(***&&)[4],type0(***&&)[4],type0(**&&)[5],
        type0**&&,type0*&&,type0*&&);
    
        ~ForceFieldEAMDMD();
        void init();
        void fin();
        
        void init_static();
        void fin_static();        
        

        static void ml_new(PyMethodDef&,PyMethodDef&,PyMethodDef&);
    };
}
#endif
 


