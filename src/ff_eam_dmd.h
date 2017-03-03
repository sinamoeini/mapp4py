#ifndef __MAPP__ff_eam_dmd__
#define __MAPP__ff_eam_dmd__
#include "ff_dmd.h"
namespace MAPP_NS
{    
    class ForceFieldEAMDMD : public ForceFieldDMD
    {
    private:
        size_t nr,nrho;
        type0 dr,drho,dr_inv,drho_inv,rho_max;
        type0(** F_arr)[5];
        type0(*** r_phi_arr)[4];
        type0(*** r_rho_arr)[4];

        
        
        Vec<type0>* dE_ptr;
        Vec<type0>* ddE_ptr;
        Vec<type0>* mu_ptr;
        Vec<type0>* cv_ptr;
        
        Vec<type0>* vec0;
        Vec<type0>* vec1;
        Vec<type0>* vec2;
        Vec<type0>* vec3;
        
        
        /*--------------------------------------------*/
        type0 kbT,beta;
        type0* c_0;
        type0* c_1;
        
        /*--------------------------------------------*/
        type0* rho_phi;
        type0* drho_phi_dr;
        type0* drho_phi_dalpha;
        size_t max_pairs;
        
        /*--------------------------------------------*/
        
        type0* xi;
        type0* wi_0;
        type0* wi_1;
        const int N;
        
        /*--------------------------------------------*/
        
        type0 calc_ent(type0);
        void calc_mu();
        


        type0* M_IJ;
        
        /*--------------------------------------------*/
        
        void calc_Q(elem_type&,type0&,type0&,type0&,type0&,type0&);
        void calc_Q(elem_type&,type0&,type0&,type0&,type0&,type0&,type0&);
        void __calc_Q(elem_type&,type0&,type0&,type0&,type0&,type0&,type0&,type0&,type0&);
        /*--------------------------------------------*/
                
    protected:
        void force_calc_static();
        void force_calc();
        void energy_calc();
        void dc();
    public:
        ForceFieldEAMDMD(class AtomsDMD*,
        type0,type0,size_t,size_t,
        type0(***&&)[4],type0(***&&)[4],type0(**&&)[5],
        type0**&&,type0*&&);
    
        ~ForceFieldEAMDMD();
        void init();
        void fin();
        void set_temp(type0);
        
        void init_static();
        void fin_static();
        void operator()(Vec<type0>*,Vec<type0>*);
        type0 update_J(type0,type0*,type0*);
        
        void update_J(type0*);
        void update_J2(type0*);
        void kernel(Vec<type0>*,Vec<type0>*);
        void kernel2(Vec<type0>*,Vec<type0>*);
        
        type0 ddc_norm();
        static void ml_new(PyMethodDef&,PyMethodDef&,PyMethodDef&);
    };
}
#endif
 


