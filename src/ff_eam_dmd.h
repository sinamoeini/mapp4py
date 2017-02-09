#ifndef __MAPP__ff_eam_dmd__
#define __MAPP__ff_eam_dmd__
#include "ff_dmd.h"
namespace MAPP_NS
{    
    class ForceFieldEAMDMD : public ForceFieldDMD
    {
    private:
        void allocate();
        void deallocate();
        
        
        size_t nr,nrho,stride;
        type0 dr,drho,dr_inv,drho_inv,rho_max;
        type0(** F_arr)[5];
        type0(*** r_phi_arr)[4];
        type0(*** r_rho_arr)[4];
        size_t** type2rho_pair_ij;
        size_t** type2rho_pair_ji;
        size_t** type2phi_pair_ij;
        
        Vec<type0>* cv_ptr;
        Vec<type0>* E_ptr;
        Vec<type0>* dE_ptr;
        Vec<type0>* mu_ptr;
        Vec<type0>* crd_ptr;
        Vec<type0>* s_ptr;
        Vec<type0>* t_ptr;
        
        
        /*--------------------------------------------*/
        type0 kbT,beta;
        type0* c_0;
        type0* c_1;
        type0* g_fac;
        
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

        
        type0* psi_IJ;
        type0* psi_JI;
        type0* phi_IJ;
        
        type0* psi_r_IJ;
        type0* psi_r_JI;
        type0* phi_r_IJ;
        type0* psi_alpha_IJ;
        type0* psi_alpha_JI;
        type0* phi_alpha_IJ;
        
        
        int* phi_psi_cmp;
        int* phi_psi_sz;
        int phi_psi_sz_sz;
        size_t n_phi_psi;
        type0* M_IJ;
        int M_N_sz_sz;
        void force_calc_static(bool);
        /*--------------------------------------------*/
        void calc_Q(int&,type0&,type0&,type0&,type0&,type0&);
        void calc_Q(int&,type0&,type0&,type0&,type0&,type0&,type0&);
        Vec<type0>* x_tmp_ptr;
        type0 alpha_tmp,alpha_inv_tmp;
        /*--------------------------------------------*/
        /*
        template<class F>
        void integrate(const type0& r_c,const type0& r,const type0& alpha,F& f)
        {
            type0 upper=(r+r_c)/alpha;
            type0 lower=(r-r_c)/alpha;
            if(lower>xi[N-1]) return;
            for(int i=0;i<N && xi[i]<upper;i++)
                if(xi[i]>lower)
                    f(r-xi[i]*alpha);
        }*/
    protected:
        void force_calc();
        void energy_calc();
        void dc();
        type0 ddc_norm();
        void ddc(type0*);
    public:
        ForceFieldEAMDMD(class AtomsDMD*&,
        type0,type0,size_t,size_t,
        type0(***&&)[4],type0(***&&)[4],type0(**&&)[5],
        type0**&&);
    
        ~ForceFieldEAMDMD();
        void init();
        void fin();
        void set_temp(type0);
        
        void operator()(Vec<type0>*,Vec<type0>*);
        void init_static();
        type0 update_J(type0,type0*,type0*);
        
        static void ml_new(PyMethodDef&,PyMethodDef&,PyMethodDef&);
    };
}
#endif
 


