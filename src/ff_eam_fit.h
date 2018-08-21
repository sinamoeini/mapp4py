#ifndef __MAPP__ff_eam_fit__
#define __MAPP__ff_eam_fit__
#include "ff_md.h"
namespace MAPP_NS
{
    class ForceFieldEAMFit: public ForceFieldMD
    {
    private:       
    protected:
        void __force_calc();
        void __energy_calc();
        void pre_xchng_energy(GCMC*);
        type0 xchng_energy(GCMC*);
        void post_xchng_energy(GCMC*);
    public:
        ForceFieldEAMFit(class AtomsMD*,type0(***&&)[2],size_t**&&,type0(***&&)[2],size_t**&&,type0(*&&)[3]);
        ForceFieldEAMFit(class AtomsMD*,type0***&,type0***&,size_t**&,type0***&,type0***&,size_t**&,type0**&);
        ~ForceFieldEAMFit();
        
        static void ml_new(PyMethodDef&);
        void init();
        void fin();
        void init_xchng();
        void fin_xchng();
        void calc_deriv(bool***&,type0 (***&)[1+__nvoigt__]
        ,bool***&,type0 (***&)[1+__nvoigt__]
        ,bool**&,type0(**&)[1+__nvoigt__]);
        size_t get_rFeH(type0*&,int*&);
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        bool alloc;
        type0 *** rho_A;
        type0 *** rho_R;
        size_t** rho_sz;
        
        type0 *** phi_A;
        type0 *** phi_R;
        size_t** phi_sz;
        
        type0** F_A;
        
        void sort_AR_ij(type0*&,type0*&,size_t);
        type0** rho_cut;
        type0** phi_cut;
        
        
        Vec<type0>* rho_ptr;
        Vec<type0>* S_ptr;
        
        void prep4deriv();
        
        void drho_A(elem_type,elem_type,type0,type0,type0(&)[1+__nvoigt__]);
        void dphi_A(elem_type,elem_type,type0,type0,type0(&)[1+__nvoigt__]);
        void dFH_A(type0(&)[1+__nvoigt__],type0(&)[1+__nvoigt__],type0(&)[1+__nvoigt__]);
        
        
        void calc_DFH(type0&,type0&,type0&,type0&);
        void calc_DdFH(type0&,type0&,type0&,type0&);
        type0 EHH(type0);
        type0 dEHH(type0);
        type0 ddEHH(type0);
        
        
        
        type0 calc_F(elem_type,type0);
        type0 calc_dF(elem_type,type0);
        type0 calc_ddF(elem_type,type0);
        
        type0 calc_rho(elem_type,elem_type,type0);
        type0 calc_drho(elem_type,elem_type,type0);
        type0 calc_phi(elem_type,elem_type,type0);
        type0 calc_dphi(elem_type,elem_type,type0);

        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
    };
    
    
}


#endif
