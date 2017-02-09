#ifdef FF_Style
    FFStyle(ForceField_eam,eam)
#else
#ifndef __MAPP__ff_eam__
#define __MAPP__ff_eam__
#include "ff_md.h"
namespace MAPP_NS
{
    class ForceFieldEAM: public ForceFieldMD
    {
    private:
        Vec<type0>* rho_ptr;
        Vec<type0>* F_ptr;
        Vec<type0>* rho_xchng_ptr;
        Vec<type0>* F_xchng_ptr;
        
        size_t nr,nrho;
        type0 dr,drho,dr_inv,drho_inv,rho_max;
        type0(** F_arr)[7];
        type0(*** r_phi_arr)[7];
        type0(*** rho_arr)[7];
        
        
        /*--------------------------------------------*/
        type0* drhoi_dr;
        type0* drhoj_dr;
        size_t max_pairs;
        /*--------------------------------------------*/
    protected:
        
        void force_calc();
        void energy_calc();
        void pre_xchng_energy(GCMC*);
        type0 xchng_energy(GCMC*);
        void post_xchng_energy(GCMC*);
    public:
        ForceFieldEAM(AtomsMD*&,
        type0,type0,size_t,size_t,
        type0(***&&)[7],type0(***&&)[7],type0(**&&)[7],
        type0**&&);
        ~ForceFieldEAM();
        
        static void ml_new(PyMethodDef&,PyMethodDef&,PyMethodDef&);
        
        void init();
        void fin();
        void init_xchng();
        void fin_xchng();
    };
    

    
}
#endif
#endif
