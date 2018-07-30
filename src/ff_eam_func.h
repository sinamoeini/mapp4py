#ifndef __MAPP__ff_eam_func__
#define __MAPP__ff_eam_func__
#include "ff_md.h"
namespace MAPP_NS
{
    class ForceFieldEAMFunc: public ForceFieldMD
    {
    private:        
        class EAMFunc* eam_func;
        Vec<type0>* rho_ptr;
        Vec<type0>* F_ptr;
        Vec<type0>* rho_xchng_ptr;
        Vec<type0>* F_xchng_ptr;
        elem_type* elem_map;
    protected:
        void force_calc();
        void energy_calc();
        void pre_xchng_energy(GCMC*);
        type0 xchng_energy(GCMC*);
        void post_xchng_energy(GCMC*);
    public:
        ForceFieldEAMFunc(class AtomsMD*,class EAMFunc*,elem_type*&&);
        ~ForceFieldEAMFunc();
        
        static void ml_new(PyMethodDef&);
        void init();
        void fin();
        void init_xchng();
        void fin_xchng();
        
        
    };
    
    
}

#endif
