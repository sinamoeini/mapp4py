#ifndef __MAPP__ff_fs__
#define __MAPP__ff_fs__
#include "ff_md.h"
namespace MAPP_NS
{
    class ForceFieldFS: public ForceFieldMD
    {
    private:        
        type0** t1;
        type0** t2;
        type0** cut_sq_phi;
        type0** cut_sq_rho;
        type0** cut_phi;
        type0** cut_rho;
        type0** k1;
        type0** k2;
        type0** k3;
        type0* A;
        
        Vec<type0>* rho_ptr;
        Vec<type0>* F_ptr;
        Vec<type0>* rho_xchng_ptr;
        Vec<type0>* F_xchng_ptr;
    protected:
        void __force_calc();
        void __energy_calc();
        void pre_xchng_energy(GCMC*);
        type0 xchng_energy(GCMC*);
        void post_xchng_energy(GCMC*);
    public:
        ForceFieldFS(class AtomsMD*,type0*&&,
        type0**&&,type0**&&,type0**&&,type0**&&,
        type0**&&,type0**&&,type0**&&);
        ~ForceFieldFS();
        
        static void ml_new(PyMethodDef&);
        void init();
        void fin();
        void init_xchng();
        void fin_xchng();
        
        
    };
}
#endif
