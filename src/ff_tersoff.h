#ifndef __MAPP__ff_tersoff__
#define __MAPP__ff_tersoff__
#include "ff_md.h"
struct _object;
typedef _object PyObject;
namespace MAPP_NS
{
    class ForceFieldTersoff: public ForceFieldMD
    {
    private:
        
        //type0*** lambda_2;
        type0 zeta_ij_k(const type0&,const type0&,
        const type0&,const type0&,const type0&,const type0&,
        const type0&,const type0&,const type0&,const type0&,
        const type0&,const type0*&,const type0*&);
        
    protected:
        void __force_calc();
        void __energy_calc();
        void pre_xchng_energy(GCMC*);
        type0 xchng_energy(GCMC*);
        void post_xchng_energy(GCMC*);
    public:
        ForceFieldTersoff(class AtomsMD*);
        ~ForceFieldTersoff();
        
        static void ml_new(PyMethodDef&);

        void init();
        void fin();
        void init_xchng();
        void fin_xchng();

    };
}
#endif
