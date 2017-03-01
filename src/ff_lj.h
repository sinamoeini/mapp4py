#ifndef __MAPP__ff_lj__
#define __MAPP__ff_lj__
#include "ff_md.h"
struct _object;
typedef _object PyObject;
namespace MAPP_NS
{
    class ForceFieldLJ: public ForceFieldMD
    {
    private:
        const bool shift;
        type0** sigma;
        type0** epsilon;
        type0** offset;
        
    protected:
        void force_calc();
        void energy_calc();
        void pre_xchng_energy(GCMC*);
        type0 xchng_energy(GCMC*);
        void post_xchng_energy(GCMC*);
    public:
        ForceFieldLJ(class AtomsMD*,type0**&&,type0**&&,type0**&&,bool);
        ~ForceFieldLJ();
        
        static void ml_new(PyMethodDef&);

        void init();
        void fin();
        void init_xchng();
        void fin_xchng();

    };
}
#endif
