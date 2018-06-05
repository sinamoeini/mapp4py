#ifndef __MAPP__ff_eam_func__
#define __MAPP__ff_eam_func__
#include "ff_md.h"
namespace MAPP_NS
{
    class ForceFieldEAMFunc: public ForceFieldMD
    {
    private:        
        
    protected:
        void force_calc();
        void energy_calc();
        void pre_xchng_energy(GCMC*);
        type0 xchng_energy(GCMC*);
        void post_xchng_energy(GCMC*);
    public:
        ForceFieldEAMFunc(class AtomsMD*);
        ~ForceFieldEAMFunc();
        
        static void ml_new(PyMethodDef&);
        void init();
        void fin();
        void init_xchng();
        void fin_xchng();
        
        
    };
    
    
}

class Foo{
public:
    void bar(){
        std::cout << "Hello" << std::endl;
    }
};
#endif
