#ifndef __MAPP__dynamic_dmd__
#define __MAPP__dynamic_dmd__

#include "dynamic.h"

namespace MAPP_NS
{
    class vec;
    template<typename> class Vec;
    class DynamicDMD: public Dynamic
    {
    private:
        const int c_dim;
        const type0 alpha_scale;
        void store_x0();
        bool decide();
        class ForceFieldDMD* ff;
        class AtomsDMD* atoms;
    protected:
        
        Vec<type0>* x0;
        Vec<type0>* alpha0;
    
    public:
        DynamicDMD(class AtomsDMD*,class ForceFieldDMD*,bool,
        std::initializer_list<vec*>,std::initializer_list<vec*>,std::initializer_list<vec*> = {});
        ~DynamicDMD();
        
        void update(vec**,int);
        void update(vec*);
        void init();
        void fin();
    };
}



#endif
