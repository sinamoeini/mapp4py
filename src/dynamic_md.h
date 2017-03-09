#ifndef __MAPP__dynamic_md__
#define __MAPP__dynamic_md__

#include "dynamic.h"

namespace MAPP_NS
{
    class vec;
    template<typename> class Vec;
    class DynamicMD: public Dynamic
    {
    private:
        void store_x0();
        bool decide();
        class ForceFieldMD* ff;
        class AtomsMD* atoms;
    protected:
        Vec<type0>* x0;
    public:        
        DynamicMD(class AtomsMD*,class ForceFieldMD*,bool,
                  std::initializer_list<vec*>,std::initializer_list<vec*>,std::initializer_list<vec*> = {});
        ~DynamicMD();
        
        void update(vec**,int);
        void update(vec*);
        void init_xchng();
        void fin_xchng();
        void init();
        void fin();
        
        
    };
}
#endif
