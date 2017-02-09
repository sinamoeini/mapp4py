#ifndef __MAPP__dynamic_dmd__
#define __MAPP__dynamic_dmd__

#include <mpi.h>
#include "global.h"
#include "exchange.h"

namespace MAPP_NS
{
    class vec;
    template<typename> class Vec;
    class DynamicDMD
    {
    private:
        const int c_dim;
        const type0 alpha_scale;
        void store_x0();
        bool decide();
        class ForceFieldDMD* ff;
        class AtomsDMD* atoms;
    protected:
        void store_arch_vecs();
        void restore_arch_vecs();
        
        
        const bool box_chng;
        
        vec** arch_vecs;
        int narch_vecs;
        int nxchng_vecs;
        int nupdt_vecs;
        Vec<type0>* x0;
        Vec<type0>* alpha0;
        Vec<unsigned int>* id_arch;
        
        
        
        MPI_Comm& world;
        Exchange* xchng;
        Update* updt;
        const type0 skin;
    
    public:
        DynamicDMD(class AtomsDMD*,class ForceFieldDMD*,
        bool,vec* const *,int,vec* const *,int,vec* const *,int);
        DynamicDMD(class AtomsDMD*,class ForceFieldDMD*,bool,
        std::initializer_list<vec*>,std::initializer_list<vec*>,std::initializer_list<vec*>);
        DynamicDMD(class AtomsDMD*,class ForceFieldDMD*,bool,
        std::initializer_list<vec*>,std::initializer_list<vec*>);
        ~DynamicDMD();
        
        void add_xchng(vec*);
        void add_updt(vec*);
        
        void update(vec**,int);
        void update(vec*);
        void init();
        void fin();
    };
}



#endif
