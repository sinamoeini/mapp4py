#ifndef __MAPP__dynamic__
#define __MAPP__dynamic__

#include <initializer_list>
#include <mpi.h>
#include "global.h"
#include "exchange.h"

namespace MAPP_NS
{
    class vec;
    template<typename> class Vec;
    class Dynamic
    {
    private:
        void store_x0();
        bool decide();
        class ForceField* ff;
        class Atoms* atoms;
    protected:
        void store_arch_vecs();
        void restore_arch_vecs();
        void create_dynamic_vecs();
        void destroy_dynamic_vecs();
        
        const bool chng_box;
        
        vec** arch_vecs;
        int narch_vecs;
        vec** xchng_comp_vecs;
        int nxchng_comp_vecs;
        vec** updt_vecs;
        int nupdt_vecs;
        int nupdt_vecs_full;
        int nxchng_vecs_full;
        vec** arch_vecs_full;
        int narch_vecs_full;
        
        
        MPI_Comm& world;
        Exchange* xchng;
#ifdef NEW_UPDATE
        Update* updt;
#else
        OldUpdate* updt;
#endif
        const type0 skin;
        
    public:
        Dynamic(class Atoms*,class ForceField*,bool,
        std::initializer_list<vec*>,std::initializer_list<vec*>,
        std::initializer_list<vec*>,std::initializer_list<vec*>,
        std::initializer_list<vec*>,std::initializer_list<vec*>);
        virtual ~Dynamic();
        
        virtual void add_xchng(vec*);
        virtual void add_updt(vec*);

        
        
    };
}

#endif 
