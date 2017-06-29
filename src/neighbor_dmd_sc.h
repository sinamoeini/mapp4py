#ifndef __MAPP__neighbor_dmd_sc__
#define __MAPP__neighbor_dmd_sc__
#include "neighbor_dmd.h"
namespace MAPP_NS
{
    class NeighborDMDSC:public NeighborDMD
    {
    private:
    protected:
    public:
        NeighborDMDSC(class AtomsDMD*,type0**&,type0*&);
        ~NeighborDMDSC();
        
        void create_list(bool);
        void init();
        void fin();
        
        void create_2nd_list(){};
        void init_static(){};
        void fin_static(){};
    
        void mark_redndnt_ph(byte*);
        void rename_atoms(int*);
    };
}
#endif
