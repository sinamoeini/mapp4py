#ifndef __MAPP__neighbor_md__
#define __MAPP__neighbor_md__
#include "neighbor.h"
namespace MAPP_NS
{
    class NeighborMD:public Neighbor
    {
    private:
        Vec<elem_type>*& elem;
        type0**& cut_sk_sq;
    protected:
    public:
        NeighborMD(class AtomsMD*,type0**&);
        ~NeighborMD();
        
        void create_list(bool);
        void init();
        void fin();
    };
}
#endif
