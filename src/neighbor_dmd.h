#ifndef __MAPP__neighbor_dmd__
#define __MAPP__neighbor_dmd__
#include "neighbor.h"
namespace MAPP_NS
{
    class NeighborDMD:public Neighbor
    {
    private:
        Vec<elem_type>*& elem;
        Vec<type0>*& c_vec;
        type0**& cut_sk;
        type0*& rsq_crd;
        const type0 scl;
        class AtomsDMD* atoms_dmd;
    protected:
    public:
        NeighborDMD(class AtomsDMD*,type0**&,type0*&);
        ~NeighborDMD();
        
        void create_list(bool);
        void create_2nd_list();
        void init();
        void fin();
        
        int** neighbor_list_2nd;
        int* neighbor_list_size_2nd;
        int neighbor_list_size_size_2nd;
        
    };
}
#endif
