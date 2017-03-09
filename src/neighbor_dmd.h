#ifndef __MAPP__neighbor_dmd__
#define __MAPP__neighbor_dmd__
#include "neighbor.h"
namespace MAPP_NS
{
    class NeighborDMD:public Neighbor
    {
    private:
        Vec<elem_type>*& elem;
        type0**& cut_sk;
        type0*& rsq_crd;
        const type0 scl;
        class AtomsDMD* atoms;
    protected:
    public:
        NeighborDMD(class AtomsDMD*,type0**&,type0*&);
        ~NeighborDMD();
        
        void create_list(bool);
        void create_2nd_list();
        void init();
        void fin();
        
        void init_static();
        void fin_static();
        
        int** neighbor_list_2nd;
        int* neighbor_list_size_2nd;
        int neighbor_list_size_size_2nd;
        
        size_t no_pairs_2nd;
        
        void mark_redndnt_ph(byte*);
        void rename_atoms(int*);
    };
}
#endif
