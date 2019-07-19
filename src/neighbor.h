#ifndef __MAPP__neighbor__
#define __MAPP__neighbor__
#include "atoms.h"
#include "xmath.h"
#include <stdio.h>
namespace MAPP_NS
{
    template<typename> class Vec;
    template<const int> class Cell;
    class Neighbor
    {
    private:
        //class Atoms* atoms;
    protected:
        int no_neigh_lists;        
    public:
        bool pair_wise;
        int** neighbor_list;
        int* neighbor_list_size;
        int neighbor_list_size_size;
        size_t no_pairs;

        
        Cell<1>* cell;
        
        Neighbor(class Atoms*);
        virtual ~Neighbor();
        
        
        void print_stats();
        
        virtual void mark_redndnt_ph(byte*)=0;
        virtual void rename_atoms(int*)=0;
    };
}
#endif









