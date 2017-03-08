#ifndef __MAPP__sgcmc__
#define __MAPP__sgcmc__
#include "gcmc.h"
namespace MAPP_NS
{
    
    /*--------------------------------------------
     allocation for this constructor has 3 levels:
     
     0. buffers that have constant size during the
     life of this object. They are allocated in
     constructor and destroyed destructor:
     
     cut_s;
     s_lo_ph;
     s_hi_ph;
     cell_size;
     ncells_per_dim;
     cell_denom;
     icell_coord;
     jcell_coord;
     nimages_per_dim; 
     *nimages_per_dim;
     ins_s_trials; (!! NOT *ins_s_trials)
     rel_neigh_lst_coord;
     rel_neigh_lst;
     
     1. buffers whose size are dependent on box 
     dimensions and domain dimensions:
     
     head_atm;
     *ins_s_trials;
     ins_cell;
     ins_cell_coord;
     ins_buff;
     ins_cell;
     del_lst;
     
     2. buffers whose sizes are decided on fly
     
     del_ids;
     
     --------------------------------------------*/
    
    
    class SGCMC:public GCMC
    {
    private:
        bool dof_empty;
        const int m;

        int ngas_before;
        
        // size dim
        type0 cell_size[__dim__];
        int N_cells[__dim__];
        int B_cells[__dim__];
        int icell_coord[__dim__];
        int jcell_coord[__dim__];

        int* rel_neigh_lst_coord;
        
        
        int* head_atm;
        int* cell_coord_buff;
        
        
        // self interaction
        int iself;
        int ineigh;
        int nneighs;
        
        
        int icell,jcell;
        int jatm_next;
        int n_cells;
        
    
        Vec<int>* cell_vec_p;
        Vec<int>* next_vec_p;
        
        void find_cell_no(type0*&,int&);
        void find_cell_coord(type0*&,int*&);
                
        
        void prep_s_x_buff();

        
        void refresh();
        
        void next_jatm_reg();
        void next_jatm_self();
        void (SGCMC::*next_jatm_p)();
        void reset_tag();
        
        type0* s_x_buff;
        void attmpt();
        type0 delta_u;
    protected:
        void ins_succ();
        void del_succ();
        void box_setup();
        void box_dismantle();

    public:
        SGCMC(class AtomsMD*&,class ForceFieldMD*&,class DynamicMD*&,int,elem_type,type0,type0,int);
        ~SGCMC();
        
        void init();
        void fin();
        
        
        void xchng(bool,int);
        
        void next_iatm();
        void next_jatm();
        void next_icomm();
        
        void reset_iatm();
        void reset_jatm();
        void reset_icomm();
    };
}
#endif

