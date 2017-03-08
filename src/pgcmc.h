#ifndef __MAPP__pgcmc__
#define __MAPP__pgcmc__
#include "gcmc.h"
namespace MAPP_NS
{
    
    class PGCMC:public GCMC
    {
    private:
#ifdef MPP_DEBUG
        type0 tot_du_test;
#endif
        bool dof_empty;
        
        int N_c[__dim__];
        int N_s[__dim__];
        int N_curr_comms[__dim__];
        int i_curr_comms[__dim__];
        bool prll_dim[__dim__];
        
        int max_n_cncrcy;
        int max_N_cncrcy[__dim__];
        
        static constexpr int comm_buff_size=__dim__*sizeof(type0)+1;
        byte comm_buff[comm_buff_size];
        
        int ip;
        int n_p;
        int op_vec[__dim__];
        int p_vec[__dim__];
        int N_p[__dim__];
        int B_p[__dim__];
        
        int n_prll;
        int N_prll[__dim__];
        int B_prll[__dim__];
        
        int n_pcomm;
        int N_pcomm[__dim__];
        int B_pcomm[__dim__];

        int n_comm;
        int N_comm[__dim__];
        int B_comm[__dim__];

        
        const int m;
        
        //static stuff
        //allocate in constructor
        //deallocate in destructor
        int N_cells[__dim__];
        int B_cells[__dim__];
        int icell_coord[__dim__];
        int jcell_coord[__dim__];
        type0 cell_size[__dim__];


        //dynamic determined by m
        //allocate in constructor
        //deallocate in destructor
        int* rel_neigh_lst_coord;
        
        //dynamic determined by box
        //allocate in setup_box
        //deallocate in dismantle_box
        int* head_atm;
        
        //dynamic determined by box & comm
        //allocate in setup_box
        //deallocate in dismantle_box
        int** cell_coord_buff;
        type0** s_x_buff;
        type0** s_buff;
        
        //dynamic determined by comm
        //allocate in setup_comm
        //deallocate in dismantle_comm
        int* ntrial_atms;
        int* roots;
        int* gcmc_mode;
        MPI_Comm* gcmc_world;
        MPI_Comm* comms;
        MPI_Comm** curr_comms;
        int** comm_id;
        type0** lcl_vars_comm;
        type0** vars_comm;
        int* success;
        int* int_buff;
        int int_buff_sz;
        
        int jatm_next;
        int iself;
        int ineigh;
        int nneighs;
        
        
        int icell,jcell;
        int n_cells;
        
    
        Vec<int>* cell_vec_p;
        Vec<int>* next_vec_p;
        
        void find_cell_no(type0*&,int&);
        void find_cell_coord(type0*&,int*&);
                
        
        void prep_s_x_buff();

                
        void next_jatm_reg();
        void next_jatm_self();
        void (PGCMC::*next_jatm_p)();
        
        void attmpt();
        
        /*---------------------------------------------------------------------*/
        

        
        class Random* lcl_random;
        int n_vars;
        int n_s;
        
        void comms_setup(int,int);
        void comms_dismantle();
        void create_comm_pattern();
        int n_curr_comms;
        void success2tag();
        void reset_tag();
        int prev_p,next_p,origin_p;
        type0 delta_u;
        int new_id;
        void decide();
        void finalize();
        void split_sing();
        void split_mult();
        /*---------------------------------------------------------------------*/

    protected:
        void ins_succ();
        void del_succ();
        void box_setup();
        void box_dismantle();

        
        
    public:
        PGCMC(class AtomsMD*&, class ForceFieldMD*&,class DynamicMD*&,int,elem_type,type0,type0,int);
        ~PGCMC();
        
        void init();
        void fin();
        
        
        void xchng(bool,int);
        
        void post_xchng();
        
        void next_iatm();
        void next_jatm();
        void next_icomm();
        
        void reset_iatm();
        void reset_jatm();
        void reset_icomm();

    };
    

    
}
#endif

