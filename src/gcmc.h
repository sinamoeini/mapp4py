/*--------------------------------------------
 Created by Sina on 06/29/16.
 Copyright (c) 2013 MIT. All rights reserved.
 --------------------------------------------*/
#ifndef __MAPP__gcmc__
#define __MAPP__gcmc__
#include "global.h"
#include <mpi.h>
namespace MAPP_NS
{
    enum{NOEX_MODE,DEL_MODE,INS_MODE};
    enum{MINE_FLAG,INTERACT_FLAG,NONEINTERACT_FLAG};
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
    template<typename> class Vec;
    class GCMC
    {
    private:
    protected:
        
                
        int igas,gas_id,ngas_lcl;
        elem_type gas_type;
        type0 vol;
        //constants
        type0 gas_mass,beta,kbT,T,mu,lambda,sigma,z_fac;
        
        int& natms_lcl;
        int& natms_ph;
        type0 cut;
        type0**& cut_sq;
        type0 (&s_lo)[__dim__];
        type0 (&s_hi)[__dim__];
        
        // size dim
        type0 s_buff[__dim__];
        type0 vel_buff[__dim__];
        type0 cut_s[__dim__];
        type0 s_lo_ph[__dim__];
        type0 s_hi_ph[__dim__];
        
        int nimages_per_dim[__dim__][2];
        type0** s_trials;

        
        int del_idx;
        
        int* del_ids;
        int del_ids_sz,del_ids_cpcty;
        int max_id;

        Vec<type0>* s_vec_p;

        
        
        class Random* random;
        
        
        int itrial_atm,ntrial_atms,max_ntrial_atms;
        
        virtual void ins_succ()=0;
        virtual void del_succ()=0;
        virtual void box_setup();
        virtual void box_dismantle();
        void add_del_id(int*,int);
        int get_new_id();

        class DynamicMD*& dynamic;
        MPI_Comm& world;
        
        
        
        /*************************************************************************/
        class AtomsMD*& atoms;
        class ForceFieldMD*& ff;
        
        type0 mvv_lcl[__nvoigt__];
#ifdef GCMCDEBUG
        type0 tot_delta_u_lcl;
#endif
    public:
        GCMC(class AtomsMD*&, class ForceFieldMD*&,class DynamicMD*&,elem_type,type0,type0,int);
        virtual ~GCMC();
        
        virtual void init();
        virtual void fin();

        virtual void xchng(bool,int)=0;

        
        virtual void next_iatm()=0;
        virtual void next_jatm()=0;
        virtual void next_icomm()=0;
        
        virtual void reset_iatm()=0;
        virtual void reset_jatm()=0;
        virtual void reset_icomm()=0;
        
        int iatm;
        int niatms;
        elem_type& ielem;
        type0* ix;
        type0* jx;
        int jatm;
        elem_type jelem;
        type0 rsq;
        int xchng_mode;
        int dof_diff;
        int ngas;
        bool im_root;

        Vec<int>* tag_vec_p;
        int icomm;
        type0* lcl_vars;
        type0* vars;
        MPI_Comm* curr_comm;
        int curr_root;
        bool root_succ;
        
        type0 mvv[__nvoigt__];
#ifdef GCMCDEBUG
        type0 tot_delta_u;
#endif
    };
}
#endif

