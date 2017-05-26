/*--------------------------------------------
 Created by Sina on 06/05/13.
 Copyright (c) 2013 MIT. All rights reserved.
 --------------------------------------------*/
#include "atoms_md.h"
#include "pgcmc.h"
#include "memory.h"
#include "random.h"
#include "elements.h"
#include "neighbor_md.h"
#include "ff_md.h"
#include "xmath.h"
#include "MAPP.h"
#include "comm.h"
#include "dynamic_md.h"
using namespace MAPP_NS;
/*--------------------------------------------
 constructor
 --------------------------------------------*/
PGCMC::PGCMC(AtomsMD*& __atoms,ForceFieldMD*&__ff,DynamicMD*& __dynamic,int __m,elem_type __gas_type,type0 __mu,type0 __T,int seed):
GCMC(__atoms,__ff,__dynamic,__gas_type,__mu,__T,seed),
m(__m)
{
    n_comm=0;
    comm_id=NULL;
    success=NULL;
    ntrial_atms=NULL;
    gcmc_mode=NULL;
    gcmc_world=NULL;
    curr_comms=NULL;
    vars_comm=NULL;
    lcl_vars_comm=NULL;
    int_buff=NULL;
    roots=NULL;
    comms=NULL;
    
    s_buff=NULL;
    s_x_buff=NULL;
    cell_coord_buff=NULL;
    
    n_cells=0;
    head_atm=NULL;
    
    max_ntrial_atms=0;
    cell_coord_buff=NULL;
    s_x_buff=NULL;
    s_buff=NULL;
    
    
    
    
    
    lcl_random=new Random(seed+atoms->comm_rank);
    
    
    
    /*--------------------------------------------------
     find the relative neighbor list for cells here we 
     figure out the number of neighboring cells and also 
     allocate the memory for rel_neigh_lst and 
     rel_neigh_lst_coord. values for rel_neigh_lst_coord
     is assigned here,but finding the values for 
     rel_neigh_lst requires knowledge of the box and 
     domain, it will be differed to create()
     --------------------------------------------------*/
    int countr[__dim__];
    for(int i=0;i<__dim__;i++)
        countr[i]=-m;
    int max_no_neighs=1;
    for(int i=0;i<__dim__;i++)
        max_no_neighs*=2*m+1;
    nneighs=0;
    rel_neigh_lst_coord=new int[max_no_neighs*__dim__];
    int* rel_neigh_lst_coord_=rel_neigh_lst_coord;
    int sum;
    int rc_sq=m*m;
    for(int i=0;i<max_no_neighs;i++)
    {
        sum=0;
        for(int j=0;j<__dim__;j++)
        {
            sum+=countr[j]*countr[j];
            if(countr[j]!=0)
                sum+=1-2*std::abs(countr[j]);
        }
        
        if(sum<rc_sq)
        {
            for(int j=0;j<__dim__;j++)
                rel_neigh_lst_coord_[j]=countr[j];
            rel_neigh_lst_coord_+=__dim__;
            nneighs++;
        }
        
        countr[0]++;
        for(int j=0;j<__dim__-1;j++)
            if(countr[j]==m+1)
            {
                countr[j]=-m;
                countr[j+1]++;
            }
    }
    rel_neigh_lst_coord_=new int[nneighs*__dim__];
    memcpy(rel_neigh_lst_coord_,rel_neigh_lst_coord,nneighs*__dim__*sizeof(int));
    delete [] rel_neigh_lst_coord;
    rel_neigh_lst_coord=rel_neigh_lst_coord_;
    
    

}
/*--------------------------------------------
 destructor
 --------------------------------------------*/
PGCMC::~PGCMC()
{
    delete [] rel_neigh_lst_coord;
    delete lcl_random;    
    

}
/*--------------------------------------------
 
 --------------------------------------------*/
void PGCMC::init()
{
    dof_empty=atoms->dof->is_empty();
    GCMC::init();
    box_setup();
}
/*--------------------------------------------
 
 --------------------------------------------*/
void PGCMC::fin()
{
    box_dismantle();
    GCMC::fin();
}
/*--------------------------------------------
 
 --------------------------------------------*/
void PGCMC::box_setup()
{
    box_dismantle();
    GCMC::box_setup();
    n_cells=1;
    for(int i=0;i<__dim__;i++)
    {
        cell_size[i]=cut_s[i]/static_cast<type0>(m);
        N_cells[i]=static_cast<int>((s_hi[i]-s_lo[i])/cell_size[i])+1;
        
        B_cells[i]=n_cells;
        n_cells*=N_cells[i];
    }
    
    head_atm=new int[n_cells];
    comms_setup(ff->gcmc_n_vars,ff->gcmc_n_cutoff);
    
    
    s_buff=new type0*[max_n_cncrcy];
    *s_buff=new type0[max_n_cncrcy*__dim__];
    for(int i=1;i<max_n_cncrcy;i++)
        s_buff[i]=s_buff[i-1]+__dim__;
    
    cell_coord_buff=new int*[max_n_cncrcy];
    *cell_coord_buff=new int[max_n_cncrcy*max_ntrial_atms*__dim__];
    for(int i=1;i<max_n_cncrcy;i++)
        cell_coord_buff[i]=cell_coord_buff[i-1]+max_ntrial_atms*__dim__;
    
    s_x_buff=new type0*[max_n_cncrcy];
    *s_x_buff=new type0[max_n_cncrcy*max_ntrial_atms*__dim__];
    for(int i=1;i<max_n_cncrcy;i++)
        s_x_buff[i]=s_x_buff[i-1]+max_ntrial_atms*__dim__;
}
/*--------------------------------------------
 release the memory allocated by box_setup
 --------------------------------------------*/
void PGCMC::box_dismantle()
{
    
    comms_dismantle();
    
    if(s_buff) delete [] *s_buff;
    delete [] s_buff;
    s_buff=NULL;
    
    if(s_x_buff) delete [] *s_x_buff;
    delete [] s_x_buff;
    s_x_buff=NULL;
    
    if(cell_coord_buff) delete [] *cell_coord_buff;
    delete [] cell_coord_buff;
    cell_coord_buff=NULL;
    
    delete [] head_atm;
    head_atm=NULL;
    
    n_cells=0;
}
/*--------------------------------------------
 construct the bin list
 --------------------------------------------*/
void PGCMC::xchng(bool box_chng,int nattmpts)
{
    dynamic->init_xchng();
    if(box_chng) box_setup();
    
    /*--------------------------------------------------
     here we allocate the memory for cell_vec & next_vec
     --------------------------------------------------*/
    if(ff->gcmc_tag_enabled) tag_vec_p=new Vec<int>(atoms,1);
    else tag_vec_p=NULL;
    cell_vec_p=new Vec<int>(atoms,1);
    next_vec_p=new Vec<int>(atoms,1);
    s_vec_p=new Vec<type0>(atoms,__dim__);
    memcpy(s_vec_p->begin(),atoms->x->begin(),sizeof(type0)*natms_lcl*__dim__);
    
    /*--------------------------------------------------
     here we reset head_atm
     --------------------------------------------------*/
    for(int i=0;i<n_cells;i++) head_atm[i]=-1;
    
    /*--------------------------------------------------
     here we assign values for cell_vec, next_vec &
     headt_atm; it is supposed that we have fractional
     coordinates at this point
     --------------------------------------------------*/
    int* next_vec=next_vec_p->begin();
    int* cell_vec=cell_vec_p->begin();
    type0* s=atoms->x->begin()+(natms_lcl-1)*__dim__;
    for(int i=natms_lcl-1;i>-1;i--,s-=__dim__)
    {
        find_cell_no(s,cell_vec[i]);
        next_vec[i]=head_atm[cell_vec[i]];
        head_atm[cell_vec[i]]=i;
    }
    
    ff->neighbor->create_list(box_chng);
    
    ff->init_xchng();
    for(int i=0;i<atoms->ndynamic_vecs;i++)
        atoms->dynamic_vecs[i]->resize(natms_lcl);
    
    ngas=0;
    elem_type* elem=atoms->elem->begin();
    for(int i=0;i<natms_lcl;i++)
        if(elem[i]==gas_type) ngas++;
    MPI_Allreduce(&ngas,&tot_ngas,1,MPI_INT,MPI_SUM,world);
    
    
    next_jatm_p=&PGCMC::next_jatm_reg;
    
    
    int nx=nattmpts/n_prll;
    if(nattmpts%n_prll) nx++;

#ifdef GCMCDEBUG
    tot_delta_u_lcl=0.0;
#endif
    dof_diff=0;
    Algebra::zero<__nvoigt__>(mvv_lcl);
    
    for(int i=0;i<nx;i++) attmpt();

    MPI_Allreduce(mvv_lcl,mvv,__nvoigt__,Vec<type0>::MPI_T,MPI_SUM,world);
#ifdef GCMCDEBUG
    MPI_Allreduce(&tot_delta_u_lcl,&tot_delta_u,1,Vec<type0>::MPI_T,MPI_SUM,world);
#endif
    
    ff->fin_xchng();
    memcpy(atoms->x->begin(),s_vec_p->begin(),sizeof(type0)*natms_lcl*__dim__);
    
    delete s_vec_p;
    delete next_vec_p;
    delete cell_vec_p;
    delete tag_vec_p;
    dynamic->fin_xchng();

}
/*--------------------------------------------
 this must be used only for local atoms
 --------------------------------------------*/
inline void PGCMC::find_cell_no(type0*& s,int& cell_no)
{
    cell_no=0;
    for(int i=0;i<__dim__;i++)
        cell_no+=B_cells[i]*MIN(static_cast<int>((s[i]-s_lo[i])/cell_size[i]),N_cells[i]-1);
}
/*--------------------------------------------
 
 --------------------------------------------*/
inline void PGCMC::find_cell_coord(type0*& s,int*& cell_coord)
{
    for(int i=0;i<__dim__;i++)
    {
        if(s[i]<s_lo[i])
            cell_coord[i]=-static_cast<int>((s_lo[i]-s[i])/cell_size[i])-1;
        else if(s_hi[i]<=s[i])
            cell_coord[i]=MAX(static_cast<int>((s[i]-s_lo[i])/cell_size[i]),N_cells[i]-1);
        else
            cell_coord[i]=MIN(static_cast<int>((s[i]-s_lo[i])/cell_size[i]),N_cells[i]-1);
    }
}
/*--------------------------------------------
 
 --------------------------------------------*/
void PGCMC::prep_s_x_buff()
{
    
    int count[__dim__];
    int n_per_dim[__dim__];
    type0* buff;
    int* cell_coord;
    //type0 (&H)[__dim__][__dim__]=atoms->H;
    type0 s;
    int no;
    for(int i=0;i<n_curr_comms;i++)
    {
        if(gcmc_mode[i]==NOEX_MODE)
        {
            ntrial_atms[i]=0;
            continue;
        }
        
        ntrial_atms[i]=1;
        
        for(int j=0;j<__dim__;j++)
        {
            no=0;
            s=s_buff[i][j];
            for(type0 s_=s;s_<s_hi_ph[j];s_++)
                if(s_lo_ph[j]<=s_ && s_<s_hi_ph[j])
                    s_trials[j][no++]=s_;
            
            for(type0 s_=s-1.0;s_lo_ph[j]<=s_;s_--)
                if(s_lo_ph[j]<=s_ && s_<s_hi_ph[j])
                    s_trials[j][no++]=s_;
            
            ntrial_atms[i]*=no;
            n_per_dim[j]=no;
        }
        
        
        if(ntrial_atms[i]==0) continue;
        
        buff=s_x_buff[i];
        cell_coord=cell_coord_buff[i];
        
        for(int j=0;j<__dim__;j++) count[j]=0;
        for(int j=0;j<ntrial_atms[i];j++)
        {
            for(int k=0;k<__dim__;k++)
                buff[k]=s_trials[k][count[k]];
            
            find_cell_coord(buff,cell_coord);
            
            //XMatrixVector::s2x<__dim__>(buff,H);
            Algebra::S2X<__dim__>(atoms->__h,buff);
            
            count[0]++;
            for(int k=0;k<__dim__-1;k++)
                if(count[k]==n_per_dim[k])
                {
                    count[k]=0;
                    count[k+1]++;
                }
            cell_coord+=__dim__;
            buff+=__dim__;
        }
    }
}
/*--------------------------------------------
 this is used for insertion trial
 --------------------------------------------*/
void PGCMC::reset_iatm()
{
    itrial_atm=-1;
    next_iatm();
}
/*--------------------------------------------
 this is used for insertion trial
 --------------------------------------------*/
void PGCMC::next_iatm()
{
    itrial_atm++;
    
    if(itrial_atm==ntrial_atms[icomm])
    {
        /*--------------------------------------------------
         if we have reached the end of the list
         --------------------------------------------------*/
        iatm=-1;
        return;
    }
    
    /*--------------------------------------------------
     give the atom a number that cannot be the
     same as the existing ones:
     
     iatm=natms_lcl+natms_ph+itrial_atm;
     
     natms_lcl+natms_ph <= iatm
     iatm < natms_lcl+natms_ph+ntrial_atms
     --------------------------------------------------*/
    if(itrial_atm==0 && im_root && xchng_mode==DEL_MODE)
        iatm=del_idx;
    else
        iatm=natms_lcl+itrial_atm;
    
    /*--------------------------------------------------
     assign the cell number and the already calculated
     cell coordinates
     --------------------------------------------------*/
    for(int i=0;i<__dim__;i++)
        icell_coord[i]=cell_coord_buff[icomm][__dim__*itrial_atm+i];
    
    /*--------------------------------------------------
     assign the position of iatm
     --------------------------------------------------*/
    ix=s_x_buff[icomm]+itrial_atm*__dim__;
}
/*--------------------------------------------
 this is used for insertion trial
 --------------------------------------------*/
void PGCMC::reset_jatm()
{
    ineigh=0;
    jatm_next=-1;
    next_jatm_reg();
}
/*--------------------------------------------
 find the next lcl atom
 --------------------------------------------*/
inline void PGCMC::next_jatm_reg()
{
    while(true)
    {
        while(jatm_next==-1 && ineigh<nneighs)
        {
            bool lcl=true;
            jcell=0;
            for(int i=0;i<__dim__ && lcl;i++)
            {
                jcell_coord[i]=icell_coord[i]+rel_neigh_lst_coord[ineigh*__dim__+i];
                jcell+=B_cells[i]*jcell_coord[i];
                
                if(jcell_coord[i]<0 || jcell_coord[i]>N_cells[i]-1)
                    lcl=false;
            }
            
            if(lcl)
                jatm_next=head_atm[jcell];
            ineigh++;
        }
        
        jatm=jatm_next;
        if(jatm==-1)
        {
            
            if(itrial_atm==0 && im_root)
            {
                iself=0;
                next_jatm_p=&PGCMC::next_jatm_self;
                return next_jatm_self();
            }
            else return;
        }
        
        jatm_next=next_vec_p->begin()[jatm];
        
        if(jatm==iatm) continue;
        
        jx=atoms->x->begin()+__dim__*jatm;
        jelem=atoms->elem->begin()[jatm];
        rsq=Algebra::RSQ<__dim__>(ix,jx);
        if(rsq>=cut_sq[ielem][jelem]) continue;
        
        if(tag_vec_p) tag_vec_p->begin()[jatm]=icomm;
        return;
    }
}
/*--------------------------------------------
 find the next interactin image atom
 --------------------------------------------*/
inline void PGCMC::next_jatm_self()
{
    while(true)
    {
        iself++;
        if(iself==ntrial_atms[icomm])
        {
            next_jatm_p=&PGCMC::next_jatm_reg;
            jatm=-1;
            return;
        }
        
        jatm=natms_lcl+iself;
        jx=s_x_buff[icomm]+iself*__dim__;
        jelem=gas_type;
        rsq=Algebra::RSQ<__dim__>(ix,jx);
        
        if(rsq>=cut_sq[ielem][jelem]) continue;
        
        return;
    }
}
/*--------------------------------------------
 find the next interactin image atom
 --------------------------------------------*/
void PGCMC::next_jatm()
{
    (this->*next_jatm_p)();
}
/*--------------------------------------------
 
 --------------------------------------------*/
void PGCMC::next_icomm()
{
    icomm++;
    if(icomm==n_curr_comms)
    {
        icomm=-1;
        return;
    }
    xchng_mode=gcmc_mode[icomm];
    lcl_vars=lcl_vars_comm[icomm];
    vars=vars_comm[icomm];
    curr_comm=curr_comms[icomm];
    curr_root=roots[icomm];
}
/*--------------------------------------------
 
 --------------------------------------------*/
void PGCMC::reset_icomm()
{
    icomm=-1;
    next_icomm();
}
/*--------------------------------------------
 
 --------------------------------------------*/
void PGCMC::success2tag()
{
    int* tag=tag_vec_p->begin();
    for(int i=0;i<natms_lcl;i++)
        if(tag[i]!=-1)
            tag[i]=success[tag[i]];
}
/*--------------------------------------------
 
 --------------------------------------------*/
void PGCMC::reset_tag()
{
    int* tag=tag_vec_p->begin();
    for(int i=0;i<natms_lcl;i++)
        tag[i]=-1;
}
/*--------------------------------------------
 attempt
 --------------------------------------------*/
void PGCMC::attmpt()
{
    create_comm_pattern();


    if(im_root)
    {
        
        if(lcl_random->uniform()<0.5 || ngas==0)
        {
            gcmc_mode[0]=INS_MODE;
            comm_buff[0]=0;
            for(int i=0;i<__dim__;i++)
                s_buff[0][i]=s_lo[i]+lcl_random->uniform()*(s_hi[i]-s_lo[i]);
            memcpy(comm_buff+1,s_buff[0],__dim__*sizeof(type0));
        }
        else
        {
            gcmc_mode[0]=DEL_MODE;
            comm_buff[0]=1;
            igas=static_cast<int>(ngas*lcl_random->uniform());
            
            int n=igas;
            int icount=-1;
            elem_type* elem=atoms->elem->begin();
            del_idx=0;
            for(;icount!=n;del_idx++)
                if(elem[del_idx]==gas_type) icount++;
            del_idx--;
            gas_id=atoms->id->begin()[del_idx];
            
            memcpy(s_buff[0],s_vec_p->begin()+del_idx*__dim__,__dim__*sizeof(type0));
            memcpy(comm_buff+1,s_buff[0],__dim__*sizeof(type0));
        }
        
        MPI_Bcast(comm_buff,comm_buff_size,MPI_BYTE,roots[0],*curr_comms[0]);

    }
    else
        for(int i=0;i<n_curr_comms;i++)
        {
            MPI_Bcast(comm_buff,comm_buff_size,MPI_BYTE,roots[i],*curr_comms[i]);
            memcpy(s_buff[i],comm_buff+1,__dim__*sizeof(type0));
            
            if(comm_buff[0]==0)
                gcmc_mode[i]=INS_MODE;
            else
                gcmc_mode[i]=DEL_MODE;
        }
    

    prep_s_x_buff();

    if(tag_vec_p) reset_tag();

    ff->pre_xchng_energy_timer(this);

    delta_u=ff->xchng_energy_timer(this);

    root_succ=false;
    if(im_root) decide();

    for(int i=0;i<n_curr_comms;i++)
        MPI_Bcast(&success[i],1,MPI_INT,roots[i],*curr_comms[i]);

    if(tag_vec_p) success2tag();

    ff->post_xchng_energy_timer(this);

    finalize();

}
/*--------------------------------------------
 things to do after a successful deletion
 --------------------------------------------*/
void PGCMC::del_succ()
{
    
    /*--------------------------------------------------
     0.0 find the link to del_idx
     --------------------------------------------------*/
    int* p_2_idx=head_atm+cell_vec_p->begin()[del_idx];
    while(*p_2_idx!=del_idx)
        p_2_idx=next_vec_p->begin()+*p_2_idx;
    /*--------------------------------------------------
     0.1 replace link to del_idx with link from del_idx
     --------------------------------------------------*/
    *p_2_idx=next_vec_p->begin()[del_idx];
    
    if(del_idx!=natms_lcl-1)
    {
        /*--------------------------------------------------
         1.0 find the first link to after del_idx
         --------------------------------------------------*/
        int* p_2_first_aft_del=head_atm+cell_vec_p->begin()[natms_lcl-1];
        while(*p_2_first_aft_del<del_idx)
            p_2_first_aft_del=next_vec_p->begin()+*p_2_first_aft_del;
        
        /*--------------------------------------------------
         1.1 find the link to natms_lcl-1
         --------------------------------------------------*/
        int* p_2_last=p_2_first_aft_del;
        while(*p_2_last!=natms_lcl-1)
            p_2_last=next_vec_p->begin()+*p_2_last;
        
        
        /*--------------------------------------------------
         1.2 remove the link to natms_lcl-1 and end it ther
         --------------------------------------------------*/
        *p_2_last=-1;
        
        /*--------------------------------------------------
         1.3 insert the new link at del_idx, but for now it
         is at natms_lcl-1 after the move operation which
         happens next it will go del_idx
         --------------------------------------------------*/
        next_vec_p->begin()[natms_lcl-1]=*p_2_first_aft_del;
        *p_2_first_aft_del=del_idx;
        
    }
    type0 __gas_mass=-gas_mass;
    Algebra::DyadicV<__dim__>(__gas_mass,atoms->x_d->begin()+del_idx*__dim__,mvv_lcl);
    atoms->del(del_idx);
    ngas--;
}
/*--------------------------------------------
 things to do after a successful insertion
 --------------------------------------------*/
void PGCMC::ins_succ()
{
    atoms->add();
    for(int i=0;i<__dim__;i++) vel_buff[i]=lcl_random->gaussian()*sigma;
    memcpy(atoms->x->begin()+(natms_lcl-1)*__dim__,s_x_buff[0],__dim__*sizeof(type0));
    memcpy(atoms->x_d->begin()+(natms_lcl-1)*__dim__,vel_buff,__dim__*sizeof(type0));
    atoms->elem->begin()[natms_lcl-1]=gas_type;
    atoms->id->begin()[natms_lcl-1]=new_id;
    if(tag_vec_p) tag_vec_p->begin()[natms_lcl-1]=-1;
    if(!dof_empty)
    {
        bool* dof=atoms->dof->begin()+(natms_lcl-1)*__dim__;
        for(int i=0;i<__dim__;i++) dof[i]=true;
    }
    
    
    Algebra::DyadicV<__dim__>(gas_mass,atoms->x_d->begin()+(natms_lcl-1)*__dim__,mvv_lcl);
    
    memcpy(s_vec_p->begin()+(natms_lcl-1)*__dim__,s_buff[0],__dim__*sizeof(type0));
    
    int cell_=0;
    for(int i=0;i<__dim__;i++) cell_+=B_cells[i]*cell_coord_buff[0][i];
    cell_vec_p->begin()[natms_lcl-1]=cell_;
    
    
    int* nxt_p=head_atm+cell_;
    while(*nxt_p!=-1)
        nxt_p=next_vec_p->begin()+*nxt_p;
    *nxt_p=natms_lcl-1;
    next_vec_p->begin()[natms_lcl-1]=-1;
    ngas++;
}

/*--------------------------------------------
 setup the communications
 --------------------------------------------*/
void PGCMC::comms_setup(int n_vars_,int n_s_)
{
    n_vars=n_vars_;
    n_s=n_s_;
    ip=atoms->comm_rank;
    n_p=atoms->comm_size;
 
    memcpy(p_vec,atoms->comm.coords,sizeof(int)*__dim__);
    
    
    /*-----------------------------------------------------------------------------------------
     n_p:           total number of processor
     n_prll:        number of parallel operations
     n_pcomm:       number of processors in each individual communication
     n_comm:        number of communication that each processor has
     max_n_cncrcy:  maximum number of possible concurrent communications that a processor has
     -----------------------------------------------------------------------------------------*/
    
    
    
    max_n_cncrcy=1;
    n_p=n_prll=n_pcomm=n_comm=1;
    int comm_id_sz=0;
    for(int i=0;i<__dim__;i++)
    {
        N_p[i]=atoms->comm.dims[i];
        N_c[i]=static_cast<int>(ceil(static_cast<type0>(N_p[i])*cut_s[i]));
        N_s[i]=static_cast<int>(ceil(static_cast<type0>(n_s*N_p[i])*cut_s[i]));
        N_prll[i]=MAX(N_p[i]/(N_s[i]+1),1);
        if(N_prll[i]>1)
        {
            N_pcomm[i]=2*N_c[i]+1;
            N_comm[i]=2*N_c[i]+1;
            prll_dim[i]=true;
            max_N_cncrcy[i]=(2*N_c[i])/(N_s[i]+1)+1;
            //1<=max_N_cncrcy[i] && max_N_cncrcy[i]<=2
        }
        else
        {
            N_pcomm[i]=N_p[i];
            N_comm[i]=1;
            prll_dim[i]=false;
            max_N_cncrcy[i]=1;
        }
        
        B_p[i]=n_p;
        B_prll[i]=n_prll;
        B_pcomm[i]=n_pcomm;
        B_comm[i]=n_comm;
        
        n_p*=N_p[i];
        n_prll*=N_prll[i];
        n_pcomm*=N_pcomm[i];
        n_comm*=N_comm[i];
    
        max_n_cncrcy*=max_N_cncrcy[i];;
        
        comm_id_sz+=max_N_cncrcy[i];
    }
    
    comm_id=new int*[__dim__];
    *comm_id=new int[comm_id_sz];
    for(int i=1;i<__dim__;i++) comm_id[i]=comm_id[i-1]+max_N_cncrcy[i-1];
    
    success=new int[max_n_cncrcy];
    ntrial_atms=new int[max_n_cncrcy];
    gcmc_mode=new int[max_n_cncrcy];
    curr_comms=new MPI_Comm*[max_n_cncrcy];
    
    vars_comm=new type0*[max_n_cncrcy];
    *vars_comm=new type0[max_n_cncrcy*n_vars];
    for(int i=1;i<max_n_cncrcy;i++) vars_comm[i]=vars_comm[i-1]+n_vars;
    
    lcl_vars_comm=new type0*[max_n_cncrcy];
    *lcl_vars_comm=new type0[max_n_cncrcy*n_vars];
    for(int i=1;i<max_n_cncrcy;i++) lcl_vars_comm[i]=lcl_vars_comm[i-1]+n_vars;
    
    int_buff=new int[2+n_prll];
    
    roots=new int[n_comm];
    comms=new MPI_Comm[n_comm];
    gcmc_world=new MPI_Comm;
    split_mult();
    
}
/*--------------------------------------------
 split communications
 --------------------------------------------*/
void PGCMC::split_sing()
{
    int ip_vec[__dim__];
    for(int i=0;i<__dim__;i++)
        ip_vec[i]=0;
    
    MPI_Comm dummy;
    int jcomm,d,jkey,jcolor,jrank,jsize;
    jkey=0;
    for(int i=0;i<__dim__;i++)
        jkey+=p_vec[i]*B_p[i];
    MPI_Comm_split(world,0,jkey,gcmc_world);
    for(int i=0;i<n_p;i++)
    {
        jcomm=0;
        jkey=0;
        jcolor=0;
        for(int j=0;j<__dim__ && jcolor!=MPI_UNDEFINED;j++)
        {
            
            if(prll_dim[j])
            {
                d=p_vec[j]-ip_vec[j];
                if(d<0) d+=N_p[j];
            
                if(d<2*N_c[j]+1)
                {
                    jcolor+=ip_vec[j]*B_p[j];
                    jcomm+=d*B_comm[j];
                    d-=N_c[j];
                    if(d<0) d+=2*N_c[i]+1;
                    jkey+=d*B_pcomm[j];
                    
                }
                else
                    jcolor=MPI_UNDEFINED;
            }
            else
                jkey+=p_vec[j]*B_pcomm[j];
        }
        
        if(jcolor!=MPI_UNDEFINED)
        {
            MPI_Comm_split(world,jcolor,jkey,&comms[jcomm]);
            //Sanity check
            MPI_Comm_rank(comms[jcomm],&jrank);
            MPI_Comm_size(comms[jcomm],&jsize);
            
            /*
            if(jrank!=jkey || n_pcomm!=jsize)
                Error::abort_sing("[%d,%d,%d] %d/%d!=%d/%d color=%d",p_vec[0],p_vec[1],p_vec[2],jkey,n_pcomm,jrank,jsize,jcolor);
             */
        }
        else
            MPI_Comm_split(world,jcolor,ip,&dummy);
        
        ip_vec[0]++;
        for(int j=0;j<__dim__-1;j++)
            if(ip_vec[j]==N_p[j])
            {
                ip_vec[j]=0;
                ip_vec[j+1]++;
            }
    }
}
/*--------------------------------------------
 split communications
 --------------------------------------------*/
void PGCMC::split_mult()
{
    int hi[__dim__];
    int mid[__dim__];
    int cntr[__dim__];
    int no=1;
    for(int i=0;i<__dim__;i++)
    {
        if(prll_dim[i])
        {
            mid[i]=2*N_c[i]+1;
            hi[i]=mid[i]+(N_p[i]%(2*N_c[i]+1));
        }
        else
            mid[i]=hi[i]=1;
        
        no*=hi[i];
        cntr[i]=0;
    }
    
    MPI_Comm dummy;
    int jcomm,d,t,r,s,l,jkey,jcolor,jrank,jsize;
    jkey=0;
    for(int i=0;i<__dim__;i++)
        jkey+=p_vec[i]*B_p[i];
    MPI_Comm_split(world,0,jkey,gcmc_world);
    
    for(int i=0;i<no;i++)
    {
        jcomm=0;
        jkey=0;
        jcolor=0;
        for(int j=0;j<__dim__ && jcolor!=MPI_UNDEFINED;j++)
        {
            
            if(prll_dim[j])
            {
                
                if(cntr[j]<mid[j])
                {
                    d=p_vec[j]-cntr[j];
                    l=N_p[j]/mid[j];
                }
                else
                {
                    d=p_vec[j]-((cntr[j]-mid[j])+(N_p[j]/mid[j])*mid[j]);
                    l=1;
                }
                
                if(d<0)
                    d+=N_p[j];
                t=d/mid[j];
                if(t<l)
                {
                    r=d%mid[j];
                    jcomm+=r*B_comm[j];
                    s=p_vec[j]+N_c[j]-r;
                    if(s<0)
                        s+=N_p[j];
                    if(s>=N_p[j])
                        s-=N_p[j];
                    jcolor+=B_p[j]*s;
                    
                    r-=N_c[j];
                    if(r<0) r+=mid[j];
                    jkey+=r*B_pcomm[j];
                    
                }
                else
                    jcolor=MPI_UNDEFINED;
            }
            else
                jkey+=p_vec[j]*B_pcomm[j];
        }
        
        if(jcolor!=MPI_UNDEFINED)
        {
            MPI_Comm_split(world,jcolor,jkey,&comms[jcomm]);
            //Sanity check
            MPI_Comm_rank(comms[jcomm],&jrank);
            MPI_Comm_size(comms[jcomm],&jsize);
            /*
            if(jrank!=jkey || n_pcomm!=jsize)
                Error::abort_sing("[%d,%d,%d] %d/%d!=%d/%d color=%d",p_vec[0],p_vec[1],p_vec[2],jkey,n_pcomm,jrank,jsize,jcolor);
             */
        }
        else
            MPI_Comm_split(world,jcolor,ip,&dummy);
        
        cntr[0]++;
        for(int j=0;j<__dim__-1;j++)
            if(cntr[j]==hi[j])
            {
                cntr[j]=0;
                cntr[j+1]++;
            }
    }
}
/*--------------------------------------------
 dismantle communications
 --------------------------------------------*/
void PGCMC::comms_dismantle()
{
    if(comm_id) delete [] *comm_id;
    delete [] comm_id;
    comm_id=NULL;
    delete [] success;
    success=NULL;
    delete [] ntrial_atms;
    ntrial_atms=NULL;
    delete [] gcmc_mode;
    gcmc_mode=NULL;
    delete [] curr_comms;
    curr_comms=NULL;
    if(vars_comm) delete *vars_comm;
    delete [] vars_comm;
    vars_comm=NULL;
    if(lcl_vars_comm) delete *lcl_vars_comm;
    delete [] lcl_vars_comm;
    lcl_vars_comm=NULL;
    delete [] int_buff;
    int_buff=NULL;
    delete [] roots;
    roots=NULL;
    if(gcmc_world) MPI_Comm_free(gcmc_world);
    delete gcmc_world;
    gcmc_world=NULL;
    for(int i=0;i<n_comm;i++) MPI_Comm_free(&comms[i]);
    delete [] comms;
    comms=NULL;
    
    n_comm=0;
}
/*--------------------------------------------
 create random communication pattern
 --------------------------------------------*/
void PGCMC::create_comm_pattern()
{
    int n,d,r,l;
    n_curr_comms=1;
    im_root=true;
    origin_p=0;
    int iprll=0;
    for(int i=0;i<__dim__;i++)
    {

        op_vec[i]=static_cast<int>(N_p[i]*random->uniform());
        origin_p+=op_vec[i]*B_p[i];

        N_curr_comms[i]=0;
        if(!prll_dim[i])
        {
            comm_id[i][0]=op_vec[i];
            if(op_vec[i]!=p_vec[i]) im_root=false;
            N_curr_comms[i]=1;
            continue;
        }
        
        n=p_vec[i]-op_vec[i];
        if(n<0) n+=N_p[i];
        d=n/(N_s[i]+1);
        l=n%(N_s[i]+1);
        
        
        if(d==N_prll[i])
        {
            l+=N_s[i]+1;
            d--;
        }
        iprll+=d*B_prll[i];
        

        
        
        if(l!=0) im_root=false;
        
        if(d==N_prll[i]-1)
            r=N_p[i]-n;
        else
            r=N_s[i]+1-l;
        
        if((d%2==0 && d!=N_prll[i]-1) || (d%2==1 && d==N_prll[i]-2))
        {
            //taking care of my left
            if(l<=N_c[i])
                comm_id[i][N_curr_comms[i]++]=N_c[i]+l;
            //taking care of my right
            if(r<=N_c[i])
                comm_id[i][N_curr_comms[i]++]=N_c[i]-r;
        }
        else
        {
            //reverse the order if odd
            //taking care of my right
            if(r<=N_c[i])
                comm_id[i][N_curr_comms[i]++]=N_c[i]-r;
            //taking care of my left
            if(l<=N_c[i])
                comm_id[i][N_curr_comms[i]++]=N_c[i]+l;
        }
        
        
        n_curr_comms*=N_curr_comms[i];
    }
    
    
    for(int i=0;i<__dim__;i++)
        i_curr_comms[i]=0;
    
    int ic,iroot;
    for(int i=0;i<n_curr_comms;i++)
    {
        ic=0;
        iroot=0;
        for(int j=0;j<__dim__;j++)
        {
            if(prll_dim[j])
                ic+=comm_id[j][i_curr_comms[j]]*B_comm[j];
            else
                iroot+=comm_id[j][i_curr_comms[j]]*B_pcomm[j];
        }

        roots[i]=iroot;
        curr_comms[i]=&comms[ic];
        
        i_curr_comms[0]++;
        for(int j=0;j<__dim__-1;j++)
            if(i_curr_comms[j]==N_curr_comms[j])
            {
                i_curr_comms[j]=0;
                i_curr_comms[j+1]++;
            }
    }
    

    if(im_root)
    {
        

        int next=iprll-1;
        int prev=iprll+1;
        
        if(next!=-1)
        {
            next_p=0;
            for(int i=__dim__-1;i>-1;i--)
            {
                next_p+=((op_vec[i]+(next/B_prll[i])*(N_s[i]+1))%N_p[i])*B_p[i];
                next%=B_prll[i];
            }
        }
        else
            next_p=-1;
        
        if(prev!=n_prll)
        {
            prev_p=0;
            for(int i=__dim__-1;i>-1;i--)
            {
                prev_p+=((op_vec[i]+(prev/B_prll[i])*(N_s[i]+1))%N_p[i])*B_p[i];
                prev%=B_prll[i];
            }
        }
        else
            prev_p=-1;
        

    }
}
/*--------------------------------------------
 
 --------------------------------------------*/
void PGCMC::decide()
{
    if(prev_p!=-1)
    {
        MPI_Recv(&int_buff_sz,1,MPI_INT,prev_p,0,*gcmc_world,MPI_STATUS_IGNORE);
        MPI_Recv(int_buff,int_buff_sz,MPI_INT,prev_p,1,*gcmc_world,MPI_STATUS_IGNORE);
    }
    else
    {
        int_buff_sz=2;
        int_buff[0]=0;
        int_buff[1]=del_ids_sz;
    }
    
    success[0]=-1;
    type0 fac;
    if(gcmc_mode[0]==INS_MODE)
    {
        fac=z_fac*vol/((static_cast<type0>(tot_ngas+int_buff[0])+1.0)*exp(beta*delta_u));
        if(lcl_random->uniform()<fac)
        {
            
#ifdef GCMCDEBUG
            tot_delta_u_lcl+=delta_u;
#endif
            root_succ=true;
            success[0]=0;
            int_buff[0]++;
            
            if(int_buff_sz>2)
            {
                new_id=int_buff[int_buff_sz-1];
                int_buff_sz--;
            }
            else
            {
                if(int_buff[1]>0)
                    new_id=del_ids[int_buff[1]-1];
                else
                    new_id=max_id+1-int_buff[1];
                
                int_buff[1]--;
            }


            
        }

    }
    else
    {
        fac=static_cast<type0>(tot_ngas+int_buff[0])*exp(beta*delta_u)/(z_fac*vol);
        
        if(lcl_random->uniform()<fac)
        {
#ifdef GCMCDEBUG
            tot_delta_u_lcl-=delta_u;
#endif
            root_succ=true;
            success[0]=0;
            int_buff[0]--;
            int_buff[int_buff_sz]=gas_id;
            int_buff_sz++;
        }
    }
    
    if(next_p!=-1)
    {
        MPI_Send(&int_buff_sz,1,MPI_INT,next_p,0,*gcmc_world);
        MPI_Send(int_buff,int_buff_sz,MPI_INT,next_p,1,*gcmc_world);
    }
    
    if(root_succ)
    {
        if(gcmc_mode[0]==INS_MODE) ins_succ();
        else del_succ();
    }
}
/*--------------------------------------------
 
 --------------------------------------------*/
void PGCMC::finalize()
{
    MPI_Bcast(&int_buff_sz,1,MPI_INT,origin_p,*gcmc_world);
    MPI_Bcast(int_buff,int_buff_sz,MPI_INT,origin_p,*gcmc_world);
    
    
    if(int_buff[1]>0)
        del_ids_sz=int_buff[1];
    else
    {
        del_ids_sz=0;
        max_id-=int_buff[1];
    }
    add_del_id(int_buff+2,int_buff_sz-2);

    tot_ngas+=int_buff[0];
    atoms->natms+=int_buff[0];
    dof_diff+=int_buff[0]*__dim__;
}




