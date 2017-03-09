/*--------------------------------------------
 Created by Sina on 06/05/13.
 Copyright (c) 2013 MIT. All rights reserved.
 --------------------------------------------*/
#include "elements.h"
#include "sgcmc.h"
#include "memory.h"
#include "random.h"
#include "neighbor_md.h"
#include "ff_md.h"
#include "xmath.h"
#include "atoms_md.h"
#include "MAPP.h"
#include "dynamic_md.h"
using namespace MAPP_NS;
/*--------------------------------------------
 constructor
 --------------------------------------------*/
SGCMC::SGCMC(AtomsMD*& __atoms,ForceFieldMD*&__ff,DynamicMD*& __dynamic,int __m,elem_type __gas_type,type0 __mu,type0 __T,int seed):
GCMC(__atoms,__ff,__dynamic,__gas_type,__mu,__T,seed),
m(__m)
{
    s_x_buff=NULL;

    
    head_atm=NULL;
    cell_coord_buff=NULL;
    n_cells=0;
    
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
SGCMC::~SGCMC()
{
    delete [] rel_neigh_lst_coord;
}
/*--------------------------------------------
 
 --------------------------------------------*/
void SGCMC::init()
{
    dof_empty=atoms->dof->is_empty();
    GCMC::init();
    MPI_Scan(&ngas,&ngas_before,1,MPI_INT,MPI_SUM,world);
    ngas_before-=ngas;
    box_setup();
    vars=new type0[ff->gcmc_n_vars];
    lcl_vars=new type0[ff->gcmc_n_vars];
}
/*--------------------------------------------
 
 --------------------------------------------*/
void SGCMC::fin()
{
    delete [] vars;
    delete [] lcl_vars;
    box_dismantle();
    GCMC::fin();
}
/*--------------------------------------------
 
 --------------------------------------------*/
void SGCMC::box_setup()
{
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
    cell_coord_buff=new int[max_ntrial_atms*__dim__];
    s_x_buff=new type0[__dim__*max_ntrial_atms];
}
/*--------------------------------------------
 
 --------------------------------------------*/
void SGCMC::box_dismantle()
{
    delete [] head_atm;
    delete [] cell_coord_buff;
    delete [] s_x_buff;
    
    head_atm=NULL;
    cell_coord_buff=NULL;
    s_x_buff=NULL;
    n_cells=0;
}
/*--------------------------------------------
 construct the bin list
 --------------------------------------------*/
void SGCMC::xchng(bool box_chng,int nattmpts)
{
    curr_comm=&world;
    curr_root=0;
        
    
    dynamic->init_xchng();
    if(box_chng)
        box_setup();
    
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
    MPI_Scan(&ngas,&ngas_before,1,MPI_INT,MPI_SUM,world);
    MPI_Allreduce(&ngas,&tot_ngas,1,MPI_INT,MPI_SUM,world);
    ngas_before-=ngas;
    

    next_jatm_p=&SGCMC::next_jatm_reg;
    dof_diff=0;
    for(int i=0;i<nattmpts;i++)
        attmpt();

    
    ff->fin_xchng();
    memcpy(atoms->x->begin(),s_vec_p->begin(),sizeof(type0)*natms_lcl*__dim__);
    
    delete tag_vec_p;
    delete s_vec_p;
    delete next_vec_p;
    delete cell_vec_p;
    dynamic->fin_xchng();


}
/*--------------------------------------------
 this must be used only for local atoms
 --------------------------------------------*/
inline void SGCMC::find_cell_no(type0*& s,int& cell_no)
{
    cell_no=0;
    for(int i=0;i<__dim__;i++)
        cell_no+=B_cells[i]*MIN(static_cast<int>((s[i]-s_lo[i])/cell_size[i]),N_cells[i]-1);
}
/*--------------------------------------------
 
 --------------------------------------------*/
inline void SGCMC::find_cell_coord(type0*& s,int*& cell_coord)
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
 attempt an insertion
 --------------------------------------------*/
void SGCMC::attmpt()
{
    if(random->uniform()<0.5)
    {
        xchng_mode=INS_MODE;
        im_root=true;
        int iproc=-1;
        for(int i=0;i<__dim__;i++)
        {
            s_buff[i]=random->uniform();
            if(s_buff[i]<s_lo[i] || s_buff[i]>=s_hi[i])
                im_root=false;
        }
        if(im_root)
            iproc=atoms->comm_rank;
        MPI_Allreduce(&iproc,&curr_root,1,MPI_INT,MPI_MAX,world);
    }
    else
    {
        if(tot_ngas==0)
            return;
        xchng_mode=DEL_MODE;
        igas=static_cast<int>(tot_ngas*random->uniform());
        int iproc=-1;
        im_root=false;
        if(ngas_before<=igas && igas<ngas_before+ngas)
        {
            im_root=true;
            iproc=atoms->comm_rank;
            int n=igas-ngas_before;
            int icount=-1;
            elem_type* elem=atoms->elem->begin();
            
            del_idx=0;
            for(;icount!=n;del_idx++)
                if(elem[del_idx]==gas_type) icount++;
            del_idx--;
            gas_id=atoms->id->begin()[del_idx];
            memcpy(s_buff,s_vec_p->begin()+del_idx*__dim__,__dim__*sizeof(type0));
        }
        MPI_Allreduce(&iproc,&curr_root,1,MPI_INT,MPI_MAX,world);
        MPI_Bcast(s_buff,__dim__,Vec<type0>::MPI_T,curr_root,world);
        MPI_Bcast(&gas_id,1,MPI_INT,curr_root,world);
    }
    
    prep_s_x_buff();
    if(tag_vec_p) reset_tag();
    ff->pre_xchng_energy_timer(this);
    
    delta_u=ff->xchng_energy_timer(this);
    MPI_Bcast(&delta_u,1,Vec<type0>::MPI_T,curr_root,world);
    
    type0 fac;
    root_succ=false;
    if(xchng_mode==INS_MODE)
    {
        fac=z_fac*vol/((static_cast<type0>(tot_ngas)+1.0)*exp(beta*delta_u));
        if(random->uniform()<fac)
        {
            root_succ=true;
            ins_succ();
            ff->post_xchng_energy_timer(this);
        }
    }
    else
    {
        fac=static_cast<type0>(tot_ngas)*exp(beta*delta_u)/(z_fac*vol);
        if(random->uniform()<fac)
        {
            root_succ=true;
            del_succ();
            ff->post_xchng_energy_timer(this);
        }
    }
    
}
/*--------------------------------------------
 things to do after a successful insertion
 --------------------------------------------*/
void SGCMC::ins_succ()
{

    int new_id=get_new_id();
    for(int i=0;i<__dim__;i++) vel_buff[i]=random->gaussian()*sigma;
    if(im_root)
    {
        atoms->add();
        
        memcpy(atoms->x->begin()+(natms_lcl-1)*__dim__,s_x_buff,__dim__*sizeof(type0));
        memcpy(atoms->x_d->begin()+(natms_lcl-1)*__dim__,vel_buff,__dim__*sizeof(type0));
        atoms->elem->begin()[natms_lcl-1]=gas_type;
        atoms->id->begin()[natms_lcl-1]=new_id;
        if(tag_vec_p) tag_vec_p->begin()[natms_lcl-1]=-1;
        if(!dof_empty)
        {
            bool* dof=atoms->dof->begin()+(natms_lcl-1)*__dim__;
            for(int i=0;i<__dim__;i++) dof[i]=true;
        }
        
        
        memcpy(s_vec_p->begin()+(natms_lcl-1)*__dim__,s_buff,__dim__*sizeof(type0));
        
        int cell_=0;
        for(int i=0;i<__dim__;i++) cell_+=B_cells[i]*cell_coord_buff[i];
        cell_vec_p->begin()[natms_lcl-1]=cell_;
        
        
        int* nxt_p=head_atm+cell_;
        while(*nxt_p!=-1)
            nxt_p=next_vec_p->begin()+*nxt_p;
        *nxt_p=natms_lcl-1;
        next_vec_p->begin()[natms_lcl-1]=-1;
        ngas++;
    }
    else
    {
        if(atoms->comm_rank>curr_root)
            ngas_before++;
    }
    
    dof_diff+=__dim__;
    tot_ngas++;
    atoms->natms++;    
}
/*--------------------------------------------
 things to do after a successful deletion
 --------------------------------------------*/
void SGCMC::del_succ()
{

    add_del_id(&gas_id,1);

    if(im_root)
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
        
        
        atoms->del(del_idx);
        ngas--;
    }
    else
    {
        if(igas<ngas_before)
            ngas_before--;
    }
    dof_diff-=__dim__;
    tot_ngas--;
    atoms->natms--;
}
/*--------------------------------------------
 
 --------------------------------------------*/
void SGCMC::prep_s_x_buff()
{
    
    im_root=true;
    for(int i=0;i<__dim__;i++)
        if(s_buff[i]<s_lo[i] || s_buff[i]>=s_hi[i])
            im_root=false;
    
    
    int n_per_dim[__dim__];
    type0 s;
    int no;
    ntrial_atms=1;
    for(int i=0;i<__dim__;i++)
    {
        no=0;
        s=s_buff[i];
        
        for(type0 s_=s;s_<s_hi_ph[i];s_++)
            if(s_lo_ph[i]<=s_ && s_<s_hi_ph[i])
                s_trials[i][no++]=s_;
        
        for(type0 s_=s-1.0;s_lo_ph[i]<=s_;s_--)
            if(s_lo_ph[i]<=s_ && s_<s_hi_ph[i])
                s_trials[i][no++]=s_;
        
        ntrial_atms*=no;
        n_per_dim[i]=no;
    }
    
    
    int count[__dim__];
    type0* buff=s_x_buff;
    int* cell_coord=cell_coord_buff;
    type0 (&H)[__dim__][__dim__]=atoms->H;
    for(int i=0;i<__dim__;i++) count[i]=0;
    for(int i=0;i<ntrial_atms;i++)
    {
        for(int j=0;j<__dim__;j++)
            buff[j]=s_trials[j][count[j]];
        find_cell_coord(buff,cell_coord);
        XMatrixVector::s2x<__dim__>(buff,H);
        
        count[0]++;
        for(int j=0;j<__dim__-1;j++)
            if(count[j]==n_per_dim[j])
            {
                count[j]=0;
                count[j+1]++;
            }
        cell_coord+=__dim__;
        buff+=__dim__;
    }
}
/*--------------------------------------------
 this is used for insertion trial
 --------------------------------------------*/
void SGCMC::reset_iatm()
{
    itrial_atm=-1;
    next_iatm();
}
/*--------------------------------------------
 this is used for insertion trial
 --------------------------------------------*/
void SGCMC::next_iatm()
{
    itrial_atm++;
    
    if(itrial_atm==ntrial_atms)
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
        icell_coord[i]=cell_coord_buff[__dim__*itrial_atm+i];
    
    /*--------------------------------------------------
     assign the position of iatm
     --------------------------------------------------*/
    ix=s_x_buff+itrial_atm*__dim__;
    
}
/*--------------------------------------------
 this is used for insertion trial
 --------------------------------------------*/
void SGCMC::reset_jatm()
{
    ineigh=0;
    jatm_next=-1;
    next_jatm_reg();
    
}
/*--------------------------------------------
 find the next lcl atom
 --------------------------------------------*/
inline void SGCMC::next_jatm_reg()
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
                next_jatm_p=&SGCMC::next_jatm_self;
                return next_jatm_self();
            }
            else return;
        }
        
        jatm_next=next_vec_p->begin()[jatm];
        
        if(jatm==iatm) continue;
        
        jx=atoms->x->begin()+__dim__*jatm;
        jelem=atoms->elem->begin()[jatm];
        
        rsq=XMatrixVector::rsq<__dim__>(ix,jx);
        if(rsq>=cut_sq[ielem][jelem]) continue;
        
        if(tag_vec_p) tag_vec_p->begin()[jatm]=icomm;
        return;
    }
}
/*--------------------------------------------
 find the next interactin image atom
 --------------------------------------------*/
inline void SGCMC::next_jatm_self()
{
    
    while(true)
    {
        iself++;
        if(iself==ntrial_atms)
        {
            next_jatm_p=&SGCMC::next_jatm_reg;
            jatm=-1;
            return;
        }
        
        jatm=natms_lcl+iself;
        jx=s_x_buff+iself*__dim__;
        jelem=gas_type;
        rsq=XMatrixVector::rsq<__dim__>(ix,jx);
        if(rsq>=cut_sq[ielem][jelem]) continue;
        
        return;
    }
}
/*--------------------------------------------
 find the next interactin image atom
 --------------------------------------------*/
void SGCMC::next_jatm()
{
    (this->*next_jatm_p)();
}
/*--------------------------------------------
 
 --------------------------------------------*/
void SGCMC::next_icomm()
{
    icomm=-1;
}
/*--------------------------------------------
 
 --------------------------------------------*/
void SGCMC::reset_icomm()
{
    icomm=0;
}
/*--------------------------------------------
 
 --------------------------------------------*/
inline void SGCMC::refresh()
{
    for(int i=0;i<n_cells;i++)
        head_atm[i]=-1;
    int* next_vec=next_vec_p->begin()+natms_lcl-1;
    int* cell_vec=cell_vec_p->begin()+natms_lcl-1;
    for(int i=natms_lcl-1;i>-1;i--,next_vec--,cell_vec--)
    {
        *next_vec=head_atm[*cell_vec];
        head_atm[*cell_vec]=i;
    }
}
/*--------------------------------------------
 
 --------------------------------------------*/
void SGCMC::reset_tag()
{
    int* tag=tag_vec_p->begin();
    for(int i=0;i<natms_lcl;i++)
        tag[i]=-1;
}






