#include "dynamic_md.h"
#include "comm.h"
#include "atoms_md.h"
#include "xmath.h"
#include "timer.h"
#include "MAPP.h"
#include "ff_styles.h"
#include "neighbor.h"
using namespace MAPP_NS;
/*--------------------------------------------
 
 --------------------------------------------*/
DynamicMD::DynamicMD(AtomsMD* __atoms,ForceFieldMD* __ff,
bool __box_chng,
vec* const * updt_vecs_,int nupdt_vecs_,
vec* const * xchng_vecs_,int nxchng_vecs_,
vec* const * arch_vecs_,int narch_vecs_):
atoms(__atoms),
world(__atoms->comm.world),
skin(__atoms->comm.skin),
ff(__ff),
box_chng(__box_chng)
{
    auto is_in=[](const vec* v,vec* const * vs,int nvs)->bool
    {
        for(int i=0;i<nvs;i++)
            if(vs[i]==v) return true;
        return false;
    };
    
    int new_nvecs=atoms->nvecs-narch_vecs_;
    vec** new_vecs=new vec*[new_nvecs];
    new_vecs[0]=atoms->x;
    memcpy(new_vecs+1,updt_vecs_,nupdt_vecs_*sizeof(vec*));
    new_vecs[nupdt_vecs_+1]=atoms->id;
    memcpy(new_vecs+nupdt_vecs_+2,xchng_vecs_,nxchng_vecs_*sizeof(vec*));
    int ivec=new_nvecs-1;
    vec** vecs=atoms->vecs;
    
    
    for(int i=0;i<atoms->nvecs;i++)
    {
        if(is_in(vecs[i],new_vecs,new_nvecs) || is_in(vecs[i],arch_vecs_,narch_vecs_))
            continue;
        new_vecs[ivec--]=vecs[i];
    }
    
    nxchng_vecs=nxchng_vecs_+nupdt_vecs_+2;
    nupdt_vecs=nupdt_vecs_+1;
    delete [] atoms->vecs;
    atoms->vecs=new_vecs;
    atoms->nvecs=new_nvecs;
    

    
    narch_vecs=narch_vecs_;
    arch_vecs=NULL;
    if(!narch_vecs) return;
    
    arch_vecs=new vec*[narch_vecs];
    memcpy(arch_vecs,arch_vecs_,narch_vecs*sizeof(vec*));
}
/*--------------------------------------------
 
 --------------------------------------------*/
DynamicMD::DynamicMD(AtomsMD* __atoms,ForceFieldMD* __ff,
bool __box_chng,std::initializer_list<vec*> v0,
std::initializer_list<vec*> v1
,std::initializer_list<vec*> v2):
DynamicMD(__atoms,__ff,__box_chng,v0.begin(),
static_cast<int>(v0.size()),v1.begin(),
static_cast<int>(v1.size()),v2.begin(),
static_cast<int>(v2.size()))
{
}
/*--------------------------------------------
 
 --------------------------------------------*/
DynamicMD::DynamicMD(AtomsMD* __atoms,ForceFieldMD* __ff,
bool __box_chng,std::initializer_list<vec*> v0,
std::initializer_list<vec*> v1):
DynamicMD(__atoms,__ff,__box_chng,v0.begin(),
static_cast<int>(v0.size()),v1.begin(),
static_cast<int>(v1.size()),NULL,0)
{
}
/*--------------------------------------------
 
 --------------------------------------------*/
void DynamicMD::add_xchng(vec* v)
{
    vec** vecs=atoms->vecs;
    int nvces=atoms->nvecs;
    int ivec=nxchng_vecs;
    for(;vecs[ivec]!=v && ivec<nvces;ivec++){}
    vecs[ivec]=vecs[nxchng_vecs];
    vecs[nxchng_vecs++]=v;
}
/*--------------------------------------------
 
 --------------------------------------------*/
void DynamicMD::add_updt(vec* v)
{
    vec** vecs=atoms->vecs;
    int nvces=atoms->nvecs;
    int ivec=nupdt_vecs;
    for(;vecs[ivec]!=v && ivec<nvces;ivec++){}
    vecs[ivec]=vecs[nupdt_vecs];
    vecs[nupdt_vecs++]=v;
    if(ivec<nxchng_vecs)
        nxchng_vecs++;
}
/*--------------------------------------------
 
 --------------------------------------------*/
DynamicMD::~DynamicMD()
{
    delete [] arch_vecs;
}
/*--------------------------------------------
 init a simulation
 --------------------------------------------*/
void DynamicMD::init()
{
    store_arch_vecs();
    x0=new Vec<type0>(atoms,__dim__);

    ff->dynamic=this;
    ff->init();
    ff->neighbor->init();
    atoms->max_cut=ff->max_cut+atoms->comm.skin;
    
    xchng=new Exchange(atoms,nxchng_vecs);
    updt=new Update(atoms,nupdt_vecs,nxchng_vecs);
    
    atoms->x2s_lcl();
    xchng->full_xchng();
    updt->reset();
    updt->list();
    ff->neighbor->create_list(true);
    store_x0();
}
/*--------------------------------------------
 update one vectors
 --------------------------------------------*/
void DynamicMD::fin()
{
    ff->neighbor->fin();
    ff->fin();
    
    delete updt;
    delete xchng;
    delete x0;
    
    restore_arch_vecs();
    for(int ivec=0;ivec<atoms->nvecs;ivec++)
    {
        atoms->vecs[ivec]->vec_sz=atoms->natms;
        atoms->vecs[ivec]->shrink_to_fit();
    }
}
/*--------------------------------------------
 
 --------------------------------------------*/
void DynamicMD::store_x0()
{
    int last_atm=atoms->natms;
    if(box_chng) last_atm+=atoms->natms_ph;
    memcpy(x0->begin(),atoms->x->begin(),last_atm*__dim__*sizeof(type0));
}
/*--------------------------------------------
 
 --------------------------------------------*/
inline bool DynamicMD::decide()
{
    type0 skin_sq=0.25*skin*skin;
    int succ,succ_lcl=1;
    type0* x_vec=atoms->x->begin();
    type0* x0_vec=x0->begin();
    int last_atm=atoms->natms;
    if(box_chng) last_atm+=atoms->natms_ph;
    for(int iatm=0;succ_lcl && iatm<last_atm;iatm++,x0_vec+=__dim__,x_vec+=__dim__)
        if(Algebra::RSQ<__dim__>(x0_vec,x_vec)>skin_sq)
            succ_lcl=0;

    MPI_Allreduce(&succ_lcl,&succ,1,MPI_INT,MPI_MIN,world);
    if(succ) return true;
    return false;
}
/*--------------------------------------------
 update one vectors
 --------------------------------------------*/
void DynamicMD::update(vec* updt_vec)
{
    update(&updt_vec,1);
}
/*--------------------------------------------
 update a number of vectors
 --------------------------------------------*/
void DynamicMD::update(vec** updt_vecs,int nupdt_vecs)
{
    bool x_xst=false;
    for(int ivec=0;x_xst==false && ivec<nupdt_vecs;ivec++)
        if(updt_vecs[ivec]==atoms->x)
            x_xst=true;
    if(x_xst==false)
    {
        if(nupdt_vecs==1)
            updt->update(updt_vecs[0],false);
        else
            updt->update(updt_vecs,nupdt_vecs,false);
        return;
    }
    
    
    if(box_chng)
    {
        if(nupdt_vecs==1)
            updt->update(atoms->x,true);
        else
            updt->update(updt_vecs,nupdt_vecs,true);

        if(decide())
            return;
        
        atoms->x2s_lcl();
        xchng->full_xchng();
        
        updt->reset();
        updt->list();
        ff->neighbor->create_list(box_chng);
        store_x0();
    }
    else
    {
        if(decide())
        {
            if(nupdt_vecs==1)
                updt->update(atoms->x,true);
            else
                updt->update(updt_vecs,nupdt_vecs,true);
            return;
        }

        atoms->x2s_lcl();
        xchng->full_xchng();
        
        updt->list();
        ff->neighbor->create_list(box_chng);
        
        store_x0();
    }
}
/*--------------------------------------------
 
 --------------------------------------------*/
void DynamicMD::init_xchng()
{
    atoms->x2s_lcl();
    xchng->full_xchng();
    updt->list();
}
/*--------------------------------------------
 
 --------------------------------------------*/
void DynamicMD::fin_xchng()
{
    updt->list();
    ff->neighbor->create_list(box_chng);
    store_x0();
}
/*--------------------------------------------
 
 --------------------------------------------*/
void DynamicMD::store_arch_vecs()
{
    if(!narch_vecs) return;
    id_arch=new Vec<unsigned int>(atoms,1);
    unsigned int* id_0=atoms->id->begin();
    unsigned int* id_1=id_arch->begin();
    memcpy(id_1,id_0,atoms->natms*sizeof(unsigned int));
    atoms->pop(id_arch);
}
/*--------------------------------------------
 
 --------------------------------------------*/
#include "memory.h"
void DynamicMD::restore_arch_vecs()
{
    if(!narch_vecs) return;
    
    const int tot_p=Communication::get_size(world);
    const int my_p=Communication::get_rank(world);
    MPI_Comm& __world=world;
    
    
    int byte_sz=0;
    for(int ivec=0;ivec<narch_vecs;ivec++)
        byte_sz+=arch_vecs[ivec]->byte_sz;
    
    int natms_old=id_arch->vec_sz;
    unsigned int* id_old=id_arch->begin();
    int* key_old=NULL;
    if(natms_old) key_old=new int[natms_old];
    for(int i=0;i<natms_old;i++) key_old[i]=i;
    XMath::quicksort(key_old,key_old+natms_old
    ,[&id_old](int* rank_i,int* rank_j){return (id_old[*rank_i]<id_old[*rank_j]);}
    ,[](int* rank_i,int* rank_j){std::swap(*rank_i,*rank_j);}
    );
    
    int natms_new=atoms->id->vec_sz;
    unsigned int* id_new=atoms->id->begin();
    int* key_new=NULL;
    if(natms_new) key_new=new int[natms_new];
    for(int i=0;i<natms_new;i++) key_new[i]=i;
    XMath::quicksort(key_new,key_new+natms_new
    ,[&id_new](int* rank_i,int* rank_j){return (id_new[*rank_i]<id_new[*rank_j]);}
    ,[](int* rank_i,int* rank_j){std::swap(*rank_i,*rank_j);}
    );
    
    
    int nsnd=0;
    int* idx_lst_snd=NULL;
    unsigned int* id_lst_snd=NULL;
    if(natms_old)
    {
        idx_lst_snd=new int[natms_old];
        id_lst_snd=new unsigned int[natms_old];
    }
    
    int nrcv=0;
    int* idx_lst_rcv=NULL;
    unsigned int* id_lst_rcv=NULL;
    if(natms_new)
    {
        idx_lst_rcv=new int[natms_new];
        id_lst_rcv=new unsigned int[natms_new];
    }
    
    int nkeep=0;
    int* idx_keep_old=NULL;
    int* idx_keep_new=NULL;
    if(MIN(natms_old,natms_new))
    {
        idx_keep_old=new int[MIN(natms_old,natms_new)];
        idx_keep_new=new int[MIN(natms_old,natms_new)];
    }
    
    int i_old=0,i_new=0;
    for(;i_old<natms_old && i_new<natms_new;)
    {
        if(id_old[key_old[i_old]]<id_new[key_new[i_new]])
        {
            id_lst_snd[nsnd]=id_old[key_old[i_old]];
            idx_lst_snd[nsnd++]=key_old[i_old++];
        }
        else if(id_old[key_old[i_old]]>id_new[key_new[i_new]])
        {
            id_lst_rcv[nrcv]=id_new[key_new[i_new]];
            idx_lst_rcv[nrcv++]=key_new[i_new++];
        }
        else
        {
            idx_keep_old[nkeep]=key_old[i_old++];
            idx_keep_new[nkeep++]=key_new[i_new++];
        }
    }
    
    for(;i_old<natms_old;)
    {
        id_lst_snd[nsnd]=id_old[key_old[i_old]];
        idx_lst_snd[nsnd++]=key_old[i_old++];
    }
    for(;i_new<natms_new;)
    {
        id_lst_rcv[nrcv]=id_new[key_new[i_new]];
        idx_lst_rcv[nrcv++]=key_new[i_new++];
    }
    
    delete [] key_old;
    delete [] key_new;
    
    Memory::shrink_to_fit(id_lst_snd,nsnd,natms_old);
    Memory::shrink_to_fit(idx_lst_snd,nsnd,natms_old);
    Memory::shrink_to_fit(id_lst_rcv,nrcv,natms_new);
    Memory::shrink_to_fit(idx_lst_rcv,nrcv,natms_new);
    Memory::shrink_to_fit(idx_keep_old,nkeep,MIN(natms_old,natms_new));
    Memory::shrink_to_fit(idx_keep_new,nkeep,MIN(natms_old,natms_new));

    
    auto sort_idx_by_p=
    [&__world,&my_p,&tot_p] (const unsigned  int* lst,int lst_sz,const unsigned int* sub_lst,int sub_lst_sz,
            int*& sub_lst_idx,int* comm_size)->void
    {
        for(int ip=0;ip<tot_p;ip++)
            comm_size[ip]=0;
        
        int max_lst_sz;
        MPI_Allreduce(&lst_sz,&max_lst_sz,1,MPI_INT,MPI_MAX,__world);
        int mother_lst_sz;
        unsigned int* mother_lst=NULL;
        if(max_lst_sz) mother_lst=new unsigned int[max_lst_sz];
        
        int fnd_sz,ufnd_sz;
        
        unsigned int* ufnd=NULL;
        int* ufnd_idx=NULL;
        int* fnd_idx=NULL;
        if(sub_lst_sz)
        {
            ufnd=new unsigned int[sub_lst_sz];
            ufnd_idx=new int[sub_lst_sz];
            fnd_idx=new int[sub_lst_sz];
        }
        
        fnd_sz=0;
        ufnd_sz=sub_lst_sz;
        memcpy(ufnd,sub_lst,sub_lst_sz*sizeof(unsigned int));
        memcpy(ufnd_idx,sub_lst_idx,sub_lst_sz*sizeof(int));
        for(int imother_lst,iufnd,ufnd_tmp_sz,ip=0;ip<tot_p;ip++)
        {
            if(ip==my_p)
            {
                mother_lst_sz=lst_sz;
                memcpy(mother_lst,lst,mother_lst_sz*sizeof(unsigned int));
            }
            MPI_Bcast(&mother_lst_sz,1,MPI_INT,ip,__world);
            MPI_Bcast(mother_lst,mother_lst_sz*sizeof(unsigned int),MPI_BYTE,ip,__world);
            
            if(ip==my_p || ufnd_sz==0 || mother_lst_sz==0) continue;
            
            ufnd_tmp_sz=0;
            for(imother_lst=0,iufnd=0;imother_lst<mother_lst_sz && iufnd<ufnd_sz;)
            {
                if(mother_lst[imother_lst]<ufnd[iufnd])
                    imother_lst++;
                else if(mother_lst[imother_lst]>ufnd[iufnd])
                {
                    ufnd_idx[ufnd_tmp_sz]=ufnd_idx[iufnd];
                    ufnd[ufnd_tmp_sz++]=ufnd[iufnd++];
                }
                else
                {
                    fnd_idx[fnd_sz++]=ufnd_idx[iufnd++];
                    imother_lst++;
                    comm_size[ip]++;
                }
            }
            
            for(;iufnd<ufnd_sz;)
            {
                ufnd_idx[ufnd_tmp_sz]=ufnd_idx[iufnd];
                ufnd[ufnd_tmp_sz++]=ufnd[iufnd++];
            }
            
            ufnd_sz=ufnd_tmp_sz;
            
        }
        
        delete [] sub_lst_idx;
        sub_lst_idx=fnd_idx;
        delete [] ufnd;
        delete [] ufnd_idx;
        delete [] mother_lst;
    };
   
    int* nsnd_comm=new int[tot_p];
    int* nrcv_comm=new int[tot_p];
    sort_idx_by_p(id_lst_rcv,nrcv,id_lst_snd,nsnd,idx_lst_snd,nsnd_comm);
    sort_idx_by_p(id_lst_snd,nsnd,id_lst_rcv,nrcv,idx_lst_rcv,nrcv_comm);
    delete [] id_lst_snd;
    delete [] id_lst_rcv;
 

    byte** snd_buff=new byte*[tot_p];
    *snd_buff=NULL;
    if(nsnd)
    {
        *snd_buff=new byte[byte_sz*nsnd];
        for(int ip=1;ip<tot_p;ip++)
            snd_buff[ip]=snd_buff[ip-1]+byte_sz*nsnd_comm[ip-1];
    }
    
    byte** rcv_buff=new byte*[tot_p];
    *rcv_buff=NULL;
    if(nrcv)
    {
        *rcv_buff=new byte[byte_sz*nrcv];
        for(int ip=1;ip<tot_p;ip++)
            rcv_buff[ip]=rcv_buff[ip-1]+byte_sz*nrcv_comm[ip-1];
    }
    
    byte* _snd_buff=*snd_buff;
    for(int ivec=0;ivec<narch_vecs;ivec++)
        for(int i=0;i<nsnd;i++)
            arch_vecs[ivec]->cpy(_snd_buff,idx_lst_snd[i]);
    delete [] idx_lst_snd;
    
    
    for(int ivec=0;ivec<narch_vecs;ivec++)
        arch_vecs[ivec]->rearrange(idx_keep_old,idx_keep_new,nkeep,natms_new);
    delete [] idx_keep_old;
    delete [] idx_keep_new;


    for(int idisp=1;idisp<tot_p;idisp++)
    {
        int rcv_p=my_p-idisp;
        if(rcv_p<0) rcv_p+=tot_p;
        int snd_p=my_p+idisp;
        if(snd_p>=tot_p) snd_p-=tot_p;
        MPI_Sendrecv(snd_buff[snd_p],nsnd_comm[snd_p]*byte_sz,MPI_BYTE,snd_p,0,
                     rcv_buff[rcv_p],nrcv_comm[rcv_p]*byte_sz,MPI_BYTE,rcv_p,0,
                     world,MPI_STATUS_IGNORE);

    }
    delete [] nsnd_comm;
    delete [] nrcv_comm;
    
    delete [] *snd_buff;
    delete [] snd_buff;
    
    byte* _rcv_buff=*rcv_buff;
    for(int ivec=0;ivec<narch_vecs;ivec++)
        for(int i=0;i<nrcv;i++)
            arch_vecs[ivec]->pst_to(_rcv_buff,idx_lst_rcv[i]);
    
    delete [] idx_lst_rcv;
    delete [] *rcv_buff;
    delete [] rcv_buff;
    
    for(int ivec=0;ivec<narch_vecs;ivec++)
        atoms->push(arch_vecs[ivec]);
    atoms->push(id_arch);
    delete id_arch;
    id_arch=NULL;
}

