#include "comm.h"
#include "dynamic_dmd.h"
#include "atoms_dmd.h"
#include "xmath.h"
#include "timer.h"
#include "MAPP.h"
#include "ff_styles.h"
#include "neighbor_dmd.h"
using namespace MAPP_NS;
/*--------------------------------------------
 
 --------------------------------------------*/
DynamicDMD::DynamicDMD(AtomsDMD* __atoms,ForceFieldDMD* __ff,
bool __box_chng,
vec* const * __updt_vecs,size_t __nupdt_vecs,
vec* const * __xchng_vecs,size_t __nxchng_vecs,
vec* const * __arch_vecs,size_t __narch_vecs):
atoms(__atoms),
world(__atoms->comm.world),
skin(__atoms->comm.skin),
ff(__ff),
box_chng(__box_chng),
c_dim(__atoms->c->dim),
alpha_scale(__atoms->xi[__atoms->N-1]),
empty_vecs(NULL),
nempty_xchng_vecs(0),
nempty_updt_vecs(0),
nempty_arch_vecs(0)
{
    
    auto is_in=[](const vec* v,vec* const * vs,size_t nvs)->bool
    {
        for(int i=0;i<nvs;i++)
            if(vs[i]==v) return true;
        return false;
    };
    
    
    
    int __nvecs=atoms->nvecs-static_cast<int>(__narch_vecs);
    vec** __vecs=new vec*[__nvecs];
    
    __vecs[0]=atoms->x;
    __vecs[1]=atoms->alpha;
    __vecs[2]=atoms->c;
    __vecs[3]=atoms->elem;
    memcpy(__vecs+4,__updt_vecs,__nupdt_vecs*sizeof(vec*));
    nupdt_vecs=static_cast<int>(__nupdt_vecs+4);
    
    __vecs[nupdt_vecs]=atoms->id;
    memcpy(__vecs+nupdt_vecs+1,__xchng_vecs,__nxchng_vecs*sizeof(vec*));
    nxchng_vecs=static_cast<int>(__nxchng_vecs+1+__nupdt_vecs+4);
    
    int ivec=__nvecs-1;
    vec** vecs=atoms->vecs;
    for(int i=0;i<atoms->nvecs;i++)
    {
        if(is_in(vecs[i],__vecs,__nvecs) || is_in(vecs[i],__arch_vecs,__narch_vecs))
            continue;
        __vecs[ivec--]=vecs[i];
    }
    
    delete [] atoms->vecs;
    atoms->vecs=__vecs;
    atoms->nvecs=__nvecs;
    

    narch_vecs=static_cast<int>(__narch_vecs);
    arch_vecs=NULL;
    if(narch_vecs) arch_vecs=new vec*[narch_vecs];
    memcpy(arch_vecs,__arch_vecs,sizeof(vec*)*narch_vecs);
}
/*--------------------------------------------
 
 --------------------------------------------*/
DynamicDMD::~DynamicDMD()
{
    delete [] arch_vecs;
    delete [] empty_vecs;
}
/*--------------------------------------------
 
 --------------------------------------------*/
DynamicDMD::DynamicDMD(AtomsDMD* __atoms,ForceFieldDMD* __ff,
bool __box_chng,std::initializer_list<vec*> v0,
std::initializer_list<vec*> v1
,std::initializer_list<vec*> v2):
DynamicDMD(__atoms,__ff,__box_chng,
v0.begin(),v0.size(),
v1.begin(),v1.size(),
v2.begin(),v2.size())
{
}
/*--------------------------------------------
 
 --------------------------------------------*/
DynamicDMD::DynamicDMD(AtomsDMD* __atoms,ForceFieldDMD* __ff,
bool __box_chng,std::initializer_list<vec*> v0,
std::initializer_list<vec*> v1):
DynamicDMD(__atoms,__ff,__box_chng,
v0.begin(),v0.size(),
v1.begin(),v1.size(),
NULL,0)
{
}
/*--------------------------------------------
 
 --------------------------------------------*/
void DynamicDMD::add_xchng(vec* v)
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
void DynamicDMD::add_updt(vec* v)
{
    vec** vecs=atoms->vecs;
    int nvces=atoms->nvecs;
    int ivec=nupdt_vecs;
    for(;vecs[ivec]!=v && ivec<nvces;ivec++){}
    vecs[ivec]=vecs[nupdt_vecs];
    vecs[nupdt_vecs++]=v;
    // I guess this is not correct
    //if(ivec<nxchng_vecs)
      //  nxchng_vecs++;
    // instead
    if(ivec>=nxchng_vecs)
        nxchng_vecs++;
}
/*--------------------------------------------
 init a simulation
 --------------------------------------------*/
void DynamicDMD::init()
{
    store_empty_vecs();
    store_arch_vecs();
    x0=new Vec<type0>(atoms,__dim__);
    alpha0=new Vec<type0>(atoms,c_dim);

    ff->dynamic=this;
    ff->init();
    ff->neighbor->init();
    atoms->max_cut=ff->max_cut+atoms->comm.skin+alpha_scale*sqrt_2*atoms->max_alpha;
    
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
void DynamicDMD::fin()
{
    ff->neighbor->fin();
    ff->fin();
    
    delete updt;
    delete xchng;
    delete alpha0;
    delete x0;
    
    restore_arch_vecs();
    for(int ivec=0;ivec<atoms->nvecs;ivec++)
    {
        atoms->vecs[ivec]->vec_sz=atoms->natms_lcl;
        atoms->vecs[ivec]->shrink_to_fit();
    }
    restore_empty_vecs();
    atoms->natms_ph=0;
}
/*--------------------------------------------
 
 --------------------------------------------*/
void DynamicDMD::store_x0()
{
    int last_atm=atoms->natms_lcl;
    if(box_chng) last_atm+=atoms->natms_ph;
    memcpy(x0->begin(),atoms->x->begin(),last_atm*__dim__*sizeof(type0));
    memcpy(alpha0->begin(),atoms->alpha->begin(),last_atm*c_dim*sizeof(type0));
}
/*--------------------------------------------
 
 --------------------------------------------*/
inline bool DynamicDMD::decide()
{
    int succ,succ_lcl=1;
    type0* x_vec=atoms->x->begin();
    type0* x0_vec=x0->begin();
    type0* alpha_vec=atoms->alpha->begin();
    type0* alpha0_vec=alpha0->begin();
    int last_atm=atoms->natms_lcl;
    if(box_chng) last_atm+=atoms->natms_ph;
    
    for(int iatm=0;succ_lcl && iatm<last_atm;iatm++,x0_vec+=__dim__,x_vec+=__dim__,alpha_vec+=c_dim,alpha0_vec+=c_dim)
    {
        type0 dr=sqrt(Algebra::RSQ<__dim__>(x0_vec,x_vec));
        type0 dalpha=alpha_vec[0]-alpha0_vec[0];
        for(int i=0;i<c_dim;i++)
            dalpha=MAX(dalpha,alpha_vec[i]-alpha0_vec[i]);
        
        if(dr+dalpha*alpha_scale>0.5*skin) succ_lcl=0;
    }

    MPI_Allreduce(&succ_lcl,&succ,1,MPI_INT,MPI_MIN,world);
    if(succ) return true;
    return false;
}
/*--------------------------------------------
 update one vectors
 --------------------------------------------*/
void DynamicDMD::update(vec* updt_vec)
{
    update(&updt_vec,1);
}
/*--------------------------------------------
 update a number of vectors
 --------------------------------------------*/
void DynamicDMD::update(vec** updt_vecs,int nupdt_vecs)
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
        atoms->max_cut=ff->max_cut+atoms->comm.skin+alpha_scale*sqrt_2*atoms->max_alpha;
        updt->reset();
        updt->list();
        ff->neighbor->create_list(true);
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
        atoms->max_cut=ff->max_cut+atoms->comm.skin+alpha_scale*sqrt_2*atoms->max_alpha;
        updt->reset();
        updt->list();
        ff->neighbor->create_list(true);
        store_x0();
    }
}
/*--------------------------------------------
 
 --------------------------------------------*/
void DynamicDMD::store_empty_vecs()
{
    for(int i=0;i<nupdt_vecs;i++)
        if(atoms->vecs[i]->is_empty())
            nempty_updt_vecs++;
    for(int i=nupdt_vecs;i<nxchng_vecs;i++)
        if(atoms->vecs[i]->is_empty())
            nempty_xchng_vecs++;
    
    for(int i=0;i<narch_vecs;i++)
        if(arch_vecs[i]->is_empty())
            nempty_arch_vecs++;
    int nempty_vecs=nempty_updt_vecs+
    nempty_xchng_vecs+nempty_arch_vecs;
    if(!nempty_vecs) return;
        
        
    empty_vecs=new vec*[nempty_vecs];
    
    int iempty_vec=0;
    if(nempty_xchng_vecs+nempty_updt_vecs)
    {
        int __nvecs=atoms->nvecs-
        nempty_xchng_vecs-nempty_updt_vecs;
        int ivec=0;
        vec** __vecs=new vec*[__nvecs];
        for(int i=0;i<nxchng_vecs;i++)
        {
            if(atoms->vecs[i]->is_empty())
                empty_vecs[iempty_vec++]=atoms->vecs[i];
            else
                __vecs[ivec++]=atoms->vecs[i];
        }
        
        memcpy(__vecs+ivec,atoms->vecs+nxchng_vecs,sizeof(vec*)*(atoms->nvecs-nxchng_vecs));
        delete [] atoms->vecs;
        atoms->vecs=__vecs;
        atoms->nvecs=__nvecs;
        nxchng_vecs-=nempty_xchng_vecs+nempty_updt_vecs;
        nupdt_vecs-=nempty_updt_vecs;
    }
    
    if(!nempty_arch_vecs) return;
    int __narch_vecs=narch_vecs-nempty_arch_vecs;
    vec** __arch_vecs=new vec*[__narch_vecs];
    int ivec=0;
    
    for(int i=0;i<narch_vecs;i++)
    {
        if(arch_vecs[i]->is_empty())
            empty_vecs[iempty_vec++]=arch_vecs[i];
        else
            __arch_vecs[ivec++]=arch_vecs[i];
    }
    
    delete [] arch_vecs;
    arch_vecs=__arch_vecs;
    narch_vecs=__narch_vecs;
}
/*--------------------------------------------
 
 --------------------------------------------*/
void DynamicDMD::restore_empty_vecs()
{
    int nempty_vecs=nempty_updt_vecs+
    nempty_xchng_vecs+nempty_arch_vecs;
    if(!nempty_vecs) return;
    
    if(nempty_xchng_vecs+nempty_updt_vecs)
    {
        int __nvecs=atoms->nvecs+
        nempty_xchng_vecs+nempty_updt_vecs;
        vec** __vecs=new vec*[__nvecs];
        vec** ___vecs=__vecs;
        memcpy(___vecs,atoms->vecs,nupdt_vecs*sizeof(vec*));
        ___vecs+=nupdt_vecs;
        memcpy(___vecs,empty_vecs,nempty_updt_vecs*sizeof(vec*));
        ___vecs+=nempty_updt_vecs;
        memcpy(___vecs,atoms->vecs+nupdt_vecs,(nxchng_vecs-nupdt_vecs)*sizeof(vec*));
        ___vecs+=nxchng_vecs-nupdt_vecs;
        memcpy(___vecs,empty_vecs+nempty_updt_vecs,nempty_xchng_vecs*sizeof(vec*));
        ___vecs+=nempty_xchng_vecs;
        memcpy(___vecs,atoms->vecs+nxchng_vecs,(atoms->nvecs-nxchng_vecs)*sizeof(vec*));
        
        
        delete [] atoms->vecs;
        atoms->vecs=__vecs;
        atoms->nvecs=__nvecs;
        nxchng_vecs+=nempty_xchng_vecs+nempty_updt_vecs;
        nupdt_vecs+=nempty_updt_vecs;
    }
    
    if(nempty_arch_vecs)
    {
        int __narch_vecs=narch_vecs+nempty_arch_vecs;
        vec** __arch_vecs=new vec*[__narch_vecs];
        memcpy(__arch_vecs,arch_vecs,sizeof(vec*)*narch_vecs);
        memcpy(__arch_vecs+narch_vecs,empty_vecs+nempty_xchng_vecs+nempty_updt_vecs,sizeof(vec*)*nempty_arch_vecs);
        delete [] arch_vecs;
        arch_vecs=__arch_vecs;
        narch_vecs=__narch_vecs;
    }
    
    delete [] empty_vecs;
    empty_vecs=NULL;
    nempty_updt_vecs=nempty_xchng_vecs=nempty_arch_vecs=0;
}
/*--------------------------------------------
 
 --------------------------------------------*/
void DynamicDMD::store_arch_vecs()
{
    if(!narch_vecs) return;
    id_arch=new Vec<unsigned int>(atoms,1);
    unsigned int* id_0=atoms->id->begin();
    unsigned int* id_1=id_arch->begin();
    memcpy(id_1,id_0,atoms->natms_lcl*sizeof(unsigned int));
    atoms->pop(id_arch);
}
/*--------------------------------------------
 
 --------------------------------------------*/
#include "memory.h"
void DynamicDMD::restore_arch_vecs()
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

