#include "atoms_styles.h"
#include "dynamic.h"
#include "ff_styles.h"
#include "xmath.h"
#include "memory.h"

/*
 form simulation comes these EXCLUSIVE lists of vectors
 0. exchange
 1. update
 2. archive
 
 there are some default vectors and the ones given by the simulation should not include them
 exchange: id
 update:   x,elem
 


  ------------------------------------ -----------------     -----------------     ---------------
 | active vectors                     | passive vectors |   | archive vectors |   | empty vectors |
  ------------------ -----------------  size: natms_lcl |   | size: natms_lcl |   | size: natms   |
 |                  |                 |      +natms_ph  |   |      +natms_ph  |   |      +natms   |
 | exchange vectors | update vectors  |                 |   |                 |   |               |        
 | size: natms_lcl  | size: natms_lcl |                 |   |                 |   |               |        
 |                  |      +natms_ph  |                 |   |                 |   |               |        
  ------------------ ----------------- -----------------     -----------------     ---------------
         |                  |                                       |                      ^        
         |                  |                                       |                      |        
          ------------------ --------------------------------------- ----------------------

 passive vectors are the ones that are used for 
 manipulation they are temporary and are used in
 spot
 
 at init()
 go through exchange, update, and archive if they 
 are any empty vectors segregate them into empty
 vectors category
 
 pop all the empty vectors from vector stack.
 pop remaining archive vectors from vector stack.
 store the archive vectors along with initial id
 vector
 
 
 now sort the vector stack as follows
   0. all exchange vectors
   1. all update vectors
   2. remaining vectors (passive)
   ** in addition the VERY FIRST exchange vector 
   and the VERY FIRST update vector should be id 
   and x, respectively.
 
 
 
 
 at fin()
   0. restore the archive vectors and push
      them into stack
   1. push the empty vectors to stack
 
 */
/*--------------------------------------------
 
 --------------------------------------------*/
Dynamic::Dynamic(Atoms* __atoms,ForceField* __ff,bool __box_chng,
std::initializer_list<vec*> __updt_vecs_def,std::initializer_list<vec*> __updt_vecs,
std::initializer_list<vec*> __xchng_comp_vecs_def,std::initializer_list<vec*> __xchng_comp_vecs,
std::initializer_list<vec*> __arch_vecs_def,std::initializer_list<vec*> __arch_vecs):
atoms(__atoms),
world(__atoms->comm.world),
skin(__atoms->comm.skin),
ff(__ff),
box_chng(__box_chng),

updt_vecs(__updt_vecs_def.size()+__updt_vecs.size()==0 ? NULL: new vec*[__updt_vecs_def.size()+__updt_vecs.size()]),
nupdt_vecs(static_cast<int>(__updt_vecs.size()+__updt_vecs_def.size())),

xchng_comp_vecs(__xchng_comp_vecs_def.size()+__xchng_comp_vecs.size()==0 ? NULL: new vec*[__xchng_comp_vecs_def.size()+__xchng_comp_vecs.size()]),
nxchng_comp_vecs(static_cast<int>(__xchng_comp_vecs.size()+__xchng_comp_vecs_def.size())),

arch_vecs(__arch_vecs_def.size()+__arch_vecs.size()==0 ? NULL: new vec*[__arch_vecs_def.size()+__arch_vecs.size()]),
narch_vecs(static_cast<int>(__arch_vecs.size()+__arch_vecs_def.size())),

nupdt_vecs_full(0),
nxchng_vecs_full(0),
arch_vecs_full(NULL),
narch_vecs_full(0)
{
    memcpy(updt_vecs,__updt_vecs_def.begin(),__updt_vecs_def.size()*sizeof(vec*));
    memcpy(updt_vecs+__updt_vecs_def.size(),__updt_vecs.begin(),__updt_vecs.size()*sizeof(vec*));
    
    memcpy(xchng_comp_vecs,__xchng_comp_vecs_def.begin(),__xchng_comp_vecs_def.size()*sizeof(vec*));
    memcpy(xchng_comp_vecs+__xchng_comp_vecs_def.size(),__xchng_comp_vecs.begin(),__xchng_comp_vecs.size()*sizeof(vec*));
    
    memcpy(arch_vecs,__arch_vecs_def.begin(),__arch_vecs_def.size()*sizeof(vec*));
    memcpy(arch_vecs+__arch_vecs_def.size(),__arch_vecs.begin(),__arch_vecs.size()*sizeof(vec*));
}
/*--------------------------------------------
 
 --------------------------------------------*/
Dynamic::~Dynamic()
{
    delete [] arch_vecs;
    delete [] xchng_comp_vecs;
    delete [] updt_vecs;
}
/*--------------------------------------------
 
 --------------------------------------------*/
void Dynamic::add_xchng(vec* v)
{
    vec** __xchng_comp_vecs=new vec*[nxchng_comp_vecs+1];
    memcpy(__xchng_comp_vecs,xchng_comp_vecs,nxchng_comp_vecs*sizeof(vec*));
    delete [] xchng_comp_vecs;
    xchng_comp_vecs=__xchng_comp_vecs;
    xchng_comp_vecs[nxchng_comp_vecs++]=v;
}
/*--------------------------------------------
 
 --------------------------------------------*/
void Dynamic::add_updt(vec* v)
{
    vec** __updt_vecs=new vec*[nupdt_vecs+1];
    memcpy(__updt_vecs,updt_vecs,nupdt_vecs*sizeof(vec*));
    delete [] updt_vecs;
    updt_vecs=__updt_vecs;
    updt_vecs[nupdt_vecs++]=v;
}
/*--------------------------------------------
 
 --------------------------------------------*/
void Dynamic::create_dynamic_vecs()
{
    
    int dynamic_vecs_cpcty=atoms->nvecs;
    vec** dynamic_vecs;
    Memory::alloc(dynamic_vecs,dynamic_vecs_cpcty);
    int ndynamic_vecs=0;
    
    nupdt_vecs_full=0;
    for(int i=0;i<nupdt_vecs;i++)
        if(!updt_vecs[i]->is_empty())
        {
            dynamic_vecs[ndynamic_vecs++]=updt_vecs[i];
            nupdt_vecs_full++;
        }
    
    nxchng_vecs_full=nupdt_vecs_full;
    for(int i=0;i<nxchng_comp_vecs;i++)
        if(!xchng_comp_vecs[i]->is_empty())
        {
            dynamic_vecs[ndynamic_vecs++]=xchng_comp_vecs[i];
            nxchng_vecs_full++;
        }
    
    
    auto is_in=[](vec*& v,vec**& vs,int& nvs)->bool
    {
        
        for(int i=0;i<nvs;i++)
            if(vs[i]==v) return true;
        return false;
    };
    
    vec** all_vecs=atoms->vecs;
    int nall_vecs=atoms->nvecs;
    for(int i=0;i<nall_vecs;i++)
    {
        if(all_vecs[i]->is_empty()) continue;
        if(is_in(all_vecs[i],dynamic_vecs,nxchng_vecs_full)) continue;
        if(is_in(all_vecs[i],arch_vecs_full,narch_vecs_full)) continue;
        dynamic_vecs[ndynamic_vecs++]=all_vecs[i];
    }
    
    
    
    Memory::shrink_to_fit(dynamic_vecs,ndynamic_vecs,dynamic_vecs_cpcty);
    delete [] atoms->dynamic_vecs;
    atoms->dynamic_vecs=dynamic_vecs;
    atoms->ndynamic_vecs=ndynamic_vecs;
}
/*--------------------------------------------
 
 --------------------------------------------*/
void Dynamic::destroy_dynamic_vecs()
{
    delete [] atoms->dynamic_vecs;
    atoms->dynamic_vecs=NULL;
    atoms->ndynamic_vecs=0;
    
    nxchng_vecs_full=nupdt_vecs_full=0;
}
/*--------------------------------------------
 
 --------------------------------------------*/
void Dynamic::store_arch_vecs()
{
    
    Memory::alloc(arch_vecs_full,narch_vecs+1);
    narch_vecs_full=0;
    for(int i=0;i<narch_vecs;i++)
        if(!arch_vecs[i]->is_empty())
            arch_vecs_full[narch_vecs_full++]=arch_vecs[i];
    if(narch_vecs_full)
        arch_vecs_full[narch_vecs_full++]=new Vec<unsigned int>(atoms,1,"id_arch");
    Memory::shrink_to_fit(arch_vecs_full,narch_vecs_full,narch_vecs+1);
    if(!narch_vecs_full) return;
    
    unsigned int* id_0=atoms->id->begin();
    unsigned int* id_1=reinterpret_cast<unsigned int *>(arch_vecs_full[narch_vecs_full-1]->begin());
    memcpy(id_1,id_0,atoms->natms_lcl*sizeof(unsigned int));
}
/*--------------------------------------------
 
 --------------------------------------------*/
void Dynamic::restore_arch_vecs()
{
    if(!narch_vecs_full) return;
    
    const int tot_p=Communication::get_size(world);
    const int my_p=Communication::get_rank(world);
    MPI_Comm& __world=world;
    
    
    int byte_sz=0;
    for(int ivec=0;ivec<narch_vecs_full-1;ivec++)
        byte_sz+=arch_vecs_full[ivec]->byte_sz;
    
    Vec<unsigned int>* id_arch=dynamic_cast<Vec<unsigned int>*>(arch_vecs_full[narch_vecs_full-1]);
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
    for(int ivec=0;ivec<narch_vecs_full-1;ivec++)
        for(int i=0;i<nsnd;i++)
            arch_vecs_full[ivec]->cpy(_snd_buff,idx_lst_snd[i]);
    delete [] idx_lst_snd;
    
    
    for(int ivec=0;ivec<narch_vecs_full-1;ivec++)
        arch_vecs_full[ivec]->rearrange(idx_keep_old,idx_keep_new,nkeep,natms_new);
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
    for(int ivec=0;ivec<narch_vecs_full-1;ivec++)
        for(int i=0;i<nrcv;i++)
            arch_vecs_full[ivec]->pst_to(_rcv_buff,idx_lst_rcv[i]);
    
    delete [] idx_lst_rcv;
    delete [] *rcv_buff;
    delete [] rcv_buff;
    

    delete arch_vecs_full[narch_vecs_full-1];
    delete [] arch_vecs_full;
    arch_vecs_full=NULL;
    narch_vecs_full=0;
}
