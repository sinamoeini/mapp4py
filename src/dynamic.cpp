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
Dynamic::Dynamic(Atoms* __atoms,ForceField* __ff,bool __chng_box,
std::initializer_list<vec*> __updt_vecs_def,std::initializer_list<vec*> __updt_vecs,
std::initializer_list<vec*> __xchng_comp_vecs_def,std::initializer_list<vec*> __xchng_comp_vecs,
std::initializer_list<vec*> __arch_vecs_def,std::initializer_list<vec*> __arch_vecs):
atoms(__atoms),
world(__atoms->comm.world),
skin(__atoms->comm.skin),
ff(__ff),
chng_box(__chng_box),

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
        arch_vecs_full[narch_vecs_full++]=new Vec<id_type>(atoms,1,"id_arch");
    Memory::shrink_to_fit(arch_vecs_full,narch_vecs_full,narch_vecs+1);
    if(!narch_vecs_full) return;
    
    id_type* id_0=atoms->id->begin();
    id_type* id_1=reinterpret_cast<id_type*>(arch_vecs_full[narch_vecs_full-1]->begin());
    memcpy(id_1,id_0,atoms->natms_lcl*sizeof(id_type));
}
/*--------------------------------------------
 
 --------------------------------------------*/
void Dynamic::restore_arch_vecs()
{
        if(!narch_vecs_full) return;
    
    Vec<id_type>* __id_old=dynamic_cast<Vec<id_type>*>(arch_vecs_full[narch_vecs_full-1]);
    id_type* id_old=__id_old->begin();
    int natms_old=__id_old->vec_sz;
    int natms_new=atoms->id->vec_sz;
    id_type* id_new=atoms->id->begin();
    
    rearrange_vecs(world,id_old,natms_old,id_new,natms_new,arch_vecs_full,narch_vecs_full-1);
    delete arch_vecs_full[narch_vecs_full-1];
    delete [] arch_vecs_full;
    arch_vecs_full=NULL;
    narch_vecs_full=0;
}
/*--------------------------------------------
 this is an amazing function
 ids_ndd: sorted ids that belong to
 me but someone else got them
 
 ids_xtr: sorted ids that I have but
 belong to someone else
 
 ndd_idxs: the indicies of the atoms
 correspeonding to ids_ndd. they will
 be resorsted according to the
 procesors they belong to
 --------------------------------------------*/
int* Dynamic::srt_idx_by_p(MPI_Comm& __world,
const id_type* xtr_ids,int xtr_ids_sz,
const id_type* ndd_ids,int ndd_ids_sz,
int*& ndd_idxs)
{
    const int tot_p=Communication::get_size(__world);
    const int my_p=Communication::get_rank(__world);
    int* comm_size=new int[tot_p];
    
    for(int ip=0;ip<tot_p;ip++)
        comm_size[ip]=0;
    
    int max_lst_sz;
    MPI_Allreduce(&xtr_ids_sz,&max_lst_sz,1,MPI_INT,MPI_MAX,__world);
    int mother_lst_sz;
    id_type* mother_lst=max_lst_sz==0 ? NULL:new id_type[max_lst_sz];
    
    
    int fnd_sz,ufnd_sz;
    id_type* ufnd=ndd_ids_sz==0 ? NULL:new id_type[ndd_ids_sz];
    int* ufnd_idx=ndd_ids_sz==0 ? NULL:new int[ndd_ids_sz];
    int* fnd_idx=ndd_ids_sz==0 ? NULL:new int[ndd_ids_sz];
    
    fnd_sz=0;
    ufnd_sz=ndd_ids_sz;
    memcpy(ufnd,ndd_ids,ndd_ids_sz*sizeof(id_type));
    memcpy(ufnd_idx,ndd_idxs,ndd_ids_sz*sizeof(int));
    for(int imother_lst,iufnd,__ufnd_sz,ip=0;ip<tot_p;ip++)
    {
        if(ip==my_p)
        {
            mother_lst_sz=xtr_ids_sz;
            memcpy(mother_lst,xtr_ids,mother_lst_sz*sizeof(id_type));
        }
        MPI_Bcast(&mother_lst_sz,1,MPI_INT,ip,__world);
        MPI_Bcast(mother_lst,mother_lst_sz*sizeof(id_type),MPI_BYTE,ip,__world);
        
        if(ip==my_p || ufnd_sz==0 || mother_lst_sz==0) continue;
        
        __ufnd_sz=0;
        for(imother_lst=0,iufnd=0;imother_lst<mother_lst_sz && iufnd<ufnd_sz;)
        {
            if(mother_lst[imother_lst]<ufnd[iufnd])
                imother_lst++;
            else if(mother_lst[imother_lst]>ufnd[iufnd])
            {
                ufnd_idx[__ufnd_sz]=ufnd_idx[iufnd];
                ufnd[__ufnd_sz++]=ufnd[iufnd++];
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
            ufnd_idx[__ufnd_sz]=ufnd_idx[iufnd];
            ufnd[__ufnd_sz++]=ufnd[iufnd++];
        }
        
        ufnd_sz=__ufnd_sz;
        
    }
    
    delete [] ndd_idxs;
    ndd_idxs=fnd_idx;
    delete [] ufnd;
    delete [] ufnd_idx;
    delete [] mother_lst;
    
    return comm_size;
}
/*--------------------------------------------
 
 --------------------------------------------*/
void Dynamic::rearrange_vecs(MPI_Comm& __world,
id_type* id_old,int natms_old,
id_type* id_new,int natms_new,
vec** vecs,int nvecs)
{

    //rank old ids
    int* key_old=NULL;
    if(natms_old) key_old=new int[natms_old];
    for(int i=0;i<natms_old;i++) key_old[i]=i;
    XMath::quicksort(key_old,key_old+natms_old
    ,[&id_old](int* rank_i,int* rank_j){return (id_old[*rank_i]<id_old[*rank_j]);}
    ,[](int* rank_i,int* rank_j){std::swap(*rank_i,*rank_j);}
    );
    
    //rank new ids
    int* key_new=NULL;
    if(natms_new) key_new=new int[natms_new];
    for(int i=0;i<natms_new;i++) key_new[i]=i;
    XMath::quicksort(key_new,key_new+natms_new
    ,[&id_new](int* rank_i,int* rank_j){return (id_new[*rank_i]<id_new[*rank_j]);}
    ,[](int* rank_i,int* rank_j){std::swap(*rank_i,*rank_j);}
    );
    
    
    /*
     allocation for the atoms to be sent
     worst case scenario all the old ones
     have to be sent
     */
    int nsnd=0;
    int* idx_lst_snd=natms_old==0 ? NULL:new int[natms_old];
    id_type* id_lst_snd=natms_old==0 ? NULL:new id_type[natms_old];


    /*
     allocation for the atoms to be received
     worst case scenario all the new ones
     have to be received
     */
    int nrcv=0;
    int* idx_lst_rcv=natms_new==0 ? NULL:new int[natms_new];
    id_type* id_lst_rcv=natms_new==0 ? NULL:new id_type[natms_new];

    
    
    /*
     allocation for the atoms that do not
     need to be sent or received.
     i.e. they intersection of old and new
     maximum possible number of this list
     is MIN(natms_old,natms_new)
     
     we keep the record of them in two
     fashion:
     idx_keep_old: idx in old
     idx_keep_new: idx in new
     */
    int nkeep=0;
    int* idx_keep_old=MIN(natms_old,natms_new)==0 ? NULL:new int[MIN(natms_old,natms_new)];
    int* idx_keep_new=MIN(natms_old,natms_new)==0 ? NULL:new int[MIN(natms_old,natms_new)];

    /*
     this is were the members of all 3
     categories are determined by
     comparing the sorted ids
     */
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
    
    /*
     done with determining the lists
     so we delete the sorted keys
     */
    
    
    /*
     before the worst case scenario was
     considered here the memory is shrunk
     to fit.
     it is not neccesaary but nice to do
     */
    Memory::shrink_to_fit(id_lst_snd,nsnd,natms_old);
    Memory::shrink_to_fit(idx_lst_snd,nsnd,natms_old);
    Memory::shrink_to_fit(id_lst_rcv,nrcv,natms_new);
    Memory::shrink_to_fit(idx_lst_rcv,nrcv,natms_new);
    Memory::shrink_to_fit(idx_keep_old,nkeep,MIN(natms_old,natms_new));
    Memory::shrink_to_fit(idx_keep_new,nkeep,MIN(natms_old,natms_new));


    int* nsnd_comm=srt_idx_by_p(__world,id_lst_rcv,nrcv,id_lst_snd,nsnd,idx_lst_snd);
    int* nrcv_comm=srt_idx_by_p(__world,id_lst_snd,nsnd,id_lst_rcv,nrcv,idx_lst_rcv);
    delete [] id_lst_snd;
    delete [] id_lst_rcv;
    
    
    
    
    
 
    const int tot_p=Communication::get_size(__world);
    const int my_p=Communication::get_rank(__world);
    
    int byte_sz=0;
    for(int ivec=0;ivec<nvecs;ivec++)
        byte_sz+=vecs[ivec]->byte_sz;
    
    
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
    
    byte* __snd_buff=*snd_buff;
    for(int ivec=0;ivec<nvecs;ivec++)
        for(int i=0;i<nsnd;i++)
            vecs[ivec]->cpy(__snd_buff,idx_lst_snd[i]);
    delete [] idx_lst_snd;
    
    
    for(int ivec=0;ivec<nvecs;ivec++)
        vecs[ivec]->rearrange(idx_keep_old,idx_keep_new,nkeep,natms_new);
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
                     __world,MPI_STATUS_IGNORE);

    }
    delete [] nsnd_comm;
    delete [] nrcv_comm;
    
    delete [] *snd_buff;
    delete [] snd_buff;
    
    byte* __rcv_buff=*rcv_buff;
    for(int ivec=0;ivec<nvecs;ivec++)
        for(int i=0;i<nrcv;i++)
            vecs[ivec]->pst_to(__rcv_buff,idx_lst_rcv[i]);
    
    delete [] idx_lst_rcv;
    delete [] *rcv_buff;
    delete [] rcv_buff;
}
/*--------------------------------------------
 
 --------------------------------------------*/
NewDynamic::NewDynamic(Atoms* __atoms,ForceField* __ff,
std::initializer_list<vec*> __updt_vecs_def,std::initializer_list<vec*> __updt_vecs,
std::initializer_list<vec*> __xchng_comp_vecs_def,std::initializer_list<vec*> __xchng_comp_vecs,
std::initializer_list<vec*> __arch_vecs_def,std::initializer_list<vec*> __arch_vecs):
atoms(__atoms),
world(__atoms->comm.world),
skin(__atoms->comm.skin),
ff(__ff),


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
NewDynamic::~NewDynamic()
{
    delete [] arch_vecs;
    delete [] xchng_comp_vecs;
    delete [] updt_vecs;
}
/*--------------------------------------------
 
 --------------------------------------------*/
void NewDynamic::add_xchng(vec* v)
{
    vec** __xchng_comp_vecs=new vec*[nxchng_comp_vecs+1];
    memcpy(__xchng_comp_vecs,xchng_comp_vecs,nxchng_comp_vecs*sizeof(vec*));
    delete [] xchng_comp_vecs;
    xchng_comp_vecs=__xchng_comp_vecs;
    xchng_comp_vecs[nxchng_comp_vecs++]=v;
}
/*--------------------------------------------
 
 --------------------------------------------*/
void NewDynamic::add_updt(vec* v)
{
    vec** __updt_vecs=new vec*[nupdt_vecs+1];
    memcpy(__updt_vecs,updt_vecs,nupdt_vecs*sizeof(vec*));
    delete [] updt_vecs;
    updt_vecs=__updt_vecs;
    updt_vecs[nupdt_vecs++]=v;
}
/*--------------------------------------------
 
 --------------------------------------------*/
void NewDynamic::create_dynamic_vecs()
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
void NewDynamic::destroy_dynamic_vecs()
{
    delete [] atoms->dynamic_vecs;
    atoms->dynamic_vecs=NULL;
    atoms->ndynamic_vecs=0;
    
    nxchng_vecs_full=nupdt_vecs_full=0;
}
/*--------------------------------------------
 
 --------------------------------------------*/
void NewDynamic::store_arch_vecs()
{
    
    Memory::alloc(arch_vecs_full,narch_vecs+1);
    narch_vecs_full=0;
    for(int i=0;i<narch_vecs;i++)
        if(!arch_vecs[i]->is_empty())
            arch_vecs_full[narch_vecs_full++]=arch_vecs[i];
    if(narch_vecs_full)
        arch_vecs_full[narch_vecs_full++]=new Vec<id_type>(atoms,1,"id_arch");
    Memory::shrink_to_fit(arch_vecs_full,narch_vecs_full,narch_vecs+1);
    if(!narch_vecs_full) return;
    
    id_type* id_0=atoms->id->begin();
    id_type* id_1=reinterpret_cast<id_type *>(arch_vecs_full[narch_vecs_full-1]->begin());
    memcpy(id_1,id_0,atoms->natms_lcl*sizeof(id_type));
}
/*--------------------------------------------
 
 --------------------------------------------*/
void NewDynamic::restore_arch_vecs()
{
    if(!narch_vecs_full) return;
    
    Vec<id_type>* __id_old=dynamic_cast<Vec<id_type>*>(arch_vecs_full[narch_vecs_full-1]);
    id_type* id_old=__id_old->begin();
    int natms_old=__id_old->vec_sz;
    int natms_new=atoms->id->vec_sz;
    id_type* id_new=atoms->id->begin();
    
    Dynamic::rearrange_vecs(world,id_old,natms_old,id_new,natms_new,arch_vecs_full,narch_vecs_full-1);
    delete arch_vecs_full[narch_vecs_full-1];
    delete [] arch_vecs_full;
    arch_vecs_full=NULL;
    narch_vecs_full=0;
}
