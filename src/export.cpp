#include "export.h"
#include "atoms_styles.h"
#include <string>
#include "memory.h"
using namespace MAPP_NS;
/*--------------------------------------------
 
 --------------------------------------------*/
Export::Export(std::initializer_list<const char*> __def_vec_names,std::string* __user_vec_names,size_t __nuser_vecs):
atoms(NULL),
nvecs(0),
vecs(NULL),
vec_names(NULL)
{

    
    size_t tot_len=0;
    const char* const* ___def_vec_names=__def_vec_names.begin();
    for(size_t i=0;i<__def_vec_names.size();i++)
        tot_len+=std::strlen(___def_vec_names[i])+1;
    
    bool* valid_user_vecs=NULL;
    Memory::alloc(valid_user_vecs,__nuser_vecs);
    for(size_t i=0;i<__nuser_vecs;i++) valid_user_vecs[i]=true;
    
    
    
    nusr_vecs=0;
    
    
    
    for(size_t i=0;i<__nuser_vecs;i++)
    {
        
        for(int j=0;j<i && valid_user_vecs[i];j++)
            if(std::strcmp(__user_vec_names[i].c_str(),__user_vec_names[j].c_str())==0)
                valid_user_vecs[i]=false;
        
        for(size_t j=0;j<__def_vec_names.size() && valid_user_vecs[i];j++)
            if(std::strcmp(__user_vec_names[i].c_str(),___def_vec_names[j])==0)
                valid_user_vecs[i]=false;
        
        
        if(valid_user_vecs[i])
        {
            tot_len+=std::strlen(__user_vec_names[i].c_str())+1;
            nusr_vecs++;
        }
    }
    ndef_vecs=static_cast<int>(__def_vec_names.size());
    nvecs=nusr_vecs+ndef_vecs;
    vecs=new vec*[nvecs];
    vec_names=new char*[nvecs];
    char* buff=new char[tot_len];
    nvecs=0;
    
    for(size_t i=0;i<__def_vec_names.size();i++)
    {
        tot_len=std::strlen(___def_vec_names[i])+1;
        memcpy(buff,___def_vec_names[i],tot_len*sizeof(char));
        vec_names[nvecs++]=buff;
        buff+=tot_len;
    }
    
    for(size_t i=0;i<__nuser_vecs;i++)
    {
        if(!valid_user_vecs[i]) continue;
        
        tot_len=__user_vec_names[i].size()+1;
        memcpy(buff,__user_vec_names[i].c_str(),tot_len*sizeof(char));
        vec_names[nvecs++]=buff;
        buff+=tot_len;
    }
    
    Memory::dealloc(valid_user_vecs);
}
/*--------------------------------------------
 
 --------------------------------------------*/
Export::~Export()
{
    delete [] *vec_names;
    delete [] vec_names;
    delete [] vecs;
}
/*--------------------------------------------
 
 --------------------------------------------*/
void Export::add_to_default(const char* def_name)
{
    ptrdiff_t tot_len=(strchr(vec_names[nvecs-1],'\0')-vec_names[0])+1;
    size_t def_name_len=std::strlen(def_name)+1;
    char* buff=new char[tot_len+def_name_len];
    memcpy(buff,def_name,def_name_len*sizeof(char));
    memcpy(buff+def_name_len,*vec_names,tot_len*sizeof(char));
    char ** __vec_names=new char*[nvecs+1];
    *__vec_names=buff;
    buff+=def_name_len;
    for(int i=0;i<nvecs;i++)
    {
        __vec_names[i+1]=buff;
        buff+=std::strlen(buff)+1;
    }
    
    nvecs++;
    ndef_vecs++;
    delete [] *vec_names;
    delete [] vec_names;
    delete [] vecs;
    vec_names=__vec_names;
    vecs=new vec*[nvecs];
}
/*--------------------------------------------
 
 --------------------------------------------*/
void Export::find_vecs()
{
    auto find_vec=[](const char* name,vec** vs,int nvs)->vec*
    {
        for(int i=0;i<nvs;i++)
            if(vs[i]->name && std::strcmp(name,vs[i]->name)==0)
                return vs[i];
        return NULL;
    };
    
    vec** dynamic_vecs=atoms->dynamic_vecs;
    int ndynamic_vecs=atoms->ndynamic_vecs;
    vec** __vecs=atoms->vecs;
    int __nvecs=atoms->nvecs;
    vec* v;
    
    ndims=0;
    for(int i=0;i<nvecs;i++)
    {
        v=find_vec(vec_names[i],dynamic_vecs,ndynamic_vecs);
        if(!v)
        {
            vec* __v=find_vec(vec_names[i],__vecs,__nvecs);
            if(!__v)
                throw "cannot print vector "+std::string(vec_names[i])+", it does not exist";
            else
            {
                if(__v->is_empty())
                    throw "cannot print vector "+std::string(vec_names[i])+", it is empty";
                else
                    throw "cannot print vector "+std::string(vec_names[i])+", it is not included in this simulation (it is archived)";
            }
        }
        ndims+=v->ndim_dump();
        
        vecs[i]=v;
    }    
}
/*--------------------------------------------
 
 --------------------------------------------*/
void Export::gather(vec** vecs,int nvecs)
{
    
    for(int i=0;i<nvecs;i++)
        vecs[i]->init_dump();
    
    int comm_size=vecs[0]->atoms->comm_size;
    int comm_rank=vecs[0]->atoms->comm_rank;
    MPI_Comm& world=vecs[0]->atoms->world;
    int natms_lcl=vecs[0]->atoms->natms_lcl;
    int natms_rcvd=natms_lcl;
    
    for(int i=1;i<comm_size;i++)
    {
        if(comm_rank==i)
        {
            
            MPI_Send(&natms_lcl,1,MPI_INT,0,i,world);
            for(int j=0;j<nvecs;j++)
                MPI_Send(vecs[j]->begin_dump(),natms_lcl*vecs[j]->byte_sz,MPI_BYTE,0,i*(j+2),world);
        }
        if(comm_rank==0)
        {
            int rcv_natms_lcl;
            MPI_Recv(&rcv_natms_lcl,1,MPI_INT,i,i,world,MPI_STATUS_IGNORE);
            
            for(int j=0;j<nvecs;j++)
            {
                int size_dump=vecs[j]->byte_sz;
                byte* dump_data=reinterpret_cast<byte*>(vecs[j]->begin_dump())+natms_rcvd*size_dump;
                MPI_Recv(dump_data,rcv_natms_lcl*size_dump,MPI_BYTE,i,i*(j+2),world,MPI_STATUS_IGNORE);
            }
            
            natms_rcvd+=rcv_natms_lcl;
        }
    }
}
/*--------------------------------------------
 
 --------------------------------------------*/
void Export::release(vec** vecs,int nvecs)
{
    for(int i=0;i<nvecs;i++)
        vecs[i]->fin_dump();
}

