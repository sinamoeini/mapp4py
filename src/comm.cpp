#include "comm.h"
#include "exchange.h"
#include "xmath.h"
using namespace MAPP_NS;
#include <structmember.h>
/*------------------------------------------------------------------------------------------------------------------------------------
 
 ------------------------------------------------------------------------------------------------------------------------------------*/
PyObject* MAPP_MPI::__new__(PyTypeObject* type,PyObject* args,PyObject* kwds)
{
    Object* __self=reinterpret_cast<Object*>(type->tp_alloc(type,0));
    PyObject* self=reinterpret_cast<PyObject*>(__self);
    return self;
}
/*--------------------------------------------
 
 --------------------------------------------*/
int MAPP_MPI::__init__(PyObject* self,PyObject* args,PyObject* kwds)
{
    FuncAPI<> func("__init__");
    if(func(args,kwds)==-1) return -1;
    Object* __self=reinterpret_cast<Object*>(self);
    __self->world=MPI_COMM_WORLD;
    MPI_Comm_rank(__self->world,&(__self->rank));
    MPI_Comm_size(__self->world,&(__self->size));
    return 0;
}
/*--------------------------------------------
 
 --------------------------------------------*/
PyObject* MAPP_MPI::__alloc__(PyTypeObject* type,Py_ssize_t)
{
    Object* __self=new Object;
    __self->ob_type=type;
    __self->ob_refcnt=1;
    __self->world=NULL;
    return reinterpret_cast<PyObject*>(__self);
}
/*--------------------------------------------
 
 --------------------------------------------*/
void MAPP_MPI::__dealloc__(PyObject* self)
{
    Object* __self=reinterpret_cast<Object*>(self);
    if(__self->world!=MPI_COMM_WORLD)
        MPI_Comm_free(&(__self->world));
}
/*--------------------------------------------*/
PyTypeObject MAPP_MPI::TypeObject ={PyObject_HEAD_INIT(NULL)};
/*--------------------------------------------*/
int MAPP_MPI::setup_tp()
{
    TypeObject.tp_name="mpi";
    TypeObject.tp_doc="just a simple container for MPI_Comm";
    
    TypeObject.tp_flags=Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE;
    TypeObject.tp_basicsize=sizeof(Object);
    
    TypeObject.tp_new=__new__;
    TypeObject.tp_init=__init__;
    TypeObject.tp_alloc=__alloc__;
    TypeObject.tp_dealloc=__dealloc__;
    setup_tp_members();
    TypeObject.tp_members=tp_members;
    
    int ichk=PyType_Ready(&TypeObject);
    if(ichk<0) return ichk;
    Py_INCREF(&TypeObject);
    return ichk;
}
/*--------------------------------------------*/
PyMemberDef MAPP_MPI::tp_members[]={[0 ... 2]={NULL,0,0,0,NULL}};
/*--------------------------------------------*/
void MAPP_MPI::setup_tp_members()
{
    tp_members[0].type=T_INT;
    tp_members[0].flags=READONLY;
    tp_members[0].name=(char*)"rank";
    tp_members[0].doc=(char*)"rank of the mpi communicator";
    tp_members[0].offset=offsetof(MAPP_MPI::Object,rank);
    
    tp_members[1].type=T_INT;
    tp_members[1].flags=READONLY;
    tp_members[1].name=(char*)"size";
    tp_members[1].doc=(char*)"size of the mpi communicator";
    tp_members[1].offset=offsetof(MAPP_MPI::Object,size);
}
/*------------------------------------------------------------------------------------------------------------------------------------
 
 ------------------------------------------------------------------------------------------------------------------------------------*/
Communication::Communication(MPI_Comm __world,int(&N)[__dim__],type0(&H)[__dim__][__dim__],type0 __skin):
world(__world),
size(get_size(__world)),
rank(get_rank(__world)),
skin(__skin),
xchng_id(0)
{
    grid(N,H);
}
/*--------------------------------------------
 
 --------------------------------------------*/
Communication::Communication(MPI_Comm __world,type0(&H)[__dim__][__dim__],type0 __skin):
world(__world),
size(get_size(__world)),
rank(get_rank(__world)),
skin(__skin),
xchng_id(0)
{
    int N[__dim__]={[0 ... __dim__-1]=0};
    grid(N,H);
}
/*--------------------------------------------
 
 --------------------------------------------*/
Communication::Communication(MPI_Comm __world,type0 __skin):
world(__world),
size(get_size(__world)),
rank(get_rank(__world)),
skin(__skin),
xchng_id(0)
{
    int N[__dim__]={[0 ... __dim__-1]=0};
    type0 H[__dim__][__dim__]={[0 ... __dim__-1]={[0 ... __dim__-1]=0.0}};
    for(int i=0;i<__dim__;i++) H[i][i]=1.0;
    grid(N,H);
}
/*--------------------------------------------
 
 --------------------------------------------*/
Communication::Communication(const Communication& r):
world(r.world),
size(r.size),
rank(r.rank),
skin(r.skin),
xchng_id(r.xchng_id)
{
    for(int i=0;i<__dim__;i++)
    {
        coords[i]=r.coords[i];
        dims[i]=r.dims[i];
        neigh[i][0]=r.neigh[i][0];
        neigh[i][1]=r.neigh[i][1];
        s_lo[i]=r.s_lo[i];
        s_hi[i]=r.s_hi[i];
    }
}
/*--------------------------------------------
 
 --------------------------------------------*/
Communication& Communication::operator=(const Communication& r)
{
    this->~Communication();
    new (this) Communication(r);
    return *this;
}
/*--------------------------------------------
 
 --------------------------------------------*/
void Communication::grid(type0(&H)[__dim__][__dim__])
{
    int N[__dim__]={[0 ... __dim__-1]=0};
    grid(N,H);
}
/*--------------------------------------------
 
 --------------------------------------------*/
void Communication::grid(int(&N)[__dim__])
{
    type0 H[__dim__][__dim__]={[0 ... __dim__-1]={[0 ... __dim__-1]=0.0}};
    for(int i=0;i<__dim__;i++) H[i][i]=1.0;
    grid(N,H);
}
/*--------------------------------------------
 
 --------------------------------------------*/
void Communication::grid(int(&N)[__dim__],type0(&H)[__dim__][__dim__])
{
    int lcl_rank,nrank,nsize;
    int* p_per_n;
    int** np_arr;
    bool uni_dist;

    network_analysis(world,lcl_rank,nrank,nsize,p_per_n,np_arr,uni_dist);
    
    if(__dim__==1)
    {
        dims[0]=size;
        coords[0]=lcl_rank;
        for(int i=0;i<nrank;i++) coords[0]+=p_per_n[i];
        neigh[0][0]=lcl_rank==0?
        np_arr[nrank==0?nsize-1:nrank-1][p_per_n[nrank==0?nsize-1:nrank-1]-1]:
        np_arr[nrank][lcl_rank-1];
        
        neigh[0][1]=lcl_rank==p_per_n[nrank]-1?
        np_arr[nrank==nsize-1?0:nrank+1][0]:np_arr[nrank][lcl_rank+1];
        
        s_lo[0]=static_cast<type0>(coords[0])/static_cast<type0>(dims[0]);
        s_hi[0]=static_cast<type0>(coords[0]+1)/static_cast<type0>(dims[0]);
        
        
        delete [] *np_arr;
        delete [] np_arr;
        delete [] p_per_n;
        return;
    }
    
    
    int ndims_det=0;
    int det_size=1;
    for(int i=0;i<__dim__;i++)
        if(N[i])
        {
            det_size*=N[i];
            ndims_det++;
        }
    
    /* if only one dimension is undetermind it can be determined */
    if(ndims_det==__dim__-1)
    {
        for(int i=0;i<__dim__;i++)
            if(!N[i]) N[i]=size/det_size;
        
        det_size=size;
        ndims_det++;
    }
    
    int prin_dim=-1;
    if(uni_dist && nsize!=1)
    {
        /* sort based on least communication required */
        int key[__dim__];
        for(int i=0;i<__dim__;i++) key[i]=i;
        XMath::quicksort(key,key+__dim__
        ,[&H](int* _i,int* _j){return (H[*_i][*_i]>H[*_j][*_j]);}
        ,[](int* _i,int* _j){std::swap(*_i,*_j);}
        );
        
        int idim=0;
        while(prin_dim==-1 && idim<__dim__)
        {
            if(N[key[idim]] && N[key[idim]]==nsize)
                prin_dim=key[idim];
            if(!N[key[idim]] && (size/det_size)%nsize==0)
            {
                N[key[idim]]=nsize;
                prin_dim=key[idim];
                det_size*=nsize;
                ndims_det++;
            }
            idim++;
        }
    }
    
    if(ndims_det!=__dim__)
    {
        int __dims[__dim__];
        type0 h[__dim__];
        int icur=0;
        for(int i=0;i<__dim__;i++)
            if(!N[i]) h[icur++]=H[i][i];
        Algebra::opt_comm_grid(__dim__-ndims_det,h,size/det_size,__dims);
        
        icur=0;
        for(int i=0;i<__dim__;i++)
            if(!N[i]) N[i]=__dims[icur++];
    }
    
    if(prin_dim!=-1)
    {
        dims[prin_dim]=nsize;
        coords[prin_dim]=nrank;
        neigh[prin_dim][0]=nrank==0? np_arr[nsize-1][lcl_rank]:np_arr[nrank-1][lcl_rank];
        neigh[prin_dim][1]=nrank==nsize-1? np_arr[0][lcl_rank]:np_arr[nrank+1][lcl_rank];
        
        int __dims[__dim__-1];
        int __coords[__dim__-1];
        int __neigh[__dim__-1][2];
        int icur=0;
        for(int i=0;i<__dim__;i++)
            if(i!=prin_dim)
                __dims[icur++]=N[i];
        
        seq_analysis(rank,np_arr[nrank],p_per_n[nrank],__dims,__coords,__neigh);
        
        icur=0;
        for(int i=0;i<__dim__;i++)
            if(i!=prin_dim)
            {
                coords[i]=__coords[icur];
                dims[i]=__dims[icur];
                neigh[i][0]=__neigh[icur][0];
                neigh[i][1]=__neigh[icur][1];
                icur++;
            }
    }
    else
    {
        for(int i=0;i<__dim__;i++) dims[i]=N[i];
        int list[__dim__];
        for(int i=0;i<__dim__;i++) list[i]=1;
        MPI_Comm cartesian;
        MPI_Cart_create(world,__dim__,dims,list,1,&cartesian);
        MPI_Cart_get(cartesian,__dim__,dims,list,coords);
        for(int i=0;i<__dim__;i++)
            MPI_Cart_shift(cartesian,i,1,&neigh[i][0],&neigh[i][1]);
        MPI_Comm_free(&cartesian);
    }
    
    for(int i=0;i<__dim__;i++)
    {
        s_lo[i]=static_cast<type0>(coords[i])/static_cast<type0>(dims[i]);
        s_hi[i]=static_cast<type0>(coords[i]+1)/static_cast<type0>(dims[i]);
    }
    
    
    delete [] *np_arr;
    delete [] np_arr;
    delete [] p_per_n;
}
/*--------------------------------------------
 
 --------------------------------------------*/
int Communication::get_rank(MPI_Comm& comm)
{
    int i;
    MPI_Comm_rank(comm,&i);
    return i;
}
/*--------------------------------------------
 
 --------------------------------------------*/
int Communication::get_size(MPI_Comm& comm)
{
    int i;
    MPI_Comm_size(comm,&i);
    return i;
}
/*--------------------------------------------
 
 --------------------------------------------*/
int Communication::get_rank()
{
    int i;
    MPI_Comm_rank(MPI_COMM_WORLD,&i);
    return i;
}
/*--------------------------------------------
 
 --------------------------------------------*/
int Communication::get_size()
{
    int i;
    MPI_Comm_size(MPI_COMM_WORLD,&i);
    return i;
}
/*--------------------------------------------
 
 --------------------------------------------*/
void Communication::network_analysis(MPI_Comm& world,
int& lcl_rank,int& nrank,int& nsize,int*& p_per_n,
int**& np_arr,bool& uni_dist)
{
    int size=get_size(world);
    int rank=get_rank(world);
    
    
    int* ids=new int[size];
    char* name=new char[MPI_MAX_PROCESSOR_NAME];
    int name_sz;
    MPI_Get_processor_name(name,&name_sz);
    name_sz++;
    int cum_name_sz;
    MPI_Allreduce(&name_sz,&cum_name_sz,1,MPI_INT,MPI_SUM,world);
    char** names=new char*[size];
    *names=NULL;
    if(cum_name_sz) *names=new char[cum_name_sz];
    
    char* __name=*names;
    int __name_sz;
    for(int ip=0;ip<size;ip++)
    {
        ids[ip]=ip;
        if(ip==rank)
        {
            __name_sz=name_sz;
            memcpy(__name,name,__name_sz*sizeof(char));
        }
        
        MPI_Bcast(&__name_sz,1,MPI_INT,ip,world);
        MPI_Bcast(__name,__name_sz,MPI_CHAR,ip,world);
        names[ip]=__name;
        __name+=__name_sz;
    }
    delete [] name;
    
    
    nsize=0;
    p_per_n=NULL;
    int* p_per_n_;
    for(int i=0,last=0;i<size;i++)
    {
        if(i==rank)
            nrank=nsize;
        for(int j=i+1;j<size;j++)
            if(strcmp(names[i],names[j])==0)
            {
                if(j==rank)
                    nrank=nsize;
                std::swap(names[i+1],names[j]);
                std::swap(ids[i+1],ids[j]);
                i++;
            }
        
        p_per_n_=new int[nsize+1];
        memcpy(p_per_n_,p_per_n,nsize*sizeof(int));
        delete [] p_per_n;
        p_per_n=p_per_n_;
        p_per_n[nsize++]=i+1-last;
        last=i+1;
    }
    
    delete [] *names;
    delete [] names;
    
    np_arr=new int*[nsize];
    *np_arr=ids;
    for(int i=1;i<nsize;i++)
        np_arr[i]=np_arr[i-1]+p_per_n[i-1];
    uni_dist=true;
    for(int i=0;i<nsize && uni_dist;i++)
        if(p_per_n[i]!=p_per_n[0])
            uni_dist=false;
    
    lcl_rank=0;
    while(np_arr[nrank][lcl_rank]!=rank) lcl_rank++;
}






