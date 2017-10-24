#include "export.h"
#include "atoms_styles.h"
#include <string>
#include "memory.h"
using namespace MAPP_NS;
/*--------------------------------------------
 
 --------------------------------------------*/
ExportMD::ExportMD(std::initializer_list<const char*> __def_vec_names,
int __nevery,std::string* __user_vec_names,size_t __nuser_vecs):
Export(__def_vec_names,__nevery,__user_vec_names,__nuser_vecs)
{
}
/*--------------------------------------------
 
 --------------------------------------------*/
ExportMD::~ExportMD()
{
}
/*--------------------------------------------
 
 --------------------------------------------*/
ExportDMD::ExportDMD(std::initializer_list<const char*> __def_vec_names,
int __nevery,std::string* __user_vec_names,size_t __nuser_vecs):
Export(__def_vec_names,__nevery,__user_vec_names,__nuser_vecs)
{
}
/*--------------------------------------------
 
 --------------------------------------------*/
ExportDMD::~ExportDMD()
{
}
/*--------------------------------------------
 
 --------------------------------------------*/
Export::Export(std::initializer_list<const char*> __def_vec_names,
int __nevery,std::string* __user_vec_names,size_t __nuser_vecs):
nevery(__nevery),
nvecs(0),
vecs(NULL),
vec_names(NULL)
{

    const char* const* ___def_vec_names=__def_vec_names.begin();
    
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
            nusr_vecs++;
    }
    
    ndef_vecs=static_cast<int>(__def_vec_names.size());
    nvecs=nusr_vecs+ndef_vecs;
    vecs=new vec*[nvecs];
    vec_names=new std::string[nvecs];
    nvecs=0;
    for(size_t i=0;i<__def_vec_names.size();i++)
        vec_names[nvecs++]=___def_vec_names[i];
    
    for(size_t i=0;i<__nuser_vecs;i++)
    {
        if(!valid_user_vecs[i]) continue;
        vec_names[nvecs++]=__user_vec_names[i];
    }
    Memory::dealloc(valid_user_vecs);
}
/*--------------------------------------------
 
 --------------------------------------------*/
Export::~Export()
{
    delete [] vec_names;
    delete [] vecs;
}
/*--------------------------------------------
 
 --------------------------------------------*/
void Export::add_to_default(const char* def_name)
{
    
    std::string* __vec_names=new std::string[nvecs+1];
    __vec_names[0]=def_name;
    for(int i=0;i<nvecs;i++)
        __vec_names[i+1]=std::move(vec_names[i]);
    nvecs++;
    ndef_vecs++;
    delete [] vec_names;
    delete [] vecs;
    vec_names=__vec_names;
    vecs=new vec*[nvecs];
}
/*--------------------------------------------
 
 --------------------------------------------*/
void Export::find_vecs(Atoms* atoms)
{
    auto find_vec=[](const std::string& name,vec** vs,int nvs)->vec*
    {
        for(int i=0;i<nvs;i++)
            if(vs[i]->name && std::strcmp(name.c_str(),vs[i]->name)==0)
                return vs[i];
        return NULL;
    };
    
    vec** dynamic_vecs=atoms->dynamic_vecs;
    int ndynamic_vecs=atoms->ndynamic_vecs;
    vec** __vecs=atoms->vecs;
    int __nvecs=atoms->nvecs;
    vec* v;
    
    ndims=0;
    if(ndynamic_vecs)
    {
        for(int i=0;i<nvecs;i++)
        {
            v=find_vec(vec_names[i],dynamic_vecs,ndynamic_vecs);
            if(!v)
            {
                vec* __v=find_vec(vec_names[i],__vecs,__nvecs);
                if(!__v)
                throw "cannot export vector '"+std::string(vec_names[i])+"' to file, it does not exist";
                else
                {
                    if(__v->is_empty())
                    throw "cannot export vector '"+std::string(vec_names[i])+"' to file, it is empty";
                    else
                    throw "cannot export vector '"+std::string(vec_names[i])+"' to file, it is not included in this simulation (it is archived)";
                }
            }
            ndims+=v->ndim_dump();
            
            vecs[i]=v;
        }
    }
    else
    {
        for(int i=0;i<nvecs;i++)
        {
            v=find_vec(vec_names[i],__vecs,__nvecs);
            if(!v)
            {
                throw "cannot export vector '"+std::string(vec_names[i])+"' to file, it does not exist";
                
            }
            else
            {
                if(v->is_empty())
                    throw "cannot export vector '"+std::string(vec_names[i])+"' to file, it is empty";
            }
            ndims+=v->ndim_dump();
            
            vecs[i]=v;
        }
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
/*--------------------------------------------
 
 --------------------------------------------*/
bool Export::open(const int rank,MPI_Comm& world, const char* path,const char* mode,FILE*& fp)
{
    int err=0;
    if(rank==0)
    {
        fp=fopen(path,mode);
        if(!fp)
            err=1;
            
    }
    else fp=NULL;
    
    MPI_Bcast(&err,1,MPI_INT,0,world);
    if(err) return false;
    return true;
}
/*--------------------------------------------
 
 --------------------------------------------*/
void Export::close(const int rank,MPI_Comm& world,FILE*& fp)
{
    if(rank==0) fclose(fp);
}
/*------------------------------------------------------------------------------------------------------------------------------------
 
 ------------------------------------------------------------------------------------------------------------------------------------*/
void Export::getset_deafult_vecs(PyGetSetDef& getset)
{
    getset.name=(char*)"default_vecs";
    getset.doc=(char*)"default vectors included in this export object";
    getset.get=[](PyObject* self,void*)->PyObject*
    {
        size_t sz=reinterpret_cast<Object*>(self)->xprt->ndef_vecs;
        std::string* vec_names=reinterpret_cast<Object*>(self)->xprt->vec_names;
        size_t* sz_ptr=&sz;
        return var<std::string*>::build(vec_names,&sz_ptr);
    };
    getset.set=[](PyObject* self,PyObject* op,void*)->int
    {
        PyErr_SetString(PyExc_TypeError,"readonly attribute");
        return -1;
    };
}
/*--------------------------------------------
 
 --------------------------------------------*/
void Export::getset_extra_vecs(PyGetSetDef& getset)
{
    getset.name=(char*)"extra_vecs";
    getset.doc=(char*)"extra vectors included in this export object";
    getset.get=[](PyObject* self,void*)->PyObject*
    {
        size_t sz=reinterpret_cast<Object*>(self)->xprt->nusr_vecs;
        std::string* vec_names=reinterpret_cast<Object*>(self)->xprt->vec_names+
        reinterpret_cast<Object*>(self)->xprt->ndef_vecs;
        size_t* sz_ptr=&sz;
        return var<std::string*>::build(vec_names,&sz_ptr);
    };
    getset.set=[](PyObject* self,PyObject* op,void*)->int
    {
        VarAPI<std::string*> extra_vecs("extra_vecs");
        int ichk=extra_vecs.set(op);
        if(ichk==-1) return -1;
        int nusr_vecs=static_cast<int> (extra_vecs.__var__.size);
        int ndef_vecs=reinterpret_cast<Object*>(self)->xprt->ndef_vecs;
        int nvecs=nusr_vecs+ndef_vecs;
        
        std::string* vec_names=reinterpret_cast<Object*>(self)->xprt->vec_names;
        std::string* __vec_names=new std::string[nvecs];
        for(int i=0;i<ndef_vecs;i++)
            __vec_names[i]=std::move(vec_names[i]);
        
        for(int i=0;i<nusr_vecs;i++)
            __vec_names[i+ndef_vecs]=std::move(extra_vecs.val[i]);
        
        
        delete [] reinterpret_cast<Object*>(self)->xprt->vecs;
        delete [] vec_names;
        reinterpret_cast<Object*>(self)->xprt->vecs=new vec*[nvecs];
        reinterpret_cast<Object*>(self)->xprt->vec_names=__vec_names;
        reinterpret_cast<Object*>(self)->xprt->nvecs=nvecs;
        reinterpret_cast<Object*>(self)->xprt->nusr_vecs=nusr_vecs;
        return 0;
    };
}
/*--------------------------------------------
 
 --------------------------------------------*/
void Export::getset_nevery(PyGetSetDef& getset)
{
    getset.name=(char*)"nevery";
    getset.doc=(char*)"export nevery step";
    getset.get=[](PyObject* self,void*)->PyObject*
    {
        
        return var<int>::build(reinterpret_cast<Object*>(self)->xprt->nevery);
    };
    getset.set=[](PyObject* self,PyObject* op,void*)->int
    {
        VarAPI<int> nevery("nevery");
        nevery.logics[0]=VLogics("gt",0);
        int ichk=nevery.set(op);
        if(ichk==-1) return -1;
        reinterpret_cast<Object*>(self)->xprt->nevery=nevery.val;
        return 0;
    };
}
/*------------------------------------------------------------------------------------------------------------------------------------
 
 ------------------------------------------------------------------------------------------------------------------------------------*/
PyObject* ExportMD::__new__(PyTypeObject* type,PyObject* args,PyObject* kwds)
{
    Object* __self=reinterpret_cast<Object*>(type->tp_alloc(type,0));
    PyObject* self=reinterpret_cast<PyObject*>(__self);
    return self;
}
/*--------------------------------------------
 
 --------------------------------------------*/
int ExportMD::__init__(PyObject* self,PyObject* args,PyObject* kwds)
{
    FuncAPI<int,std::string*> f("__init__",{"nevery","extra_vecs"});
    f.noptionals=2;
    f.val<0>()=10000;
    f.logics<0>()[0]=VLogics("gt",0);
    
    
    
    if(f(args,kwds)==-1) return -1;
    Object* __self=reinterpret_cast<Object*>(self);
    __self->xprt=new ExportMD({"elem","x"},f.val<0>(),f.val<1>(),f.v<1>().size);
    return 0;
}
/*--------------------------------------------
 
 --------------------------------------------*/
PyObject* ExportMD::__alloc__(PyTypeObject* type,Py_ssize_t)
{
    Object* __self=new Object;
    __self->ob_type=type;
    __self->ob_refcnt=1;
    __self->xprt=NULL;
    return reinterpret_cast<PyObject*>(__self);
}
/*--------------------------------------------
 
 --------------------------------------------*/
void ExportMD::__dealloc__(PyObject* self)
{
    Object* __self=reinterpret_cast<Object*>(self);
    delete __self->xprt;
    __self->xprt=NULL;
    delete __self;
}
/*--------------------------------------------
 
 --------------------------------------------*/
PyObject* ExportMD::__call__(PyObject* self,PyObject* args,PyObject* kwds)
{
    FuncAPI<OP<AtomsMD>,std::string> f("__call__",{"atoms","file_name"});
    
    if(f(args,kwds)==-1) return NULL;
    
    reinterpret_cast<ExportMD::Object*>(self)->xprt->atoms=reinterpret_cast<AtomsMD::Object*>(f.val<0>().ob)->atoms;
    try
    {
        reinterpret_cast<ExportMD::Object*>(self)->xprt->write(f.val<1>().c_str());
    }
    catch(std::string& err_msg)
    {
        reinterpret_cast<ExportMD::Object*>(self)->xprt->atoms=NULL;
        PyErr_SetString(PyExc_TypeError,err_msg.c_str());
        return NULL;
    }
    
    reinterpret_cast<ExportMD::Object*>(self)->xprt->atoms=NULL;
    Py_RETURN_NONE;
}
/*--------------------------------------------*/
PyTypeObject ExportMD::TypeObject={PyObject_HEAD_INIT(NULL)};
/*--------------------------------------------*/
int ExportMD::setup_tp()
{
    TypeObject.tp_name="mapp.md.export";
    TypeObject.tp_doc="export";
    
    TypeObject.tp_flags=Py_TPFLAGS_DEFAULT;
    TypeObject.tp_basicsize=sizeof(Object);
    
    TypeObject.tp_new=__new__;
    TypeObject.tp_init=__init__;
    TypeObject.tp_alloc=__alloc__;
    TypeObject.tp_dealloc=__dealloc__;
    TypeObject.tp_call=__call__;
    setup_tp_methods();
    TypeObject.tp_methods=methods;
    setup_tp_getset();
    TypeObject.tp_getset=getset;
    
    int ichk=PyType_Ready(&TypeObject);
    if(ichk<0) return ichk;
    Py_INCREF(&TypeObject);
    return ichk;
}
/*--------------------------------------------*/
PyGetSetDef ExportMD::getset[]=EmptyPyGetSetDef(1);
/*--------------------------------------------*/
void ExportMD::setup_tp_getset()
{
}
/*--------------------------------------------*/
PyMethodDef ExportMD::methods[]=EmptyPyMethodDef(1);
/*--------------------------------------------*/
void ExportMD::setup_tp_methods()
{
}
/*------------------------------------------------------------------------------------------------------------------------------------
 
 ------------------------------------------------------------------------------------------------------------------------------------*/
PyObject* ExportDMD::__new__(PyTypeObject* type,PyObject* args,PyObject* kwds)
{
    Object* __self=reinterpret_cast<Object*>(type->tp_alloc(type,0));
    PyObject* self=reinterpret_cast<PyObject*>(__self);
    return self;
}
/*--------------------------------------------
 
 --------------------------------------------*/
int ExportDMD::__init__(PyObject* self,PyObject* args,PyObject* kwds)
{
    FuncAPI<int,std::string*> f("__init__",{"nevery","extra_vecs"});
    f.noptionals=2;
    f.val<0>()=10000;
    f.logics<0>()[0]=VLogics("gt",0);
    
    
    
    if(f(args,kwds)==-1) return -1;
    Object* __self=reinterpret_cast<Object*>(self);
    __self->xprt=new ExportDMD({"x","alpha","c"},f.val<0>(),f.val<1>(),f.v<1>().size);
    return 0;
}
/*--------------------------------------------
 
 --------------------------------------------*/
PyObject* ExportDMD::__alloc__(PyTypeObject* type,Py_ssize_t)
{
    Object* __self=new Object;
    __self->ob_type=type;
    __self->ob_refcnt=1;
    __self->xprt=NULL;
    return reinterpret_cast<PyObject*>(__self);
}
/*--------------------------------------------
 
 --------------------------------------------*/
void ExportDMD::__dealloc__(PyObject* self)
{
    Object* __self=reinterpret_cast<Object*>(self);
    delete __self->xprt;
    __self->xprt=NULL;
    delete __self;
}
/*--------------------------------------------
 
 --------------------------------------------*/
PyObject* ExportDMD::__call__(PyObject* self,PyObject* args,PyObject* kwds)
{
    FuncAPI<OP<AtomsDMD>,std::string> f("__call__",{"atoms","file_name"});
    
    if(f(args,kwds)==-1) return NULL;
    
    reinterpret_cast<ExportDMD::Object*>(self)->xprt->atoms=reinterpret_cast<AtomsDMD::Object*>(f.val<0>().ob)->atoms;
    try
    {
        reinterpret_cast<ExportDMD::Object*>(self)->xprt->write(f.val<1>().c_str());
    }
    catch(std::string& err_msg)
    {
        reinterpret_cast<ExportDMD::Object*>(self)->xprt->atoms=NULL;
        PyErr_SetString(PyExc_TypeError,err_msg.c_str());
        return NULL;
    }
    
    reinterpret_cast<ExportDMD::Object*>(self)->xprt->atoms=NULL;
    Py_RETURN_NONE;
}
/*--------------------------------------------*/
PyTypeObject ExportDMD::TypeObject={PyObject_HEAD_INIT(NULL)};
/*--------------------------------------------*/
int ExportDMD::setup_tp()
{
    TypeObject.tp_name="mapp.dmd.export";
    TypeObject.tp_doc="export";
    
    TypeObject.tp_flags=Py_TPFLAGS_DEFAULT;
    TypeObject.tp_basicsize=sizeof(Object);
    
    TypeObject.tp_new=__new__;
    TypeObject.tp_init=__init__;
    TypeObject.tp_alloc=__alloc__;
    TypeObject.tp_dealloc=__dealloc__;
    TypeObject.tp_call=__call__;
    
    setup_tp_methods();
    TypeObject.tp_methods=methods;
    setup_tp_getset();
    TypeObject.tp_getset=getset;
    
    int ichk=PyType_Ready(&TypeObject);
    if(ichk<0) return ichk;
    Py_INCREF(&TypeObject);
    return ichk;
}
/*--------------------------------------------*/
PyGetSetDef ExportDMD::getset[]=EmptyPyGetSetDef(1);
/*--------------------------------------------*/
void ExportDMD::setup_tp_getset()
{
}
/*--------------------------------------------*/
PyMethodDef ExportDMD::methods[]=EmptyPyMethodDef(1);
/*--------------------------------------------*/
void ExportDMD::setup_tp_methods()
{
}

