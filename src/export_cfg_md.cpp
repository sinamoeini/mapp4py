#include "export_cfg_md.h"
#include "atoms_md.h"
#include "elements.h"
#include "print.h"
using namespace MAPP_NS;
/*--------------------------------------------
 
 --------------------------------------------*/
ExportCFGMD::ExportCFGMD(int __nevery,const std::string& __pattern,
std::string* user_vec_names,size_t nuservecs,bool __sort):
Export({"elem","x"},user_vec_names,nuservecs),
pattern(__pattern+".%09d.cfg"),
sort(__sort)
{
    if(sort) add_to_default("id");
}
/*--------------------------------------------
 
 --------------------------------------------*/
ExportCFGMD::~ExportCFGMD()
{
}
/*--------------------------------------------
 
 --------------------------------------------*/
void ExportCFGMD::write_header(FILE* fp)
{
    if(atoms->comm_rank) return;
    
    fprintf(fp,"Number of particles = %d\n",atoms->natms);
    fprintf(fp,"A = %lf Angstrom (basic length-scale)\n",1.0);
    
    for(int i=0;i<__dim__;i++)
        for(int j=0;j<__dim__;j++)
            fprintf(fp,"H0(%d,%d) = %lf A\n",i+1,j+1,atoms->H[i][j]);
    
    if(x_d_inc)
    {
        fprintf(fp,"R = 1.0");
        fprintf(fp,"entry_count = %d\n",ndims-__dim__);
        
        vec** usr_vecs=vecs+ndef_vecs;
        int icmp=0;
        for(int i=0;i<nusr_vecs;i++)
        {
            if(usr_vecs[i]==atoms->x_d) continue;
            int d=usr_vecs[i]->ndim_dump();
            for(int j=0;j<d;j++)
                fprintf(fp,"auxiliary[%d] = %s_%d [reduced unit]\n",icmp++,usr_vecs[i]->name,j);
        }
        
    }
    else
    {
        fprintf(fp,".NO_VELOCITY.\n");
        fprintf(fp,"entry_count = %d\n",ndims);
        
        vec** usr_vecs=vecs+ndef_vecs;
        int icmp=0;
        for(int i=0;i<nusr_vecs;i++)
        {
            int d=usr_vecs[i]->ndim_dump();
            for(int j=0;j<d;j++)
                fprintf(fp,"auxiliary[%d] = %s_%d [reduced unit]\n",icmp++,usr_vecs[i]->name,j);
        }
        
    }
}
/*--------------------------------------------
 
 --------------------------------------------*/
void ExportCFGMD::write_body_sort(FILE* fp)
{
    gather(vecs,nvecs);
    
    if(atoms->comm_rank==0)
    {
        atoms->x2s_dump();
        if(x_d_inc) atoms->x_d2s_d_dump();
        int natms=atoms->natms;
        unsigned int* id=atoms->id->begin_dump();
        unsigned int* id_map=natms==0 ? NULL:new unsigned int[natms];
        for(int i=0;i<natms;i++) id_map[id[i]]=i;
        
        char** elem_names=atoms->elements.names;
        type0* masses=atoms->elements.masses;
        elem_type* elem=atoms->elem->begin_dump();
        int curr_elem=-1;
        int __curr_elem;
        unsigned int iatm;
        
        vec** usr_vecs=vecs+ndef_vecs;
        for(int i=0;i<natms;i++)
        {
            iatm=id_map[i];
            
            __curr_elem=elem[iatm];
            if(__curr_elem!=curr_elem)
            {
                curr_elem=__curr_elem;
                fprintf(fp,"%lf\n",masses[curr_elem]);
                fprintf(fp,"%s\n",elem_names[curr_elem]);
            }
            
            atoms->x->print(fp,iatm);
            
            for(int j=0;j<nusr_vecs;j++)
                usr_vecs[j]->print(fp,iatm);
            
            fprintf(fp,"\n");
        }
        delete [] id_map;
    }
    
    release(vecs,nvecs);
}
/*--------------------------------------------
 
 --------------------------------------------*/
void ExportCFGMD::write_body(FILE* fp)
{
    gather(vecs,nvecs);
    
    if(atoms->comm_rank==0)
    {
        atoms->x2s_dump();
        if(x_d_inc) atoms->x_d2s_d_dump();
        int natms=atoms->natms;
        
        char** elem_names=atoms->elements.names;
        type0* masses=atoms->elements.masses;
        elem_type* elem=atoms->elem->begin_dump();
        int curr_elem=-1;
        int __curr_elem;
        
        vec** usr_vecs=vecs+ndef_vecs;
        for(int i=0;i<natms;i++)
        {
            
            __curr_elem=elem[i];
            if(__curr_elem!=curr_elem)
            {
                curr_elem=__curr_elem;
                fprintf(fp,"%lf\n",masses[curr_elem]);
                fprintf(fp,"%s\n",elem_names[curr_elem]);
            }
            
            atoms->x->print(fp,i);
            
            for(int j=0;j<nusr_vecs;j++)
                usr_vecs[j]->print(fp,i);
            
            fprintf(fp,"\n");
        }
    }
    
    release(vecs,nvecs);
}
/*--------------------------------------------
 
 --------------------------------------------*/
void ExportCFGMD::init()
{
    try
    {
        find_vecs();
    }
    catch(std::string& err_msg)
    {
        throw err_msg;
    }
    
    x_d_inc=false;
    for(int i=0;i<nvecs && !x_d_inc;i++)
        if(std::strcmp(vec_names[i],"x_d")==0)
            x_d_inc=true;
}
/*--------------------------------------------
 
 --------------------------------------------*/
void ExportCFGMD::fin()
{
    
}
/*------------------------------------------------------------------------------------------------------------------------------------
 
 ------------------------------------------------------------------------------------------------------------------------------------*/
PyObject* ExportCFGMD::__new__(PyTypeObject* type,PyObject* args,PyObject* kwds)
{
    Object* __self=reinterpret_cast<Object*>(type->tp_alloc(type,0));
    PyObject* self=reinterpret_cast<PyObject*>(__self);
    return self;
}
/*--------------------------------------------
 
 --------------------------------------------*/
int ExportCFGMD::__init__(PyObject* self,PyObject* args,PyObject* kwds)
{
    //FuncAPI<int,std::string,std::string*,bool> f("__init__",{"nevery","file_pattern","vector_names","sort"});
    
    FuncAPI<std::string*>f("__init__",{"names"});
    //f.noptionals=2;
    //f.val<3>()=false;
    //f.logics<0>()[0]=VLogics("gt",0);
    
    if(f(args,kwds)==-1) return -1;
    //Object* __self=reinterpret_cast<Object*>(self);
    //__self->xprt=new ExportCFGMD(f.val<0>(),f.val<1>(),f.val<2>(),f.v<2>().size,f.val<3>());
    return 0;
}
/*--------------------------------------------
 
 --------------------------------------------*/
PyObject* ExportCFGMD::__alloc__(PyTypeObject* type,Py_ssize_t)
{
    Object* __self=new Object;
    __self->ob_type=type;
    __self->ob_refcnt=1;
    __self->xprt=NULL;
    return reinterpret_cast<PyObject*>(__self);
}
/*--------------------------------------------
 
 --------------------------------------------*/
void ExportCFGMD::__dealloc__(PyObject* self)
{
    Object* __self=reinterpret_cast<Object*>(self);
    delete __self->xprt;
    __self->xprt=NULL;
    delete __self;
}
/*--------------------------------------------*/
PyTypeObject ExportCFGMD::TypeObject={PyObject_HEAD_INIT(NULL)};
/*--------------------------------------------*/
void ExportCFGMD::setup_tp()
{
    TypeObject.tp_name="mapp.md.export_cfg";
    TypeObject.tp_doc="export atomeye";
    
    TypeObject.tp_flags=Py_TPFLAGS_DEFAULT;
    TypeObject.tp_basicsize=sizeof(Object);
    
    TypeObject.tp_new=__new__;
    TypeObject.tp_init=__init__;
    TypeObject.tp_alloc=__alloc__;
    TypeObject.tp_dealloc=__dealloc__;
    setup_tp_methods();
    TypeObject.tp_methods=methods;
    setup_tp_getset();
    TypeObject.tp_getset=getset;
}
/*--------------------------------------------*/
PyGetSetDef ExportCFGMD::getset[]={[0 ... 0]={NULL,NULL,NULL,NULL,NULL}};
/*--------------------------------------------*/
void ExportCFGMD::setup_tp_getset()
{
}
/*--------------------------------------------*/
PyMethodDef ExportCFGMD::methods[]={[0 ... 0]={NULL}};
/*--------------------------------------------*/
void ExportCFGMD::setup_tp_methods()
{
}

