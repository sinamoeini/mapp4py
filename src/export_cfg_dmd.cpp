#include "export_cfg_dmd.h"
#include "atoms_dmd.h"
#include "elements.h"
#include "print.h"
/*--------------------------------------------
 
 --------------------------------------------*/
ExportCFGDMD::ExportCFGDMD(const std::string& __pattern,int __nevery,
std::string* user_vec_names,size_t nuservecs,bool __sort):
ExportDMD({"x","alpha","c"},__nevery,user_vec_names,nuservecs),
pattern(__pattern+".%09d.cfg"),
sort(__sort)
{
    if(sort) add_to_default("id");
}
/*--------------------------------------------
 
 --------------------------------------------*/
ExportCFGDMD::~ExportCFGDMD()
{
}
/*--------------------------------------------
 
 --------------------------------------------*/
void ExportCFGDMD::write_header(FILE* fp)
{
    if(atoms->comm_rank) return;
    
    fprintf(fp,"Number of particles = %d\n",atoms->natms);
    fprintf(fp,"A = %lf Angstrom (basic length-scale)\n",1.0);
    
    for(int i=0;i<__dim__;i++)
        for(int j=0;j<__dim__;j++)
            fprintf(fp,"H0(%d,%d) = %lf A\n",i+1,j+1,atoms->H[i][j]);
    
    fprintf(fp,".NO_VELOCITY.\n");
    int nelems=static_cast<int>(atoms->elements.nelems);
    
    if(sort)
        fprintf(fp,"entry_count = %d\n",ndims-1);
    else
        fprintf(fp,"entry_count = %d\n",ndims);
    
    int icmp=0;
    for(int i=0;i<nelems;i++)
        fprintf(fp,"auxiliary[%d] = %s_%d [reduced unit]\n",icmp++,"alpha",i);
    for(int i=0;i<nelems;i++)
        fprintf(fp,"auxiliary[%d] = %s_%d [reduced unit]\n",icmp++,"c",i);
    
    vec** usr_vecs=vecs+ndef_vecs;
        
    for(int i=0;i<nusr_vecs;i++)
    {
        int d=usr_vecs[i]->ndim_dump();
        for(int j=0;j<d;j++)
            fprintf(fp,"auxiliary[%d] = %s_%d [reduced unit]\n",icmp++,usr_vecs[i]->name,j);
    }
}
/*--------------------------------------------
 
 --------------------------------------------*/
void ExportCFGDMD::write_body_sort(FILE* fp)
{
    gather(vecs,nvecs);
    
    if(atoms->comm_rank==0)
    {
        atoms->x2s_dump();
        int natms=atoms->natms;
        unsigned int* id=atoms->id->begin_dump();
        unsigned int* id_map=natms==0 ? NULL:new unsigned int[natms];
        for(int i=0;i<natms;i++) id_map[id[i]]=i;
        
        std::string* elem_names=atoms->elements.names;
        type0* masses=atoms->elements.masses;
        int nelems=static_cast<int>(atoms->elements.nelems);
        bool* elem_printed=nelems==0 ? NULL:new bool[nelems];
        for(int i=0;i<nelems;i++) elem_printed[i]=false;
        int curr_elem=-1;
        int __curr_elem;
        type0* c=atoms->c->begin_dump();
        type0* __c=c;
        type0 max_c;
        unsigned int iatm;
        
        vec** usr_vecs=vecs+ndef_vecs;
        for(int i=0;i<natms;i++)
        {
            iatm=id_map[i];
            
            __c=c+iatm*nelems;
            max_c=*__c;
            __curr_elem=0;
            for(int j=1;j<nelems;j++)
                if(__c[j]>max_c)
                {
                    max_c=__c[j];
                    __curr_elem=j;
                }
            
            if(__curr_elem!=curr_elem)
            {
                curr_elem=__curr_elem;
                elem_printed[curr_elem]=true;
                fprintf(fp,"%lf\n",masses[curr_elem]);
                fprintf(fp,"%s\n",elem_names[curr_elem].c_str());
            }
            
            atoms->x->print(fp,iatm);
            atoms->alpha->print(fp,iatm);
            atoms->c->print(fp,iatm);
            
            for(int j=0;j<nusr_vecs;j++)
                usr_vecs[j]->print(fp,iatm);
            
            fprintf(fp,"\n");
        }
        
        for(int i=0;i<nelems;i++)
            if(!elem_printed[i])
            {
                fprintf(fp,"%lf\n",masses[i]);
                fprintf(fp,"%s\n",elem_names[i].c_str());
            }
        
        delete [] elem_printed;
        delete [] id_map;
    }
    
    
    release(vecs,nvecs);
}
/*--------------------------------------------
 
 --------------------------------------------*/
void ExportCFGDMD::write_body(FILE* fp)
{
    gather(vecs,nvecs);
    
    if(atoms->comm_rank==0)
    {
        atoms->x2s_dump();
        int natms=atoms->natms;
        std::string* elem_names=atoms->elements.names;
        type0* masses=atoms->elements.masses;
        int nelems=static_cast<int>(atoms->elements.nelems);
        bool* elem_printed=nelems==0 ? NULL:new bool[nelems];
        for(int i=0;i<nelems;i++) elem_printed[i]=false;
        int curr_elem=-1;
        int __curr_elem;
        type0* c=atoms->c->begin_dump();
        type0* __c;
        type0 max_c;
        vec** usr_vecs=vecs+ndef_vecs;
        for(int i=0;i<natms;i++)
        {
            
            
            
            __c=c+i*nelems;
            max_c=*__c;
            __curr_elem=0;
            for(int j=1;j<nelems;j++)
                if(__c[j]>max_c)
                {
                    max_c=__c[j];
                    __curr_elem=j;
                }
            
            
            if(__curr_elem!=curr_elem)
            {
                curr_elem=__curr_elem;
                elem_printed[curr_elem]=true;
                fprintf(fp,"%lf\n",masses[curr_elem]);
                fprintf(fp,"%s\n",elem_names[curr_elem].c_str());
            }
            
            atoms->x->print(fp,i);
            atoms->alpha->print(fp,i);
            atoms->c->print(fp,i);
            for(int j=0;j<nusr_vecs;j++)
                usr_vecs[j]->print(fp,i);
            
            fprintf(fp,"\n");
        }
        
        for(int i=0;i<nelems;i++)
            if(!elem_printed[i])
            {
                fprintf(fp,"%lf\n",masses[i]);
                fprintf(fp,"%s\n",elem_names[i].c_str());
            }
        
        delete [] elem_printed;
    }

    release(vecs,nvecs);
}
/*--------------------------------------------
 
 --------------------------------------------*/
void ExportCFGDMD::write(int stps)
{
    FILE* fp;
    char* file_name=Print::vprintf(pattern.c_str(),stps);
    bool file_chk=open(atoms->comm_rank,atoms->world,file_name,"w",fp);
    delete [] file_name;
    if(!file_chk) return;
    
    
    write_header(fp);
    if(sort) write_body_sort(fp);
    else write_body(fp);
    
    close(atoms->comm_rank,atoms->world,fp);
}
/*--------------------------------------------
 
 --------------------------------------------*/
void ExportCFGDMD::init()
{
    try
    {
        find_vecs(atoms);
    }
    catch(std::string& err_msg)
    {
        throw err_msg;
    }
}
/*--------------------------------------------
 
 --------------------------------------------*/
void ExportCFGDMD::fin()
{
    
}
/*------------------------------------------------------------------------------------------------------------------------------------
 
 ------------------------------------------------------------------------------------------------------------------------------------*/
PyObject* ExportCFGDMD::__new__(PyTypeObject* type,PyObject* args,PyObject* kwds)
{
    Object* __self=reinterpret_cast<Object*>(type->tp_alloc(type,0));
    PyObject* self=reinterpret_cast<PyObject*>(__self);
    return self;
}
/*--------------------------------------------
 
 --------------------------------------------*/
int ExportCFGDMD::__init__(PyObject* self,PyObject* args,PyObject* kwds)
{
    FuncAPI<std::string,int,std::string*,bool> f("__init__",{"file_pattern","nevery","extra_vecs","sort"});
    f.noptionals=3;
    f.logics<1>()[0]=VLogics("gt",0);
    
    f.val<3>()=false;
    f.val<1>()=10000;
    
    if(f(args,kwds)==-1) return -1;
    Object* __self=reinterpret_cast<Object*>(self);
    __self->xprt=new ExportCFGDMD(f.val<0>(),f.val<1>(),f.val<2>(),f.v<2>().size,f.val<3>());
    return 0;
}
/*--------------------------------------------
 
 --------------------------------------------*/
PyObject* ExportCFGDMD::__alloc__(PyTypeObject* type,Py_ssize_t)
{
    Object* __self=new Object;
    __self->ob_type=type;
    __self->ob_refcnt=1;
    __self->xprt=NULL;
    return reinterpret_cast<PyObject*>(__self);
}
/*--------------------------------------------
 
 --------------------------------------------*/
void ExportCFGDMD::__dealloc__(PyObject* self)
{
    Object* __self=reinterpret_cast<Object*>(self);
    delete __self->xprt;
    __self->xprt=NULL;
    delete __self;
}
/*--------------------------------------------*/
PyTypeObject ExportCFGDMD::TypeObject={PyObject_HEAD_INIT(NULL)};
/*--------------------------------------------*/
int ExportCFGDMD::setup_tp()
{
    TypeObject.tp_name="mapp.dmd.export_cfg";
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
    TypeObject.tp_base=&ExportDMD::TypeObject;
    
    int ichk=PyType_Ready(&TypeObject);
    if(ichk<0) return ichk;
    Py_INCREF(&TypeObject);
    return ichk;
}
/*--------------------------------------------*/
PyGetSetDef ExportCFGDMD::getset[]={[0 ... 0]={NULL,NULL,NULL,NULL,NULL}};
/*--------------------------------------------*/
void ExportCFGDMD::setup_tp_getset()
{
}
/*--------------------------------------------*/
PyMethodDef ExportCFGDMD::methods[]={[0 ... 0]={NULL}};
/*--------------------------------------------*/
void ExportCFGDMD::setup_tp_methods()
{
}


