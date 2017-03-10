#include "export_cfg.h"
#include "atoms_dmd.h"
#include "elements.h"
#include "print.h"
using namespace MAPP_NS;
/*--------------------------------------------
 
 --------------------------------------------*/
ExportCFGDMD::ExportCFGDMD(const std::string& __pattern,
std::string* user_vec_names,size_t nuservecs,bool __sort):
Export({"x","alpha","c"},user_vec_names,nuservecs),
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
        
        char** elem_names=atoms->elements.names;
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
                fprintf(fp,"%s\n",elem_names[curr_elem]);
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
                fprintf(fp,"%s\n",elem_names[i]);
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
        char** elem_names=atoms->elements.names;
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
                fprintf(fp,"%s\n",elem_names[curr_elem]);
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
                fprintf(fp,"%s\n",elem_names[i]);
            }
        
        delete [] elem_printed;
    }

    release(vecs,nvecs);
}
/*--------------------------------------------
 
 --------------------------------------------*/
void ExportCFGDMD::write(int stps)
{
    /*
     we have a list of vectors
        defaults 
        user defined ones 
     
     
     some vectors will just send their regular 
     others will need preparing 
     
     
     */
    
    
    char* file_name=Print::vprintf(pattern.c_str(),stps);
    
    FILE* fp=NULL;
    if(atoms->comm_rank==0) fp=fopen(file_name,"w");
    delete [] file_name;
    
    write_header(fp);
    if(sort) write_body_sort(fp);
    else write_body(fp);
    
    
    if(atoms->comm_rank==0)
        fclose(fp);
}


