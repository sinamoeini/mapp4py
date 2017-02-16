/*--------------------------------------------
 Created by Sina on 07/23/13.
 Copyright (c) 2013 MIT. All rights reserved.
 --------------------------------------------*/
#include <stdlib.h>
#include "read_cfg.h"
#include "elements.h"
#include "xmath.h"
#include "memory.h"
#include "atoms_md.h"
#include "atoms_dmd.h"
#include "comm.h"
#include "print.h"
#include "api.h"
using namespace MAPP_NS;
/*--------------------------------------------
 constructor
 --------------------------------------------*/
ReadCFG::ReadCFG(MPI_Comm& world):
Read(world)
{
}
/*--------------------------------------------
 constructor
 --------------------------------------------*/
ReadCFG::~ReadCFG()
{
}
/*--------------------------------------------
 reads the header of the cfg file
 --------------------------------------------*/
void ReadCFG::read_header(Atoms* atoms,FileReader& freader,char*& line,size_t& line_cpcty)
{
    type0 basic_length=1.0;
    type0 H0[__dim__][__dim__];
    type0 eta[__dim__][__dim__];
    type0 trns[__dim__][__dim__];
    for(int i=0;i<__dim__;i++)
        for(int j=0;j<__dim__;j++)
            H0[i][j]=eta[i][j]=trns[i][j]=0;
    
    for(int i=0;i<__dim__;i++)
        trns[i][i]=H0[i][i]=1.0;
    
    size_t nargs=0;

    R=1.0;
    entry_count=6;
    vel_xst=true;
    is_ext=false;
    
    bool empty=false;
    int icmp,jcmp,tmpno;
    bool header_cmplt=false;
    type0 tmp;
    
    
    while(!header_cmplt && !freader(line,line_cpcty))
    {
        if(freader.finished) continue;
        nargs=freader.hash_remover(line);
        
        if(nargs==0)
            continue;
        
        
        if(sscanf(line," Number of particles = %d ",&tmpno)==1)
        {
            if(tmpno<0)
                throw Print::vprintf("Number of particles in %s file should be greater than 0",freader.file);
            if(tmpno==0)
                empty=true;
            atoms->tot_natms=tmpno;
            continue;
        }
        
        if(sscanf(line," A = %lf ",&tmp)==1)
        {
            if(tmp<=0.0)
                throw Print::vprintf("A in %s file should be greater than 0.0",freader.file);
            basic_length=tmp;
            continue;
        }
        if(strcmp(line,".NO_VELOCITY.")==0)
        {
            if(!is_ext)
            {
                vel_xst=false;
                is_ext=true;
            }
            continue;
        }
        
        if(sscanf(line," entry_count = %d ",&tmpno)==1)
        {
            int mincomp=__dim__+(vel_xst?__dim__:0);
            if(tmpno<mincomp)
                throw Print::vprintf("entry_count in %s should at least be equal to %d",freader.file,mincomp);
            entry_count=tmpno;
            is_ext=true;
            continue;
        }
        
        if(sscanf(line," R = %lf %*s ",&tmp)==1)
        {
            R=tmp;
            continue;
        }
        
        if(sscanf(line," auxiliary [ %d ] = %*s [ %*s ]",&icmp)==1)
        {
            if(!is_ext)
                throw Print::vprintf("auxiliary [ ... ] = ... must come after entry_count = ... in %s",freader.file);
            
            int mincomp=__dim__+(vel_xst?__dim__:0);
            if(icmp+mincomp+1>entry_count)
                throw Print::vprintf("wrong component in %s file for auxiliary[%d], %d+%d+1 > entry_count",freader.file,icmp,mincomp,icmp);
            continue;
        }
        
        if(sscanf(line," Transform ( %d , %d ) = %lf ",&icmp,&jcmp,&tmp)==3)
        {
            if(icmp<1 || icmp>__dim__ || jcmp<1 || jcmp>__dim__)
                throw Print::vprintf("wrong component in %s file for Transform(%d,%d)",freader.file,icmp,jcmp);
            trns[icmp-1][jcmp-1]=tmp;
            continue;
        }
        if(sscanf(line," H0 ( %d , %d ) = %lf  %*s",&icmp,&jcmp,&tmp)==3)
        {
            if(icmp<1 || icmp>__dim__ || jcmp<1 || jcmp>__dim__)
                throw Print::vprintf("wrong component in %s file for H0(%d,%d)",freader.file,icmp,jcmp);
            H0[icmp-1][jcmp-1]=tmp;
            continue;
        }
        
        if(sscanf(line," eta ( %d , %d ) = %lf ",&icmp,&jcmp,&tmp)==3)
        {
            if(icmp<1 || icmp>__dim__ || jcmp<1 || jcmp>__dim__)
                throw Print::vprintf("wrong component in %s file for eta(%d,%d)",freader.file,icmp,jcmp);
            eta[icmp-1][jcmp-1]=tmp;
            continue;
        }
        
        if((nargs==8 && !is_ext) || (nargs==1 && is_ext))
        {
            header_cmplt=true;
            continue;
        }
        
        throw Print::vprintf("invalid line in %s file: %s",freader.file,line);
    }
    
    
    if(freader.finished && empty)
        header_cmplt=true;
    
    if(!header_cmplt)
        throw Print::vprintf("file %s ended unexpectedly",freader.file);
    
    
    
    type0 H_x[__dim__][__dim__];
    type0 H_x_d[__dim__][__dim__];
    for(int i=0;i<__dim__;i++)
        for(int j=0;j<__dim__;j++)
            H_x_d[i][j]=0.0;
    
    for(int i=0;i<__dim__;i++)
        for(int j=0;j<__dim__;j++)
        {
            H_x[i][j]=0.0;
            for(int k=0;k<__dim__;k++)
                H_x[i][j]+=H0[i][k]*trns[k][j];
        }
    
    bool chk=true;
    for(int i=0;i<__dim__ && chk;i++)
        for(int j=0;j<__dim__ && chk;j++)
            if(eta[i][j]!=0.0)
                chk=false;
    if(!chk)
    {
        for(int i=0;i<__dim__;i++)
            for(int j=0;j<__dim__;j++)
                eta[i][j]*=2.0;
        for(int i=0;i<__dim__;i++)
            eta[i][i]++;
        
        type0 eta_sq[__dim__][__dim__];
        
        if(XMath::Msqrt(eta,eta_sq)==0)
            throw Print::vprintf("eta in %s should be positive definite",freader.file);
        
        for(int i=0;i<__dim__;i++)
            for(int j=0;j<__dim__;j++)
                H0[i][j]=H_x[i][j];
        
        for(int i=0;i<__dim__;i++)
            for(int j=0;j<__dim__;j++)
            {
                H_x[i][j]=0.0;
                for(int k=0;k<__dim__;k++)
                    H_x[i][j]+=H0[i][k]*eta_sq[k][j];
            }
    }
    
    XMatrixVector::M_2_Mlt(H_x,H0);
    atoms->comm.grid(H0);
    for(int i=0;i<__dim__;i++)
        for(int j=0;j<__dim__;j++)
            atoms->H[i][j]=H0[i][j]*basic_length;

    atoms->update_H();
}
/*--------------------------------------------
 reads the atom section of the cfg file
 --------------------------------------------*/
void ReadCFG::read_body_ext(Atoms* atoms,FileReader& freader,char*& line,size_t& line_cpcty)
{
    type0 (&s_lo)[__dim__]=atoms->comm.s_lo;
    type0 (&s_hi)[__dim__]=atoms->comm.s_hi;
    auto belong=
    [&s_lo,s_hi](type0 (&s)[__dim__])->bool
    {
        for(int i=0;i<__dim__;i++)
            if(s[i]<s_lo[i] || s[i] >=s_hi[i])
                return false;
        return true;
    };
    
    byte* buff=new byte[sizeof(unsigned int)+sizeof(elem_type)+sizeof(type0)*entry_count];
    vec_list[0]=atoms->id;
    vec_list[1]=new Vec<elem_type>(atoms,1);
    vec_list[2]=atoms->x;
    vec_list[3]=NULL;
    nvecs=3;
    if(entry_count-__dim__)
        vec_list[nvecs++]=new Vec<type0>(atoms,entry_count-__dim__);
    
    char** args=NULL;
    size_t args_cpcty=0;
    size_t nargs;
    
    type0 s[__dim__];
    byte* __buff;
    elem_type ielem=0;
    type0 mass=0.0;
    unsigned int curr_id=0;
    
    bool mass_flag=false;
    bool elem_init=false;
    
    
    while(!freader.finished)
    {
        nargs=freader.parse_line(line,args,args_cpcty);
        if(nargs==0)
        {
            freader(line,line_cpcty);
            continue;
        }
        if(nargs!=1 && nargs!=entry_count)
        {
            delete [] buff;
            delete [] args;
            throw Print::vprintf("invalid line in %s file: %s",freader.file,line);
        }
        
        if(mass_flag && nargs!=1)
        {
            delete [] buff;
            delete [] args;
            throw Print::vprintf("expected chemical symbol of element after mass in %s",freader.file);
        }
        
        
        if(nargs==1)
        {
            if(!mass_flag)
            {
                mass=atof(args[0]);
                if(mass<=0.0)
                {
                    delete [] buff;
                    delete [] args;
                    throw Print::vprintf("mass of %s %s file (%lf) should be greater than 0.0",args[0],freader.file,line,mass);
                }
                mass_flag=true;
            }
            else
            {
                try
                {
                    ielem=atoms->elements->add_type(mass,args[0]);
                }
                catch(char* err_msg)
                {
                    delete [] buff;
                    delete [] args;
                    throw err_msg;
                }
                mass_flag=false;
                elem_init=true;
            }
            
            freader(line,line_cpcty);
            continue;
        }
        
        if(!elem_init)
        {
            delete [] buff;
            delete [] args;
            throw Print::vprintf("line %s in file %s comes before any element was defined",line,freader.file);
        }
        
        for(int i=0;i<__dim__;i++)
        {
            s[i]=atof(args[i]);
            
            while(s[i]>=1.0)
                s[i]--;
            while(s[i]<0.0)
                s[i]++;
        }
        
        if(!belong(s))
        {
            curr_id++;
            freader(line,line_cpcty);
            continue;
        }
        
        
        __buff=buff;
        memcpy(__buff,&curr_id,sizeof(unsigned int));
        __buff+=sizeof(unsigned int);
        memcpy(__buff,&ielem,sizeof(elem_type));
        __buff+=sizeof(elem_type);
        
        memcpy(__buff,s,sizeof(type0)*__dim__);
        __buff+=sizeof(type0)*__dim__;
        
        for(int i=__dim__; i<entry_count;i++)
        {
            type0 v=atof(args[i]);
            memcpy(__buff,&v,sizeof(type0));
            __buff+=sizeof(type0);
        }
        
        atoms->insert(buff,vec_list,nvecs,1);
        
        curr_id++;
        freader(line,line_cpcty);
    }
    
    delete [] buff;
    delete [] args;
    
    int tot_natms;
    MPI_Allreduce(&(atoms->natms),&tot_natms,1,MPI_INT,MPI_SUM,world);
    if(tot_natms!=atoms->tot_natms)
        throw Print::vprintf("expected %d natms in %s but read %d",atoms->tot_natms,freader.file,tot_natms);
}
/*--------------------------------------------
 reads the atom section of the cfg file
 --------------------------------------------*/
void ReadCFG::read_body_std(Atoms* atoms,FileReader& freader,char*& line,size_t& line_cpcty)
{
    type0 (&s_lo)[__dim__]=atoms->comm.s_lo;
    type0 (&s_hi)[__dim__]=atoms->comm.s_hi;
    auto belong=
    [&s_lo,s_hi](type0 (&s)[__dim__])->bool
    {
        for(int i=0;i<__dim__;i++)
            if(s[i]<s_lo[i] || s[i] >=s_hi[i])
                return false;
        return true;
    };
    
    byte* buff=new byte[sizeof(unsigned int)+sizeof(elem_type)+sizeof(type0)*entry_count];
    vec_list[0]=atoms->id;
    vec_list[1]=new Vec<elem_type>(atoms,1);
    vec_list[2]=atoms->x;
    vec_list[3]=NULL;
    nvecs=3;
    if(entry_count-__dim__)
        vec_list[nvecs++]=new Vec<type0>(atoms,entry_count-__dim__);
    
    char** args=NULL;
    size_t args_cpcty=0;
    size_t nargs;
    
    type0 s[__dim__];
    byte* __buff;
    elem_type ielem=0;
    type0 mass=0.0;
    unsigned int curr_id=0;
    

    while(!freader.finished)
    {
        nargs=freader.parse_line(line,args,args_cpcty);
        if(nargs==0)
        {
            freader(line,line_cpcty);
            continue;
        }
        
        if(nargs!=8)
        {
            delete [] buff;
            delete [] args;
            throw Print::vprintf("invalid line in %s file: %s",freader.file,line);
        }
        
        mass=atof(args[0]);
        if(mass<=0.0)
        {
            delete [] buff;
            delete [] args;
            throw Print::vprintf("mass of %s %s file (%lf) should be greater than 0.0",args[0],freader.file,line,mass);
        }
        
        try
        {
            ielem=atoms->elements->add_type(mass,args[1]);
        }
        catch(char* err_msg)
        {
            delete [] buff;
            delete [] args;
            throw err_msg;
        }
        
        
        for(int i=0;i<__dim__;i++)
        {
            s[i]=atof(args[i+2]);
            
            while(s[i]>=1.0)
                s[i]--;
            while(s[i]<0.0)
                s[i]++;
        }
        
        if(!belong(s))
        {
            curr_id++;
            freader(line,line_cpcty);
            continue;
        }
        
        __buff=buff;
        memcpy(__buff,&curr_id,sizeof(unsigned int));
        __buff+=sizeof(unsigned int);
        memcpy(__buff,&ielem,sizeof(elem_type));
        __buff+=sizeof(elem_type);
        
        memcpy(__buff,s,sizeof(type0)*__dim__);
        __buff+=sizeof(type0)*__dim__;
        
        for(int i=2+__dim__; i<8;i++)
        {
            type0 v=atof(args[i]);
            memcpy(__buff,&v,sizeof(type0));
            __buff+=sizeof(type0);
        }
        
        atoms->insert(buff,vec_list,nvecs,1);
        
        curr_id++;
        freader(line,line_cpcty);
    }
    
    delete [] buff;
    delete [] args;
    
    int tot_natms;
    MPI_Allreduce(&(atoms->natms),&tot_natms,1,MPI_INT,MPI_SUM,world);
    if(tot_natms!=atoms->tot_natms)
        throw Print::vprintf("expected %d natms in %s but read %d",atoms->tot_natms,freader.file,tot_natms);
}
/*--------------------------------------------
 constructor
 --------------------------------------------*/
ReadCFGMD::ReadCFGMD(MPI_Comm& world):
ReadCFG(world)
{
}
/*--------------------------------------------
 constructor
 --------------------------------------------*/
ReadCFGMD::~ReadCFGMD()
{
}
/*--------------------------------------------
 reads the atom section of the cfg file
 --------------------------------------------*/
AtomsMD* ReadCFGMD::operator()(const char* file)
{
    FileReader* pfreader;
    try
    {
        pfreader=new FileReader(world,file);
    }
    catch (char* err_msg)
    {
        throw err_msg;
    }

    AtomsMD* atoms=new AtomsMD(world);
    char* line=NULL;
    size_t line_cpcty=MAXCHAR;
    if(line_cpcty) line=new char[line_cpcty];
    
    try
    {
        read_header(atoms,*pfreader,line,line_cpcty);
    }
    catch(char* err_msg)
    {
        delete [] line;
        delete pfreader;
        delete atoms;
        throw err_msg;
    }
    
    if(pfreader->finished)
    {
        delete pfreader;
        delete [] line;
        return atoms;
    }
    
    try
    {
        if(is_ext)
            read_body_ext(atoms,*pfreader,line,line_cpcty);
        else
            read_body_std(atoms,*pfreader,line,line_cpcty);
        
    }
    catch(char* err_msg)
    {
        delete [] line;
        delete pfreader;
        delete atoms;
        throw err_msg;
    }
    
    delete pfreader;
    delete [] line;
    
    atoms->elem=dynamic_cast<Vec<elem_type>*>(vec_list[1]);
    vec_list[1]=NULL;
    if(vel_xst)
    {
        atoms->x_d=dynamic_cast<Vec<type0>*>(vec_list[3]);
        atoms->x_d->change_dim(__dim__);
        vec_list[3]=NULL;
        
        type0 M[__dim__][__dim__];
        for(int i=0;i<__dim__;i++)
            for(int j=0;j<__dim__;j++)
                M[i][j]=R*atoms->H[i][j];
            
        type0* x_d=atoms->x_d->begin();
        int natms=atoms->natms;
        for(int i=0;i<natms;i++,x_d+=__dim__)
            XMatrixVector::s2x(x_d,M);
        
    }
    else
        delete vec_list[3];

    atoms->s2x_lcl();
    
    return atoms;
}
/*--------------------------------------------
 
 --------------------------------------------*/
void ReadCFGMD::ml_cfg(PyMethodDef& tp_methods)
{
    tp_methods.ml_flags=METH_VARARGS | METH_KEYWORDS;
    tp_methods.ml_name="cfg";
    tp_methods.ml_doc="this function reads cfg file and returns an atoms_md object";
    
    tp_methods.ml_meth=(PyCFunction)(PyCFunctionWithKeywords)
    [](PyObject* self,PyObject* args,PyObject* kwds)->PyObject*
    {
        FuncAPI<std::string,OP<MAPP_MPI>> f("cfg",{"file","mpi"});
        f.noptionals=1;
        if(f(args,kwds)==-1) return NULL;
        
        MPI_Comm world=MPI_COMM_WORLD;
        if(f.val<1>().ob)
            world=reinterpret_cast<MAPP_MPI::Object*>(f.val<1>().ob)->world;
        
        ReadCFGMD read(world);
        AtomsMD* atoms=NULL;
        try
        {
            atoms=read(f.val<0>().c_str());
        }
        catch (char* err_msg)
        {
            PyErr_SetString(PyExc_TypeError,err_msg);
            delete [] err_msg;
            return NULL;
        }
        PyObject* op=AtomsMD::TypeObject.tp_alloc(&AtomsMD::TypeObject,0);
        reinterpret_cast<AtomsMD::Object*>(op)->atoms=atoms;
        return op;
    };
}
/*--------------------------------------------
 constructor
 --------------------------------------------*/
ReadCFGDMD::ReadCFGDMD(MPI_Comm& world):
ReadCFG(world)
{
}
/*--------------------------------------------
 constructor
 --------------------------------------------*/
ReadCFGDMD::~ReadCFGDMD()
{
}
/*--------------------------------------------
 reads the atom section of the cfg file
 --------------------------------------------*/
AtomsDMD* ReadCFGDMD::operator()(int N,const char* file)
{
    FileReader* pfreader;
    try
    {
        pfreader=new FileReader(world,file);
    }
    catch (char* err_msg)
    {
        throw err_msg;
    }

    AtomsDMD* atoms=new AtomsDMD(world,N);
    char* line=NULL;
    size_t line_cpcty=MAXCHAR;
    if(line_cpcty) line=new char[line_cpcty];
    
    try
    {
        read_header(atoms,*pfreader,line,line_cpcty);
    }
    catch(char* err_msg)
    {
        delete [] line;
        delete pfreader;
        delete atoms;
        throw err_msg;
    }
    
    if(pfreader->finished)
    {
        delete pfreader;
        delete [] line;
        return atoms;
    }
    
    if(!is_ext)
    {
        delete [] line;
        delete pfreader;
        delete atoms;
        throw Print::vprintf("cannot read cfg standard file for dmd mode");
    }
    
    try
    {
        read_body_ext(atoms,*pfreader,line,line_cpcty);
    }
    catch(char* err_msg)
    {
        delete [] line;
        delete pfreader;
        delete atoms;
        throw err_msg;
    }
    
    delete pfreader;
    delete [] line;

    
    delete vec_list[1];
    vec_list[1]=NULL;
    int nelems=static_cast<int>(atoms->elements->nelems);
    if(!vec_list[3] || vec_list[3]->dim<2*nelems)
    {
        delete atoms;
        throw Print::vprintf("values for c and alpha of all elements should be provided for every atom");
    }
    
    int c_dim=nelems;
    const int vec_dim=vec_list[3]->dim;
    
    
    atoms->alpha=new Vec<type0>(atoms,c_dim);
    atoms->c=new Vec<type0>(atoms,c_dim);
    atoms->elem=new Vec<elem_type>(atoms,c_dim);
    
    type0* alpha=dynamic_cast<Vec<type0>*>(vec_list[3])->begin();
    type0* c=alpha+nelems;
    type0* __alpha=atoms->alpha->begin();
    type0* __c=atoms->c->begin();
    elem_type* __elem=atoms->elem->begin();
    
    
    int lcl_err=0;
    int max_ncmp_lcl=0,__max_ncmp=0;
    int jcmp;
    elem_type __nelems=static_cast<elem_type>(nelems);
    type0 max_alpha_lcl=0.0;
    for(int i=0;i<atoms->natms && lcl_err==0;i++,alpha+=vec_dim,c+=vec_dim,__alpha+=c_dim,__c+=c_dim,__elem+=c_dim)
    {
        __max_ncmp=0;
        jcmp=0;
        for(elem_type j=0;j<__nelems && lcl_err==0;j++)
        {
            if(c[j]>1.0 || (c[j]<0.0 && c[j]!=-1.0))
            {
                lcl_err=1;
                continue;
            }
            
            if(c[j]!=-1.0 && alpha[j]<0.0)
            {
                lcl_err=2;
                continue;
            }
            
            if(c[j]==-1.0 && alpha[j]<0.0)
                alpha[j]=0.0;
            if(c[j]!=-1.0)
            {
                __elem[__max_ncmp]=j;
                __c[__max_ncmp]=c[j];
                __alpha[__max_ncmp]=alpha[j];
                max_alpha_lcl=MAX(alpha[j],max_alpha_lcl);
                __max_ncmp++;
            }
        }
        
        for(int j=__max_ncmp;j<c_dim;j++)
        {
            __elem[j]=0;
            __c[j]=-1.0;
            __alpha[j]=0.0;
        }
        
        max_ncmp_lcl=MAX(max_ncmp_lcl,__max_ncmp);
    }
    
    int lcl_err0=lcl_err==1?1:0;
    int lcl_err1=lcl_err==2?1:0;
    int err0,err1;
    MPI_Allreduce(&lcl_err0,&err0,1,MPI_INT,MPI_MAX,world);
    if(err0)
    {
        delete atoms;
        throw Print::vprintf("all values of c should be greater or equal to 0.0 and less than or equal to 1.0, or equal to -1.0");
    }
    
    MPI_Allreduce(&lcl_err1,&err1,1,MPI_INT,MPI_MAX,world);
    
    if(err1)
    {
        delete atoms;
        throw Print::vprintf("all values of alpha corresponding to present elements (i.e. c!=-1.0) should be greater than 0.0");
    }
    
    
    int max_ncmp;
    MPI_Allreduce(&max_ncmp_lcl,&max_ncmp,1,MPI_INT,MPI_MAX,world);
    MPI_Allreduce(&max_alpha_lcl,&atoms->max_alpha,1,MPI_DOUBLE,MPI_MAX,world);
    atoms->alpha->change_dim(max_ncmp);
    atoms->elem->change_dim(max_ncmp);
    atoms->c->change_dim(max_ncmp);
    atoms->c_dim=max_ncmp;
    delete vec_list[3];
    atoms->s2x_lcl();
    
    return atoms;
}
/*--------------------------------------------
 
 --------------------------------------------*/
void ReadCFGDMD::ml_cfg(PyMethodDef& tp_methods)
{
    tp_methods.ml_flags=METH_VARARGS | METH_KEYWORDS;
    tp_methods.ml_name="cfg";
    tp_methods.ml_doc="this function reads cfg file and returns an atoms_dmd object";
    
    tp_methods.ml_meth=(PyCFunction)(PyCFunctionWithKeywords)
    [](PyObject* self,PyObject* args,PyObject* kwds)->PyObject*
    {
        FuncAPI<int,std::string,OP<MAPP_MPI>> f("cfg",{"N","file","mpi"});
        f.noptionals=1;
        f.logics<0>()[0]=VLogics("gt",0);
        if(f(args,kwds)==-1) return NULL;
        
        MPI_Comm world=MPI_COMM_WORLD;
        if(f.val<2>().ob)
            world=reinterpret_cast<MAPP_MPI::Object*>(f.val<2>().ob)->world;
        
        ReadCFGDMD read(world);
        AtomsDMD* atoms=NULL;
        try
        {
            atoms=read(f.val<0>(),f.val<1>().c_str());
        }
        catch (char* err_msg)
        {
            PyErr_SetString(PyExc_TypeError,err_msg);
            delete [] err_msg;
            return NULL;
        }
        PyObject* op=AtomsDMD::TypeObject.tp_alloc(&AtomsDMD::TypeObject,0);
        reinterpret_cast<AtomsDMD::Object*>(op)->atoms=atoms;
        return op;
    };
}


























