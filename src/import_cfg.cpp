/*--------------------------------------------
 Created by Sina on 07/23/13.
 Copyright (c) 2013 MIT. All rights reserved.
 --------------------------------------------*/
#include "elements.h"
#include <stdlib.h>
#include "import_cfg.h"
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
ImportCFG::ImportCFG(MPI_Comm& world):
Import(world)
{
}
/*--------------------------------------------
 constructor
 --------------------------------------------*/
ImportCFG::~ImportCFG()
{
}
/*--------------------------------------------
 reads the header of the cfg file
 --------------------------------------------*/
void ImportCFG::read_header(Atoms* atoms,FileReader& freader,char*& line,size_t& line_cpcty)
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
            atoms->natms=tmpno;
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
    
    Algebra::MSQ_2_MLT(H_x,H0);
    atoms->comm.grid(H0);
    for(int i=0;i<__dim__;i++)
        for(int j=0;j<__dim__;j++)
            atoms->H[i][j]=H0[i][j]*basic_length;

    atoms->update_H();
}
/*--------------------------------------------
 reads the atom section of the cfg file
 --------------------------------------------*/
void ImportCFG::read_body_ext(Atoms* atoms,FileReader& freader,char*& line,size_t& line_cpcty)
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
                    ielem=atoms->elements.add_type(mass,args[0]);
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
    
    int natms;
    MPI_Allreduce(&(atoms->natms_lcl),&natms,1,MPI_INT,MPI_SUM,world);
    if(natms!=atoms->natms)
        throw Print::vprintf("expected %d natms_lcl in %s but read %d",atoms->natms,freader.file,natms);
}
/*--------------------------------------------
 reads the atom section of the cfg file
 --------------------------------------------*/
void ImportCFG::read_body_std(Atoms* atoms,FileReader& freader,char*& line,size_t& line_cpcty)
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
            ielem=atoms->elements.add_type(mass,args[1]);
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
    
    int natms;
    MPI_Allreduce(&(atoms->natms_lcl),&natms,1,MPI_INT,MPI_SUM,world);
    if(natms!=atoms->natms)
        throw Print::vprintf("expected %d natms_lcl in %s but read %d",atoms->natms,freader.file,natms);
}
/*--------------------------------------------
 constructor
 --------------------------------------------*/
ImportCFGMD::ImportCFGMD(MPI_Comm& world):
ImportCFG(world)
{
}
/*--------------------------------------------
 constructor
 --------------------------------------------*/
ImportCFGMD::~ImportCFGMD()
{
}
/*--------------------------------------------
 reads the atom section of the cfg file
 --------------------------------------------*/
AtomsMD* ImportCFGMD::operator()(const char* file)
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

    Atoms* atoms=new Atoms(world);
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
        AtomsMD* atoms_md=new AtomsMD(world);
        *atoms_md=*atoms;
        delete atoms;
        return atoms_md;
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
    
    
    AtomsMD* atoms_md=new AtomsMD(world);
    *atoms_md=*atoms;
    memcpy(atoms_md->elem->begin(),vec_list[1]->begin(),atoms->natms_lcl*sizeof(elem_type));
    delete vec_list[1];
    vec_list[1]=NULL;
    if(vel_xst)
    {
        atoms_md->x_d->fill();
        memcpy(atoms_md->x_d->begin(),vec_list[3]->begin(),atoms->natms_lcl*__dim__*sizeof(type0));
        type0* __x_d=atoms_md->x_d->begin();
        type0* ___x_d=reinterpret_cast<type0*>(vec_list[3]->begin());
        int vec_dim=vec_list[3]->dim;
        for(int i=0;i<atoms->natms_lcl;i++)
        {
            memcpy(__x_d,___x_d,__dim__*sizeof(type0));
            __x_d+=__dim__;
            ___x_d+=vec_dim;
        }
        
        
        
        type0 M[__dim__][__dim__];
        for(int i=0;i<__dim__;i++)
            for(int j=0;j<__dim__;j++)
                M[i][j]=R*atoms->H[i][j];
            
        type0* x_d=atoms_md->x_d->begin();
        int natms_lcl=atoms->natms_lcl;
        for(int i=0;i<natms_lcl;i++,x_d+=__dim__)
            Algebra::V_mul_MLT(x_d,M,x_d);
    }

    delete vec_list[3];
    vec_list[3]=NULL;

    delete atoms;
    atoms_md->s2x_lcl();
    
    return atoms_md;
}
/*--------------------------------------------
 
 --------------------------------------------*/
void ImportCFGMD::ml_import(PyMethodDef& tp_methods)
{
    tp_methods.ml_flags=METH_VARARGS | METH_KEYWORDS | METH_CLASS;
    tp_methods.ml_name="import_cfg";
    tp_methods.ml_doc=R"---(
    import(cfg_file)
    
    Import cfg file.
        
    Parameters
    ----------
    cfg_file : string
       path to cfg_file
    
    Returns
    -------
    None
   
    )---";
    
    tp_methods.ml_meth=(PyCFunction)(PyCFunctionWithKeywords)
    [](PyObject* type,PyObject* args,PyObject* kwds)->PyObject*
    {
        FuncAPI<std::string,OP<MAPP_MPI>> f("cfg",{"cfg_file","mpi"});
        f.noptionals=1;
        if(f(args,kwds)==-1) return NULL;
        
        MPI_Comm world=MPI_COMM_WORLD;
        if(f.val<1>().ob)
            world=reinterpret_cast<MAPP_MPI::Object*>(f.val<1>().ob)->world;
        
        ImportCFGMD read(world);
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
ImportCFGDMD::ImportCFGDMD(MPI_Comm& world):
ImportCFG(world)
{
}
/*--------------------------------------------
 constructor
 --------------------------------------------*/
ImportCFGDMD::~ImportCFGDMD()
{
}
/*--------------------------------------------
 reads the atom section of the cfg file
 --------------------------------------------*/
AtomsDMD* ImportCFGDMD::operator()(int N,const char* file)
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
    
    Atoms* atoms=new Atoms(world);
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
        AtomsDMD* atoms_dmd=new AtomsDMD(world,0,N);
        *atoms_dmd=*atoms;
        delete atoms;
        return atoms_dmd;
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
    int nelems=static_cast<int>(atoms->elements.nelems);
    if(!vec_list[3] || vec_list[3]->dim<2*nelems)
    {
        delete atoms;
        throw Print::vprintf("values for c and alpha of all elements should be provided for every atom");
    }
    
    const int vec_dim=vec_list[3]->dim;
    
    

    
    type0* alpha=dynamic_cast<Vec<type0>*>(vec_list[3])->begin();
    type0* c=alpha+nelems;
    
    
    int lcl_err=0;
    int max_ncmp_lcl=0,__max_ncmp=0;
    int jcmp;
    elem_type __nelems=static_cast<elem_type>(nelems);
    type0 __c_sum;
    for(int i=0;i<atoms->natms_lcl && lcl_err==0;i++,alpha+=vec_dim,c+=vec_dim)
    {
        __max_ncmp=0;
        jcmp=0;
        __c_sum=0.0;
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
                __c_sum+=c[j];
                __max_ncmp++;
            }
        }
        
        max_ncmp_lcl=MAX(max_ncmp_lcl,__max_ncmp);
        if(__c_sum>1.0) lcl_err=3;
    }
    
    int lcl_err0=lcl_err==1?1:0;
    int lcl_err1=lcl_err==2?1:0;
    int lcl_err2=lcl_err==3?1:0;
    int err0,err1,err2;
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
    
    MPI_Allreduce(&lcl_err2,&err2,1,MPI_INT,MPI_MAX,world);
    if(err2)
    {
        delete atoms;
        throw Print::vprintf("sum of present elements of one site should not exceed 1.0");
    }
    
    int c_dim;
    MPI_Allreduce(&max_ncmp_lcl,&c_dim,1,MPI_INT,MPI_MAX,world);
    
    AtomsDMD* atoms_dmd=new AtomsDMD(world,c_dim,N);
    *atoms_dmd=*atoms;
    
    alpha=dynamic_cast<Vec<type0>*>(vec_list[3])->begin();
    c=alpha+nelems;
    type0* __alpha=atoms_dmd->alpha->begin();
    type0* __c=atoms_dmd->c->begin();
    elem_type* __elem=atoms_dmd->elem->begin();
    
    type0 max_alpha_lcl=0.0;
    for(int i=0;i<atoms->natms_lcl;i++,alpha+=vec_dim,c+=vec_dim,__alpha+=c_dim,__c+=c_dim,__elem+=c_dim)
    {
        __max_ncmp=0;
        jcmp=0;
        for(elem_type j=0;j<__nelems;j++)
        {
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
    }
    MPI_Allreduce(&max_alpha_lcl,&atoms_dmd->max_alpha,1,MPI_DOUBLE,MPI_MAX,world);
    
    delete vec_list[3];
    delete atoms;
    
    atoms_dmd->s2x_lcl();
    return atoms_dmd;
}
/*--------------------------------------------
 
 --------------------------------------------*/
void ImportCFGDMD::ml_import(PyMethodDef& tp_methods)
{
    tp_methods.ml_flags=METH_VARARGS | METH_KEYWORDS | METH_CLASS;
    tp_methods.ml_name="import_cfg";
    tp_methods.ml_doc=R"---(
    import(N,cfg_file,mpi=None)
    
    Imports cfg file :py:class:`mapp.dmd.atoms`.
    
    This is a static function that is used to import a desired system presented in `cfg (Atomeye) <http://li.mit.edu/Archive/Graphics/A>`_ format. Please see explanations below.
    
    
        
    Parameters
    ----------
    N : int
       number of gaussian quadraure abscissas
    cfg_file : string
       path to cfg_file
    
    Returns
    -------
    mapp.dmd.atoms
       object that holds the configuration of the system
    
    Notes
    -----
    
    The cfg format for a DMD simulation is different from a regular cfg file structure. This is because in addition to giving the position of every atom one needs to specify the pertinent colors (:math:`c`) and :math:`\alpha` values. Suppose that our simulation is consist of n types of atoms. Therefore, for each site/atom we need to specify :math:`3 + 2n` values (:math:`3` for position, :math:`n` for :math:`\alpha` and :math:`n` for :math:`c`). For the same reason only `extended cfg <http://li.mit.edu/Archive/Graphics/A/#extended_CFG>`_ format file can be used. Ignoring the header part of cfg file each line (except for the lines defining mass and element) should look like this:
    
    .. math::
        \underbrace{s_x\quad s_y \quad s_z}_{\text{fractional coordinates}} \quad \alpha_0 \quad \alpha_1 \quad \cdots \quad \alpha_{n-1} \quad c_0 \quad c_1 \quad \cdots \quad c_{n-1}
    
    This will take care of per atom properties. It remains to determine the element to which each of :math:`c` and :math:`\alpha` components refer to. This is the same as the elements' in the file: :math:`\alpha_0` and :math:`c_0` refer to first element appearing in the file, :math:`\alpha_1` and :math:`c_1` to second and so on.
    
    Examples
    --------
    
    ::
    
       Number of particles = 24
       A = 1.0 Angstrom (basic length-scale)
       H0(1,1) =   2.472772 A
       H0(1,2) =   0.000000 A
       H0(1,3) =   0.000000 A
       H0(2,1) =   0.000000 A
       H0(2,2) =   4.038020 A
       H0(2,3) =   0.000000 A
       H0(3,1) =   0.000000 A
       H0(3,2) =   0.000000 A
       H0(3,3) =   6.994057 A
       .NO_VELOCITY.
       entry_count = 7
       auxiliary[0] = alpha_0 [reduced unit]
       auxiliary[1] = alpha_1 [reduced unit]
       auxiliary[2] = c_0 [reduced unit]
       auxiliary[3] = c_1 [reduced unit]
       55.845000 
       Fe 
       0.000000 0.50000 0.50000 0.1 0.0 1.0 -1.0
       0.000000 0.00000 0.00000 0.1 0.0 1.0 -1.0
       0.666667 0.50000 0.16666 0.1 0.0 1.0 -1.0
       0.666667 0.00000 0.66666 0.1 0.0 1.0 -1.0
       0.333333 0.50000 0.83333 0.1 0.0 1.0 -1.0
       0.333333 0.00000 0.33333 0.1 0.0 1.0 -1.0
       1.00794 
       H 
       0.000000 0.50000 0.00000 0.0 0.1 -1.0 0.00000001
       0.000000 0.00000 0.50000 0.0 0.1 -1.0 0.00000001
       0.666666 0.50000 0.66666 0.0 0.1 -1.0 0.00000001
       0.666666 0.00000 0.16666 0.0 0.1 -1.0 0.00000001
       0.000000 0.75000 0.25000 0.0 0.1 -1.0 0.00000001
       0.000000 0.25000 0.75000 0.0 0.1 -1.0 0.00000001
       0.000000 0.75000 0.75000 0.0 0.1 -1.0 0.00000001
       0.000000 0.25000 0.25000 0.0 0.1 -1.0 0.00000001
       0.666666 0.75000 0.41666 0.0 0.1 -1.0 0.00000001
       0.333333 0.50000 0.33333 0.0 0.1 -1.0 0.00000001
       0.666666 0.25000 0.91666 0.0 0.1 -1.0 0.00000001
       0.666666 0.75000 0.91666 0.0 0.1 -1.0 0.00000001
       0.333333 0.75000 0.08333 0.0 0.1 -1.0 0.00000001
       0.333333 0.25000 0.58333 0.0 0.1 -1.0 0.00000001
       0.333333 0.00000 0.83333 0.0 0.1 -1.0 0.00000001
       0.666666 0.25000 0.41666 0.0 0.1 -1.0 0.00000001
       0.333333 0.75000 0.58333 0.0 0.1 -1.0 0.00000001
       0.333333 0.25000 0.08333 0.0 0.1 -1.0 0.00000001
    )---";
    tp_methods.ml_meth=(PyCFunction)(PyCFunctionWithKeywords)
    [](PyObject* type,PyObject* args,PyObject* kwds)->PyObject*
    {
        FuncAPI<int,std::string,OP<MAPP_MPI>> f("cfg",{"N","cfg_file","mpi"});
        f.noptionals=1;
        f.logics<0>()[0]=VLogics("gt",1);
        if(f(args,kwds)==-1) return NULL;
        
        MPI_Comm world=MPI_COMM_WORLD;
        if(f.val<2>().ob)
            world=reinterpret_cast<MAPP_MPI::Object*>(f.val<2>().ob)->world;
        
        ImportCFGDMD read(world);
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

























