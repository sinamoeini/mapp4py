#ifndef __MAPP__read_eam__
#define __MAPP__read_eam__
#include "global.h"
#include "import.h"
namespace MAPP_NS
{
    class ReadEAM
    {
    private:
        static void read_double(type0*,size_t,FileReader&,char*&,size_t&,char**&,size_t&);
        static void skip(size_t,FileReader&,char*&,size_t&);
        template<size_t N0>
        static void set_cuttoff(size_t,type0,size_t,type0(***)[N0],type0(***)[N0],type0**);
        
    public:
        
        
        
        template<size_t N0,size_t N1>
        static void funcfl(size_t,std::string*,
        type0&,type0&,size_t&,size_t&,
        type0(***&)[N0],type0(***&)[N0],type0(**&)[N1],
        type0**&,MPI_Comm=MPI_COMM_WORLD);
        
        template<size_t N0,size_t N1>
        static void setfl(size_t,char**,std::string,
        type0&,type0&,size_t&,size_t&,
        type0(***&)[N0],type0(***&)[N0],type0(**&)[N1],
        type0**&,MPI_Comm=MPI_COMM_WORLD);
        
        template<size_t N0,size_t N1>
        static void fs(size_t,char**,std::string,
        type0&,type0&,size_t&,size_t&,
        type0(***&)[N0],type0(***&)[N0],type0(**&)[N1],
        type0**&,MPI_Comm=MPI_COMM_WORLD);
        
        static type0 interpolate(type0*,size_t,type0,size_t);
        static void interpolate(size_t,type0,type0(*)[4]);
        static void interpolate(size_t,type0,type0(*)[5]);
        static void interpolate(size_t,type0,type0(*)[7]);
    };
}
using namespace MAPP_NS;
#include "print.h"
/*--------------------------------------------
 
 --------------------------------------------*/
template<size_t N0>
void ReadEAM::set_cuttoff(size_t nelems,type0 dr,size_t nr,type0(*** r_phi)[N0],type0(*** rho)[N0],type0** r_c)
{
    
    for(size_t i=0;i<nelems;i++)
        for(size_t j=0;j<i+1;j++)
        {
            size_t __nr=nr-1;
            while((r_phi[i][j][__nr][0]==0.0 && rho[i][j][__nr][0]==0.0  && rho[j][i][__nr][0]==0.0) && __nr>=0)
                __nr--;
            r_c[i][j]=r_c[j][i]=static_cast<type0>(__nr)*dr;
        }
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<size_t N0,size_t N1>
void ReadEAM::funcfl(size_t nfiles,std::string* files,
type0& dr,type0& drho,size_t& nr,size_t& nrho,
type0(***& r_phi)[N0],type0(***& rho)[N0],type0(**& F)[N1],
type0**& r_c,MPI_Comm world)
{
    FileReader* frs=new FileReader[nfiles];
    for(size_t i=0;i<nfiles;i++)
    {
        try
        {
            new (frs+i) FileReader(world,files[i]);
        }
        catch(char* err_msg)
        {
            delete [] frs;
            throw err_msg;
        }
    }
    
    size_t nelems=nfiles;
    size_t tot_nrs=0,tot_nrhos=0;
    size_t max_nr=0,max_nrho=0;
    size_t line_cpcty=0,args_cpcty=0,nargs=0;
    
    size_t* nrs=new size_t[nelems];
    size_t* nrhos=new size_t[nelems];
    type0* drs=new type0[nelems];
    type0* drhos=new type0[nelems];
    char* line=NULL;
    type0* buff=NULL;
    type0** __F=NULL;
    type0** __rho=NULL;
    type0** __Z=NULL;
    char** args=NULL;
    
    auto dealloc=[&nrs,&nrhos,&drs,&drhos,&line,&frs,&__F,&__rho,&__Z,&buff,&args]()->void
    {
        delete [] frs;
        
        delete [] nrs;
        delete [] nrhos;
        delete [] drs;
        delete [] drhos;
        
        delete [] line;
        
        delete [] buff;
        if(__F) delete [] *__F;
        delete [] __F;
        if(__rho) delete [] *__rho;
        delete [] __rho;
        if(__Z) delete [] *__Z;
        delete [] __Z;
        delete [] args;
    };
    
    
    
    for(size_t i=0;i<nelems;i++)
    {
        
        if(frs[i].skip(2))
        {
            dealloc();
            throw Print::vprintf("%s file ended immaturely",files[i].c_str());
        }
        
        
        if(frs[i](line,line_cpcty))
        {
            dealloc();
            throw Print::vprintf("%s file ended immaturely",files[i].c_str());
        }
        nargs=frs[i].hash_remover(line);
        
        
        if(sscanf(line," %zu %lf %zu %lf ",nrhos+i,drhos+i,nrs+i,drs+i)!=4)
        {
            dealloc();
            throw Print::vprintf("invalid line in %s file: %s",frs[i].file,line);
        }
        
        
        if(nrhos[i]<5)
        {
            dealloc();
            throw Print::vprintf("nrho in %s file should be larger than 5",frs[i].file);
        }
        if(nrs[i]<5)
        {
            dealloc();
            throw Print::vprintf("nr in %s file should be larger than 5",frs[i].file);
        }
        if(drhos[i]<=0.0)
        {
            dealloc();
            throw Print::vprintf("drho in %s file should be larger than 0.0",frs[i].file);
        }
        if(drs[i]<=0.0)
        {
            dealloc();
            throw Print::vprintf("dr in %s file should be larger than 0.0",frs[i].file);
        }
        
        tot_nrs+=nrs[i];
        tot_nrhos+=nrhos[i];
        max_nr=MAX(max_nr,nrs[i]);
        max_nrho=MAX(max_nrho,nrhos[i]);
    }
    
    
    buff=new type0[2*max_nr+max_nrho];
    __F=new type0*[nelems];
    *__F=new type0[tot_nrhos];
    __rho=new type0*[nelems];
    *__rho=new type0[tot_nrs];
    __Z=new type0*[nelems];
    *__Z=new type0[tot_nrs];
    
    for(size_t i=1;i<nelems;i++)
    {
        __F[i]=__F[i-1]+nrhos[i-1];
        __rho[i]=__rho[i-1]+nrs[i-1];
        __Z[i]=__Z[i-1]+nrs[i-1];
    }
    
    
    
    
    try
    {
        for(size_t i=0;i<nelems;i++)
        {
            read_double(buff,nrhos[i]+2*nrs[i],frs[i],line,line_cpcty,args,args_cpcty);
            memcpy(__F[i],buff,nrhos[i]*sizeof(type0));
            memcpy(__Z[i],buff+nrhos[i],nrs[i]*sizeof(type0));
            memcpy(__rho[i],buff+nrhos[i]+nrs[i],nrs[i]*sizeof(type0));
        }
    }
    catch (char* err_msg)
    {
        dealloc();
        throw err_msg;
    }

    type0 maxr=0.0;
    type0 maxrho=0.0;
    type0 maxdr=0.0;
    type0 maxdrho=0.0;
    
    for(size_t i=0;i<nelems;i++)
    {
        maxr=MAX(maxr,static_cast<type0>(nrs[i]-1)*drs[i]);
        maxrho=MAX(maxrho,static_cast<type0>(nrhos[i]-1)*drhos[i]);
        maxdr=MAX(maxdr,drs[i]);
        maxdrho=MAX(maxdrho,drhos[i]);
    }
    
    
    nr=static_cast<size_t>(maxr/maxdr+0.5);
    nrho=static_cast<size_t>(maxrho/maxdrho+0.5);
    dr=maxdr;
    drho=maxdrho;
    
    
    type0 (*arr0)[N0];
    using PARR0=type0 (*)[N0];
    using PPARR0=type0 (**)[N0];
    

    rho=new PPARR0[nelems];
    *rho=new PARR0[nelems*nelems];
    for(size_t i=1;i<nelems;i++) rho[i]=rho[i-1]+nelems;
    arr0=new type0[nr*nelems][N0];
    for(size_t i=0;i<nelems;i++)
    {
        for(size_t j=0;j<nelems;j++)
            rho[i][j]=arr0;
        arr0+=nr;
    }
    
    r_phi=new PPARR0[nelems];
    *r_phi=new PARR0[nelems*nelems];
    for(size_t i=1;i<nelems;i++) r_phi[i]=r_phi[i-1]+nelems;
    arr0=new type0[nr*nelems*(nelems+1)/2][N0];
    for(size_t i=0;i<nelems;i++)
        for(size_t j=0;j<i+1;j++)
        {
            r_phi[i][j]=arr0;
            r_phi[j][i]=arr0;
            arr0+=nr;
        }
    
    
    type0 (*arr1)[N1];
    using PARR1=type0 (*)[N1];
    F=new PARR1[nelems];
    arr1=new type0[nrho*nelems][N1];
    for(size_t i=0;i<nelems;i++)
    {
        F[i]=arr1;
        arr1+=nrho;
    }
    
    r_c=new type0*[nelems];
    *r_c=new type0[nelems*nelems];
    for(size_t i=1;i<nelems;i++) r_c[i]=r_c[i-1]+nelems;
    
    type0 r,rh,p,tmp0,tmp1;
    size_t l;
    for(size_t i=0;i<nelems;i++)
        for(size_t j=0;j<i+1;j++)
            for(size_t k=0;k<nr;k++)
            {
                r=static_cast<type0>(k)*dr;
                
                p=r/drs[i];
                l=static_cast<size_t>(p);
                l=MIN(l,nrs[i]-2);
                p-=l;
                p=MIN(p,1.0);
                tmp0=interpolate(__Z[i],nrs[i],p,k);
                
                
                p=r/drs[j];
                l=static_cast<size_t>(p);
                l=MIN(l,nrs[j]-2);
                p-=l;
                p=MIN(p,1.0);
                tmp1=interpolate(__Z[j],nrs[j],p,k);
                
                r_phi[i][j][k][0]=27.2*0.529*tmp0*tmp1;
            }
    
    for(size_t i=0;i<nelems;i++)
        for(size_t j=0;j<nr;j++)
        {
            r=static_cast<type0>(j)*dr;
            
            p=r/drs[i];
            l=static_cast<size_t>(p);
            l=MIN(l,nrs[i]-2);
            p-=l;
            p=MIN(p,1.0);
            
            rho[i][0][j][0]=interpolate(__rho[i],nrs[i],p,l);
        }
    
    for(size_t i=0;i<nelems;i++)
        for(size_t j=0;j<nrho;j++)
        {
            rh=static_cast<type0>(j)*drho;
            
            p=rh/drhos[i];
            l=static_cast<int> (p);
            l=MIN(l,nrhos[i]-2);
            p-=l;
            p=MIN(p,1.0);
            
            F[i][j][0]=interpolate(__F[i],nrhos[i],p,l);
        }

    
    dealloc();
    
    
    for(size_t i=0;i<nelems;i++)
    {
        interpolate(nrho,drho,F[i]);
        interpolate(nr,dr,rho[i][0]);
    }
    for(size_t i=0;i<nelems;i++)
        for(size_t j=0;j<i+1;j++)
            interpolate(nr,dr,r_phi[i][j]);
    set_cuttoff(nelems,dr,nr,r_phi,rho,r_c);
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<size_t N0,size_t N1>
void ReadEAM::setfl(size_t nelems,char** elems,std::string file,
type0& dr,type0& drho,size_t& nr,size_t& nrho,
type0(***& r_phi)[N0],type0(***& rho)[N0],type0(**& F)[N1],
type0**& r_c,MPI_Comm world)
{
    FileReader fr;
    try
    {
        new (&fr)FileReader(world,file);
    }
    catch(char* err_msg)
    {
        throw err_msg;
    }
    
    
    if(fr.skip(3))
        throw Print::vprintf("%s file ended immaturely",fr.file);
    

    size_t line_cpcty=0,args_cpcty=0,nargs=0;
    
    type0* buff=NULL;
    char* line=NULL;
    char** args=NULL;
    int* elem_ref=NULL;
    bool* found=NULL;
    
    F=NULL;
    rho=NULL;
    r_phi=NULL;
    r_c=NULL;
    
    auto dealloc=[&elem_ref,&line,&buff,&args,&found,&F,&rho,&r_phi,&r_c]()->void
    {
        if(r_phi)
        {
            if(*r_phi) delete [] **r_phi;
            delete [] *r_phi;
        }
        delete [] r_phi;
        
        if(rho)
        {
            if(*rho) delete [] **rho;
            delete [] *rho;
        }
        delete [] rho;
        
        if(F) delete [] *F;
        delete [] F;
        
        if(r_c) delete [] *r_c;
        delete [] r_c;
        
        delete [] line;
        delete [] elem_ref;
        delete [] buff;
        delete [] args;
        delete [] found;
    };
    

    if(fr(line,line_cpcty)) throw Print::vprintf("%s file ended immaturely",fr.file);
    nargs=fr.parse_line(line,args,args_cpcty);
    if(nargs<2)
    {
        dealloc();
        throw Print::vprintf("invalid line in %s file: %s",fr.file,line);
    }
    
    size_t nelem_refs=nargs-1;
    elem_ref=new int[nelem_refs];
    found=new bool[nelems];
    for(int i=0;i<nelems;i++) found[i]=false;
    for(int i=0;i<nelem_refs;i++)
    {
        int ielem=-1;
        for(int j=0;j<nelems&& ielem==-1 ;j++)
            if(!strcmp(args[1+i],elems[j]))
            {
                found[j]=true;
                ielem=j;
            }
        elem_ref[i]=ielem;
    }
    
    for(int i=0;i<nelems;i++)
    {
        if(!found[i])
        {
            dealloc();
            throw Print::vprintf("%s file does not contain parameters for element %s",
            fr.file,elems[i]);
        }
    }
    
    
    type0 __r_c;
    if(fr(line,line_cpcty)) throw Print::vprintf("%s file ended immaturely",fr.file);
    if(sscanf(line," %zu %lf %zu %lf %lf",&nrho,&drho,&nr,&dr,&__r_c)!=5)
    {
        dealloc();
        throw Print::vprintf("invalid line in %s file: %s",fr.file,line);
    }
    
    if(nrho<5)
    {
        dealloc();
        throw Print::vprintf("nrho in %s file should be larger than 5",fr.file);
    }
    if(nr<5)
    {
        dealloc();
        throw Print::vprintf("nr in %s file should be larger than 5",fr.file);
    }
    if(drho<=0.0)
    {
        dealloc();
        throw Print::vprintf("drho in %s file should be larger than 0.0",fr.file);
    }
    if(dr<=0.0)
    {
        dealloc();
        throw Print::vprintf("dr in %s file should be larger than 0.0",fr.file);
    }
    
    
    type0 (*arr0)[N0];
    using PARR0=type0 (*)[N0];
    using PPARR0=type0 (**)[N0];
    
    
    rho=new PPARR0[nelems];
    *rho=new PARR0[nelems*nelems];
    for(size_t i=1;i<nelems;i++) rho[i]=rho[i-1]+nelems;
    arr0=new type0[nr*nelems][N0];
    for(size_t i=0;i<nelems;i++)
    {
        for(size_t j=0;j<nelems;j++)
            rho[i][j]=arr0;
        arr0+=nr;
    }
    
    r_phi=new PPARR0[nelems];
    *r_phi=new PARR0[nelems*nelems];
    for(size_t i=1;i<nelems;i++) r_phi[i]=r_phi[i-1]+nelems;
    arr0=new type0[nr*nelems*(nelems+1)/2][N0];
    for(size_t i=0;i<nelems;i++)
        for(size_t j=0;j<i+1;j++)
        {
            r_phi[i][j]=arr0;
            r_phi[j][i]=arr0;
            arr0+=nr;
        }
    
    
    type0 (*arr1)[N1];
    using PARR1=type0 (*)[N1];
    F=new PARR1[nelems];
    arr1=new type0[nrho*nelems][N1];
    for(size_t i=0;i<nelems;i++)
    {
        F[i]=arr1;
        arr1+=nrho;
    }
    
    r_c=new type0*[nelems];
    *r_c=new type0[nelems*nelems];
    for(size_t i=1;i<nelems;i++) r_c[i]=r_c[i-1]+nelems;
    buff=new type0[MAX(nr+nrho,nelem_refs*(nelem_refs+1)/2*nr)];
    for(size_t i=0;i<nelem_refs;i++)
    {
        if(fr.skip(1)) throw Print::vprintf("%s file ended immaturely",fr.file);
        int ielem=elem_ref[i];
        if(ielem==-1)
        {
            try
            {
                skip(nrho+nr,fr,line,line_cpcty);
            }
            catch (char* err_msg)
            {
                dealloc();
                throw err_msg;
            }
        }
        else
        {
            try
            {
                read_double(buff,nrho+nr,fr,line,line_cpcty,args,args_cpcty);
            }
            catch (char* err_msg)
            {
                dealloc();
                throw err_msg;
            }
            
            
            
            for(size_t j=0;j<nrho;j++)
                F[ielem][j][0]=buff[j];
            for(size_t j=0;j<nr;j++)
                rho[ielem][0][j][0]=buff[j+nrho];
        }
    }
    
    try
    {
        read_double(buff,nelem_refs*(nelem_refs+1)/2*nr,fr,line,line_cpcty,args,args_cpcty);
    }
    catch (char* err_msg)
    {
        dealloc();
        throw err_msg;
    }
    
    size_t ipos=0;
    for(size_t i=0;i<nelem_refs;i++)
        for(size_t j=0;j<i+1;j++)
        {
            int ielem=elem_ref[i];
            int jelem=elem_ref[j];
            
            if(ielem==-1 || jelem==-1)
            {
                ipos+=nr;
                continue;
            }
            
            for(size_t k=0;k<nr;k++)
                r_phi[ielem][jelem][k][0]=buff[ipos++];
        }

    delete [] line;
    delete [] elem_ref;
    delete [] buff;
    delete [] args;
    delete [] found;
    
    for(size_t i=0;i<nelems;i++)
    {
        interpolate(nrho,drho,F[i]);
        interpolate(nr,dr,rho[i][0]);
    }
    for(size_t i=0;i<nelems;i++)
        for(size_t j=0;j<i+1;j++)
            interpolate(nr,dr,r_phi[i][j]);
    set_cuttoff(nelems,dr,nr,r_phi,rho,r_c);

}
/*--------------------------------------------
 
 --------------------------------------------*/
template<size_t N0,size_t N1>
void ReadEAM::fs(size_t nelems,char** elems,std::string file,
type0& dr,type0& drho,size_t& nr,size_t& nrho,
type0(***& r_phi)[N0],type0(***& rho)[N0],type0(**& F)[N1],
type0**& r_c,MPI_Comm world)
{
    FileReader fr;
    try
    {
        new (&fr)FileReader(world,file);
    }
    catch(char* err_msg)
    {
        throw err_msg;
    }
    
    
    if(fr.skip(3))
        throw Print::vprintf("%s file ended immaturely",fr.file);
    

    size_t line_cpcty=0,args_cpcty=0,nargs=0;
    
    type0* buff=NULL;
    char* line=NULL;
    char** args=NULL;
    int* elem_ref=NULL;
    bool* found=NULL;
    
    F=NULL;
    rho=NULL;
    r_phi=NULL;
    r_c=NULL;
    
    auto dealloc=[&elem_ref,&line,&buff,&args,&found,&F,&rho,&r_phi,&r_c]()->void
    {
        if(r_phi)
        {
            if(*r_phi) delete [] **r_phi;
            delete [] *r_phi;
        }
        delete [] r_phi;
        
        if(rho)
        {
            if(*rho) delete [] **rho;
            delete [] *rho;
        }
        delete [] rho;
        
        if(F) delete [] *F;
        delete [] F;
        
        if(r_c) delete [] *r_c;
        delete [] r_c;
        
        delete [] line;
        delete [] elem_ref;
        delete [] buff;
        delete [] args;
        delete [] found;
    };
    
    
    if(fr(line,line_cpcty)) throw Print::vprintf("%s file ended immaturely",fr.file);
    nargs=fr.parse_line(line,args,args_cpcty);
    if(nargs<2)
    {
        dealloc();
        throw Print::vprintf("invalid line in %s file: %s",fr.file,line);
    }
    
    size_t nelem_refs=nargs-1;
    elem_ref=new int[nelem_refs];
    found=new bool[nelems];
    for(int i=0;i<nelems;i++) found[i]=false;
    for(int i=0;i<nelem_refs;i++)
    {
        int ielem=-1;
        for(int j=0;j<nelems&& ielem==-1 ;j++)
            if(!strcmp(args[1+i],elems[j]))
            {
                found[j]=true;
                ielem=j;
            }
        elem_ref[i]=ielem;
    }
    
    for(int i=0;i<nelems;i++)
    {
        if(!found[i])
        {
            dealloc();
            throw Print::vprintf("%s file does not contain parameters for element %s",
            fr.file,elems[i]);
        }
    }
    
    
    type0 __r_c;
    
    if(fr(line,line_cpcty)) throw Print::vprintf("%s file ended immaturely",fr.file);
    if(sscanf(line," %zu %lf %zu %lf %lf",&nrho,&drho,&nr,&dr,&__r_c)!=5)
    {
        dealloc();
        throw Print::vprintf("invalid line in %s file: %s",fr.file,line);
    }
    
    if(nrho<5)
    {
        dealloc();
        throw Print::vprintf("nrho in %s file should be larger than 5",fr.file);
    }
    if(nr<5)
    {
        dealloc();
        throw Print::vprintf("nr in %s file should be larger than 5",fr.file);
    }
    if(drho<=0.0)
    {
        dealloc();
        throw Print::vprintf("drho in %s file should be larger than 0.0",fr.file);
    }
    if(dr<=0.0)
    {
        dealloc();
        throw Print::vprintf("dr in %s file should be larger than 0.0",fr.file);
    }
    
    
    type0 (*arr0)[N0];
    using PARR0=type0 (*)[N0];
    using PPARR0=type0 (**)[N0];
    
    
    rho=new PPARR0[nelems];
    *rho=new PARR0[nelems*nelems];
    for(size_t i=1;i<nelems;i++) rho[i]=rho[i-1]+nelems;
    arr0=new type0[nr*nelems*nelems][N0];
    for(size_t i=0;i<nelems;i++)
    {
        for(size_t j=0;j<nelems;j++)
        {
            rho[i][j]=arr0;
            arr0+=nr;
        }
    }
    
    r_phi=new PPARR0[nelems];
    *r_phi=new PARR0[nelems*nelems];
    for(size_t i=1;i<nelems;i++) r_phi[i]=r_phi[i-1]+nelems;
    arr0=new type0[nr*nelems*(nelems+1)/2][N0];
    for(size_t i=0;i<nelems;i++)
        for(size_t j=0;j<i+1;j++)
        {
            r_phi[i][j]=arr0;
            r_phi[j][i]=arr0;
            arr0+=nr;
        }
    
    
    type0 (*arr1)[N1];
    using PARR1=type0 (*)[N1];
    F=new PARR1[nelems];
    arr1=new type0[nrho*nelems][N1];
    for(size_t i=0;i<nelems;i++)
    {
        F[i]=arr1;
        arr1+=nrho;
    }
    
    r_c=new type0*[nelems];
    *r_c=new type0[nelems*nelems];
    for(size_t i=1;i<nelems;i++) r_c[i]=r_c[i-1]+nelems;
    
    buff=new type0[MAX(nrho+nr*nelem_refs,nelem_refs*(nelem_refs+1)/2*nr)];
    for(size_t i=0;i<nelem_refs;i++)
    {
        if(fr.skip(1)) throw Print::vprintf("%s file ended immaturely",fr.file);
        int ielem=elem_ref[i];
        if(ielem==-1)
        {
            try
            {
                skip(nrho+nr*nelem_refs,fr,line,line_cpcty);
            }
            catch (char* err_msg)
            {
                dealloc();
                throw err_msg;
            }
        }
        else
        {
            try
            {
                read_double(buff,nrho+nr*nelem_refs,fr,line,line_cpcty,args,args_cpcty);
            }
            catch (char* err_msg)
            {
                dealloc();
                throw err_msg;
            }
            
            
            
            for(size_t j=0;j<nrho;j++)
                F[ielem][j][0]=buff[j];
            size_t ipos=nrho;
            
            for(size_t j=0;j<nelem_refs;j++)
            {
                int jelem=elem_ref[j];
                if(jelem==-1)
                {
                    ipos+=nr;
                    continue;
                }
                
                for(size_t k=0;k<nr;k++)
                    rho[ielem][jelem][k][0]=buff[ipos++];
            }
                
            
        }
    }
    
    try
    {
        read_double(buff,nelem_refs*(nelem_refs+1)/2*nr,fr,line,line_cpcty,args,args_cpcty);
    }
    catch (char* err_msg)
    {
        dealloc();
        throw err_msg;
    }
    
    size_t ipos=0;
    for(size_t i=0;i<nelem_refs;i++)
        for(size_t j=0;j<i+1;j++)
        {
            int ielem=elem_ref[i];
            int jelem=elem_ref[j];
            
            if(ielem==-1 || jelem==-1)
            {
                ipos+=nr;
                continue;
            }
            
            for(size_t k=0;k<nr;k++)
                r_phi[ielem][jelem][k][0]=buff[ipos++];
        }

    delete [] line;
    delete [] elem_ref;
    delete [] buff;
    delete [] args;
    delete [] found;
    
    for(size_t i=0;i<nelems;i++)
        interpolate(nrho,drho,F[i]);
    
    for(size_t i=0;i<nelems;i++)
        for(size_t j=0;j<nelems;j++)
            interpolate(nr,dr,rho[i][j]);
            
    for(size_t i=0;i<nelems;i++)
        for(size_t j=0;j<i+1;j++)
            interpolate(nr,dr,r_phi[i][j]);
    set_cuttoff(nelems,dr,nr,r_phi,rho,r_c);
}



#endif
