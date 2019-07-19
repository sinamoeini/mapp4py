/*--------------------------------------------
 Created by Sina on 07/23/13.
 Copyright (c) 2013 MIT. All rights reserved.
 --------------------------------------------*/
#include "comm.h"
#include "import.h"
#include "print.h"
using namespace MAPP_NS;
/*--------------------------------------------
 constructor
 --------------------------------------------*/
Import::Import(MPI_Comm& __world):
world(__world)
{
}
/*--------------------------------------------
 destructor
 --------------------------------------------*/
Import::~Import()
{

}
/*--------------------------------------------
 
 --------------------------------------------*/
FileReader::FileReader():
world(MPI_COMM_NULL),
rank(0),
file(NULL),
finished(false)
{
}
/*--------------------------------------------
 
 --------------------------------------------*/
FileReader::~FileReader()
{
    if(!rank && file) fst.close();
    delete [] file;
}
/*--------------------------------------------
 
 --------------------------------------------*/
FileReader::FileReader(Communication*& comm,const char* __file):
FileReader(comm->world,__file)
{
}
/*--------------------------------------------
 
 --------------------------------------------*/
FileReader::FileReader(Communication*& comm,const std::string __file):
FileReader(comm->world,__file.c_str())
{
}
/*--------------------------------------------
 
 --------------------------------------------*/
FileReader::FileReader(MPI_Comm __world,const std::string __file):
FileReader(__world,__file.c_str())
{
}
/*--------------------------------------------
 
 --------------------------------------------*/
FileReader::FileReader(MPI_Comm __world,const char* __file):
fst(Communication::get_rank(__world)==0 ? std::ifstream(__file,std::ios_base::in): std::ifstream()),
world(__world),
rank(Communication::get_rank(__world)),
line_no(0),
file(NULL),
finished(false)
{
    bool file_xst=true;
    if(!rank)
        if(fst.fail()) file_xst=false;

    
    MPI_Bcast(&file_xst,1,MPI_BYTE,0,world);
    
    if(!file_xst)
        throw Print::vprintf("failed to open %s",__file);
    
    if(!rank)
        str_line.reserve(1024);
    
    file=new char[strlen(__file)+1];
    memcpy(file,__file,(strlen(__file)+1)*sizeof(char));
}
/*--------------------------------------------
 read line and broadcast
 --------------------------------------------*/
bool FileReader::operator()(char*& line,size_t& line_cpcty)
{
    if(finished) return true;
    size_t size=0;
    
    if(!rank)
    {
        std::getline(fst,str_line);
        finished=fst.eof();
        if(!finished) line_no++;
        size=str_line.size();
    }
    
    MPI_Bcast(&finished,1,MPI_BYTE,0,world);
    if(finished) return true;

    MPI_Bcast(&size,sizeof(size_t),MPI_BYTE,0,world);
    if(size+1>line_cpcty)
    {
        delete [] line;
        line =new char[size+1];
        line_cpcty=size+1;
    }
    
    if(!rank) memcpy(line,str_line.c_str(),(size+1)*sizeof(char));
    MPI_Bcast(line,static_cast<int>(size)+1,MPI_CHAR,0,world);
    return false;
}
/*--------------------------------------------
 
 --------------------------------------------*/
size_t FileReader::get_line_no()
{
    
    MPI_Bcast(&line_no,sizeof(size_t),MPI_BYTE,0,world);
    return line_no;
}
#include <limits>
/*--------------------------------------------
 
 --------------------------------------------*/
bool FileReader::skip(const int n)
{
    if(!rank)
    {
        for(int i=0;i<n && !fst.eof();i++)
            fst.ignore(std::numeric_limits<std::streamsize>::max(),'\n');
        finished=fst.eof();
    }
    MPI_Bcast(&finished,1,MPI_BYTE,0,world);
    return finished;
}
/*--------------------------------------------
 hash sign remover:
 it removes the statement after hash sign,
 replaces the withe spaces with space;
 for example:
 line=
 "  this  is a    test, #here is comment";
 will be converted to:
 newline=
 "this is a test,";
 nargs=4;
 --------------------------------------------*/
size_t FileReader::hash_remover(char*& line)
{
    char* last;
    last=line;
    while(*last!='#' && *last!='\0' && *last!='\n')
        last++;
    *last='\0';
    
    char* ipos=line;
    char* jpos=line;
    char* kpos=line;
    size_t nargs=0;
    while(ipos!=last)
    {
        while(isspace(*ipos) && ipos!=last)
            ipos++;
        if(ipos==last)
        {
            if(kpos!=line)
                kpos--;
            *kpos='\0';
            continue;
        }
        
        jpos=ipos;
        nargs++;
        
        while(!isspace(*ipos) && ipos!=last)
            ipos++;
        memmove(kpos,jpos,ipos-jpos);
        if(ipos==last)
            continue;
        
        kpos+=ipos-jpos;
        *kpos=' ';
        kpos++;
    }
    
    return nargs;
}
/*--------------------------------------------
 hash sign remover:
 it removes the statement after hash sign,
 replaces the withe spaces with space;
 for example:
 line=
 "  this  is a    test, #here is comment";
 will be converted to:
 newline=
 "this is a test,";
 nargs=4;
 --------------------------------------------*/
size_t FileReader::parse_line(char* ipos,char**& args,size_t& args_cpcty)
{
    char* last;
    last=ipos;
    while(*last!='#' && *last!='\0' && *last!='\n')
        last++;
    *last='\0';
    
    
    size_t nargs=0;
    while(ipos!=last)
    {
        while(isspace(*ipos) && ipos!=last)
            ipos++;
        if(ipos==last)
            continue;
        
        if(nargs+1>args_cpcty)
        {
            char** args_=new char*[nargs+1];
            memcpy(args_,args,nargs*sizeof(char*));
            delete [] args;
            args=args_;
            args_cpcty++;
        }
        args[nargs]=ipos;
        nargs++;
        
        while(!isspace(*ipos) && ipos!=last)
            ipos++;
        if(ipos==last)
            continue;
        *ipos='\0';
        ipos++;
    }
    
    return nargs;
}
/*--------------------------------------------
 hash sign remover:
 it removes the statement after hash sign,
 replaces the withe spaces with space;
 for example:
 line=
 "  this  is a    test, #here is comment";
 will be converted to:
 newline=
 "this is a test,";
 nargs=4;
 --------------------------------------------*/
size_t FileReader::concatenate(size_t nargs,char** args
,char*& line)
{
    size_t lngth=0;
    size_t pos=0;
    for(size_t i=0;i<nargs;i++)
        lngth+=strlen(args[i]);
    
    lngth+=nargs;
    line=new char[lngth];
    for(size_t i=0;i<nargs;i++)
    {
        lngth=strlen(args[i]);
        for(size_t j=0;j<lngth;j++)
        {
            line[pos]=args[i][j];
            pos++;
        }
        
        if(i!=nargs-1)
        {
            line[pos]=' ';
            pos++;
        }
        else
        {
            line[pos]='\0';
            pos++;
        }
    }
    return lngth;
}

