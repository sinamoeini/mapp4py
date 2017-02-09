#include "read_eam.h"
using namespace MAPP_NS;
/*--------------------------------------------
 
 --------------------------------------------*/
type0 ReadEAM::interpolate(type0* arr,size_t n,type0 p,size_t k)
{
    type0 coef0,coef1,coef2,coef3,tmp;
    coef0=arr[k];
    
    if(k==0)
    {
        coef1=arr[1]-arr[0];
        tmp=0.5*(arr[2]-arr[0]);
        coef2=(-arr[0]+2.0*arr[1]-arr[2])/2.0;
        coef3=(arr[0]-2.0*arr[1]+arr[2])/2.0;
        return ((coef3*p+coef2)*p+coef1)*p+coef0;
    }
    else if(k==1)
    {
        coef1=0.5*(arr[2]-arr[0]);
        tmp=((arr[0]-arr[4])+8.0*(arr[3]-arr[1]))/12.0;
        coef2=(11.0*arr[0]-28.0*arr[1]+24.0*arr[2]-8.0*arr[3]+arr[4])/12.0;
        coef3=(-5.0*arr[0]+16.0*arr[1]-18.0*arr[2]+8.0*arr[3]-arr[4])/12.0;
        return ((coef3*p+coef2)*p+coef1)*p+coef0;
    }
    else if(k==n-2)
    {
        coef1=0.5*(arr[n-1]-arr[n-3]);
        tmp=arr[n-1]-arr[n-2];
        coef2=arr[n-3]-3.0*arr[n-2]+2.0*arr[n-1];
        coef3=0.5*(-arr[n-3]+3.0*arr[n-2]-2.0*arr[n-1]);
        return ((coef3*p+coef2)*p+coef1)*p+coef0;
    }
    else if(k==n-1)
    {
        coef1=arr[n-1]-arr[n-2];
        return coef0+coef1*p;
    }
    else
    {
        coef1=((arr[k-2]-arr[k+2])+
               8.0*(arr[k+1]-arr[k-1]))/12.0;
        tmp=((arr[k-1]-arr[k+3])+
             8.0*(arr[k+2]-arr[k]))/12.0;
        
        coef2=(-2.0*arr[k-2]+15.0*arr[k-1]-28.0*arr[k]+21.0*arr[k+1]-6.0*arr[k+2])/12.0;
        coef3=(arr[k-2]-7.0*arr[k-1]+16.0*arr[k]-17.0*arr[k+1]+7.0*arr[k+2])/12.0;
        return ((coef3*p+coef2)*p+coef1)*p+coef0;
    }
}
/*--------------------------------------------

 --------------------------------------------*/
void ReadEAM::interpolate(size_t n,type0 delta,type0(*spline)[4])
{
    spline[0][1]=spline[1][0]-spline[0][0];
    spline[1][1]=0.5*(spline[2][0]-spline[0][0]);
    spline[n-2][1]=0.5*(spline[n-1][0]-spline[n-3][0]);
    spline[n-1][1]=spline[n-1][0]-spline[n-2][0];
    
    for(size_t i=2;i<n-2;i++)
        spline[i][1]=((spline[i-2][0]-spline[i+2][0])+
        8.0*(spline[i+1][0]-spline[i-1][0]))/12.0;
    
    for(size_t i=0;i<n-1;i++)
    {
        spline[i][2]=3.0*(spline[i+1][0]-spline[i][0])-
        2.0*spline[i][1]-spline[i+1][1];
        spline[i][3]=spline[i][1]+spline[i+1][1]-
        2.0*(spline[i+1][0]-spline[i][0]);
    }
    
    spline[n-1][2]=0.0;
    spline[n-1][3]=0.0;
}
/*--------------------------------------------

 --------------------------------------------*/
void ReadEAM::interpolate(size_t n,type0 delta,type0(*spline)[5])
{
    spline[0][1]=spline[1][0]-spline[0][0];
    spline[1][1]=0.5*(spline[2][0]-spline[0][0]);
    spline[n-2][1]=0.5*(spline[n-1][0]-spline[n-3][0]);
    spline[n-1][1]=spline[n-1][0]-spline[n-2][0];
    
    for(size_t i=2;i<n-2;i++)
        spline[i][1]=((spline[i-2][0]-spline[i+2][0])+
        8.0*(spline[i+1][0]-spline[i-1][0]))/12.0;
    
    for(size_t i=0;i<n-1;i++)
    {
        spline[i][2]=3.0*(spline[i+1][0]-spline[i][0])-
        2.0*spline[i][1]-spline[i+1][1];
        spline[i][3]=spline[i][1]+spline[i+1][1]-
        2.0*(spline[i+1][0]-spline[i][0]);
    }
    
    spline[n-1][2]=0.0;
    spline[n-1][3]=0.0;
    
    for(size_t i=0;i<n-1;i++)
    {
        spline[i][4]=(spline[i+1][2]-spline[i][2])/6.0-0.5*spline[i][3];
    }
    spline[n-1][4]=0.0;
    
}
/*--------------------------------------------

 --------------------------------------------*/
void ReadEAM::interpolate(size_t n,type0 delta,type0(*spline)[7])
{
    spline[0][1]=spline[1][0]-spline[0][0];
    spline[1][1]=0.5*(spline[2][0]-spline[0][0]);
    spline[n-2][1]=0.5*(spline[n-1][0]-spline[n-3][0]);
    spline[n-1][1]=spline[n-1][0]-spline[n-2][0];
    
    for(int i=2;i<n-2;i++)
        spline[i][1]=((spline[i-2][0]-spline[i+2][0])+
                      8.0*(spline[i+1][0]-spline[i-1][0]))/12.0;
    
    for(size_t i=0;i<n-1;i++)
    {
        spline[i][2]=3.0*(spline[i+1][0]-spline[i][0])-
        2.0*spline[i][1]-spline[i+1][1];
        spline[i][3]=spline[i][1]+spline[i+1][1]-
        2.0*(spline[i+1][0]-spline[i][0]);
    }
    
    spline[n-1][2]=0.0;
    spline[n-1][3]=0.0;
    for(size_t i=0;i<n;i++)
    {
        spline[i][4]=spline[i][1]/delta;
        spline[i][5]=2.0*spline[i][2]/delta;
        spline[i][6]=3.0*spline[i][3]/delta;
    }
}
/*--------------------------------------------
 
 --------------------------------------------*/
void ReadEAM::read_double(type0* buff,size_t n,FileReader& fr,char*& line,size_t& line_cpcty,char**& args,size_t& args_cpcty)
{
    size_t i=0;
    size_t nargs;
    
    while(i<n)
    {
        if(fr(line,line_cpcty))
            throw Print::vprintf("%s file ended immaturely",fr.file);
        nargs=fr.parse_line(line,args,args_cpcty);
        if(nargs+i>n)
            throw Print::vprintf("%s file has extra arguments",fr.file);
        
        for(int j=0;j<nargs;j++)
            buff[i++]=atof(args[j]);
    }
}
/*--------------------------------------------
 
 --------------------------------------------*/
void ReadEAM::skip(size_t n,FileReader& fr,char*& line,size_t& line_cpcty)
{
    size_t i=0;
    size_t nargs;
    
    while(i<n)
    {
        fr(line,line_cpcty);
        if(fr.finished) throw Print::vprintf("%s file ended immaturely",fr.file);
        nargs=fr.hash_remover(line);
        if(nargs+i>n)
            throw Print::vprintf("%s file has extra arguments",fr.file);
        i+=nargs;
    }
}
