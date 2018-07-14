#include "thermo_dynamics.h"
#include "MAPP.h"
#include "memory.h"
#include "print.h"
#include <math.h>



using namespace MAPP_NS;
/*--------------------------------------------
 constructor
 --------------------------------------------*/
ThermoDynamics::ThermoDynamics(const int __precision):
precision(__precision),
l_int(Print::max_length<int>()+1),
l_double(Print::max_length<double>()+__precision),
qs(NULL),
nqs(0)
{
}
/*--------------------------------------------
 destructor
 --------------------------------------------*/
ThermoDynamics::~ThermoDynamics()
{
    delete [] qs;
}
/*--------------------------------------------
 initiated before a run
 --------------------------------------------*/
void ThermoDynamics::init()
{
    int n=nqs*(l_double+1)+l_int+1+nqs;
    fprintf(MAPP::mapp_out," ");
    for(int i=0;i<n;i++)
        fprintf(MAPP::mapp_out,THERMO_LINE);
    fprintf(MAPP::mapp_out," \n");
    
    
    fprintf(MAPP::mapp_out,THERMO_DLMTR);
    print_header(MAPP::mapp_out,"step",l_int+1);
    
    for(int i=0;i<nqs;i++)
    {
        fprintf(MAPP::mapp_out,THERMO_DLMTR);
        qs[i].print_header(MAPP::mapp_out,l_double+1);
    }
    fprintf(MAPP::mapp_out,THERMO_DLMTR"\n");
    
    
    fprintf(MAPP::mapp_out,THERMO_DLMTR);
    for(int i=0;i<l_int+1;i++)
        fprintf(MAPP::mapp_out,THERMO_LINE);
    for(int i=0;i<nqs;i++)
    {
        fprintf(MAPP::mapp_out,THERMO_DLMTR);
        for(int j=0;j<l_double+1;j++)
            fprintf(MAPP::mapp_out,THERMO_LINE);
    }
    
    fprintf(MAPP::mapp_out,THERMO_DLMTR"\n");
    
}
/*--------------------------------------------
 print output
 --------------------------------------------*/
void ThermoDynamics::print(int step_no)
{
    fprintf(MAPP::mapp_out,THERMO_DLMTR" %0*d ",l_int-1,step_no);
    for(int i=0;i<nqs;i++)
        fprintf(MAPP::mapp_out,THERMO_DLMTR"%+*.*e ",l_double,precision,*qs[i].ptr);
    
    fprintf(MAPP::mapp_out,THERMO_DLMTR"\n");
}
/*--------------------------------------------
 finish after a run
 --------------------------------------------*/
void ThermoDynamics::fin()
{
    int n=nqs*(l_double+1)+l_int+1+nqs;
    fprintf(MAPP::mapp_out," ");
    for(int i=0;i<n;i++)
        fprintf(MAPP::mapp_out,THERMO_LINE);
    fprintf(MAPP::mapp_out," \n");
}
/*--------------------------------------------
 
 --------------------------------------------*/
ThermoQuantity::ThermoQuantity():
name(NULL),
name_lngth(0),
ptr(NULL)
{
}
/*--------------------------------------------
 
 --------------------------------------------*/
ThermoQuantity::~ThermoQuantity()
{
}
/*--------------------------------------------
 
 --------------------------------------------*/
ThermoQuantity::ThermoQuantity(const char* __name,type0& val):
name(__name),ptr(&val),
name_lngth(static_cast<int>(strlen(__name)))
{
}
/*--------------------------------------------
 
 --------------------------------------------*/
ThermoQuantity::ThermoQuantity(const ThermoQuantity& r):
name(r.name),ptr(r.ptr),
name_lngth(r.name_lngth)
{
}
/*--------------------------------------------
 
 --------------------------------------------*/
ThermoQuantity::ThermoQuantity(ThermoQuantity&& r):
name(r.name),ptr(r.ptr),
name_lngth(r.name_lngth)
{
}
/*--------------------------------------------
 
 --------------------------------------------*/
ThermoQuantity& ThermoQuantity::operator =(const ThermoQuantity& r)
{
    this->~ThermoQuantity();
    new (this) ThermoQuantity(r);
    return *this;
}
/*--------------------------------------------
 
 --------------------------------------------*/
ThermoQuantity& ThermoQuantity::operator =(ThermoQuantity&& r)
{
    this->~ThermoQuantity();
    new(this) ThermoQuantity(std::move(r));
    return *this;
}
/*--------------------------------------------
 
 --------------------------------------------*/
void ThermoQuantity::print_header(FILE* __out,int lngth)
{
    
    if(name_lngth<=lngth)
        fprintf(__out,"%*s%*s",(lngth-lngth/2)+name_lngth/2,name,lngth/2-name_lngth/2,"");
    else
        for(int i=0;i<lngth;i++) fprintf(__out,"%c",name[i]);
}

