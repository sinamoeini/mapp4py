/*--------------------------------------------
 Created by Sina on 06/27/14.
 Copyright (c) 2014 MIT. All rights reserved.
 --------------------------------------------*/
#include "elements.h"
#include "memory.h"
#include <limits>
using namespace MAPP_NS;
/*--------------------------------------------
 constructor
 --------------------------------------------*/
Elements::Elements():
nelems(0),
__nelems(0),
names(NULL),
masses(NULL),
__nelems__(nelems,"nelems")
{
}
/*--------------------------------------------
 copy constructor
 --------------------------------------------*/
Elements::Elements(const Elements& r):
nelems(r.nelems),
__nelems(r.__nelems),
names(NULL),
masses(NULL),
__nelems__(nelems,"nelems")
{
    if(nelems==0) return;
    ptrdiff_t len=0;
    if(names)
        len=(strchr(names[nelems-1],'\0')-names[0])+1;
    
    char** names=new char*[nelems];
    *names=new char[len];
    memcpy(*names,*(r.names),len*sizeof(char));
    masses=new type0[nelems];
    memcpy(masses,r.masses,nelems*sizeof(type0));
}
/*--------------------------------------------
 copy constructor
 --------------------------------------------*/
Elements::Elements(Elements&& r):
nelems(r.nelems),
__nelems(r.__nelems),
names(r.names),
masses(r.masses),
__nelems__(std::move(r.__nelems__))
{
    r.names=NULL;
    r.masses=NULL;
}
/*--------------------------------------------
 destructor
 --------------------------------------------*/
Elements::~Elements()
{
    Memory::dealloc(masses);
    Memory::dealloc(names);
}
/*--------------------------------------------
 copy assignment
 --------------------------------------------*/
Elements& Elements::operator=(const Elements& r)
{
    this->~Elements();
    new (this) Elements(r);
    return *this;
}
/*--------------------------------------------
 move assignment
 --------------------------------------------*/
Elements& Elements::operator=(Elements&& r)
{
    this->~Elements();
    new (this) Elements(std::move(r));
    return *this;
}
/*--------------------------------------------
 comparison 
 --------------------------------------------*/
bool Elements::operator==(const Elements& r)
{
    if(nelems!=r.nelems) return false;
    for(size_t i=0;i<nelems;i++)
        if(strcmp(names[i],r.names[i])!=0)
            return false;

    return true;
}
/*--------------------------------------------
 add a new element with given mass
 --------------------------------------------*/
elem_type Elements::add_type(const type0 mass,const char* name)
{
    for(elem_type i=0;i<__nelems;i++)
        if(!strcmp(name,names[i]))
            return i;
    if(__nelems==std::numeric_limits<elem_type>::max())
        throw "cannot have more than 256 elements";
    
    
    Memory::grow(masses,nelems,nelems+1);
    masses[nelems]=mass;
    
    
    size_t new_len=strlen(name)+1;
    ptrdiff_t len=0;
    if(names)
        len=(strchr(names[nelems-1],'\0')-names[0])+1;
    
    char** __names=new char*[nelems+1];
    *__names=new char[len+new_len];
    
    if(names) memcpy(*__names,*names,len*sizeof(char));
    memcpy(*__names+len,name,new_len*sizeof(char));
    for(size_t i=1;i<nelems;i++)
        __names[i]=__names[i-1]+(names[i]-names[i-1]);
    __names[nelems]=*__names+len;
    Memory::dealloc(names);
    names=__names;
    
    nelems++;
    return __nelems++;
}
/*--------------------------------------------
 find an elemrnt number given its name
 --------------------------------------------*/
elem_type Elements::find(const char* name)
{
    for(elem_type i=0;i<__nelems;i++)
        if(!strcmp(name,names[i]))
            return i;
    
    throw 0;
}
/*--------------------------------------------
 find an elemrnt number given its name
 --------------------------------------------*/
PyObject* Elements::get_dict()
{
    if(nelems==0) Py_RETURN_NONE;
    PyObject* dict= PyDict_New();
    for(size_t i=0;i<nelems;i++)
        PyDict_SetItemString(dict,names[i],PyInt_FromSize_t(i));
    return dict;
}
/*--------------------------------------------
 find a type without error
 --------------------------------------------*/
void Elements::assign_color_rad(const char* name,type0 (&mat)[4])
{
    if(0){}
    else if(!strcmp(name,"H"))
    {
        mat[0]=0.800000;
        mat[1]=0.800000;
        mat[2]=0.800000;
        mat[3]=0.435000;
    }
    else if(!strcmp(name,"He"))
    {
        mat[0]=0.643137;
        mat[1]=0.666667;
        mat[2]=0.678431;
        mat[3]=0.700000;
    }
    else if(!strcmp(name,"Li"))
    {
        mat[0]=0.700000;
        mat[1]=0.700000;
        mat[2]=0.700000;
        mat[3]=1.519900;
    }
    else if(!strcmp(name,"Be"))
    {
        mat[0]=0.643137;
        mat[1]=0.666667;
        mat[2]=0.678431;
        mat[3]=1.143000;
    }
    else if(!strcmp(name,"B"))
    {
        mat[0]=0.900000;
        mat[1]=0.400000;
        mat[2]=0.000000;
        mat[3]=0.975000;
    }
    else if(!strcmp(name,"C"))
    {
        mat[0]=0.350000;
        mat[1]=0.350000;
        mat[2]=0.350000;
        mat[3]=0.655000;
    }
    else if(!strcmp(name,"N"))
    {
        mat[0]=0.200000;
        mat[1]=0.200000;
        mat[2]=0.800000;
        mat[3]=0.750000;
    }
    else if(!strcmp(name,"O"))
    {
        mat[0]=0.800000;
        mat[1]=0.200000;
        mat[2]=0.200000;
        mat[3]=0.730000;
    }
    else if(!strcmp(name,"F"))
    {
        mat[0]=0.700000;
        mat[1]=0.850000;
        mat[2]=0.450000;
        mat[3]=0.720000;
    }
    else if(!strcmp(name,"Ne"))
    {
        mat[0]=0.643137;
        mat[1]=0.666667;
        mat[2]=0.678431;
        mat[3]=1.600000;
    }
    else if(!strcmp(name,"Na"))
    {
        mat[0]=0.600000;
        mat[1]=0.600000;
        mat[2]=0.600000;
        mat[3]=1.857900;
    }
    else if(!strcmp(name,"Mg"))
    {
        mat[0]=0.600000;
        mat[1]=0.600000;
        mat[2]=0.700000;
        mat[3]=1.604700;
    }
    else if(!strcmp(name,"Al"))
    {
        mat[0]=0.643137;
        mat[1]=0.666667;
        mat[2]=0.678431;
        mat[3]=1.431800;
    }
    else if(!strcmp(name,"Si"))
    {
        mat[0]=0.690196;
        mat[1]=0.768627;
        mat[2]=0.870588;
        mat[3]=1.175800;
    }
    else if(!strcmp(name,"P"))
    {
        mat[0]=0.100000;
        mat[1]=0.700000;
        mat[2]=0.300000;
        mat[3]=1.060000;
    }
    else if(!strcmp(name,"S"))
    {
        mat[0]=0.950000;
        mat[1]=0.900000;
        mat[2]=0.200000;
        mat[3]=1.020000;
    }
    else if(!strcmp(name,"Cl"))
    {
        mat[0]=0.150000;
        mat[1]=0.500000;
        mat[2]=0.100000;
        mat[3]=0.990000;
    }
    else if(!strcmp(name,"Ar"))
    {
        mat[0]=0.643137;
        mat[1]=0.666667;
        mat[2]=0.678431;
        mat[3]=1.900000;
    }
    else if(!strcmp(name,"K"))
    {
        mat[0]=0.576471;
        mat[1]=0.439216;
        mat[2]=0.858824;
        mat[3]=2.262000;
    }
    else if(!strcmp(name,"Ca"))
    {
        mat[0]=0.800000;
        mat[1]=0.800000;
        mat[2]=0.700000;
        mat[3]=1.975800;
    }
    else if(!strcmp(name,"Sc"))
    {
        mat[0]=0.643137;
        mat[1]=0.666667;
        mat[2]=0.678431;
        mat[3]=1.654500;
    }
    else if(!strcmp(name,"Ti"))
    {
        mat[0]=0.643137;
        mat[1]=0.666667;
        mat[2]=0.678431;
        mat[3]=1.475500;
    }
    else if(!strcmp(name,"V"))
    {
        mat[0]=0.643137;
        mat[1]=0.666667;
        mat[2]=0.678431;
        mat[3]=1.309000;
    }
    else if(!strcmp(name,"Cr"))
    {
        mat[0]=0.000000;
        mat[1]=0.800000;
        mat[2]=0.000000;
        mat[3]=1.249000;
    }
    else if(!strcmp(name,"Mn"))
    {
        mat[0]=0.643137;
        mat[1]=0.666667;
        mat[2]=0.678431;
        mat[3]=1.350000;
    }
    else if(!strcmp(name,"Fe"))
    {
        mat[0]=0.517647;
        mat[1]=0.576471;
        mat[2]=0.652941;
        mat[3]=1.241100;
    }
    else if(!strcmp(name,"Co"))
    {
        mat[0]=0.643137;
        mat[1]=0.666667;
        mat[2]=0.678431;
        mat[3]=1.253500;
    }
    else if(!strcmp(name,"Ni"))
    {
        mat[0]=0.257255;
        mat[1]=0.266667;
        mat[2]=0.271373;
        mat[3]=1.246000;
    }
    else if(!strcmp(name,"Cu"))
    {
        mat[0]=0.950000;
        mat[1]=0.790074;
        mat[2]=0.013859;
        mat[3]=1.278000;
    }
    else if(!strcmp(name,"Zn"))
    {
        mat[0]=0.643137;
        mat[1]=0.666667;
        mat[2]=0.678431;
        mat[3]=1.332500;
    }
    else if(!strcmp(name,"Ga"))
    {
        mat[0]=0.900000;
        mat[1]=0.000000;
        mat[2]=1.000000;
        mat[3]=1.350100;
    }
    else if(!strcmp(name,"Ge"))
    {
        mat[0]=0.643137;
        mat[1]=0.666667;
        mat[2]=0.678431;
        mat[3]=1.224800;
    }
    else if(!strcmp(name,"As"))
    {
        mat[0]=1.000000;
        mat[1]=1.000000;
        mat[2]=0.300000;
        mat[3]=1.200000;
    }
    else if(!strcmp(name,"Se"))
    {
        mat[0]=0.643137;
        mat[1]=0.666667;
        mat[2]=0.678431;
        mat[3]=1.160000;
    }
    else if(!strcmp(name,"Br"))
    {
        mat[0]=0.500000;
        mat[1]=0.080000;
        mat[2]=0.120000;
        mat[3]=1.140000;
    }
    else if(!strcmp(name,"Kr"))
    {
        mat[0]=0.643137;
        mat[1]=0.666667;
        mat[2]=0.678431;
        mat[3]=2.000000;
    }
    else if(!strcmp(name,"Rb"))
    {
        mat[0]=0.643137;
        mat[1]=0.666667;
        mat[2]=0.678431;
        mat[3]=2.470000;
    }
    else if(!strcmp(name,"Sr"))
    {
        mat[0]=0.643137;
        mat[1]=0.666667;
        mat[2]=0.678431;
        mat[3]=2.151300;
    }
    else if(!strcmp(name,"Y"))
    {
        mat[0]=0.643137;
        mat[1]=0.666667;
        mat[2]=0.678431;
        mat[3]=1.823700;
    }
    else if(!strcmp(name,"Zr"))
    {
        mat[0]=0.643137;
        mat[1]=0.666667;
        mat[2]=0.678431;
        mat[3]=1.615600;
    }
    else if(!strcmp(name,"Nb"))
    {
        mat[0]=0.643137;
        mat[1]=0.666667;
        mat[2]=0.678431;
        mat[3]=1.431800;
    }
    else if(!strcmp(name,"Mo"))
    {
        mat[0]=0.643137;
        mat[1]=0.666667;
        mat[2]=0.678431;
        mat[3]=1.362600;
    }
    else if(!strcmp(name,"Tc"))
    {
        mat[0]=0.643137;
        mat[1]=0.666667;
        mat[2]=0.678431;
        mat[3]=1.367500;
    }
    else if(!strcmp(name,"Ru"))
    {
        mat[0]=0.643137;
        mat[1]=0.666667;
        mat[2]=0.678431;
        mat[3]=1.352900;
    }
    else if(!strcmp(name,"Rh"))
    {
        mat[0]=0.643137;
        mat[1]=0.666667;
        mat[2]=0.678431;
        mat[3]=1.345000;
    }
    else if(!strcmp(name,"Pd"))
    {
        mat[0]=0.643137;
        mat[1]=0.666667;
        mat[2]=0.678431;
        mat[3]=1.375500;
    }
    else if(!strcmp(name,"Ag"))
    {
        mat[0]=0.643137;
        mat[1]=0.666667;
        mat[2]=0.678431;
        mat[3]=1.444700;
    }
    else if(!strcmp(name,"Cd"))
    {
        mat[0]=0.643137;
        mat[1]=0.666667;
        mat[2]=0.678431;
        mat[3]=1.489400;
    }
    else if(!strcmp(name,"In"))
    {
        mat[0]=0.643137;
        mat[1]=0.666667;
        mat[2]=0.678431;
        mat[3]=1.666200;
    }
    else if(!strcmp(name,"Sn"))
    {
        mat[0]=0.643137;
        mat[1]=0.666667;
        mat[2]=0.678431;
        mat[3]=1.537500;
    }
    else if(!strcmp(name,"Sb"))
    {
        mat[0]=0.643137;
        mat[1]=0.666667;
        mat[2]=0.678431;
        mat[3]=1.400000;
    }
    else if(!strcmp(name,"Te"))
    {
        mat[0]=0.643137;
        mat[1]=0.666667;
        mat[2]=0.678431;
        mat[3]=1.360000;
    }
    else if(!strcmp(name,"I"))
    {
        mat[0]=0.500000;
        mat[1]=0.100000;
        mat[2]=0.500000;
        mat[3]=1.330000;
    }
    else if(!strcmp(name,"Xe"))
    {
        mat[0]=0.643137;
        mat[1]=0.666667;
        mat[2]=0.678431;
        mat[3]=2.200000;
    }
    else if(!strcmp(name,"Cs"))
    {
        mat[0]=0.643137;
        mat[1]=0.666667;
        mat[2]=0.678431;
        mat[3]=2.632500;
    }
    else if(!strcmp(name,"Ba"))
    {
        mat[0]=0.643137;
        mat[1]=0.666667;
        mat[2]=0.678431;
        mat[3]=2.170500;
    }
    else if(!strcmp(name,"La"))
    {
        mat[0]=0.643137;
        mat[1]=0.666667;
        mat[2]=0.678431;
        mat[3]=1.872500;
    }
    else if(!strcmp(name,"Ce"))
    {
        mat[0]=0.800000;
        mat[1]=0.800000;
        mat[2]=0.000000;
        mat[3]=1.824300;
    }
    else if(!strcmp(name,"Pr"))
    {
        mat[0]=0.643137;
        mat[1]=0.666667;
        mat[2]=0.678431;
        mat[3]=1.836200;
    }
    else if(!strcmp(name,"Nd"))
    {
        mat[0]=0.643137;
        mat[1]=0.666667;
        mat[2]=0.678431;
        mat[3]=1.829500;
    }
    else if(!strcmp(name,"Pm"))
    {
        mat[0]=0.643137;
        mat[1]=0.666667;
        mat[2]=0.678431;
        mat[3]=1.809000;
    }
    else if(!strcmp(name,"Sm"))
    {
        mat[0]=0.643137;
        mat[1]=0.666667;
        mat[2]=0.678431;
        mat[3]=1.804000;
    }
    else if(!strcmp(name,"Eu"))
    {
        mat[0]=0.643137;
        mat[1]=0.666667;
        mat[2]=0.678431;
        mat[3]=1.984000;
    }
    else if(!strcmp(name,"Gd"))
    {
        mat[0]=1.000000;
        mat[1]=0.843137;
        mat[2]=0.000000;
        mat[3]=1.818000;
    }
    else if(!strcmp(name,"Tb"))
    {
        mat[0]=0.643137;
        mat[1]=0.666667;
        mat[2]=0.678431;
        mat[3]=1.800500;
    }
    else if(!strcmp(name,"Dy"))
    {
        mat[0]=0.643137;
        mat[1]=0.666667;
        mat[2]=0.678431;
        mat[3]=1.795100;
    }
    else if(!strcmp(name,"Ho"))
    {
        mat[0]=0.643137;
        mat[1]=0.666667;
        mat[2]=0.678431;
        mat[3]=1.788600;
    }
    else if(!strcmp(name,"Er"))
    {
        mat[0]=0.643137;
        mat[1]=0.666667;
        mat[2]=0.678431;
        mat[3]=1.779400;
    }
    else if(!strcmp(name,"Tm"))
    {
        mat[0]=0.643137;
        mat[1]=0.666667;
        mat[2]=0.678431;
        mat[3]=1.768700;
    }
    else if(!strcmp(name,"Yb"))
    {
        mat[0]=0.643137;
        mat[1]=0.666667;
        mat[2]=0.678431;
        mat[3]=1.939600;
    }
    else if(!strcmp(name,"Lu"))
    {
        mat[0]=0.643137;
        mat[1]=0.666667;
        mat[2]=0.678431;
        mat[3]=1.751500;
    }
    else if(!strcmp(name,"Hf"))
    {
        mat[0]=0.643137;
        mat[1]=0.666667;
        mat[2]=0.678431;
        mat[3]=1.597300;
    }
    else if(!strcmp(name,"Ta"))
    {
        mat[0]=0.643137;
        mat[1]=0.666667;
        mat[2]=0.678431;
        mat[3]=1.428000;
    }
    else if(!strcmp(name,"W"))
    {
        mat[0]=0.643137;
        mat[1]=0.666667;
        mat[2]=0.678431;
        mat[3]=1.370500;
    }
    else if(!strcmp(name,"Re"))
    {
        mat[0]=0.643137;
        mat[1]=0.666667;
        mat[2]=0.678431;
        mat[3]=1.380000;
    }
    else if(!strcmp(name,"Os"))
    {
        mat[0]=0.643137;
        mat[1]=0.666667;
        mat[2]=0.678431;
        mat[3]=1.367600;
    }
    else if(!strcmp(name,"Ir"))
    {
        mat[0]=0.643137;
        mat[1]=0.666667;
        mat[2]=0.678431;
        mat[3]=1.357300;
    }
    else if(!strcmp(name,"Pt"))
    {
        mat[0]=0.643137;
        mat[1]=0.666667;
        mat[2]=0.678431;
        mat[3]=1.387300;
    }
    else if(!strcmp(name,"Au"))
    {
        mat[0]=0.900000;
        mat[1]=0.800000;
        mat[2]=0.000000;
        mat[3]=1.441900;
    }
    else if(!strcmp(name,"Hg"))
    {
        mat[0]=0.643137;
        mat[1]=0.666667;
        mat[2]=0.678431;
        mat[3]=1.502500;
    }
    else if(!strcmp(name,"Tl"))
    {
        mat[0]=0.643137;
        mat[1]=0.666667;
        mat[2]=0.678431;
        mat[3]=1.728300;
    }
    else if(!strcmp(name,"Pb"))
    {
        mat[0]=0.643137;
        mat[1]=0.666667;
        mat[2]=0.678431;
        mat[3]=1.750100;
    }
    else if(!strcmp(name,"Bi"))
    {
        mat[0]=0.643137;
        mat[1]=0.666667;
        mat[2]=0.678431;
        mat[3]=1.460000;
    }
    else if(!strcmp(name,"Po"))
    {
        mat[0]=0.643137;
        mat[1]=0.666667;
        mat[2]=0.678431;
        mat[3]=1.460000;
    }
    else if(!strcmp(name,"At"))
    {
        mat[0]=0.800000;
        mat[1]=0.200000;
        mat[2]=0.200000;
        mat[3]=4.350000;
    }
    else if(!strcmp(name,"Rn"))
    {
        mat[0]=0.643137;
        mat[1]=0.666667;
        mat[2]=0.678431;
        mat[3]=1.430000;
    }
    else if(!strcmp(name,"Fr"))
    {
        mat[0]=0.643137;
        mat[1]=0.666667;
        mat[2]=0.678431;
        mat[3]=2.500000;
    }
    else if(!strcmp(name,"Ra"))
    {
        mat[0]=0.643137;
        mat[1]=0.666667;
        mat[2]=0.678431;
        mat[3]=2.140000;
    }
    else if(!strcmp(name,"Ac"))
    {
        mat[0]=0.643137;
        mat[1]=0.666667;
        mat[2]=0.678431;
        mat[3]=1.877500;
    }
    else if(!strcmp(name,"Th"))
    {
        mat[0]=0.643137;
        mat[1]=0.666667;
        mat[2]=0.678431;
        mat[3]=1.797500;
    }
    else if(!strcmp(name,"Pa"))
    {
        mat[0]=0.643137;
        mat[1]=0.666667;
        mat[2]=0.678431;
        mat[3]=1.608600;
    }
    else if(!strcmp(name,"U"))
    {
        mat[0]=0.643137;
        mat[1]=0.666667;
        mat[2]=0.678431;
        mat[3]=1.568300;
    }
    else if(!strcmp(name,"Np"))
    {
        mat[0]=0.643137;
        mat[1]=0.666667;
        mat[2]=0.678431;
        mat[3]=1.000000;
    }
    else if(!strcmp(name,"Pu"))
    {
        mat[0]=0.643137;
        mat[1]=0.666667;
        mat[2]=0.678431;
        mat[3]=1.000000;
    }
    else if(!strcmp(name,"Am"))
    {
        mat[0]=0.643137;
        mat[1]=0.666667;
        mat[2]=0.678431;
        mat[3]=1.000000;
    }
    else if(!strcmp(name,"Cm"))
    {
        mat[0]=0.643137;
        mat[1]=0.666667;
        mat[2]=0.678431;
        mat[3]=1.000000;
    }
    else if(!strcmp(name,"Bk"))
    {
        mat[0]=0.643137;
        mat[1]=0.666667;
        mat[2]=0.678431;
        mat[3]=1.000000;
    }
    else if(!strcmp(name,"Cf"))
    {
        mat[0]=0.100000;
        mat[1]=0.700000;
        mat[2]=0.300000;
        mat[3]=1.460000;
    }
    else if(!strcmp(name,"Es"))
    {
        mat[0]=0.100000;
        mat[1]=0.300000;
        mat[2]=0.700000;
        mat[3]=1.752000;
    }
    else if(!strcmp(name,"Fm"))
    {
        mat[0]=0.643137;
        mat[1]=0.666667;
        mat[2]=0.678431;
        mat[3]=1.000000;
    }
    else if(!strcmp(name,"Md"))
    {
        mat[0]=0.643137;
        mat[1]=0.666667;
        mat[2]=0.678431;
        mat[3]=1.000000;
    }
    else if(!strcmp(name,"No"))
    {
        mat[0]=0.643137;
        mat[1]=0.666667;
        mat[2]=0.678431;
        mat[3]=1.000000;
    }
    else if(!strcmp(name,"Lr"))
    {
        mat[0]=0.643137;
        mat[1]=0.666667;
        mat[2]=0.678431;
        mat[3]=1.000000;
    }
    else if(!strcmp(name,"Rf"))
    {
        mat[0]=0.643137;
        mat[1]=0.666667;
        mat[2]=0.678431;
        mat[3]=1.000000;
    }
    else if(!strcmp(name,"Db"))
    {
        mat[0]=0.900000;
        mat[1]=0.800000;
        mat[2]=0.000000;
        mat[3]=1.460000;
    }
    else if(!strcmp(name,"Sg"))
    {
        mat[0]=0.643137;
        mat[1]=0.666667;
        mat[2]=0.678431;
        mat[3]=1.000000;
    }
    else if(!strcmp(name,"Bh"))
    {
        mat[0]=1.000000;
        mat[1]=1.000000;
        mat[2]=0.000000;
        mat[3]=1.900000;
    }
    else if(!strcmp(name,"Hs"))
    {
        mat[0]=0.643137;
        mat[1]=0.666667;
        mat[2]=0.678431;
        mat[3]=1.000000;
    }
    else if(!strcmp(name,"Mt"))
    {
        mat[0]=0.643137;
        mat[1]=0.666667;
        mat[2]=0.678431;
        mat[3]=1.000000;
    }
    else
    {
        mat[0]=0.500000;
        mat[1]=0.500000;
        mat[2]=0.500000;
        mat[3]=0.435000;
    }
}
