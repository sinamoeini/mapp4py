#include "api.h"
#include "import_eam.h"
//#include <structmember.h>
#include "memory.h"
using namespace MAPP_NS;
/*--------------------------------------------
 
 --------------------------------------------*/
type0 ImportEAM::interpolate(type0* arr,size_t n,type0 p,size_t k)
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
void ImportEAM::interpolate(size_t n,type0 delta,type0(*spline)[1])
{}
/*--------------------------------------------

 --------------------------------------------*/
void ImportEAM::interpolate(size_t n,type0 delta,type0(*spline)[4])
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
void ImportEAM::interpolate(size_t n,type0 delta,type0(*spline)[5])
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
void ImportEAM::__interpolate(size_t n,type0 delta,type0(*spline)[7])
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
    for(size_t i=0;i<n;i++)
    {
        spline[i][4]=spline[i][1]/delta;
        spline[i][5]=2.0*spline[i][2]/delta;
        spline[i][6]=3.0*spline[i][3]/delta;
    }
}
/*--------------------------------------------
 
 --------------------------------------------*/
void ImportEAM::read_double(type0* buff,size_t n,FileReader& fr,char*& line,size_t& line_cpcty,char**& args,size_t& args_cpcty)
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
        
        for(size_t j=0;j<nargs;j++)
            buff[i++]=atof(args[j]);
    }
}
/*--------------------------------------------
 
 --------------------------------------------*/
void ImportEAM::skip(size_t n,FileReader& fr,char*& line,size_t& line_cpcty)
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

/*--------------------------------------------
 
 --------------------------------------------*/
void ImportEAM::ml_read_eam(PyMethodDef& method_0,PyMethodDef& method_1,PyMethodDef& method_2)
{
    method_0.ml_flags=METH_VARARGS | METH_KEYWORDS;
    method_0.ml_name="import_funcfl";
    method_0.ml_meth=(PyCFunction)(PyCFunctionWithKeywords)(
    [](PyObject* self,PyObject* args,PyObject* kwds)->PyObject*
    {
        FuncAPI<std::string*> f("import_funcfl",{"funcfl_files"});
        if(f(args,kwds)) return NULL;
        
        size_t nelems=f.v<0>().size;
        size_t nr,nrho;
        type0 dr,drho;
        type0** r_c;
        type0(** F)[1]=NULL;
        type0(*** r_phi)[1]=NULL;
        type0(*** rho)[1]=NULL;
        try
        {
            ImportEAM::funcfl(nelems,f.val<0>(),dr,drho,nr,nrho,r_phi,rho,F,r_c);
        }
        catch(char* err_msg)
        {
            PyErr_SetString(PyExc_TypeError,err_msg);
            delete [] err_msg;
            return NULL;
        }
        
        type0* data=NULL;
        Memory::alloc(data,nelems*nrho+nelems*nr+nelems*(nelems+1)*nr/2);
        type0* __data=data;
        for(size_t i=0;i<nelems;i++)
            for(size_t in=0;in<nrho;in++,++__data)
                *__data=F[i][in][0];
        
        
        
        
        for(size_t i=0;i<nelems;i++)
            for(size_t in=0;in<nr;in++,++__data)
                *__data=rho[i][0][in][0];
        
        for(size_t i=0;i<nelems;i++)
            for(size_t j=i;j<nelems;j++)
            {
                *__data=0.0;
                ++__data;
                for(size_t in=1;in<nr;in++,++__data)
                    *__data=r_phi[i][j][in][0]/(static_cast<type0>(in)*dr);
            }
        Memory::dealloc(F);
        Memory::dealloc(rho);
        Memory::dealloc(r_phi);
        
        
        
        size_t* nrhop=&nrho;
        size_t* nrp=&nr;
        PyObject* op;
        
        __data=data;
        PyObject* F_obj=PyList_New(nelems);
        for(size_t i=0;i<nelems;i++)
        {
            op=var<type0*>::build(__data,&nrhop);
            PyList_SET_ITEM(F_obj,i,op);
            __data+=nrho;
        }
        
        
        PyObject* rho_obj=PyList_New(nelems);
        for(size_t i=0;i<nelems;i++)
        {
            op=var<type0*>::build(__data,&nrp);
            PyList_SET_ITEM(rho_obj,i,op);
            __data+=nr;
        }
        PyObject* phi_obj=PyList_New(nelems);
        for(size_t i=0;i<nelems;i++)
        {
            PyList_SET_ITEM(phi_obj,i,PyList_New(nelems));
            for(size_t j=0;j<nelems;j++)
            {
                if(j<i)
                {
                    op=PyList_GET_ITEM(PyList_GET_ITEM(phi_obj,j),i);
                    Py_INCREF(op);
                }
                else
                {
                    op=var<type0*>::build(__data,&nrp);
                    __data+=nr;
                }
                
                PyList_SET_ITEM(PyList_GET_ITEM(phi_obj,i),j,op);
            }
            
            
        }
        
        
        
        Memory::dealloc(data);
        
        PyObject* ans=PyList_New(5);
        PyList_SET_ITEM(ans,0,F_obj);
        PyList_SET_ITEM(ans,1,rho_obj);
        PyList_SET_ITEM(ans,2,phi_obj);
        PyList_SET_ITEM(ans,3,var<type0>::build(drho));
        PyList_SET_ITEM(ans,4,var<type0>::build(dr));
        return ans;
    });
    method_0.ml_doc=(char*)R"---(
    import_eam_funcfl(funcfl_files)
   
    Importing FuncFL file/s
    
    Imports FuncFL file/s
    
    Parameters
    ----------
    funcfl_files : string[nelems]
        list of relative paths to DYNAMO files with FuncFL format
    
    Returns
    -------
    None
   
    Notes
    -----
    This is tabulated form of Embedded Atom Method (EAM) potential
    
    
    Examples
    --------
    Ni
    
    ::
    
        >>> from mapp import md
        >>> sim=md.cfg("configs/Ni.cfg")
        >>> sim.ff_eam_funcfl("potentials/Ni_u3.eam")
    
    

    )---";
    
    method_1.ml_flags=METH_VARARGS | METH_KEYWORDS;
    method_1.ml_name="import_setfl";
    method_1.ml_meth=(PyCFunction)(PyCFunctionWithKeywords)(
    [](PyObject* self,PyObject* args,PyObject* kwds)->PyObject*
    {

        FuncAPI<std::string> f("import_setfl",{"setfl_file"});
        if(f(args,kwds)) return NULL;
        
        size_t nelems=0;
        size_t nr,nrho;
        type0 dr,drho;
        type0** r_c;
        type0(** F)[1]=NULL;
        type0(*** r_phi)[1]=NULL;
        type0(*** rho)[1]=NULL;
        try
        {
            ImportEAM::setfl(nelems,NULL,f.val<0>(),dr,drho,nr,nrho,r_phi,rho,F,r_c);
        }
        catch(char* err_msg)
        {
            PyErr_SetString(PyExc_TypeError,err_msg);
            delete [] err_msg;
            return NULL;
        }
        
        type0* data=NULL;
        Memory::alloc(data,nelems*nrho+nelems*nr+nelems*(nelems+1)*nr/2);
        type0* __data=data;
        for(size_t i=0;i<nelems;i++)
            for(size_t in=0;in<nrho;in++,++__data)
                *__data=F[i][in][0];
        
        
        
        
        for(size_t i=0;i<nelems;i++)
            for(size_t in=0;in<nr;in++,++__data)
                *__data=rho[i][0][in][0];
        
        for(size_t i=0;i<nelems;i++)
            for(size_t j=i;j<nelems;j++)
            {
                *__data=0.0;
                ++__data;
                for(size_t in=1;in<nr;in++,++__data)
                    *__data=r_phi[i][j][in][0]/(static_cast<type0>(in)*dr);
            }
        Memory::dealloc(F);
        Memory::dealloc(rho);
        Memory::dealloc(r_phi);
        
        
        
        size_t* nrhop=&nrho;
        size_t* nrp=&nr;
        PyObject* op;
        
        __data=data;
        PyObject* F_obj=PyList_New(nelems);
        for(size_t i=0;i<nelems;i++)
        {
            op=var<type0*>::build(__data,&nrhop);
            PyList_SET_ITEM(F_obj,i,op);
            __data+=nrho;
        }
        

        PyObject* rho_obj=PyList_New(nelems);
        for(size_t i=0;i<nelems;i++)
        {
            op=var<type0*>::build(__data,&nrp);
            PyList_SET_ITEM(rho_obj,i,op);
            __data+=nr;
        }
        PyObject* phi_obj=PyList_New(nelems);
        for(size_t i=0;i<nelems;i++)
        {
            PyList_SET_ITEM(phi_obj,i,PyList_New(nelems));
            for(size_t j=0;j<nelems;j++)
            {
                if(j<i)
                {
                    op=PyList_GET_ITEM(PyList_GET_ITEM(phi_obj,j),i);
                    Py_INCREF(op);
                }
                else
                {
                    op=var<type0*>::build(__data,&nrp);
                    __data+=nr;
                }
                
                PyList_SET_ITEM(PyList_GET_ITEM(phi_obj,i),j,op);
            }
            
            
        }
        
        
        
        Memory::dealloc(data);
        
        PyObject* ans=PyList_New(5);
        PyList_SET_ITEM(ans,0,F_obj);
        PyList_SET_ITEM(ans,1,rho_obj);
        PyList_SET_ITEM(ans,2,phi_obj);
        PyList_SET_ITEM(ans,3,var<type0>::build(drho));
        PyList_SET_ITEM(ans,4,var<type0>::build(dr));
        return ans;
    });
    method_1.ml_doc=(char*)R"---(
    ff_eam_setfl(setfl_file)
   
    Tabulated EAM force field given by a single SetFL file
    
    Assigns EAM force field to system
    
    Parameters
    ----------
    setfl_file : string
        relative path to DYNAMO file with SetFL format
    
    Returns
    -------
    None
   
    Notes
    -----
    This is tabulated form of Embedded Atom Method (EAM) potential
    
    
    Examples
    --------
    Cu
    
    ::
    
        >>> from mapp import md
        >>> sim=md.cfg("configs/Cu.cfg")
        >>> sim.ff_eam_setfl("potentials/Cu_mishin.eam.alloy")

    
    
    )---";
    
    
    method_2.ml_flags=METH_VARARGS | METH_KEYWORDS;
    method_2.ml_name="import_fs";
    method_2.ml_meth=(PyCFunction)(PyCFunctionWithKeywords)(
    [](PyObject* self,PyObject* args,PyObject* kwds)->PyObject*
    {

        FuncAPI<std::string> f("import_fs",{"fs_file"});
        if(f(args,kwds)) return NULL;
        
        size_t nelems=0;
        size_t nr,nrho;
        type0 dr,drho;
        type0** r_c;
        type0(** F)[1]=NULL;
        type0(*** r_phi)[1]=NULL;
        type0(*** rho)[1]=NULL;
        try
        {
            ImportEAM::fs(nelems,NULL,f.val<0>(),dr,drho,nr,nrho,r_phi,rho,F,r_c);
        }
        catch(char* err_msg)
        {
            PyErr_SetString(PyExc_TypeError,err_msg);
            delete [] err_msg;
            return NULL;
        }
        
        type0* data=NULL;
        Memory::alloc(data,nelems*nrho+nelems*nelems*nr+nelems*(nelems+1)*nr/2);
        type0* __data=data;
        for(size_t i=0;i<nelems;i++)
            for(size_t in=0;in<nrho;in++,++__data)
                *__data=F[i][in][0];
        
        
        
        
        for(size_t i=0;i<nelems;i++)
            for(size_t j=0;j<nelems;j++)
                for(size_t in=0;in<nr;in++,++__data)
                    *__data=rho[i][j][in][0];
        
        for(size_t i=0;i<nelems;i++)
            for(size_t j=i;j<nelems;j++)
            {
                *__data=0.0;
                ++__data;
                for(size_t in=1;in<nr;in++,++__data)
                    *__data=r_phi[i][j][in][0]/(static_cast<type0>(in)*dr);
            }
        Memory::dealloc(F);
        Memory::dealloc(rho);
        Memory::dealloc(r_phi);
        
        
        
        size_t* nrhop=&nrho;
        size_t* nrp=&nr;
        PyObject* op;
        
        __data=data;
        PyObject* F_obj=PyList_New(nelems);
        for(size_t i=0;i<nelems;i++)
        {
            op=var<type0*>::build(__data,&nrhop);
            PyList_SET_ITEM(F_obj,i,op);
            __data+=nrho;
        }
        
        
        PyObject* rho_obj=PyList_New(nelems);
        for(size_t i=0;i<nelems;i++)
        {
            PyList_SET_ITEM(rho_obj,i,PyList_New(nelems));
            for(size_t j=0;j<nelems;j++)
            {
                op=var<type0*>::build(__data,&nrp);
                PyList_SET_ITEM(PyList_GET_ITEM(rho_obj,i),j,op);
                __data+=nr;
            }
            
        }
        
        
        PyObject* phi_obj=PyList_New(nelems);
        for(size_t i=0;i<nelems;i++)
        {
            PyList_SET_ITEM(phi_obj,i,PyList_New(nelems));
            for(size_t j=0;j<nelems;j++)
            {
                if(j<i)
                {
                    op=PyList_GET_ITEM(PyList_GET_ITEM(phi_obj,j),i);
                    Py_INCREF(op);
                }
                else
                {
                    op=var<type0*>::build(__data,&nrp);
                    __data+=nr;
                }
                
                PyList_SET_ITEM(PyList_GET_ITEM(phi_obj,i),j,op);
            }
            
            
        }
        
        
        
        Memory::dealloc(data);
        
        PyObject* ans=PyList_New(5);
        PyList_SET_ITEM(ans,0,F_obj);
        PyList_SET_ITEM(ans,1,rho_obj);
        PyList_SET_ITEM(ans,2,phi_obj);
        PyList_SET_ITEM(ans,3,var<type0>::build(drho));
        PyList_SET_ITEM(ans,4,var<type0>::build(dr));
        return ans;
    });
    method_2.ml_doc=(char*)R"---(
    ff_eam_fs(fs_file)
   
    Tabulated Finnis-Sinclair EAM
    
    Assigns Finnis-Sinclair EAM force field to system. For explanation of the parameter see the Notes section.
    
    Parameters
    ----------
    fs_file : string
        relative path to DYNAMO file with fs format
    
    Returns
    -------
    None
   
    Notes
    -----
    This is tabulated form of Finnis-Sinclair Embedded Atom Method (EAM) potential
    
    
    Examples
    --------
    Iron Hydrogrn mixture
    ::
    
        >>> from mapp import md
        >>> sim=md.cfg("configs/FeH.cfg")
        >>> sim.ff_eam_fs("potentials/FeH.eam.fs")
    
    

    )---";
    
}
 

