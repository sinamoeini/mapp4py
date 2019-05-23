#include "ff_eam_potfit.h"
using namespace MAPP_NS;
/*--------------------------------------------
 python constructor
 --------------------------------------------*/
ForceFieldEAMPotFit<ForceFieldEAMPotFitAckOgata::nelems>* ForceFieldEAMPotFitAckOgata::get_new_ff(AtomsMD* atoms,PyObject* args,PyObject* kwds)
{
    std::string __names[nelems]={"Fe","H"};
    for(size_t i=0;i<atoms->elements.nelems;i++)
    {
        bool found=false;
        for(size_t j=0;j<nelems && !found;j++)
            if(strcmp(__names[j].c_str(),atoms->elements.names[i].c_str())==0) found=true;
        if(!found)
        {
            PyErr_Format(PyExc_TypeError,"this forcefield does not support elemnt %s",atoms->elements.names[i].c_str());
            return NULL;
        }
    }
    
    
    FuncAPI<OB<PyListObject,PyList_Type>,OB<PyListObject,PyList_Type>,OB<PyListObject,PyList_Type>,
    OB<PyListObject,PyList_Type>,OB<PyListObject,PyList_Type>,
    OB<PyListObject,PyList_Type>,OB<PyListObject,PyList_Type>>
    f("ff_eam_fit_ao",{"A_phi_FeFe","A_phi_FeH","A_phi_HH","A_rho_Fe","A_rho_H","A_F_Fe","A_F_H"});
    if(f(args,kwds)) return NULL;
    
    
    
    PotFitPairFunc* __phi_ptr[nelems][nelems];
    PotFitPairFunc* __rho_ptr[nelems][nelems];
    PotFitEmbFunc* __F_ptr[nelems];
    __phi_ptr[0][0]=new PotFitPhiSpl("phi_FeFe");
    __phi_ptr[1][0]=new PotFitPhiSpl("phi_FeH");
    __phi_ptr[0][1]=__phi_ptr[1][0];
    __phi_ptr[1][1]=new PotFitPhiSpl("phi_HH");
    
    __rho_ptr[0][0]=new PotFitRhoSpl("rho_Fe");
    __rho_ptr[1][0]=new PotFitRhoAng("rho_H");
    __rho_ptr[0][1]=__rho_ptr[0][0];
    __rho_ptr[1][1]=__rho_ptr[1][0];
    
    __F_ptr[0]=new PotFitEmbAck("F_Fe");
    __F_ptr[1]=new PotFitEmbAng("F_H");
    
    type0* data=NULL;
    size_t tot_nvars=0;
    int set=0;
    
    if(set==0) set=__phi_ptr[0][0]->set_init(f.val<0>(),data,tot_nvars);
    if(set==0) set=__phi_ptr[0][1]->set_init(f.val<1>(),data,tot_nvars);
    if(set==0) set=__phi_ptr[1][1]->set_init(f.val<2>(),data,tot_nvars);
    
    if(set==0) set=__rho_ptr[0][0]->set_init(f.val<3>(),data,tot_nvars);
    if(set==0) set=__rho_ptr[1][0]->set_init(f.val<4>(),data,tot_nvars);
    
    if(set==0) set=__F_ptr[0]->set_init(f.val<5>(),data,tot_nvars);
    if(set==0) set=__F_ptr[1]->set_init(f.val<6>(),data,tot_nvars);
    
    if(set)
    {
        delete __F_ptr[1];
        delete __F_ptr[0];
        delete __phi_ptr[1][1];
        delete __phi_ptr[0][1];
        delete __phi_ptr[0][0];
        delete __rho_ptr[1][0];
        delete __rho_ptr[0][0];
        return NULL;
    }
    
    type0* __data=data;
    
    __phi_ptr[0][0]->vars=__data;
    __data+=__phi_ptr[0][0]->nvars;
    __phi_ptr[0][1]->vars=__data;
    __data+=__phi_ptr[0][1]->nvars;
    __phi_ptr[1][1]->vars=__data;
    __data+=__phi_ptr[1][1]->nvars;
    
    
    __rho_ptr[0][0]->vars=__data;
    __data+=__rho_ptr[0][0]->nvars;
    __rho_ptr[1][0]->vars=__data;
    __data+=__rho_ptr[1][0]->nvars;
    
    
    __F_ptr[0]->vars=__data;
    __data+=__F_ptr[0]->nvars;
    __F_ptr[1]->vars=__data;
    __data+=__F_ptr[1]->nvars;
    
    return new ForceFieldEAMPotFit<nelems>(atoms,__names,data,tot_nvars,__phi_ptr,__rho_ptr,__F_ptr);
}
/*--------------------------------------------
 python constructor
 --------------------------------------------*/
void ForceFieldEAMPotFitAckOgata::ml_new(PyMethodDef& tp_methods)
{
    tp_methods.ml_flags=METH_VARARGS | METH_KEYWORDS;
    tp_methods.ml_name="ff_eam_ao";
    tp_methods.ml_meth=(PyCFunction)(PyCFunctionWithKeywords)(
    [](PyObject* self,PyObject* args,PyObject* kwds)->PyObject*
    {
        AtomsMD::Object* __self=reinterpret_cast<AtomsMD::Object*>(self);
        ForceFieldMD* new_ff=get_new_ff(__self->atoms,args,kwds);
        if(new_ff==NULL) return NULL;
        delete __self->ff;
        __self->ff=new_ff;
        Py_RETURN_NONE;
    });

    tp_methods.ml_doc="";
}

/*--------------------------------------------
 python constructor
 --------------------------------------------*/
ForceFieldEAMPotFit<ForceFieldEAMPotFitAckJP::nelems>* ForceFieldEAMPotFitAckJP::get_new_ff(AtomsMD* atoms,PyObject* args,PyObject* kwds)
{
    std::string __names[nelems]={"Fe","H"};
    for(size_t i=0;i<atoms->elements.nelems;i++)
    {
        bool found=false;
        for(size_t j=0;j<nelems && !found;j++)
            if(strcmp(__names[j].c_str(),atoms->elements.names[i].c_str())==0) found=true;
        if(!found)
        {
            PyErr_Format(PyExc_TypeError,"this forcefield does not support elemnt %s",atoms->elements.names[i].c_str());
            return NULL;
        }
    }
    
    
    FuncAPI<OB<PyListObject,PyList_Type>,OB<PyListObject,PyList_Type>,OB<PyListObject,PyList_Type>,
    OB<PyListObject,PyList_Type>,OB<PyListObject,PyList_Type>,
    OB<PyListObject,PyList_Type>,OB<PyListObject,PyList_Type>>
    f("ff_eam_fit_ajp",{"A_phi_FeFe","A_phi_FeH","A_phi_HH","A_rho_Fe","A_rho_H","A_F_Fe","A_F_H"});
    if(f(args,kwds)) return NULL;
    
    
    
    PotFitPairFunc* __phi_ptr[nelems][nelems];
    PotFitPairFunc* __rho_ptr[nelems][nelems];
    PotFitEmbFunc* __F_ptr[nelems];
    __phi_ptr[0][0]=new PotFitPhiSpl("phi_FeFe");
    __phi_ptr[1][0]=new PotFitPhiSpl("phi_FeH");
    __phi_ptr[0][1]=__phi_ptr[1][0];
    __phi_ptr[1][1]=new PotFitPhiSpl("phi_HH");
    
    __rho_ptr[0][0]=new PotFitRhoSpl("rho_Fe");
    __rho_ptr[1][0]=new PotFitRhoJP("rho_H");
    __rho_ptr[0][1]=__rho_ptr[0][0];
    __rho_ptr[1][1]=__rho_ptr[1][0];
    
    __F_ptr[0]=new PotFitEmbAck("F_Fe");
    __F_ptr[1]=new PotFitEmbJP("F_H");
    
    type0* data=NULL;
    size_t tot_nvars=0;
    int set=0;
    
    if(set==0) set=__phi_ptr[0][0]->set_init(f.val<0>(),data,tot_nvars);
    if(set==0) set=__phi_ptr[0][1]->set_init(f.val<1>(),data,tot_nvars);
    if(set==0) set=__phi_ptr[1][1]->set_init(f.val<2>(),data,tot_nvars);
    
    if(set==0) set=__rho_ptr[0][0]->set_init(f.val<3>(),data,tot_nvars);
    if(set==0) set=__rho_ptr[1][0]->set_init(f.val<4>(),data,tot_nvars);
    
    if(set==0) set=__F_ptr[0]->set_init(f.val<5>(),data,tot_nvars);
    if(set==0) set=__F_ptr[1]->set_init(f.val<6>(),data,tot_nvars);
    
    if(set)
    {
        delete __F_ptr[1];
        delete __F_ptr[0];
        delete __phi_ptr[1][1];
        delete __phi_ptr[0][1];
        delete __phi_ptr[0][0];
        delete __rho_ptr[1][0];
        delete __rho_ptr[0][0];
        return NULL;
    }
    
    type0* __data=data;
    
    __phi_ptr[0][0]->vars=__data;
    __data+=__phi_ptr[0][0]->nvars;
    __phi_ptr[0][1]->vars=__data;
    __data+=__phi_ptr[0][1]->nvars;
    __phi_ptr[1][1]->vars=__data;
    __data+=__phi_ptr[1][1]->nvars;
    
    
    __rho_ptr[0][0]->vars=__data;
    __data+=__rho_ptr[0][0]->nvars;
    __rho_ptr[1][0]->vars=__data;
    __data+=__rho_ptr[1][0]->nvars;
    
    
    __F_ptr[0]->vars=__data;
    __data+=__F_ptr[0]->nvars;
    __F_ptr[1]->vars=__data;
    __data+=__F_ptr[1]->nvars;
    
    return new ForceFieldEAMPotFit<nelems>(atoms,__names,data,tot_nvars,__phi_ptr,__rho_ptr,__F_ptr);
}
/*--------------------------------------------
 python constructor
 --------------------------------------------*/
void ForceFieldEAMPotFitAckJP::ml_new(PyMethodDef& tp_methods)
{
    tp_methods.ml_flags=METH_VARARGS | METH_KEYWORDS;
    tp_methods.ml_name="ff_eam_ajp";
    tp_methods.ml_meth=(PyCFunction)(PyCFunctionWithKeywords)(
    [](PyObject* self,PyObject* args,PyObject* kwds)->PyObject*
    {
        AtomsMD::Object* __self=reinterpret_cast<AtomsMD::Object*>(self);
        ForceFieldMD* new_ff=get_new_ff(__self->atoms,args,kwds);
        if(new_ff==NULL) return NULL;
        delete __self->ff;
        __self->ff=new_ff;
        Py_RETURN_NONE;
    });

    tp_methods.ml_doc="";
}
