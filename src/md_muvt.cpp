#include "md_muvt.h"
#include "thermo_dynamics.h"
#include "dynamic_md.h"
#include "elements.h"
#include "ff_md.h"
#include "pgcmc.h"
#include "xmath.h"
/*--------------------------------------------
 
 --------------------------------------------*/
MDMuVT::MDMuVT(type0 __mu,type0 __T,type0 __dt,std::string __gas_elem_name,int __seed):
MDNVT(__T,__dt),
seed(__seed),
mu(__mu),
nevery(1000),
nattempts(1000),
gas_elem_name(__gas_elem_name)
{
}
/*--------------------------------------------
 
 --------------------------------------------*/
MDMuVT::~MDMuVT()
{
}
/*--------------------------------------------
 pre run check it throw excepctions
 --------------------------------------------*/
void MDMuVT::pre_run_chk(AtomsMD* atoms,ForceFieldMD* ff)
{
    //check if configuration is loaded
    if(!atoms)
        throw std::string("cannot start md without initial conditions");
    
    //check if force field is loaded
    if(!ff)
        throw std::string("cannot start md without governing equations (force field)");
    
    if(std::isnan(atoms->kB))
        throw std::string("boltzmann constant should be set prior to MD");
    
    if(std::isnan(atoms->hP))
        throw std::string("planck constant should be set prior to GCMC");
    
    elem_type ielem;
    try
    {
        ielem=atoms->elements.find(gas_elem_name.c_str());
    }
    catch(int)
    {
        throw "atom "+gas_elem_name+" is not assigned to configuration";
    }
    
    
    if(!atoms->x_dof->is_empty())
    {
        bool* dof=atoms->x_dof->begin();
        elem_type* elem=atoms->elem->begin();
        int err_lcl=0;
        for(int i=0;i<atoms->natms_lcl;i++,dof+=__dim__)
        {
            if(elem[i]==ielem)
                Algebra::Do<__dim__>::func([&dof,&err_lcl](int i){ if(!dof[i]) err_lcl=1;});
        }
        
        int err=0;
        MPI_Allreduce(&err_lcl,&err,1,MPI_INT,MPI_MAX,atoms->world);
        if(err)
            throw std::string("cannot fix any degrees of freedom of any of the gas atoms (")+gas_elem_name+std::string(")");
    }
}
/*--------------------------------------------
 
 --------------------------------------------*/
void MDMuVT::pre_init()
{
    MDNVT::pre_init();
}
/*--------------------------------------------
 
 --------------------------------------------*/
void MDMuVT::init()
{
    pre_init();
    
    dynamic=new DynamicMD(atoms,ff,false,{},{atoms->x_d,atoms->x_dof},{});
    dynamic->init();
    
    if(xprt)
    {
        try
        {
            xprt->atoms=atoms;
            xprt->init();
        }
        catch(std::string& err_msg)
        {
            fin();
            throw err_msg;
        }
    }
}
/*--------------------------------------------
 
 --------------------------------------------*/
void MDMuVT::fin()
{
    if(xprt)
    {
        xprt->fin();
        xprt->atoms=NULL;
    }
    
    dynamic->fin();
    delete dynamic;
    dynamic=NULL;
}
/*--------------------------------------------
 
 --------------------------------------------*/
void MDMuVT::run(int nsteps)
{
    
    int step=atoms->step;
    
    PGCMC gcmc(atoms,ff,dynamic,1,atoms->elements.find(gas_elem_name.c_str()),mu,T,seed);
    gcmc.init();
    type0 gas_frac=static_cast<type0>(gcmc.ngas)/static_cast<type0>(atoms->natms);
#ifdef GCMCDEBUG
    type0 delta_u=0.0;
    FILE* fp_debug=NULL;
    if(atoms->comm_rank==0)
    {
        fp_debug=fopen("gcmc_debug","w");
        fprintf(fp_debug,"step\tacc_du\tgcmc_du\trel_error\terror\n");
    }
#endif
    
    ff->force_calc();
    
    int nevery_xprt=xprt==NULL ? 0:xprt->nevery;
    if(nevery_xprt) xprt->write(step);
    
    ThermoDynamics thermo(6,"T",T_part,"PE",atoms->pe,
    "S[0][0]",S_part[0][0],
    "S[1][1]",S_part[1][1],
    "S[2][2]",S_part[2][2],
    "S[1][2]",S_part[2][1],
    "S[2][0]",S_part[2][0],
    "S[0][1]",S_part[1][0],
    "GAS_FRAC",gas_frac);
    
    if(ntally) thermo.init();
    Algebra::DoLT<__dim__>::func([this](const int i,const int j)
    {
        S_part[i][j]=atoms->S_pe[i][j]-mvv[i+j*__dim__-j*(j+1)/2]/atoms->vol;
    });
    
    
    if(ntally) thermo.print(step);
    
    type0 fac,fac_x_d=1.0;
    if(dof_empty)
    {
        for(int istep=0;istep<nsteps;istep++)
        {
            // particle thermostat
            fac_x_d*=fac=thermo_part(T_part/T,ndof_part);
            fac*=fac;
            Algebra::Do<__nvoigt__>::func([this,&fac](int i){mvv[i]*=fac;});
            T_part*=fac;
            
            update_x_d__x(fac_x_d);
            
            if((istep+1)%nevery)
                dynamic->update<true>();
            else
            {
    #ifdef GCMCDEBUG
                dynamic->update<true>();

                ff->force_calc();
                delta_u=atoms->pe;
    #endif
                gcmc.xchng(false,nattempts);
                gas_frac=static_cast<type0>(gcmc.ngas)/static_cast<type0>(atoms->natms);
                ndof_part+=static_cast<type0>(gcmc.dof_diff);
                Algebra::Do<__nvoigt__>::func([this,&gcmc](const int i){mvv[i]+=gcmc.mvv[i];});
                T_part=Algebra::Tr_DyadicV<__dim__>(mvv)/(ndof_part*kB);
            }

            ff->force_calc();
            
    #ifdef GCMCDEBUG
            if((istep+1)%nevery==0)
            {
                delta_u-=atoms->pe;
                if(atoms->comm_rank==0)
                    fprintf(fp_debug,"%d\t%e\t%e\t%e\t%e\n",istep,-delta_u,gcmc.tot_delta_u,fabs((delta_u+gcmc.tot_delta_u)/delta_u),fabs(delta_u+gcmc.tot_delta_u));
            }
    #endif
            
            
            update_x_d();
            
            // particle thermostat
            fac_x_d=fac=thermo_part(T_part/T,ndof_part);
            fac*=fac;
            Algebra::Do<__nvoigt__>::func([this,&fac](int i){mvv[i]*=fac;});
            T_part*=fac;
            
            
            Algebra::DoLT<__dim__>::func([this](const int i,const int j)
            {
                S_part[i][j]=atoms->S_pe[i][j]-mvv[i+j*__dim__-j*(j+1)/2]/atoms->vol;
            });
            
            if(ntally && (istep+1)%ntally==0) thermo.print(step+istep+1);
            if(nevery_xprt && (istep+1)%nevery_xprt==0) xprt->write(step+istep+1);
        }
    }
    else
    {
        for(int istep=0;istep<nsteps;istep++)
        {
            // particle thermostat
            fac_x_d*=fac=thermo_part(T_part/T,ndof_part);
            fac*=fac;
            Algebra::Do<__nvoigt__>::func([this,&fac](int i){mvv[i]*=fac;});
            T_part*=fac;
            
            update_x_d__x_w_dof(fac_x_d);
            
            if((istep+1)%nevery)
                dynamic->update<true>();
            else
            {

                gcmc.xchng(false,nattempts);
                gas_frac=static_cast<type0>(gcmc.ngas)/static_cast<type0>(atoms->natms);
                ndof_part+=static_cast<type0>(gcmc.dof_diff);
                Algebra::Do<__nvoigt__>::func([this,&gcmc](const int i){mvv[i]+=gcmc.mvv[i];});
                T_part=Algebra::Tr_DyadicV<__dim__>(mvv)/(ndof_part*kB);
            }

            ff->force_calc();
            
            update_x_d_w_dof();
            
            // particle thermostat
            fac_x_d=fac=thermo_part(T_part/T,ndof_part);
            fac*=fac;
            Algebra::Do<__nvoigt__>::func([this,&fac](int i){mvv[i]*=fac;});
            T_part*=fac;
            
            
            Algebra::DoLT<__dim__>::func([this](const int i,const int j)
            {
                S_part[i][j]=atoms->S_pe[i][j]-mvv[i+j*__dim__-j*(j+1)/2]/atoms->vol;
            });
            
            if(ntally && (istep+1)%ntally==0) thermo.print(step+istep+1);
            if(nevery_xprt && (istep+1)%nevery_xprt==0) xprt->write(step+istep+1);
        }
    }
    if(ntally && nsteps%ntally) thermo.print(step+nsteps);
    if(nevery_xprt && nsteps%nevery_xprt) xprt->write(step+nsteps);
    
    if(dof_empty)
        update_x_d_final(fac_x_d);
    else
        update_x_d_final_w_dof(fac_x_d);
    
    if(ntally) thermo.fin();
    
    gcmc.fin();
    
    atoms->step+=nsteps;
    
#ifdef GCMCDEBUG
    if(atoms->comm_rank==0) fclose(fp_debug);
#endif
}
/*------------------------------------------------------------------------------------------------------------------------------------
 
 ------------------------------------------------------------------------------------------------------------------------------------*/
PyObject* MDMuVT::__new__(PyTypeObject* type,PyObject* args,PyObject* kwds)
{
    Object* __self=reinterpret_cast<Object*>(type->tp_alloc(type,0));
    PyObject* self=reinterpret_cast<PyObject*>(__self);
    return self;
}
/*--------------------------------------------
 
 --------------------------------------------*/
int MDMuVT::__init__(PyObject* self,PyObject* args,PyObject* kwds)
{
    FuncAPI<type0,type0,type0,std::string,int> f("__init__",{"mu","T","dt","gas_element","seed"});
    
    f.logics<1>()[0]=VLogics("gt",0.0);
    f.logics<2>()[0]=VLogics("gt",0.0);
    f.logics<4>()[0]=VLogics("gt",0);
    if(f(args,kwds)==-1) return -1;
    
    Object* __self=reinterpret_cast<Object*>(self);
    __self->md=new MDMuVT(f.val<0>(),f.val<1>(),f.val<2>(),f.val<3>(),f.val<4>());
    __self->xprt=NULL;
    
    return 0;
}
/*--------------------------------------------
 
 --------------------------------------------*/
PyObject* MDMuVT::__alloc__(PyTypeObject* type,Py_ssize_t)
{
    Object* __self=new Object;
    Py_TYPE(__self)=type;
    Py_REFCNT(__self)=1;
    __self->md=NULL;
        __self->xprt=NULL;
    return reinterpret_cast<PyObject*>(__self);
}
/*--------------------------------------------
 
 --------------------------------------------*/
void MDMuVT::__dealloc__(PyObject* self)
{
    Object* __self=reinterpret_cast<Object*>(self);
    delete __self->md;
    __self->md=NULL;
    if(__self->xprt) Py_DECREF(__self->xprt);
    __self->xprt=NULL;
    delete __self;
}
/*--------------------------------------------*/
PyTypeObject MDMuVT::TypeObject ={PyObject_HEAD_INIT(NULL)};
/*--------------------------------------------*/
int MDMuVT::setup_tp()
{
    TypeObject.tp_name="mapp.md.muvt";
    TypeObject.tp_doc=R"---(
    __init__(mu,T,dt,gas_element,seed)
    
    :math:`\mu VT` ensemble
    
    Molecular dynamics of grand canonical ensemble
        
    Parameters
    ----------
    mu : double
       Chemical potential
    T : double
       Temperature of the ensemble
    dt : double
       Time step for simulation
    gas_element : string
       Gas element
    seed : int
       Random seed for Grand Canonical Monte Carlo
    
    Notes
    -----
       * Thermostat Details
          Nose-Hoover chain

    
    )---";
    
    TypeObject.tp_flags=Py_TPFLAGS_DEFAULT;
    TypeObject.tp_basicsize=sizeof(Object);
    
    TypeObject.tp_new=__new__;
    TypeObject.tp_init=__init__;
    TypeObject.tp_alloc=__alloc__;
    TypeObject.tp_dealloc=__dealloc__;
    setup_tp_getset();
    TypeObject.tp_getset=getset;
    
    TypeObject.tp_base=&MDNVT::TypeObject;
    
    int ichk=PyType_Ready(&TypeObject);
    if(ichk<0) return ichk;
    Py_INCREF(&TypeObject);
    GET_WRAPPER_DOC(TypeObject,__init__)=NULL;
    return ichk;
}
/*--------------------------------------------*/
PyGetSetDef MDMuVT::getset[]=EmptyPyGetSetDef(5);
/*--------------------------------------------*/
void MDMuVT::setup_tp_getset()
{
    getset_nevery(getset[0]);
    getset_nattempts(getset[1]);
    getset_seed(getset[2]);
    getset_gas_element(getset[3]);
}
/*--------------------------------------------
 
 --------------------------------------------*/
void MDMuVT::getset_nevery(PyGetSetDef& getset)
{
    getset.name=(char*)"nevery";
    getset.doc=(char*)R"---(
    (int) gcmc period
    
    Perform GCMC every nevery steps.
    )---";
    getset.get=[](PyObject* self,void*)->PyObject*
    {
        return var<int>::build(reinterpret_cast<Object*>(self)->md->nevery);
    };
    getset.set=[](PyObject* self,PyObject* op,void*)->int
    {
        VarAPI<int> nevery("nevery");
        nevery.logics[0]=VLogics("gt",0);
        int ichk=nevery.set(op);
        if(ichk==-1) return -1;
        
        reinterpret_cast<Object*>(self)->md->nevery=nevery.val;
        return 0;
    };
}
/*--------------------------------------------
 
 --------------------------------------------*/
void MDMuVT::getset_nattempts(PyGetSetDef& getset)
{
    getset.name=(char*)"nattempts";
    getset.doc=(char*)R"---(
    (int) number of gcmc trials
    
    Number of deletion/insertion attempts per GCMC.
    )---";
    getset.get=[](PyObject* self,void*)->PyObject*
    {
        return var<int>::build(reinterpret_cast<Object*>(self)->md->nattempts);
    };
    getset.set=[](PyObject* self,PyObject* op,void*)->int
    {
        VarAPI<int> nattempts("nattempts");
        nattempts.logics[0]=VLogics("gt",0);
        int ichk=nattempts.set(op);
        if(ichk==-1) return -1;
        
        reinterpret_cast<Object*>(self)->md->nattempts=nattempts.val;
        return 0;
    };
}
/*--------------------------------------------
 
 --------------------------------------------*/
void MDMuVT::getset_seed(PyGetSetDef& getset)
{
    getset.name=(char*)"seed";
    getset.doc=(char*)R"---(
    (int) random seed
    
    Random seed for Grand Canonical Monte Carlo
    )---";
    getset.get=[](PyObject* self,void*)->PyObject*
    {
        return var<int>::build(reinterpret_cast<Object*>(self)->md->seed);
    };
    getset.set=[](PyObject* self,PyObject* op,void*)->int
    {
        VarAPI<int> seed("seed");
        seed.logics[0]=VLogics("gt",0);
        int ichk=seed.set(op);
        if(ichk==-1) return -1;
        
        reinterpret_cast<Object*>(self)->md->seed=seed.val;
        return 0;
    };
}
/*--------------------------------------------
 
 --------------------------------------------*/
void MDMuVT::getset_gas_element(PyGetSetDef& getset)
{
    getset.name=(char*)"gas_element";
    getset.doc=(char*)R"---(
    (int) gas element
    
    Gas element for gcmc
    )---";
    getset.get=[](PyObject* self,void*)->PyObject*
    {
        return var<std::string>::build(reinterpret_cast<Object*>(self)->md->gas_elem_name);
    };
    getset.set=[](PyObject* self,PyObject* op,void*)->int
    {
        VarAPI<std::string> gas_element("gas_element");
        int ichk=gas_element.set(op);
        if(ichk==-1) return -1;
        reinterpret_cast<Object*>(self)->md->gas_elem_name=gas_element.val;
        return 0;
    };
}







