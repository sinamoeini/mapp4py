#include "dae.h"
#include "atoms_dmd.h"
#include "dynamic_dmd.h"
#include "ff_dmd.h"
#include "memory.h"
using namespace MAPP_NS;
/*--------------------------------------------
 
 --------------------------------------------*/
DAE::DAE():
max_nsteps(1000),
a_tol(sqrt(std::numeric_limits<type0>::epsilon())),
min_dt(std::numeric_limits<type0>::epsilon()),
c(NULL),
c_d(NULL),
xprt(NULL),
ncs(0),
ntally(1000),
nreset(0),
chng_box(false),
S_dof{[0 ... __dim__-1]={[0 ... __dim__-1]=false}},
S{[0 ... __dim__-1]={[0 ... __dim__-1]=NAN}}
{
}
/*--------------------------------------------
 
 --------------------------------------------*/
DAE::~DAE()
{
}
/*--------------------------------------------
 
 --------------------------------------------*/
void DAE::pre_run_chk(AtomsDMD* __atoms, ForceFieldDMD* __ff)
{
    //check if configuration is loaded
    if(!__atoms)
        throw std::string("cannot start dae run without initial conditions");
    
    //check if force field is loaded
    if(!__ff)
        throw std::string("cannot start dae run without governing equations (force field)");
    
    if(std::isnan(__atoms->kB))
        throw std::string("boltzmann constant should be set prior to dae run");
    
    if(std::isnan(__atoms->hP))
        throw std::string("planck constant should be set prior to dae run");
    
    //check to see if the H_dof components are consistent with stoms->dof
    
    if(!__atoms->dof->is_empty())
    {
        bool* dof=__atoms->dof->begin();
        int __dof_lcl[__dim__]{[0 ... __dim__-1]=0};
        for(int i=0;i<__atoms->natms_lcl;i++,dof+=__dim__)
            Algebra::Do<__dim__>::func([&dof,&__dof_lcl](int i){ if(!dof[i]) __dof_lcl[i]=1;});
        
        int __dof[__dim__]{[0 ... __dim__-1]=0};
        MPI_Allreduce(__dof_lcl,__dof,__dim__,MPI_INT,MPI_MAX,atoms->world);
        std::string err_msg=std::string();
        for(int i=0;i<__dim__;i++)
            for(int j=i;j<__dim__;j++)
                if(S_dof[i][j] && __dof[i])
                {
                    if(!err_msg.empty()) err_msg+="\n";
                    err_msg+="cannot impose stress component ["+Print::to_string(i)+"]["+Print::to_string(j)
                    +"] while any of the atoms do not have degree freedom in "+Print::to_string(i)
                    +" direction";
                }
        
        if(!err_msg.empty()) throw err_msg;
    }
}
/*--------------------------------------------
 
 --------------------------------------------*/
int DAE::calc_ndofs(AtomsDMD* __atoms)
{
    int ndof_lcl=__atoms->c_dim*__atoms->natms_lcl;
    int n=ndof_lcl;
    type0* c=__atoms->c->begin();
    if(!__atoms->dof_c->is_empty())
    {
        bool* c_dof=__atoms->dof_c->begin();
        for(int i=0;i<n;i++)
            if(!c_dof[i] || c[i]<0.0) ndof_lcl--;
    }
    else
        for(int i=0;i<n;i++)
            if(c[i]<0.0) ndof_lcl--;
    int ndof;
    MPI_Allreduce(&ndof_lcl,&ndof,1,MPI_INT,MPI_SUM,__atoms->world);
    return ndof;
}
/*--------------------------------------------
 
 --------------------------------------------*/
void DAE::init_static()
{
    //static related
    ncs=atoms->natms_lcl*c_dim;
    //a_tol_sqrt_nc_dofs=a_tol*sqrt(static_cast<type0>(calc_ndofs(atoms)));
    c=atoms->c->begin();
    c_d=ff->c_d->begin();
}
/*--------------------------------------------
 
 --------------------------------------------*/
void DAE::fin_static()
{
    //static related
    c=c_d=NULL;
}
/*--------------------------------------------
 
 --------------------------------------------*/
void DAE::init()
{
    c_dim=atoms->c_dim;
    ff->c_d->fill();
    dynamic=new DynamicDMD(atoms,ff,chng_box,{},{},{});
    dynamic->init();
    a_tol_sqrt_nc_dofs=a_tol*sqrt(static_cast<type0>(calc_ndofs(atoms)));
}
/*--------------------------------------------
 
 --------------------------------------------*/
void DAE::fin()
{
    dynamic->fin();
    delete dynamic;
    dynamic=NULL;
}
/*--------------------------------------------
 
 --------------------------------------------*/
void DAE::min_error()
{
    VecTens<type0,2> x(atoms,chng_box,atoms->H,atoms->x,atoms->alpha);
    VecTens<type0,2> f(atoms,chng_box,ff->f,ff->f_alpha);
    VecTens<type0,2> h(atoms,chng_box,__dim__,c_dim);
    

    vec* uvecs[2];
    uvecs[0]=atoms->x;
    uvecs[1]=atoms->alpha;
    type0 norm,res;

    
    __GMRES__<VecTens<type0,2>> gmres(10,atoms,chng_box,__dim__,c_dim);
    auto J=[this](VecTens<type0,2>& x,VecTens<type0,2>& Jx)->void
    {
        ff->J_timer(x,Jx);
        
        if(chng_box)
        {
            type0 dlog_vol=0.0;
            type0 (&H)[__dim__][__dim__]=atoms->H;
        
            Algebra::Do<__dim__>::func([&H,&x,&dlog_vol](int i)
            {
                dlog_vol+=x.A[i][i]/H[i][i];
            });
    
            type0 tmp=dlog_vol*atoms->vol;
            Algebra::DoLT<__dim__>::func([&Jx,&tmp,this](int i,int j)
            {
                if(!S_dof[i][j]) Jx.A[i][j]=0.0;
                else Jx.A[i][j]-=S[i][j]*tmp;
            });
        }
    };
    
    
    res=chng_box ? ff->prep_timer(f,S):ff->prep_timer(f);
    //printf("res %e %e\n",res,a_tol_sqrt_nc_dofs);
    
    int istep=0;
    for(;istep<1000 && res/a_tol_sqrt_nc_dofs>1.0;istep++)
    {
        /*
        printf("%d  %e  %e %e %e %e %e %e\n",istep,res/a_tol_sqrt_nc_dofs
        ,atoms->S_fe[0][0]
        ,atoms->S_fe[1][1]
        ,atoms->S_fe[2][2]
        ,atoms->S_fe[1][2]
        ,atoms->S_fe[2][0]
        ,atoms->S_fe[0][1]);*/
        
        gmres.solve(J,f,0.005*a_tol_sqrt_nc_dofs,norm,h);
        x+=h;

        type0 max_alpha_lcl=0.0;
        const int n=atoms->natms_lcl*atoms->alpha->dim;
        type0* alpha_vec=atoms->alpha->begin();
        type0* c_vec=atoms->c->begin();
        for(int i=0;i<n;i++)
            if(c_vec[i]>0.0) max_alpha_lcl=MAX(max_alpha_lcl,alpha_vec[i]);
        MPI_Allreduce(&max_alpha_lcl,&atoms->max_alpha,1,Vec<type0>::MPI_T,MPI_MAX,atoms->world);
        
        if(chng_box)
            atoms->update_H();
        
        dynamic->update(uvecs,2);
        
        res=chng_box ? ff->prep_timer(f,S):ff->prep_timer(f);
        
    }
    
    //printf("%d res %e %e\n",istep,a0,a1);
    //printf("res %e\n",res);
}
/*------------------------------------------------------------------------------------------------------------------------------------
 
 ------------------------------------------------------------------------------------------------------------------------------------*/
PyObject* DAE::__new__(PyTypeObject* type,PyObject* args,PyObject* kwds)
{
    Object* __self=reinterpret_cast<Object*>(type->tp_alloc(type,0));
    PyObject* self=reinterpret_cast<PyObject*>(__self);
    return self;
}
/*--------------------------------------------
 
 --------------------------------------------*/
int DAE::__init__(PyObject* self,PyObject* args,PyObject* kwds)
{
    FuncAPI<> f("__init__");
    
    if(f(args,kwds)==-1) return -1;
    Object* __self=reinterpret_cast<Object*>(self);
    __self->dae=new DAE();
    __self->xprt=NULL;
    
    return 0;
}
/*--------------------------------------------
 
 --------------------------------------------*/
PyObject* DAE::__alloc__(PyTypeObject* type,Py_ssize_t)
{
    Object* __self=new Object;
    __self->ob_type=type;
    __self->ob_refcnt=1;
    __self->dae=NULL;
    __self->xprt=NULL;
    return reinterpret_cast<PyObject*>(__self);
}
/*--------------------------------------------
 
 --------------------------------------------*/
void DAE::__dealloc__(PyObject* self)
{
    Object* __self=reinterpret_cast<Object*>(self);
    delete __self->dae;
    __self->dae=NULL;
    if(__self->xprt) Py_DECREF(__self->xprt);
    __self->xprt=NULL;
    delete __self;
}
/*--------------------------------------------*/
PyTypeObject DAE::TypeObject={PyObject_HEAD_INIT(NULL)};
/*--------------------------------------------*/
int DAE::setup_tp()
{
    TypeObject.tp_name="mapp.dmd.dae";
    TypeObject.tp_doc="chemical integration";
    
    TypeObject.tp_flags=Py_TPFLAGS_DEFAULT;
    TypeObject.tp_basicsize=sizeof(Object);
    
    TypeObject.tp_new=__new__;
    TypeObject.tp_init=__init__;
    TypeObject.tp_alloc=__alloc__;
    TypeObject.tp_dealloc=__dealloc__;
    setup_tp_methods();
    TypeObject.tp_methods=methods;
    setup_tp_getset();
    TypeObject.tp_getset=getset;
    
    int ichk=PyType_Ready(&TypeObject);
    if(ichk<0) return ichk;
    Py_INCREF(&TypeObject);
    return ichk;
    
}
/*--------------------------------------------*/
PyGetSetDef DAE::getset[]={[0 ... 7]={NULL,NULL,NULL,NULL,NULL}};
/*--------------------------------------------*/
void DAE::setup_tp_getset()
{
    getset_a_tol(getset[0]);
    getset_max_nsteps(getset[1]);
    getset_min_dt(getset[2]);
    getset_nreset(getset[3]);
    getset_ntally(getset[4]);
    getset_S(getset[5]);
    getset_export(getset[6]);
}
/*--------------------------------------------*/
PyMethodDef DAE::methods[]={[0 ... 0]={NULL}};
/*--------------------------------------------*/
void DAE::setup_tp_methods()
{
}
/*--------------------------------------------
 
 --------------------------------------------*/
void DAE::getset_a_tol(PyGetSetDef& getset)
{
    getset.name=(char*)"a_tol";
    getset.doc=(char*)R"---(
    (double) LTE tolerance
    
    Absolute error tolerence in local trucation error
    )---";
    getset.get=[](PyObject* self,void*)->PyObject*
    {
        return var<type0>::build(reinterpret_cast<Object*>(self)->dae->a_tol);
    };
    getset.set=[](PyObject* self,PyObject* op,void*)->int
    {
        VarAPI<type0> a_tol("a_tol");
        a_tol.logics[0]=VLogics("gt",0.0)*VLogics("lt",1.0);
        int ichk=a_tol.set(op);
        if(ichk==-1) return -1;
        reinterpret_cast<Object*>(self)->dae->a_tol=a_tol.val;
        return 0;
    };
}
/*--------------------------------------------
 
 --------------------------------------------*/
void DAE::getset_max_nsteps(PyGetSetDef& getset)
{
    getset.name=(char*)"max_nsteps";
    getset.doc=(char*)R"---(
    (int) maximum number of steps
    
    Maximum number of steps to achieve energy minimization
    )---";
    getset.get=[](PyObject* self,void*)->PyObject*
    {
        return var<int>::build(reinterpret_cast<Object*>(self)->dae->max_nsteps);
    };
    getset.set=[](PyObject* self,PyObject* op,void*)->int
    {
        VarAPI<int> max_nsteps("max_nsteps");
        max_nsteps.logics[0]=VLogics("ge",0);
        int ichk=max_nsteps.set(op);
        if(ichk==-1) return -1;
        reinterpret_cast<Object*>(self)->dae->max_nsteps=max_nsteps.val;
        return 0;
    };
}
/*--------------------------------------------
 
 --------------------------------------------*/
void DAE::getset_min_dt(PyGetSetDef& getset)
{
    getset.name=(char*)"min_dt";
    getset.doc=(char*)R"---(
    (double) minimum time step
    
    Minimum time step
    )---";
    getset.get=[](PyObject* self,void*)->PyObject*
    {
        return var<type0>::build(reinterpret_cast<Object*>(self)->dae->min_dt);
    };
    getset.set=[](PyObject* self,PyObject* op,void*)->int
    {
        VarAPI<type0> min_dt("min_dt");
        min_dt.logics[0]=VLogics("gt",0.0);
        int ichk=min_dt.set(op);
        if(ichk==-1) return -1;
        reinterpret_cast<Object*>(self)->dae->min_dt=min_dt.val;
        return 0;
    };
}
/*--------------------------------------------
 
 --------------------------------------------*/
void DAE::getset_nreset(PyGetSetDef& getset)
{
    getset.name=(char*)"nreset";
    getset.doc=(char*)R"---(
    (int) configuration reset period
    
    Number of steps to reset and readjust. If set to 0 no readjustment will occur
    )---";
    getset.get=[](PyObject* self,void*)->PyObject*
    {
        int niters=reinterpret_cast<Object*>(self)->dae->nreset;
        return var<int>::build(niters,NULL);
    };
    getset.set=[](PyObject* self,PyObject* op,void*)->int
    {
        VarAPI<int> nreset("nreset");
        nreset.logics[0]=VLogics("ge",0);
        int ichk=nreset.set(op);
        if(ichk==-1) return -1;
        reinterpret_cast<Object*>(self)->dae->nreset=nreset.val;
        return 0;
    };
}
/*--------------------------------------------
 
 --------------------------------------------*/
void DAE::getset_ntally(PyGetSetDef& getset)
{
    getset.name=(char*)"ntally";
    getset.doc=(char*)R"---(
    (int) thermodynamic tallying period
    
    Number of steps to be taken from one thermodynamics output to the next.
    )---";
    getset.get=[](PyObject* self,void*)->PyObject*
    {
        return var<int>::build(reinterpret_cast<Object*>(self)->dae->ntally);
    };
    getset.set=[](PyObject* self,PyObject* op,void*)->int
    {
        VarAPI<int> ntally("ntally");
        ntally.logics[0]=VLogics("ge",0);
        int ichk=ntally.set(op);
        if(ichk==-1) return -1;
        reinterpret_cast<Object*>(self)->dae->ntally=ntally.val;
        return 0;
    };
}
/*--------------------------------------------
 
 --------------------------------------------*/
void DAE::getset_S(PyGetSetDef& getset)
{
    getset.name=(char*)"S";
    getset.doc=(char*)R"---(
    (symm<double[dim][dim]>) external stress tensor
    
    External stress imposed on system, here dim is the dimension of simulation
    )---";
    getset.get=[](PyObject* self,void*)->PyObject*
    {
        return var<symm<type0[__dim__][__dim__]>>::build(reinterpret_cast<Object*>(self)->dae->S);
    };
    getset.set=[](PyObject* self,PyObject* op,void*)->int
    {
        VarAPI<symm<type0[__dim__][__dim__]>> S("S");
        int ichk=S.set(op);
        if(ichk==-1) return -1;
        
        bool (&__S_dof)[__dim__][__dim__]=reinterpret_cast<Object*>(self)->dae->S_dof;
        type0 (&__S)[__dim__][__dim__]=reinterpret_cast<Object*>(self)->dae->S;
        bool& __chng_box=reinterpret_cast<Object*>(self)->dae->chng_box;
        __chng_box=false;
        Algebra::DoLT<__dim__>::func([&__S_dof,&__S,&__chng_box,&S](int i,int j)
        {
            
            if(std::isnan(S.val[i][j]))
            {
                __S[i][j]=__S[j][i]=NAN;
                __S_dof[i][j]=__S_dof[j][i]=false;
            }
            else
            {
                __S[i][j]=__S[j][i]=S.val[i][j];
                __S_dof[i][j]=__S_dof[j][i]=true;
                __chng_box=true;
            }
        });
        
        return 0;
    };
}
/*--------------------------------------------
 
 --------------------------------------------*/
void DAE::getset_export(PyGetSetDef& getset)
{
    getset.name=(char*)"export";
    getset.doc=(char*)R"---(
    (mapp.dmd.export) export object
    
    Export object to record the snapshots of the system while minimizing
    )---";
    getset.get=[](PyObject* self,void*)->PyObject*
    {
        ExportDMD::Object* xprt=reinterpret_cast<Object*>(self)->xprt;
        if(!xprt) Py_RETURN_NONE;
        Py_INCREF(xprt);
        return reinterpret_cast<PyObject*>(xprt);
    };
    getset.set=[](PyObject* self,PyObject* op,void*)->int
    {
        VarAPI<OP<ExportDMD>> xprt("export");
        int ichk=xprt.set(op);
        if(ichk==-1) return -1;
        if(reinterpret_cast<Object*>(self)->xprt) Py_DECREF(reinterpret_cast<Object*>(self)->xprt);
        Py_INCREF(xprt.val.ob);
        reinterpret_cast<Object*>(self)->xprt=reinterpret_cast<ExportDMD::Object*>(xprt.val.ob);
        return 0;
    };
}
/*--------------------------------------------
 
 --------------------------------------------*/
void DAE::ml_run(PyMethodDef& tp_methods)
{
    tp_methods.ml_flags=METH_VARARGS | METH_KEYWORDS;
    tp_methods.ml_name="run";    
    tp_methods.ml_meth=(PyCFunction)(PyCFunctionWithKeywords)
    [](PyObject* self,PyObject* args,PyObject* kwds)->PyObject*
    {
        Object* __self=reinterpret_cast<Object*>(self);
        FuncAPI<OP<AtomsDMD>,type0> f("run",{"atoms","t"});
        f.logics<1>()[0]=VLogics("ge",0.0);
        if(f(args,kwds)) return NULL;
        
        AtomsDMD* __atoms=reinterpret_cast<AtomsDMD::Object*>(f.val<0>().ob)->atoms;
        ForceFieldDMD* __ff=reinterpret_cast<AtomsDMD::Object*>(f.val<0>().ob)->ff;
        ExportDMD* __xprt=__self->xprt==NULL ? NULL:__self->xprt->xprt;
        try
        {
            __self->dae->pre_run_chk(__atoms,__ff);
        }
        catch(std::string& err_msg)
        {
            PyErr_SetString(PyExc_TypeError,err_msg.c_str());
            return NULL;
        }
        
        __self->dae->atoms=__atoms;
        __self->dae->ff=__ff;
        __self->dae->xprt=__xprt;
        
        try
        {
            __self->dae->init();
        }
        catch(std::string& err_msg)
        {
            __self->dae->xprt=NULL;
            __self->dae->ff=NULL;
            __self->dae->atoms=NULL;
            PyErr_SetString(PyExc_TypeError,err_msg.c_str());
            return NULL;
        }
        

        
        __self->dae->run(f.val<1>());
        

        __self->dae->fin();
        
        __self->dae->xprt=NULL;
        __self->dae->ff=NULL;
        __self->dae->atoms=NULL;
        
        Py_RETURN_NONE;
    };
    
    tp_methods.ml_doc=(char*)R"---(
    run(atoms,t)
   
    Execute DEA
    
    This method starts differential-algebraic equations for a given atoms object and number of time.
    
    Parameters
    ----------
    atoms : mapp.dmd.atoms
        System of interest
    t : double
        Desired time
        
    Returns
    -------
    None

    )---";
}


