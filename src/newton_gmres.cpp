#include "newton_gmres.h"
#include "xmath.h"
#include "atoms_dmd.h"
#include "ff_styles.h"
#include "dynamic_dmd.h"
#include "thermo_dynamics.h"
#include <cmath>
/*--------------------------------------------
 
 --------------------------------------------*/
NewtonGMRES::NewtonGMRES():
max_nsteps(1000),
a_tol(sqrt(std::numeric_limits<type0>::epsilon())),
xprt(NULL),
ntally(1000),
m(0),
chng_box(false),
S_dof{DESIG2(__dim__,__dim__,false)},
S{DESIG2(__dim__,__dim__,NAN)}
{
}
/*--------------------------------------------
 
 --------------------------------------------*/
NewtonGMRES::~NewtonGMRES()
{
}
/*--------------------------------------------
 
 --------------------------------------------*/
void NewtonGMRES::pre_run_chk(AtomsDMD* __atoms, ForceFieldDMD* __ff)
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
        int __dof_lcl[__dim__]{DESIG(__dim__,0)};
        for(int i=0;i<__atoms->natms_lcl;i++,dof+=__dim__)
            Algebra::Do<__dim__>::func([&dof,&__dof_lcl](int i){ if(!dof[i]) __dof_lcl[i]=1;});
        
        int __dof[__dim__]{DESIG(__dim__,0)};
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
void NewtonGMRES::init()
{
    dynamic=new DynamicDMD(atoms,ff,chng_box,{},{},{});
    dynamic->init();
    c_dim=atoms->c_dim;
    a_tol_sqrt_nc_dofs=a_tol*sqrt(static_cast<type0>(atoms->natms*(__dim__+c_dim)));
    
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
void NewtonGMRES::fin()
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
void NewtonGMRES::run(int nsteps)
{
    
    
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
    
    
    int step=atoms->step;
    int nevery_xprt=xprt==NULL ? 0:xprt->nevery;
    
    
    
    __GMRES__<VecTens<type0,2>> gmres(m,max_ngmres_iters,atoms,chng_box,__dim__,c_dim);
    VecTens<type0,2> x(atoms,chng_box,atoms->H,atoms->x,atoms->alpha);
    VecTens<type0,2> f(atoms,chng_box,ff->f,ff->f_alpha);
    
    
    VecTens<type0,1> f_x(atoms,false,ff->f);
    VecTens<type0,1> f_alpha(atoms,false,ff->f_alpha);
    VecTens<type0,2> h(atoms,chng_box,__dim__,c_dim);

    
    vec* uvecs[2];
    uvecs[0]=atoms->x;
    uvecs[1]=atoms->alpha;
    type0 norm,err,err_x,err_alpha;
    
    err=(chng_box ? ff->prep_timer(f,S):ff->prep_timer(f))/a_tol_sqrt_nc_dofs;
    err_alpha=sqrt(f_alpha*f_alpha)/a_tol_sqrt_nc_dofs;
    err_x=sqrt(err*err-err_alpha*err_alpha);
    ThermoDynamics thermo(12,
    "REL_ERR",err,
    "REL_ERR_alpha",err_x,
    "REL_ERR_x",err_alpha,
    "FE",atoms->fe/*,
    "S[0][0]",atoms->S_fe[0][0],
    "S[1][1]",atoms->S_fe[1][1],
    "S[2][2]",atoms->S_fe[2][2],
    "S[1][2]",atoms->S_fe[2][1],
    "S[2][0]",atoms->S_fe[2][0],
    "S[0][1]",atoms->S_fe[1][0]*/
    );
    thermo.init();
    int istep=0;
    if(ntally) thermo.print(step);
    if(nevery_xprt) xprt->write(step);
    
    type0 eta_max=0.999,gamma=0.9,err_prev;
    type0 eta=eta_max,eta_prev,eta_A;
    
    for(;istep<nsteps && err>1.0;istep++)
    {
        
        gmres.solve_restart(J,f,a_tol_sqrt_nc_dofs,norm,h);
        
        //gmres.solve_restart(J,f,eta*a_tol_sqrt_nc_dofs,norm,h);
        //gmres.solve_restart(J,f,eta*err,norm,h);
        
        //printf("%e %d\n",eta,gmres.iter);
        
        
        x+=h;
        
        
        const int n=atoms->natms_lcl*c_dim;
        type0* alpha_vec=atoms->alpha->begin();
        type0 max_alpha_lcl=0.0;
        type0* c_vec=atoms->c->begin();
        for(int i=0;i<n;i++)
            if(c_vec[i]>=0.0) max_alpha_lcl=MAX(max_alpha_lcl,alpha_vec[i]);
        MPI_Allreduce(&max_alpha_lcl,&atoms->max_alpha,1,Vec<type0>::MPI_T,MPI_MAX,atoms->world);
        
        if(chng_box) atoms->update_H();
        
        dynamic->update(uvecs,2);
        
        err_prev=err;
        err=(chng_box ? ff->prep_timer(f,S):ff->prep_timer(f))/a_tol_sqrt_nc_dofs;
        err_alpha=sqrt(f_alpha*f_alpha)/a_tol_sqrt_nc_dofs;
        err_x=sqrt(err*err-err_alpha*err_alpha);
        if(ntally && (istep+1)%ntally==0) thermo.print(step+istep+1);
        if(nevery_xprt && (istep+1)%nevery_xprt==0) xprt->write(step+istep+1);
        
        
        
        eta_prev=eta;
        eta_A=gamma*(err/err_prev)*(err/err_prev);
        
        if(gamma*eta_prev*eta_prev<=0.1)
        {
            eta=MIN(eta_max,eta_A);
        }
        else
        {
            eta=MIN(eta_max,MAX(eta_A,gamma*eta_prev*eta_prev));
        }
        
        
    }
    
    if(ntally && istep%ntally) thermo.print(step+istep);
    if(nevery_xprt && istep%nevery_xprt) xprt->write(step+istep);
    
    if(ntally) thermo.fin();

    atoms->step+=istep;
}
/*------------------------------------------------------------------------------------------------------------------------------------
 
 ------------------------------------------------------------------------------------------------------------------------------------*/
PyObject* NewtonGMRES::__new__(PyTypeObject* type,PyObject* args,PyObject* kwds)
{
    Object* __self=reinterpret_cast<Object*>(type->tp_alloc(type,0));
    PyObject* self=reinterpret_cast<PyObject*>(__self);
    return self;
}
/*--------------------------------------------
 
 --------------------------------------------*/
int NewtonGMRES::__init__(PyObject* self,PyObject* args,PyObject* kwds)
{
    FuncAPI<> f("__init__");
    
    if(f(args,kwds)==-1) return -1;
    Object* __self=reinterpret_cast<Object*>(self);
    __self->ngmres=new NewtonGMRES();
    __self->xprt=NULL;
    
    return 0;
}
/*--------------------------------------------
 
 --------------------------------------------*/
PyObject* NewtonGMRES::__alloc__(PyTypeObject* type,Py_ssize_t)
{
    Object* __self=new Object;
    __self->ob_type=type;
    __self->ob_refcnt=1;
    __self->ngmres=NULL;
    __self->xprt=NULL;
    return reinterpret_cast<PyObject*>(__self);
}
/*--------------------------------------------
 
 --------------------------------------------*/
void NewtonGMRES::__dealloc__(PyObject* self)
{
    Object* __self=reinterpret_cast<Object*>(self);
    delete __self->ngmres;
    __self->ngmres=NULL;
    if(__self->xprt) Py_DECREF(__self->xprt);
    __self->xprt=NULL;
    delete __self;
}
/*--------------------------------------------*/
PyTypeObject NewtonGMRES::TypeObject={PyObject_HEAD_INIT(NULL)};
/*--------------------------------------------*/
int NewtonGMRES::setup_tp()
{
    TypeObject.tp_name="mapp.dmd.newton_gmres";
    TypeObject.tp_doc="Newton GMRES";
    
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
PyGetSetDef NewtonGMRES::getset[]=EmptyPyGetSetDef(7);
/*--------------------------------------------*/
void NewtonGMRES::setup_tp_getset()
{
    getset_S(getset[0]);
    getset_max_ngmres_iters(getset[1]);
    getset_a_tol(getset[2]);
    getset_ntally(getset[3]);
    getset_export(getset[4]);
    getset_m(getset[5]);
}
/*--------------------------------------------*/
PyMethodDef NewtonGMRES::methods[]=EmptyPyMethodDef(2);
/*--------------------------------------------*/
void NewtonGMRES::setup_tp_methods()
{
    ml_run(methods[0]);
}
/*--------------------------------------------
 
 --------------------------------------------*/
void NewtonGMRES::getset_a_tol(PyGetSetDef& getset)
{
    getset.name=(char*)"a_tol";
    getset.doc=(char*)R"---(
    (double) LTE tolerance
    
    Absolute error tolerence in local trucation error
    )---";
    getset.get=[](PyObject* self,void*)->PyObject*
    {
        return var<type0>::build(reinterpret_cast<Object*>(self)->ngmres->a_tol);
    };
    getset.set=[](PyObject* self,PyObject* op,void*)->int
    {
        VarAPI<type0> a_tol("a_tol");
        a_tol.logics[0]=VLogics("gt",0.0)*VLogics("lt",1.0);
        int ichk=a_tol.set(op);
        if(ichk==-1) return -1;
        reinterpret_cast<Object*>(self)->ngmres->a_tol=a_tol.val;
        return 0;
    };
}
/*--------------------------------------------
 
 --------------------------------------------*/
void NewtonGMRES::getset_ntally(PyGetSetDef& getset)
{
    getset.name=(char*)"ntally";
    getset.doc=(char*)R"---(
    (int) thermodynamic tallying period
    
    Number of steps to be taken from one thermodynamics output to the next.
    )---";
    getset.get=[](PyObject* self,void*)->PyObject*
    {
        return var<int>::build(reinterpret_cast<Object*>(self)->ngmres->ntally);
    };
    getset.set=[](PyObject* self,PyObject* op,void*)->int
    {
        VarAPI<int> ntally("ntally");
        ntally.logics[0]=VLogics("ge",0);
        int ichk=ntally.set(op);
        if(ichk==-1) return -1;
        reinterpret_cast<Object*>(self)->ngmres->ntally=ntally.val;
        return 0;
    };
}
/*--------------------------------------------
 
 --------------------------------------------*/
void NewtonGMRES::getset_S(PyGetSetDef& getset)
{
    getset.name=(char*)"S";
    getset.doc=(char*)R"---(
    (symm<double[dim][dim]>) external stress tensor
    
    External stress imposed on system, here dim is the dimension of simulation
    )---";
    getset.get=[](PyObject* self,void*)->PyObject*
    {
        return var<symm<type0[__dim__][__dim__]>>::build(reinterpret_cast<Object*>(self)->ngmres->S);
    };
    getset.set=[](PyObject* self,PyObject* op,void*)->int
    {
        VarAPI<symm<type0[__dim__][__dim__]>> S("S");
        int ichk=S.set(op);
        if(ichk==-1) return -1;
        
        bool (&__S_dof)[__dim__][__dim__]=reinterpret_cast<Object*>(self)->ngmres->S_dof;
        type0 (&__S)[__dim__][__dim__]=reinterpret_cast<Object*>(self)->ngmres->S;
        bool& __chng_box=reinterpret_cast<Object*>(self)->ngmres->chng_box;
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
void NewtonGMRES::getset_export(PyGetSetDef& getset)
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
void NewtonGMRES::getset_max_ngmres_iters(PyGetSetDef& getset)
{
    getset.name=(char*)"max_ngmres_iters";
    getset.doc=(char*)R"---(
    (int) maximim number of gmres iterations
    
    Maximum number of iterations of linear solver (GMRES)
    )---";
    getset.get=[](PyObject* self,void*)->PyObject*
    {
        return var<int>::build(reinterpret_cast<Object*>(self)->ngmres->max_ngmres_iters);
    };
    getset.set=[](PyObject* self,PyObject* op,void*)->int
    {
        VarAPI<int> max_ngmres_iters("max_ngmres_iters");
        max_ngmres_iters.logics[0]=VLogics("gt",0);
        int ichk=max_ngmres_iters.set(op);
        if(ichk==-1) return -1;
        reinterpret_cast<Object*>(self)->ngmres->max_ngmres_iters=max_ngmres_iters.val;
        return 0;
    };
}
/*--------------------------------------------
 
 --------------------------------------------*/
void NewtonGMRES::getset_m(PyGetSetDef& getset)
{
    getset.name=(char*)"m";
    getset.doc=(char*)R"---(
    (int) maximim number of gmres iterations
    
    Maximum number of restart for GMRES
    )---";
    getset.get=[](PyObject* self,void*)->PyObject*
    {
        return var<int>::build(reinterpret_cast<Object*>(self)->ngmres->m);
    };
    getset.set=[](PyObject* self,PyObject* op,void*)->int
    {
        VarAPI<int> max_ngmres_iters("m");
        max_ngmres_iters.logics[0]=VLogics("ge",0);
        int ichk=max_ngmres_iters.set(op);
        if(ichk==-1) return -1;
        reinterpret_cast<Object*>(self)->ngmres->max_ngmres_iters=max_ngmres_iters.val;
        return 0;
    };
}
/*--------------------------------------------
 
 --------------------------------------------*/
void NewtonGMRES::ml_run(PyMethodDef& tp_methods)
{
    tp_methods.ml_flags=METH_VARARGS | METH_KEYWORDS;
    tp_methods.ml_name="run";    
    tp_methods.ml_meth=(PyCFunction)(PyCFunctionWithKeywords)(
    [](PyObject* self,PyObject* args,PyObject* kwds)->PyObject*
    {
        Object* __self=reinterpret_cast<Object*>(self);
        FuncAPI<OP<AtomsDMD>,int> f("run",{"atoms","max_nsteps"});
        f.logics<1>()[0]=VLogics("ge",0);        
        if(f(args,kwds)) return NULL;
        
        AtomsDMD* __atoms=reinterpret_cast<AtomsDMD::Object*>(f.val<0>().ob)->atoms;
        ForceFieldDMD* __ff=reinterpret_cast<AtomsDMD::Object*>(f.val<0>().ob)->ff;
        ExportDMD* __xprt=__self->xprt==NULL ? NULL:__self->xprt->xprt;
        try
        {
            __self->ngmres->pre_run_chk(__atoms,__ff);
        }
        catch(std::string& err_msg)
        {
            PyErr_SetString(PyExc_TypeError,err_msg.c_str());
            return NULL;
        }
        
        __self->ngmres->atoms=__atoms;
        __self->ngmres->ff=__ff;
        __self->ngmres->xprt=__xprt;
        
        try
        {
            __self->ngmres->init();
        }
        catch(std::string& err_msg)
        {
            __self->ngmres->xprt=NULL;
            __self->ngmres->ff=NULL;
            __self->ngmres->atoms=NULL;
            PyErr_SetString(PyExc_TypeError,err_msg.c_str());
            return NULL;
        }
        

        
        __self->ngmres->run(f.val<1>());
        

        __self->ngmres->fin();
        
        __self->ngmres->xprt=NULL;
        __self->ngmres->ff=NULL;
        __self->ngmres->atoms=NULL;
        
        Py_RETURN_NONE;
    });
    
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
    
    
    
    
    
    
