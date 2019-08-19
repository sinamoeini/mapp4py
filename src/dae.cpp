#include "dae.h"
#include "MAPP.h"
#include "atoms_dmd.h"
#include "dynamic_dmd.h"
#include "ff_dmd.h"
#include "memory.h"
#include "min_cg_dmd.h"
using namespace MAPP_NS;
/*--------------------------------------------
 
 --------------------------------------------*/
DAE::DAE():
x_err_tol(sqrt(std::numeric_limits<type0>::epsilon())),
alpha_err_tol(sqrt(std::numeric_limits<type0>::epsilon())),
ncs_lcl(0),
chng_box(false),
xprt(NULL),
max_nnewton_iters(5),
max_ngmres_iters(5),
ntally(1000),
nreset(0),
max_nsteps(1000),
a_tol(sqrt(std::numeric_limits<type0>::epsilon())),
min_dt(std::numeric_limits<type0>::epsilon()),
c(NULL),
c_d(NULL)
{
    Algebra::set<__dim__*__dim__>(S_dof[0],false);
    Algebra::set<__dim__*__dim__>(S[0],std::numeric_limits<type0>::quiet_NaN());
    Algebra::set<__dim__*__dim__>(S_err_tol[0],std::numeric_limits<type0>::epsilon());
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
    
    if(chng_box && !__atoms->x_dof->is_empty())
    {
        bool* dof=__atoms->x_dof->begin();
        int __dof_lcl[__dim__];
        Algebra::set<__dim__>(__dof_lcl,0);
        for(int i=0;i<__atoms->natms_lcl;i++,dof+=__dim__)
            Algebra::Do<__dim__>::func([&dof,&__dof_lcl](int i){ if(!dof[i]) __dof_lcl[i]=1;});
        
        int __dof[__dim__];
        Algebra::set<__dim__>(__dof,0);
        MPI_Allreduce(__dof_lcl,__dof,__dim__,MPI_INT,MPI_MAX,__atoms->world);
        std::string err_msg=std::string();
        for(int i=0;i<__dim__;i++)
            for(int j=i;j<__dim__;j++)
                if(S_dof[i][j] && __dof[i])
                {
                    err_msg+="\nyou have atoms that do not have degree freedom in "+Print::to_string(i)
                    +" direction they will be deformed affinely due to defined degrees of freedom H_dof["+Print::to_string(i)+"]["+Print::to_string(j)+"]";
                }
        
        //if(!err_msg.empty()) throw err_msg;
        if(!err_msg.empty() && ntally)
        {
            err_msg="Warning:"+err_msg+"\n";
            fprintf(MAPP::mapp_out,"%s",err_msg.c_str());
        }
    }
    
}
/*--------------------------------------------
 
 --------------------------------------------*/
void DAE::init_static()
{
    //static related
    ncs_lcl=atoms->natms_lcl*c_dim;
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
    nerr_mins=0;
    c_dim=atoms->c_dim;
    ff->c_d->fill();
    ls=new LineSearchBrent();
    bool __H_dof[__dim__][__dim__];
    Algebra::set<__dim__*__dim__>(__H_dof[0],false);
    min=new MinCGDMD(1.0e-9,__H_dof,false,1.0,0.1,ls);
    min->atoms=atoms;
    min->ff=ff;
    min->init();
    min->ntally=0;
    min->__init<false,true,true,false>();
    updt=ff->updt;
    
    ff->calc_ndof();
    a_tol_sqrt_nc_dof=a_tol*sqrt(static_cast<type0>(ff->nc_dof));
    int n=ff->nx_dof+ff->nalpha_dof;
    Algebra::DoLT<__dim__>::func([this,&n](int i,int j)
    {if(!std::isnan(S[i][j])) ++n;});
    sqrt_nx_nalpha_nS_dof=sqrt(static_cast<type0>(n));
    a_tol_sqrt_nx_nalpha_nS_dof=a_tol*sqrt_nx_nalpha_nS_dof;
    x_err_tol_sqrt_ndof=x_err_tol*sqrt(static_cast<type0>(ff->nx_dof));
    alpha_err_tol_sqrt_ndof=alpha_err_tol*sqrt(static_cast<type0>(ff->nalpha_dof));
}
/*--------------------------------------------
 
 --------------------------------------------*/
void DAE::fin()
{
    updt=NULL;
    min->__fin<false,true,true,false>();
    min->ff=NULL;
    min->atoms=NULL;
    delete min;
    min=NULL;
    delete ls;
    ls=NULL;
}
/*--------------------------------------------
 
 --------------------------------------------*/
type0 DAE::calc_err()
{
    if(!chng_box) return sqrt(ff->err_sq_x+ff->err_sq_alpha)/sqrt_nx_nalpha_nS_dof;
    type0 err_sq=0.0;
    type0 (&S_fe)[__dim__][__dim__]=atoms->S_fe;
    Algebra::DoLT<__dim__>::func([&S_fe,&err_sq,this](int i,int j)
    {if(!std::isnan(S[i][j])) err_sq+=Algebra::pow<2>(S[i][j]-S_fe[i][j]);});
    err_sq*=atoms->vol*atoms->vol;
    err_sq+=ff->err_sq_x+ff->err_sq_alpha;
    return sqrt(err_sq)/sqrt_nx_nalpha_nS_dof;
}
/*--------------------------------------------
 
 --------------------------------------------*/
void DAE::min_error_true()
{
    min->__run<false,true,true,false,false,false>(ls,1000);
    MinDMDHandler<false,true,true,false> &handler=*reinterpret_cast<MinDMDHandler<false,true,true,false>*>(&min->handler_buff);
    typedef typename MinDMDHandler<false,true,true,false>::VECTENS0 VECTENS0;
    typedef typename MinDMDHandler<false,true,true,false>::VECTENS1 VECTENS1;
    VECTENS1& f=handler.f;
    VECTENS1& h=handler.h;
    VECTENS0& x=handler.x;
    
    /*
     dynamic_sub;
     ff;
     atoms;
     c_dim;
     alpha_scale;
     x0;
     alpha0;
     some vile hackery so I do not need to create another dynamic object or cast
     */   
    static_assert(sizeof(NewDynamicDMD<true,true,true>)==sizeof(NewDynamicDMD<false,true,true>),"DynamicDMD size mismatch");
    static_assert(alignof(NewDynamicDMD<true,true,true>)==alignof(NewDynamicDMD<false,true,true>),"DynamicDMD align mismatch");
    NewDynamicDMD<true,true,true>* dynamic_ptr=reinterpret_cast<NewDynamicDMD<true,true,true>*>(&handler.dynamic);
    dynamic_ptr->reset();
    
    
    
    GMRES<VECTENS1> gmres(max_ngmres_iters,atoms,__dim__,true,c_dim,true,c_dim,false);
    auto J=[this](VECTENS1& dx,VECTENS1& Jdx)->void
    {
        ff->Jnew(dx.A,dx.vecs[0],dx.vecs[1],Jdx.A,Jdx.vecs[0],Jdx.vecs[1]);
        
        type0 dlog_vol=0.0;
        Algebra::Do<__dim__>::func([&dx,&dlog_vol](int i)
        {
            dlog_vol+=dx.A[i][i];
        });

        type0 tmp=dlog_vol*atoms->vol;
        Algebra::DoLT<__dim__>::func([&Jdx,&tmp,this](int i,int j)
        {
            if(!S_dof[i][j]) Jdx.A[i][j]=0.0;
            else Jdx.A[i][j]-=S[i][j]*tmp;
        });
    };
    
    
    auto get_res=[this,&f]()->type0
    {
        ff->prepJ_n_res(f.vecs[0],f.vecs[1]);
        type0 res_sq=ff->err_sq_x+ff->err_sq_alpha;
        type0 vol_neg=-atoms->vol;
        Algebra::DoLT<__dim__>::func([&res_sq,&f,this,&vol_neg](int i,int j)
         {
             if(S_dof[i][j])
             {
                 f.A[i][j]=ff->F_H[i][j]-S[i][j]*vol_neg;
                 res_sq+=f.A[i][j]*f.A[i][j];
             }
             else
                 f.A[i][j]=0.0;
         });
        return sqrt(res_sq);
    };
    
    type0 r,norm,res=get_res();
    int istep=0;
    for(;istep<max_nnewton_iters && res/a_tol_sqrt_nx_nalpha_nS_dof>1.0;istep++)
    {
        gmres.solve(J,f,0.005*a_tol_sqrt_nx_nalpha_nS_dof,norm,h);
        
        //similar to MinDMDHandler<BC,X,ALPHA,C>::prep_x_d()
        const int natms_lcl=atoms->natms_lcl;
        type0* xvec=atoms->x->begin();
        type0* dxvec=h.vecs[0]->begin();
        type0 __A[__dim__][__dim__];
        Algebra::V_eq<__dim__*__dim__>(&h.A[0][0],&__A[0][0]);
        Algebra::MLT_mul_MLT(atoms->H,__A,h.A);
        for(int i=0;i<natms_lcl;i++,dxvec+=__dim__,xvec+=__dim__)
            Algebra::V_mul_MLT_add_in(xvec,h.A,dxvec);
        
        handler.max_dalpha=std::numeric_limits<type0>::infinity();
        type0 r_lcl=1.0;
        handler.max_alpha_lcl_alpha(r_lcl,h.vecs[1]->begin());
        MPI_Allreduce(&r_lcl,&r,1,Vec<type0>::MPI_T,MPI_MIN,atoms->world);
        
        if(r==1.0)
            x+=h;
        else
        {
            r*=0.5;
            x+=r*h;
        }
        
        dynamic_ptr->update();
        res=get_res();
        
    }
    if(istep) nerr_mins++;
    
    handler.dynamic.reset();
}
/*--------------------------------------------
 
 --------------------------------------------*/
void DAE::min_error_false()
{
    min->__run<false,true,true,false,false,false>(ls,1000);
    MinDMDHandler<false,true,true,false> &handler=*reinterpret_cast<MinDMDHandler<false,true,true,false>*>(&min->handler_buff);
    typedef typename MinDMDHandler<false,true,true,false>::VECTENS0 VECTENS0;
    typedef typename MinDMDHandler<false,true,true,false>::VECTENS1 VECTENS1;
    VECTENS1& f=handler.f;
    VECTENS1& h=handler.h;
    VECTENS0& x=handler.x;
    
    
    NewDynamicDMD<false,true,true>* dynamic_ptr=&handler.dynamic;
    
    
    GMRES<VECTENS1> gmres(max_ngmres_iters,atoms,__dim__,true,c_dim,true,c_dim,false);
    auto J=[this](VECTENS1& dx,VECTENS1& Jdx)->void
    {
        ff->Jnew(dx.vecs[0],dx.vecs[1],Jdx.vecs[0],Jdx.vecs[1]);
    };
    
    
    auto get_res=[this,&f]()->type0
    {
        ff->prepJ_n_res(f.vecs[0],f.vecs[1]);
        type0 res_sq=ff->err_sq_x+ff->err_sq_alpha;
        return sqrt(res_sq);
    };
    
    
    type0 r,norm,res=get_res();
    int istep=0;
    for(;istep<max_nnewton_iters && res/a_tol_sqrt_nx_nalpha_nS_dof>1.0;istep++)
    {
        gmres.solve(J,f,0.005*a_tol_sqrt_nx_nalpha_nS_dof,norm,h);
                
        handler.max_dalpha=std::numeric_limits<type0>::infinity();
        type0 r_lcl=1.0;
        handler.max_alpha_lcl_alpha(r_lcl,h.vecs[1]->begin());
        MPI_Allreduce(&r_lcl,&r,1,Vec<type0>::MPI_T,MPI_MIN,atoms->world);
        
        if(r==1.0)
            x+=h;
        else
        {
            r*=0.5;
            x+=r*h;
        }
        
        dynamic_ptr->update();
        res=get_res();
        
    }
    if(istep) nerr_mins++;
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
    Py_TYPE(__self)=type;
    Py_REFCNT(__self)=1;
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
    TypeObject.tp_name="mapp4py.dmd.dae";
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
PyGetSetDef DAE::getset[]=EmptyPyGetSetDef(8);
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
PyMethodDef DAE::methods[]=EmptyPyMethodDef(1);
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
    (mapp4py.dmd.export) export object
    
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
    tp_methods.ml_meth=(PyCFunction)(PyCFunctionWithKeywords)(
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
    });
    
    tp_methods.ml_doc=(char*)R"---(
    run(atoms,t)
   
    Execute DEA
    
    This method starts differential-algebraic equations for a given atoms object and number of time.
    
    Parameters
    ----------
    atoms : mapp4py.dmd.atoms
        System of interest
    t : double
        Desired time
        
    Returns
    -------
    None

    )---";
}


