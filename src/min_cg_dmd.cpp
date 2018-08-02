#include "min_cg_dmd.h"
#include <stdlib.h>
using namespace MAPP_NS;
/*--------------------------------------------
 constructor
 --------------------------------------------*/
MinCGDMD::MinCGDMD(type0 __e_tol,
bool(&__H_dof)[__dim__][__dim__],bool __affine,type0 __max_dx,type0 __max_dalpha,LineSearch* __ls):
Min(__e_tol,__H_dof,__affine,__max_dx,__ls),
atoms(NULL),
ff(NULL),
max_dalpha(__max_dalpha),
xprt(NULL)
{
}
/*--------------------------------------------
 destructor
 --------------------------------------------*/
MinCGDMD::~MinCGDMD()
{
    atoms=NULL;
    ff=NULL;
}
/*--------------------------------------------
 pre run check it throw excepctions
 --------------------------------------------*/
void MinCGDMD::pre_run_chk(AtomsDMD* atoms,ForceFieldDMD* ff)
{
    try
    {
        Min::pre_run_chk(atoms,ff);
    }
    catch (std::string& err_msg)
    {
        throw err_msg;
    }
    
    if(std::isnan(atoms->kB))
        throw std::string("boltzmann constant should be set prior to minimizatiom");
    
    if(std::isnan(atoms->hP))
        throw std::string("planck constant should be set prior to minimizatiom");

    
    if(std::isnan(atoms->temp))
        throw std::string("temperature should be set prior to minimizatiom");
}
/*--------------------------------------------
 
 --------------------------------------------*/
void MinCGDMD::force_calc()
{
    ff->derivative_timer();

    if(chng_box)
    {
        type0 (&S_fe)[__dim__][__dim__]=atoms->S_fe;
        type0 neg_v=-atoms->vol;
        Algebra::DoLT<__dim__>::func([this,&S_fe,&neg_v](int i,int j)
        {f.A[i][j]=H_dof[i][j] ? S_fe[i][j]*neg_v:0.0;});
        
    }
}
/*--------------------------------------------
 
 --------------------------------------------*/
void MinCGDMD::prep()
{
    x_d=h;
    if(!chng_box) return;
    
    const int natms_lcl=atoms->natms_lcl;
    type0* xvec=x0.vecs[0]->begin();
    type0* hvec=h.vecs[0]->begin();
    type0* x_dvec=x_d.vecs[0]->begin();
    Algebra::MLT_mul_MLT(atoms->H,h.A,x_d.A);
    
    
    if(affine)
    {
        for(int iatm=0;iatm<natms_lcl;iatm++,xvec+=__dim__,x_dvec+=__dim__)
            Algebra::V_mul_MLT(xvec,h.A,x_dvec);
    }
    else
    {
        for(int iatm=0;iatm<natms_lcl;iatm++,xvec+=__dim__,x_dvec+=__dim__,hvec+=__dim__)
            Algebra::V_mul_MLT_add_in(xvec,h.A,x_dvec);
    }
}
/*--------------------------------------------
 
 --------------------------------------------*/
type0 MinCGDMD::calc_ndofs()
{
    return 0.0;
}
/*--------------------------------------------
 init before a run
 --------------------------------------------*/
void MinCGDMD::init()
{
    const int c_dim=atoms->c->dim;
    x.~VecTens();
    new (&x) VecTens<type0,2>(atoms,chng_box,atoms->H,atoms->x,atoms->alpha);
    f.~VecTens();
    new (&f) VecTens<type0,2>(atoms,chng_box,ff->f,ff->f_alpha);
    h.~VecTens();
    new (&h) VecTens<type0,2>(atoms,chng_box,__dim__,c_dim);
    x0.~VecTens();
    new (&x0) VecTens<type0,2>(atoms,chng_box,__dim__,c_dim);
    x_d.~VecTens();
    new (&x_d) VecTens<type0,2>(atoms,chng_box,__dim__,c_dim);
    f0.~VecTens();
    new (&f0) VecTens<type0,2>(atoms,chng_box,__dim__,c_dim);
    
    dynamic=new DynamicDMD(atoms,ff,chng_box,{},
    {atoms->x_dof,atoms->alpha_dof,atoms->c_dof,h.vecs[0],h.vecs[1],x0.vecs[0],x0.vecs[1],x_d.vecs[0],x_d.vecs[1],f0.vecs[0],f0.vecs[1]},{});
    
    dynamic->init();
    
    uvecs[0]=atoms->x;
    uvecs[1]=atoms->alpha;
    
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
 finishing minimization
 --------------------------------------------*/
void MinCGDMD::fin()
{
    if(xprt)
    {
        xprt->fin();
        xprt->atoms=NULL;
    }
    
    uvecs[1]=NULL;
    uvecs[0]=NULL;
    
    dynamic->fin();
    delete dynamic;
    dynamic=NULL;
    
    f0.~VecTens();
    x_d.~VecTens();
    x0.~VecTens();
    h.~VecTens();
    f.~VecTens();
    x.~VecTens();
}
/*--------------------------------------------
 min
 --------------------------------------------*/
void MinCGDMD::run(int nsteps)
{
    if(dynamic_cast<LineSearchGoldenSection*>(ls))
        return run(dynamic_cast<LineSearchGoldenSection*>(ls),nsteps);
    
    if(dynamic_cast<LineSearchBrent*>(ls))
        return run(dynamic_cast<LineSearchBrent*>(ls),nsteps);
    
    if(dynamic_cast<LineSearchBackTrack*>(ls))
        return run(dynamic_cast<LineSearchBackTrack*>(ls),nsteps);
}

/*--------------------------------------------
 
 --------------------------------------------*/
void MinCGDMD::refine(int n,int nsteps)
{
    const int c_dim=atoms->c->dim;
    
    
    x.~VecTens();
    new (&x) VecTens<type0,2>(atoms,false,atoms->H,atoms->x,atoms->alpha);
    f.~VecTens();
    new (&f) VecTens<type0,2>(atoms,false,ff->f,ff->f_alpha);
    h.~VecTens();
    new (&h) VecTens<type0,2>(atoms,false,__dim__,c_dim);
    
    dynamic=new DynamicDMD(atoms,ff,false,{atoms->x_dof,atoms->alpha_dof,atoms->c_dof},{},{});
    dynamic->init();
    uvecs[0]=atoms->x;
    uvecs[1]=atoms->alpha;
    type0 norm,res;
    
    
    
    __GMRES__<VecTens<type0,2>> gmres(n,atoms,false,__dim__,c_dim);
    auto J=[this](VecTens<type0,2>& x,VecTens<type0,2>& Jx)->void
    {
        ff->J(x,Jx);
    };
    
    
    res=ff->prep_timer(f);
    for(int istep=0;istep<nsteps && res>1.0e-8;istep++)
    {
        printf("res %e\n",res);
        gmres.solve(J,f,e_tol,norm,h);
        
        
        //printf("-h.f %e\n",-(h*f));
        
        x+=h;
        
        type0 max_alpha_lcl=0.0;
        const int n=atoms->natms_lcl*atoms->alpha->dim;
        type0* alpha_vec=atoms->alpha->begin();
        type0* c_vec=atoms->c->begin();
        for(int i=0;i<n;i++)
            if(c_vec[i]>=0.0) max_alpha_lcl=MAX(max_alpha_lcl,alpha_vec[i]);
        MPI_Allreduce(&max_alpha_lcl,&atoms->max_alpha,1,Vec<type0>::MPI_T,MPI_MAX,atoms->world);
        
        dynamic->update(uvecs,2);
        res=ff->prep_timer(f);
        
    }
    
    
    uvecs[1]=NULL;
    uvecs[0]=NULL;
    
    dynamic->fin();
    delete dynamic;
    dynamic=NULL;
    
    h.~VecTens();
    f.~VecTens();
    x.~VecTens();
}
/*--------------------------------------------
 
 --------------------------------------------*/
type0 MinCGDMD::F(type0 alpha)
{
    x=x0+alpha*x_d;
    type0 max_alpha_lcl=0.0;
    const int n=atoms->natms_lcl*atoms->alpha->dim;
    type0* alpha_vec=atoms->alpha->begin();
    type0* c_vec=atoms->c->begin();
    for(int i=0;i<n;i++)
        if(c_vec[i]>=0.0) max_alpha_lcl=MAX(max_alpha_lcl,alpha_vec[i]);
    MPI_Allreduce(&max_alpha_lcl,&atoms->max_alpha,1,Vec<type0>::MPI_T,MPI_MAX,atoms->world);
    
    if(chng_box)
        atoms->update_H();

    dynamic->update(uvecs,2);
    return ff->value_timer();
}
/*--------------------------------------------
 inner product of f and h
 --------------------------------------------*/
type0 MinCGDMD::dF(type0 alpha,type0& drev)
{
    x=x0+alpha*h;
    type0 max_alpha_lcl=0.0;
    const int n=atoms->natms_lcl*atoms->alpha->dim;
    type0* alpha_vec=atoms->alpha->begin();
    type0* c_vec=atoms->c->begin();
    for(int i=0;i<n;i++)
        if(c_vec[i]>=0.0) max_alpha_lcl=MAX(max_alpha_lcl,alpha_vec[i]);
    MPI_Allreduce(&max_alpha_lcl,&atoms->max_alpha,1,Vec<type0>::MPI_T,MPI_MAX,atoms->world);
    
    if(chng_box)
        atoms->update_H();
    
    dynamic->update(uvecs,2);
    force_calc();
    
    drev=-(f*h);
    return atoms->fe;
}
/*--------------------------------------------
 find maximum h
 --------------------------------------------*/
void MinCGDMD::ls_prep(type0& dfa,type0& h_norm,type0& max_a)
{
    
    h_norm=h*h;
    
    dfa=-f_h;
    
    if(h_norm==0.0)
    {
        max_a=0.0;
        dfa=0.0;
        return;
    }
    
    if(dfa>=0.0)
    {
        max_a=0.0;
        dfa=1.0;
        return;
    }
    
    h_norm=sqrt(h_norm);
    
    //type0 max_a_lcl=max_dx;
    type0 max_x_d_lcl=0.0;
    type0* x_dvec=x_d.vecs[0]->begin();
    const int natms_lcl=atoms->natms_lcl;
    int n=natms_lcl*__dim__;
    for(int i=0;i<n;i++)
        max_x_d_lcl=MAX(max_x_d_lcl,fabs(x_dvec[i]));

    
    type0 max_alpha_ratio=std::numeric_limits<type0>::infinity();
    type0* x_d_alphavec=x_d.vecs[1]->begin();
    type0* alphavec=x.vecs[1]->begin();
    n=natms_lcl*atoms->alpha->dim;
    type0 max_x_d_alpha_lcl=0.0;
    for(int i=0;i<n;i++)
    {
        if(x_d_alphavec[i]<0.0)
            max_alpha_ratio=MIN(-alphavec[i]/x_d_alphavec[i],max_alpha_ratio);
        max_x_d_alpha_lcl=MAX(max_x_d_alpha_lcl,fabs(x_d_alphavec[i]));
    }
    type0 max_a_lcl=MIN(fabs(max_dx/max_x_d_lcl),max_alpha_ratio);
    max_a_lcl=MIN(max_a_lcl,fabs(max_dalpha/max_x_d_alpha_lcl));

    MPI_Allreduce(&max_a_lcl,&max_a,1,Vec<type0>::MPI_T,MPI_MIN,atoms->world);
    
}
/*--------------------------------------------
 reset to initial position
 --------------------------------------------*/
void MinCGDMD::F_reset()
{
    x=x0;
    const int n=atoms->natms_lcl*atoms->alpha->dim;
    type0 max_alpha_lcl=0.0;
    type0* alpha_vec=atoms->alpha->begin();
    type0* c_vec=atoms->c->begin();
    for(int i=0;i<n;i++)
        if(c_vec[i]>=0.0) max_alpha_lcl=MAX(max_alpha_lcl,alpha_vec[i]);
    MPI_Allreduce(&max_alpha_lcl,&atoms->max_alpha,1,Vec<type0>::MPI_T,MPI_MAX,atoms->world);
    if(chng_box) atoms->update_H();
    dynamic->update(atoms->x);
}
/*------------------------------------------------------------------------------------------------------------------------------------
 
 ------------------------------------------------------------------------------------------------------------------------------------*/
PyObject* MinCGDMD::__new__(PyTypeObject* type,PyObject* args,PyObject* kwds)
{
    Object* __self=reinterpret_cast<Object*>(type->tp_alloc(type,0));
    PyObject* self=reinterpret_cast<PyObject*>(__self);
    return self;
}
/*--------------------------------------------
 
 --------------------------------------------*/
int MinCGDMD::__init__(PyObject* self,PyObject* args,PyObject* kwds)
{
    FuncAPI<type0,symm<bool[__dim__][__dim__]>,bool,type0,type0,OP<LineSearch>> f("__init__",{"e_tol","H_dof","affine","max_dx","max_dalpha","ls"});
    f.noptionals=6;
    f.logics<0>()[0]=VLogics("ge",0.0);
    f.logics<3>()[0]=VLogics("gt",0.0);
    f.logics<4>()[0]=VLogics("gt",0.0);
    
    //set the defualts
    f.val<0>()=sqrt(std::numeric_limits<type0>::epsilon());
    for(int i=0;i<__dim__;i++) for(int j=0;j<__dim__;j++)f.val<1>()[i][j]=false;
    f.val<2>()=false;
    f.val<3>()=1.0;
    f.val<4>()=0.1;
    PyObject* empty_tuple=PyTuple_New(0);
    PyObject* empty_dict=PyDict_New();
    PyObject* __ls=LineSearchBackTrack::__new__(&LineSearchBackTrack::TypeObject,empty_tuple,empty_dict);
    LineSearchBackTrack::__init__(__ls,empty_tuple,empty_dict);
    Py_DECREF(empty_dict);
    Py_DECREF(empty_tuple);
    f.val<5>().ob=__ls;
    
    
    if(f(args,kwds)==-1) return -1;
    
    
    
    Object* __self=reinterpret_cast<Object*>(self);
    Py_INCREF(f.val<5>().ob);
    __self->ls=reinterpret_cast<LineSearch::Object*>(f.val<5>().ob);
    __self->min=new MinCGDMD(f.val<0>(),f.val<1>(),f.val<2>(),f.val<3>(),f.val<4>(),&(__self->ls->ls));
    __self->xprt=NULL;
    
    return 0;
}
/*--------------------------------------------
 
 --------------------------------------------*/
PyObject* MinCGDMD::__alloc__(PyTypeObject* type,Py_ssize_t)
{
    Object* __self=new Object;
    Py_TYPE(__self)=type;
    Py_REFCNT(__self)=1;
    __self->min=NULL;
    __self->ls=NULL;
    __self->xprt=NULL;
    return reinterpret_cast<PyObject*>(__self);
}
/*--------------------------------------------
 
 --------------------------------------------*/
void MinCGDMD::__dealloc__(PyObject* self)
{
    Object* __self=reinterpret_cast<Object*>(self);
    delete __self->min;
    __self->min=NULL;
    if(__self->ls) Py_DECREF(__self->ls);
    __self->ls=NULL;
    if(__self->xprt) Py_DECREF(__self->xprt);
    __self->xprt=NULL;
    delete __self;
}
/*--------------------------------------------*/
PyTypeObject MinCGDMD::TypeObject={PyObject_HEAD_INIT(NULL)};
/*--------------------------------------------*/
int MinCGDMD::setup_tp()
{
    TypeObject.tp_name="mapp.dmd.min_cg";
    TypeObject.tp_doc=R"---(
    __init__(e_tol=1.0e-8,H_dof=[[False],[False,False],[False,False,False]],affine=False,max_dx=1.0,max_dalpha=0.1,ls=mapp.dmd.ls_bt())
    
    CG minimization algorithm
        
    Parameters
    ----------
    e_tol : double
       Energy tolerance criterion for stopping minimization
    H_dof : symm<bool[dim][dim]>
       Unitcell degrees of freedom during minimization, here dim is the dimension of simulation
    affine : bool
       If set to True atomic displacements would be affine
    max_dx : double
       Maximum displacement of any atom in one step of minimization
    max_dalpha : double
       Maximum change in alpha component of any atom in one step of minimization
    ls : mapp.ls
       Line search method

    Notes
    -----
    Cojugate Gradient (CG) algorithm for minimization.    

    )---";
    
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
    GET_WRAPPER_DOC(TypeObject,__init__)=NULL;
    return ichk;
}
/*--------------------------------------------*/
PyGetSetDef MinCGDMD::getset[]=EmptyPyGetSetDef(9);
/*--------------------------------------------*/
void MinCGDMD::setup_tp_getset()
{
    getset_e_tol(getset[0]);
    getset_H_dof(getset[1]);
    getset_affine(getset[2]);
    getset_max_dx(getset[3]);
    getset_max_dalpha(getset[4]);
    getset_ls(getset[5]);
    getset_ntally(getset[6]);
    getset_export(getset[7]);
}
/*--------------------------------------------*/
PyMethodDef MinCGDMD::methods[]=EmptyPyMethodDef(3);
/*--------------------------------------------*/
void MinCGDMD::setup_tp_methods()
{
    ml_run(methods[0]);
    ml_refine(methods[1]);
}
/*--------------------------------------------
 
 --------------------------------------------*/
void MinCGDMD::getset_max_dalpha(PyGetSetDef& getset)
{
    getset.name=(char*)"max_dalpha";
    getset.doc=(char*)R"---(
    (double) mximum alpha change
    
    Maximum change in alpha component of any atom in one step of minimization
    )---";
    getset.get=[](PyObject* self,void*)->PyObject*
    {
        return var<type0>::build(reinterpret_cast<Object*>(self)->min->max_dx);
    };
    getset.set=[](PyObject* self,PyObject* op,void*)->int
    {
        VarAPI<type0> max_dalpha("max_dalpha");
        max_dalpha.logics[0]=VLogics("gt",0.0);
        int ichk=max_dalpha.set(op);
        if(ichk==-1) return -1;
        reinterpret_cast<Object*>(self)->min->max_dalpha=max_dalpha.val;
        return 0;
    };
}
/*--------------------------------------------
 
 --------------------------------------------*/
void MinCGDMD::getset_export(PyGetSetDef& getset)
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
void MinCGDMD::ml_run(PyMethodDef& tp_methods)
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
            __self->min->pre_run_chk(__atoms,__ff);
        }
        catch(std::string& err_msg)
        {
            PyErr_SetString(PyExc_TypeError,err_msg.c_str());
            return NULL;
        }
        
        __self->min->atoms=__atoms;
        __self->min->ff=__ff;
        __self->min->xprt=__xprt;
        
        try
        {
            __self->min->init();
        }
        catch(std::string& err_msg)
        {
            __self->min->xprt=NULL;
            __self->min->ff=NULL;
            __self->min->atoms=NULL;
            PyErr_SetString(PyExc_TypeError,err_msg.c_str());
            return NULL;
        }
        
        __self->min->run(f.val<1>());
        
        __self->min->fin();
        
        __self->min->xprt=NULL;
        __self->min->ff=NULL;
        __self->min->atoms=NULL;
        
        Py_RETURN_NONE;
    });
    
    tp_methods.ml_doc=(char*)R"---(
    run(atoms,max_nsteps)
   
    Execute minimization
    
    This method starts the energy minimization for a given atoms object and maximum number of steps.
    
    Parameters
    ----------
    atoms : mapp.md.atoms
        System of interest
    max_nsteps : int
        Maximum number of steps to achieve energy minimization
        
    Returns
    -------
    None

    )---";
}



/*--------------------------------------------
 
 --------------------------------------------*/
void MinCGDMD::ml_refine(PyMethodDef& tp_methods)
{
    tp_methods.ml_flags=METH_VARARGS | METH_KEYWORDS;
    tp_methods.ml_name="refine";
    tp_methods.ml_meth=(PyCFunction)(PyCFunctionWithKeywords)(
    [](PyObject* self,PyObject* args,PyObject* kwds)->PyObject*
    {
        Object* __self=reinterpret_cast<Object*>(self);
        FuncAPI<OP<AtomsDMD>,int,int> f("refine",{"atoms","m","max_nsteps"});
        f.logics<1>()[0]=VLogics("gt",0);
        f.logics<2>()[0]=VLogics("ge",0);
        if(f(args,kwds)) return NULL;
        
        AtomsDMD* __atoms=reinterpret_cast<AtomsDMD::Object*>(f.val<0>().ob)->atoms;
        ForceFieldDMD* __ff=reinterpret_cast<AtomsDMD::Object*>(f.val<0>().ob)->ff;
        ExportDMD* __xprt=__self->xprt==NULL ? NULL:__self->xprt->xprt;
        
        __self->min->atoms=__atoms;
        __self->min->ff=__ff;
        __self->min->xprt=__xprt;
        
        
        
        __self->min->refine(f.val<1>(),f.val<2>());
        
        
        __self->min->xprt=NULL;
        __self->min->ff=NULL;
        __self->min->atoms=NULL;
        
        Py_RETURN_NONE;
    });
    
    tp_methods.ml_doc=(char*)R"---(
    run(atoms,max_nsteps)
   
    Execute minimization
    
    This method starts the energy minimization for a given atoms object and maximum number of steps.
    
    Parameters
    ----------
    atoms : mapp.md.atoms
        System of interest
    max_nsteps : int
        Maximum number of steps to achieve energy minimization
        
    Returns
    -------
    None

    )---";
}
