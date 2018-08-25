#include "dae.h"
#include "MAPP.h"
#include "atoms_dmd.h"
#include "dynamic_dmd.h"
#include "ff_dmd.h"
#include "memory.h"
#ifdef MINCG_W_NEWTON
#include "min_cg_dmd.h"
#endif
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
S_dof{DESIG2(__dim__,__dim__,false)},
S{DESIG2(__dim__,__dim__,NAN)},
max_nnewton_iters(5),
max_ngmres_iters(5)
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
    
    if(chng_box && !__atoms->x_dof->is_empty())
    {
        bool* dof=__atoms->x_dof->begin();
        int __dof_lcl[__dim__]{DESIG(__dim__,0)};
        for(int i=0;i<__atoms->natms_lcl;i++,dof+=__dim__)
            Algebra::Do<__dim__>::func([&dof,&__dof_lcl](int i){ if(!dof[i]) __dof_lcl[i]=1;});
        
        int __dof[__dim__]{DESIG(__dim__,0)};
        MPI_Allreduce(__dof_lcl,__dof,__dim__,MPI_INT,MPI_MAX,__atoms->world);
        std::string err_msg=std::string();
        for(int i=0;i<__dim__;i++)
            for(int j=i;j<__dim__;j++)
                if(S_dof[i][j] && __dof[i])
                {
                    /*
                     if(!err_msg.empty()) err_msg+="\n";
                     err_msg+="cannot impose stress component ["+Print::to_string(i)+"]["+Print::to_string(j)
                     +"] while any of the atoms do not have degree freedom in "+Print::to_string(i)
                     +" direction";
                     */
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
    nerr_mins=0;
    c_dim=atoms->c_dim;
    ff->c_d->fill();
#ifdef MINCG_W_NEWTON
    ls=new LineSearchBrent();
    bool __H_dof[__dim__][__dim__];
    for(int i=0;i<__dim__;i++) for(int j=0;j<__dim__;j++) __H_dof[i][j]=false;
    min= new MinCGDMD(1.0e-9,__H_dof,false,1.0,0.1,ls);
    min->atoms=atoms;
    min->ff=ff;
    min->init();
    min->ntally=0;
    dynamic=min->dynamic;
#else
    dynamic=new DynamicDMD(atoms,ff,chng_box,{atoms->c_dof},{atoms->x_dof,atoms->alpha_dof},{});
    dynamic->init();
#endif
    
    ff->calc_ndof();
    a_tol_sqrt_nc_dof=a_tol*sqrt(static_cast<type0>(ff->nc_dof));
    int n=ff->nx_dof+ff->nalpha_dof;
    Algebra::DoLT<__dim__>::func([this,&n](int i,int j)
    {if(!std::isnan(S[i][j])) ++n;});
    sqrt_nx_nalpha_nS_dof=sqrt(static_cast<type0>(n));
    a_tol_sqrt_nx_nalpha_nS_dof=a_tol*sqrt_nx_nalpha_nS_dof;
}
/*--------------------------------------------
 
 --------------------------------------------*/
void DAE::fin()
{
#ifdef MINCG_W_NEWTON
    min->fin();
    min->ff=NULL;
    min->atoms=NULL;
    delete min;
    min=NULL;
    delete ls;
    ls=NULL;
#else
    dynamic->fin();
    delete dynamic;
#endif
    dynamic=NULL;
}
/*--------------------------------------------
 
 --------------------------------------------*/
type0 DAE::calc_err()
{
    if(!chng_box) return ff->err/sqrt_nx_nalpha_nS_dof;
    type0 err_sq=0.0;
    type0 (&S_fe)[__dim__][__dim__]=atoms->S_fe;
    Algebra::DoLT<__dim__>::func([&S_fe,&err_sq,this](int i,int j)
    {if(!std::isnan(S[i][j])) err_sq+=(S[i][j]-S_fe[i][j])*(S[i][j]-S_fe[i][j]);});
    err_sq*=(atoms->vol)*(atoms->vol);
    err_sq+=(ff->err)*(ff->err);
    return sqrt(err_sq)/sqrt_nx_nalpha_nS_dof;
}
/*--------------------------------------------
 
 --------------------------------------------*/
void DAE::min_error_true()
{
#ifdef MINCG_W_NEWTON
    int __step=atoms->step;
    min->run(10000);
    atoms->step=__step;
#endif
    
    VecTens<type0,2> x(atoms,true,atoms->H,atoms->x,atoms->alpha);
    VecTens<type0,2> f(atoms,true,ff->F_H,ff->f,ff->f_alpha);
    VecTens<type0,2> h(atoms,true,__dim__,c_dim);

#ifndef NEW_UPDTAE
    vec* uvecs[2];
    uvecs[0]=atoms->x;
    uvecs[1]=atoms->alpha;
#endif
    type0 norm,res,res_sq;
    
    __GMRES__<VecTens<type0,2>> gmres(max_ngmres_iters,atoms,true,__dim__,c_dim);
/*
    auto J=[this](VecTens<type0,2>& dx,VecTens<type0,2>& Jdx)->void
    {
#ifdef OLD_UPDATE
        dynamic->update(dx.vecs[0],dx.A);
        dynamic->update(dx.vecs[1]);
#else
        dynamic->update(dx.A,dx.vecs[0],dx.vecs[1]);
#endif
        
        type0* __vec=ff->J(dx.vecs[0],dx.vecs[1],Jdx.vecs[0],Jdx.vecs[1]);
        Algebra::DyadicV_2_MLT(__vec,Jdx.A);
        
        type0 dlog_vol=0.0;
        type0 (&H)[__dim__][__dim__]=atoms->H;
    
        Algebra::Do<__dim__>::func([&H,&dx,&dlog_vol](int i)
        {
            dlog_vol+=dx.A[i][i]/H[i][i];
        });

        type0 tmp=dlog_vol*atoms->vol;
        Algebra::DoLT<__dim__>::func([&Jdx,&tmp,this](int i,int j)
        {
            if(!S_dof[i][j]) Jdx.A[i][j]=0.0;
            else Jdx.A[i][j]-=S[i][j]*tmp;
        });
    };
*/
    auto J=[this](VecTens<type0,2>& dx,VecTens<type0,2>& Jdx)->void
    {
#ifdef OLD_UPDATE
        dynamic->update(dx.vecs[0]);
        dynamic->update(dx.vecs[1]);
#else
        dynamic->update(dx.vecs[0],dx.vecs[1]);
#endif
        const int n=atoms->natms_lcl+atoms->natms_ph;
        type0* __dx=dx.vecs[0]->begin();
        type0* __x=atoms->x->begin();
        for(int i=0;i<n;i++,__dx+=__dim__,__x+=__dim__)
            Algebra::V_mul_MLT_add_in(__x,dx.A,__dx);
        
        
        type0* __vec=ff->J(dx.vecs[0],dx.vecs[1],Jdx.vecs[0],Jdx.vecs[1]);
        Algebra::DyadicV_2_MLT(__vec,Jdx.A);
        
        type0 dlog_vol=0.0;
        type0 (&H)[__dim__][__dim__]=atoms->H;
    
        Algebra::Do<__dim__>::func([&H,&dx,&dlog_vol](int i)
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
    
    res_sq=ff->prepJ_n_res(f.vecs[0],f.vecs[1]);
    type0 vol_neg=-atoms->vol;
    Algebra::DoLT<__dim__>::func([&res_sq,&f,this,&vol_neg](int i,int j)
     {
         if(!std::isnan(S[i][j]))
         {
             f.A[i][j]=ff->F_H[i][j]-S[i][j]*vol_neg;
             res_sq+=f.A[i][j]*f.A[i][j];
         }
         else
             f.A[i][j]=0.0;
     });
    
    
    res=sqrt(res_sq);
    
    type0 r;
    int istep=0;
    for(;istep<max_nnewton_iters && res/a_tol_sqrt_nx_nalpha_nS_dof>1.0;istep++)
    {
        gmres.solve(J,f,0.005*a_tol_sqrt_nx_nalpha_nS_dof,norm,h);
        
        const int n0=atoms->natms_lcl;
        type0* __x=atoms->x->begin();
        type0* __dx=h.vecs[0]->begin();
        for(int i=0;i<n0;i++,__dx+=__dim__,__x+=__dim__)
            Algebra::V_mul_MLT_add_in(__x,h.A,__dx);
        Algebra::MLT_mul_MLT(atoms->H,h.A,h.A);
        
        const int n=atoms->natms_lcl*c_dim;
        type0* alpha_vec=atoms->alpha->begin();
        type0*  halpha_vec=h.vecs[1]->begin();
        type0 r_lcl=1.0,tmp;
        for(int i=0;i<n;i++)
        {
            tmp=alpha_vec[i]+r_lcl*halpha_vec[i];
            if(tmp<0.0)
            {
                r_lcl=-alpha_vec[i]/halpha_vec[i];
                while(alpha_vec[i]+r_lcl*halpha_vec[i]<=0.0)
                    r_lcl=nextafter(r_lcl,0.0);
            }
        }
        MPI_Allreduce(&r_lcl,&r,1,Vec<type0>::MPI_T,MPI_MIN,atoms->world);
        
        if(r==1.0)
            x+=h;
        else
        {
            r*=0.5;
            x+=r*h;
        }
        type0 max_alpha_lcl=0.0;
        
        
        type0* c_vec=atoms->c->begin();
        for(int i=0;i<n;i++)
            if(c_vec[i]>=0.0) max_alpha_lcl=MAX(max_alpha_lcl,alpha_vec[i]);
        MPI_Allreduce(&max_alpha_lcl,&atoms->max_alpha,1,Vec<type0>::MPI_T,MPI_MAX,atoms->world);
        
        atoms->update_H();
#ifdef OLD_UPDATE
        dynamic->update(uvecs,2);
#else
        dynamic->update<true,true>();
#endif
        
        
        res_sq=ff->prepJ_n_res(f.vecs[0],f.vecs[1]);
        type0 vol_neg=-atoms->vol;
        
        Algebra::DoLT<__dim__>::func([&res_sq,&f,this,&vol_neg](int i,int j)
         {
             if(!std::isnan(S[i][j]))
             {
                 f.A[i][j]=ff->F_H[i][j]-S[i][j]*vol_neg;
                 res_sq+=f.A[i][j]*f.A[i][j];
             }
             else
                 f.A[i][j]=0.0;
         });

        res=sqrt(res_sq);
        
    }
    if(istep) nerr_mins++;
}
/*--------------------------------------------
 
 --------------------------------------------*/
void DAE::min_error_false()
{
#ifdef MINCG_W_NEWTON
    int __step=atoms->step;
    min->run(10000);
    atoms->step=__step;
#endif
    
    VecTens<type0,2> x(atoms,false,atoms->H,atoms->x,atoms->alpha);
    VecTens<type0,2> f(atoms,false,ff->F_H,ff->f,ff->f_alpha);
    VecTens<type0,2> h(atoms,false,__dim__,c_dim);

#ifndef NEW_UPDTAE
    vec* uvecs[2];
    uvecs[0]=atoms->x;
    uvecs[1]=atoms->alpha;
#endif
    type0 norm,res,res_sq;
    
    __GMRES__<VecTens<type0,2>> gmres(max_ngmres_iters,atoms,false,__dim__,c_dim);
    auto J=[this](VecTens<type0,2>& x,VecTens<type0,2>& Jx)->void
    {
#ifdef OLD_UPDATE
        dynamic->update(x.vecs[0]);
        dynamic->update(x.vecs[1]);
#else
        dynamic->update(x.vecs[0],x.vecs[1]);
#endif
        
        ff->J(x.vecs[0],x.vecs[1],Jx.vecs[0],Jx.vecs[1]);

    };
    
    res_sq=ff->prepJ_n_res(f.vecs[0],f.vecs[1]);
    res=sqrt(res_sq);
    
    type0 r;
    int istep=0;
    for(;istep<max_nnewton_iters && res/a_tol_sqrt_nx_nalpha_nS_dof>1.0;istep++)
    {
        gmres.solve(J,f,0.005*a_tol_sqrt_nx_nalpha_nS_dof,norm,h);
        
        
        const int n=atoms->natms_lcl*c_dim;
        type0* alpha_vec=atoms->alpha->begin();
        type0*  halpha_vec=h.vecs[1]->begin();
        type0 r_lcl=1.0,tmp;
        for(int i=0;i<n;i++)
        {
            tmp=alpha_vec[i]+r_lcl*halpha_vec[i];
            if(tmp<0.0)
            {
                r_lcl=-alpha_vec[i]/halpha_vec[i];
                while(alpha_vec[i]+r_lcl*halpha_vec[i]<=0.0)
                    r_lcl=nextafter(r_lcl,0.0);
            }
        }
        MPI_Allreduce(&r_lcl,&r,1,Vec<type0>::MPI_T,MPI_MIN,atoms->world);

        if(r==1.0)
            x+=h;
        else
        {
            r*=0.5;
            x+=r*h;
        }
        type0 max_alpha_lcl=0.0;
        
        
        type0* c_vec=atoms->c->begin();
        for(int i=0;i<n;i++)
            if(c_vec[i]>=0.0) max_alpha_lcl=MAX(max_alpha_lcl,alpha_vec[i]);
        MPI_Allreduce(&max_alpha_lcl,&atoms->max_alpha,1,Vec<type0>::MPI_T,MPI_MAX,atoms->world);

#ifdef OLD_UPDATE
        dynamic->update(uvecs,2);
#else
        dynamic->update<true,true>();
#endif
        
        
        res_sq=ff->prepJ_n_res(f.vecs[0],f.vecs[1]);
        res=sqrt(res_sq);
        
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
    atoms : mapp.dmd.atoms
        System of interest
    t : double
        Desired time
        
    Returns
    -------
    None

    )---";
}
/*--------------------------------------------
     
 --------------------------------------------*/
#include <iostream>
#include "random.h"
#include "ff_eam_dmd.h"
/*--------------------------------------------
 natms = 250;
 No = 100;
 SetDirectory[NotebookDirectory[]];
 A = Import["data.txt", "Data"];
 Manipulate[
  ListLinePlot[
   {
    Table[{A[[i*No + j, 1]], A[[i*No + j, 2]]}, {j, 1, No}],
    Table[{A[[i*No + j, 1]], A[[i*No + j, 3]]}, {j, 1, No}]
    },
   PlotRange -> All,
   Epilog -> Inset[Graphics[Text[Style[ToString[i], Large]]]]
   ]
  , {i, 0, natms - 1, 1}]
 --------------------------------------------*/
/*
void DAE::ml_Jtest(PyMethodDef& tp_methods)
{
    tp_methods.ml_flags=METH_VARARGS | METH_KEYWORDS;
    tp_methods.ml_name="Jtest";
    tp_methods.ml_meth=(PyCFunction)(PyCFunctionWithKeywords)(
    [](PyObject* self,PyObject* args,PyObject* kwds)->PyObject*
    {
        
        FuncAPI<OP<AtomsDMD>> fff("Jtest",{"atoms"});
        if(fff(args,kwds)) return NULL;
        
        AtomsDMD* atoms=reinterpret_cast<AtomsDMD::Object*>(fff.val<0>().ob)->atoms;
        ForceFieldDMD* ff=reinterpret_cast<AtomsDMD::Object*>(fff.val<0>().ob)->ff;
        
        
        bool chng_box=false;
        int c_dim=atoms->c_dim;
        auto J=[&ff](VecTens<type0,2>& x,VecTens<type0,2>& Jx)->void
        {
            ff->J_timer(x,Jx);
            
        };
        
        Random rand(3541684);
        constexpr int nvecs=100;
        type0 delta=1.0e-9;
        
        VecTens<type0,2> x(atoms,chng_box,atoms->H,atoms->x,atoms->alpha);
        VecTens<type0,2> f0(atoms,chng_box,__dim__,c_dim);
        VecTens<type0,2> f(atoms,chng_box,__dim__,c_dim);
        VecTens<type0,2> h(atoms,chng_box,__dim__,c_dim);
        VecTens<type0,2> Jh(atoms,chng_box,__dim__,c_dim);
        VecTens<type0,2> x0(atoms,chng_box,__dim__,c_dim);
        
        VecTens<type0,2>* dFs=new VecTens<type0,2>[nvecs];
        for(int i=0;i<nvecs;i++)
        {
            dFs[i].~VecTens();
            new (dFs+i) VecTens<type0,2>(atoms,chng_box,__dim__,c_dim);
        }
        
        DynamicDMD* dynamic=new DynamicDMD(atoms,ff,chng_box,{h.vecs[0],h.vecs[1]},
        {x0.vecs[0],x0.vecs[1],f0.vecs[0],f0.vecs[1],Jh.vecs[0],Jh.vecs[1]},{});
        
        for(int i=0;i<nvecs;i++)
        {
            dynamic->add_xchng(dFs[i].vecs[0]);
            dynamic->add_xchng(dFs[i].vecs[1]);
            
        }
        
        
        dynamic->init();
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        x0=x;
        
        vec* uvecs[2];
        uvecs[0]=atoms->x;
        uvecs[1]=atoms->alpha;
        int natms_lcl=atoms->natms_lcl;
        
        
        
        
        
        type0* __h=h.vecs[0]->begin();
        for(int i=0;i<natms_lcl*__dim__;i++)
        {
            
            if(rand.uniform()>0.5)
                __h[i]=rand.uniform()*0.01/(delta*nvecs);
            else
                __h[i]=-rand.uniform()*0.01/(delta*nvecs);

        }
        
        __h=h.vecs[1]->begin();
        for(int i=0;i<natms_lcl*c_dim;i++)
        {
            __h[i]=0.0;
            
            if(rand.uniform()>0.5)
                __h[i]=rand.uniform()*0.01/(delta*nvecs);
            else
                __h[i]=-rand.uniform()*0.01/(delta*nvecs);
            
        }
        
        
        
        
        
        
        ff->prep_timer(f0);
        J(h,Jh);
        
        
        
        for(int ivec=0;ivec<nvecs;ivec++)
        {
            type0 __delta=delta*ivec;
            x=x0+__delta*h;
            
            
            
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
            
            printf("%d\n",ivec);
            
            ff->prep_timer(f.vecs[0],f.vecs[0]);
            
            
            dFs[ivec]=f-f0;
            
        }
        
        
        x=x0;
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
        
        
        
        
        
        
        
        
        
        
        
        natms_lcl=atoms->natms_lcl;
        
        FILE* fp;
        fp=fopen("/Users/sina/Desktop/data_x.txt","w");
        for(int i=0;i<natms_lcl*__dim__;i++)
            for(int ivec=0;ivec<nvecs;ivec++)
                fprintf(fp,"%e\t%e\t%e\n",ivec*delta,dFs[ivec].vecs[0]->begin()[i],-Jh.vecs[0]->begin()[i]*ivec*delta);
        fclose(fp);
        
        
        fp=fopen("/Users/sina/Desktop/data_alpha.txt","w");
        for(int i=0;i<natms_lcl*c_dim;i++)
            for(int ivec=0;ivec<nvecs;ivec++)
                fprintf(fp,"%e\t%e\t%e\n",ivec*delta,dFs[ivec].vecs[1]->begin()[i],-Jh.vecs[1]->begin()[i]*ivec*delta);
        fclose(fp);
        
        dynamic->fin();
        delete dynamic;
        
        Py_RETURN_NONE;
    });
    
    tp_methods.ml_doc=(char*)R"---(


    )---";
}
*/

