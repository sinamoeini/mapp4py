#include "md_nst.h"
#include "xmath.h"
#include "atoms_md.h"
#include "elements.h"
#include "ff_styles.h"
#include "thermo_dynamics.h"
#include "dynamic_md.h"

using namespace MAPP_NS;

/*--------------------------------------------
 
 --------------------------------------------*/
MDNST::MDNST(type0(&__S)[__dim__][__dim__],type0 __T,type0 __dt):
MDNVT(__T,__dt),
MLT0{DESIG2(__dim__,__dim__,0.0)},
MLT1{DESIG2(__dim__,__dim__,0.0)},
MLT2{DESIG2(__dim__,__dim__,0.0)},
S_dev{DESIG2(__dim__,__dim__,0.0)},
B_ref{DESIG2(__dim__,__dim__,0.0)},
V_H{DESIG2(__dim__,__dim__,0.0)},
thermo_baro(__dt/2.0,1000.0*__dt,3,1),
S_dof{DESIG2(__dim__,__dim__,false)},
S{DESIG2(__dim__,__dim__,0.0)},
tau{DESIG2(__dim__,__dim__,0.0)},
nreset(0)
{
    
    Algebra::DoLT<__dim__>::func([this,&__S,&__dt]
    (int i,int j)
    {
        if(std::isnan(__S[i][j]))
        {
            S_dof[i][j]=S_dof[j][i]=false;
            S[i][j]=S[j][i]=NAN;
        }
        else
        {
            S_dof[i][j]=S_dof[j][i]=true;
            S[i][j]=S[j][i]=__S[i][j];
        }
        
        tau[i][j]=tau[j][i]=1000.0*__dt;
        
    });
}
/*--------------------------------------------
 
 --------------------------------------------*/
MDNST::~MDNST()
{
}
/*--------------------------------------------
 
 --------------------------------------------*/
void MDNST::change_dt(type0 __dt)
{
    if(__dt==dt) return;
    dt=__dt;
    dt2=__dt/2;
    
    int L=thermo_part.L;
    int niters=thermo_part.niters;
    type0 t_relax=thermo_part.t_relax;
    
    thermo_part.~ThermostatNHC();
    new (&thermo_part) ThermostatNHC(__dt/2.0,t_relax,L,niters);
    
    
    L=thermo_baro.L;
    niters=thermo_baro.niters;
    t_relax=thermo_baro.t_relax;
    
    thermo_baro.~ThermostatNHC();
    new (&thermo_baro) ThermostatNHC(__dt/2.0,t_relax,L,niters);
}
/*--------------------------------------------
 
 --------------------------------------------*/
void MDNST::update_x_d__x__x_d(type0 xi)
{

    MDMath::calc(Algebra::Tr(V_H)/ndof_part,dt,dt2,V_H,MLT0,MLT1,MLT2);
    type0* RESTRICT f=ff->f->begin();
    type0* RESTRICT x=atoms->x->begin();
    type0* RESTRICT x_d=atoms->x_d->begin();
    elem_type* RESTRICT elem=atoms->elem->begin();
    type0* m=atoms->elements.masses;
    type0 m_i,m_inv;
    const int natms0=atoms->natms_lcl;
    type0 dx_lcl[__dim__]={DESIG(__dim__,0.0)};
    for(int i=0;i<natms0;++i)
    {
        m_inv=1.0/m[*elem];
        MDMath::____NONAME0<__dim__,__dim__>::func(xi,x_d,&MLT0[0][0],m_inv,f,&MLT1[0][0]);
        MDMath::____NONAME1<__dim__,__dim__>::func(x,&V_H[0][0],x_d,v0);
        MDMath::____NONAME2<__dim__,__dim__>::func(v0,&MLT2[0][0]);
        
        Algebra::V_add<__dim__>(v0,x);
        Algebra::V_add<__dim__>(v0,dx_lcl);
        
        f+=__dim__;
        x_d+=__dim__;
        x+=__dim__;
        ++elem;
    }
    
    type0 dx[__dim__]={DESIG(__dim__,0.0)};
    MPI_Allreduce(dx_lcl,dx,__dim__,Vec<type0>::MPI_T,MPI_SUM,atoms->world);
    type0 natms=static_cast<type0>(atoms->natms);
    Algebra::Do<__dim__>::func([&dx,natms](const int i){dx[i]/=natms;});
    x=atoms->x->begin();
    for(int i=0;i<natms0;++i,x+=__dim__)
        Algebra::Do<__dim__>::func([&dx,&x](const int j){x[j]-=dx[j];});
    
    
    Algebra::MLT_mul_MLT(MLT2,V_H,MLT2);
    Algebra::Do<__dim__>::func([this](int i){MLT2[i][i]++;});
    Algebra::MLT_mul_MLT(atoms->H,MLT2,atoms->H);
    atoms->update_H();
    
#ifdef OLD_UPDATE
    dynamic->update(atoms->x);
#else
    dynamic->update<true>();
#endif
    ff->force_calc();
    
    
    f=ff->f->begin();
    x_d=atoms->x_d->begin();
    elem=atoms->elem->begin();
    Algebra::zero<__nvoigt__>(__vec_lcl);
    const int natms1=atoms->natms_lcl;
    for(int i=0;i<natms1;++i)
    {
        m_i=m[*elem];
        m_inv=1.0/m_i;
        
        MDMath::____NONAME0<__dim__,__dim__>::func(x_d,&MLT0[0][0],m_inv,f,&MLT1[0][0]);
        Algebra::DyadicV<__dim__>(m_i,x_d,__vec_lcl);
        
        f+=__dim__;
        x_d+=__dim__;
        ++elem;
    }
    MPI_Allreduce(__vec_lcl,mvv,__nvoigt__,Vec<type0>::MPI_T,MPI_SUM,atoms->world);
    T_part=Algebra::Tr_DyadicV<__dim__>(mvv)/(ndof_part*kB);
}
/*--------------------------------------------
 
 --------------------------------------------*/
void MDNST::update_x_d__x__x_d_w_dof(type0 xi)
{

    MDMath::calc(Algebra::Tr(V_H)/ndof_part,dt,dt2,V_H,MLT0,MLT1,MLT2);
    type0* RESTRICT f=ff->f->begin();
    type0* RESTRICT x=atoms->x->begin();
    type0* RESTRICT x_d=atoms->x_d->begin();
    elem_type* RESTRICT elem=atoms->elem->begin();
    bool* RESTRICT dof=atoms->x_dof->begin();
    
    type0* m=atoms->elements.masses;
    type0 m_i,m_inv;
    const int natms0=atoms->natms_lcl;
    type0 dx_lcl[__dim__]={DESIG(__dim__,0.0)};
    type0 __x_d[__dim__]={DESIG(__dim__,0.0)};
    for(int i=0;i<natms0;++i)
    {
        m_inv=1.0/m[*elem];
        Algebra::V_eq<__dim__>(x_d,__x_d);
        MDMath::____NONAME0<__dim__,__dim__>::func(xi,x_d,&MLT0[0][0],m_inv,f,&MLT1[0][0]);
        Algebra::Do<__dim__>::func([&__x_d,&x_d,&dof](const int j){if(!dof[j]) x_d[j]=__x_d[j];});
        MDMath::____NONAME1<__dim__,__dim__>::func(x,&V_H[0][0],x_d,v0);
        MDMath::____NONAME2<__dim__,__dim__>::func(v0,&MLT2[0][0]);
        
        Algebra::V_add<__dim__>(v0,x);
        Algebra::V_add<__dim__>(v0,dx_lcl);
        
        f+=__dim__;
        x_d+=__dim__;
        x+=__dim__;
        dof+=__dim__;
        ++elem;
    }
    
    type0 dx[__dim__]={DESIG(__dim__,0.0)};
    MPI_Allreduce(dx_lcl,dx,__dim__,Vec<type0>::MPI_T,MPI_SUM,atoms->world);
    type0 natms=static_cast<type0>(atoms->natms);
    Algebra::Do<__dim__>::func([&dx,natms](const int i){dx[i]/=natms;});
    x=atoms->x->begin();
    for(int i=0;i<natms0;++i,x+=__dim__)
        Algebra::Do<__dim__>::func([&dx,&x](const int j){x[j]-=dx[j];});
    
    
    Algebra::MLT_mul_MLT(MLT2,V_H,MLT2);
    Algebra::Do<__dim__>::func([this](int i){MLT2[i][i]++;});
    Algebra::MLT_mul_MLT(atoms->H,MLT2,atoms->H);
    atoms->update_H();

#ifdef OLD_UPDATE
    dynamic->update(atoms->x);
#else
    dynamic->update<true>();
#endif
    ff->force_calc();
    
    
    f=ff->f->begin();
    x_d=atoms->x_d->begin();
    elem=atoms->elem->begin();
    dof=atoms->x_dof->begin();
    Algebra::zero<__nvoigt__>(__vec_lcl);
    const int natms1=atoms->natms_lcl;
    for(int i=0;i<natms1;++i)
    {
        m_i=m[*elem];
        m_inv=1.0/m_i;
        
        Algebra::V_eq<__dim__>(x_d,__x_d);
        MDMath::____NONAME0<__dim__,__dim__>::func(x_d,&MLT0[0][0],m_inv,f,&MLT1[0][0]);
        Algebra::Do<__dim__>::func([&__x_d,&x_d,&dof](const int j){if(!dof[j]) x_d[j]=__x_d[j];});
        Algebra::DyadicV<__dim__>(m_i,x_d,__vec_lcl);
        
        f+=__dim__;
        x_d+=__dim__;
        dof+=__dim__;
        ++elem;
    }
    MPI_Allreduce(__vec_lcl,mvv,__nvoigt__,Vec<type0>::MPI_T,MPI_SUM,atoms->world);
    T_part=Algebra::Tr_DyadicV<__dim__>(mvv)/(ndof_part*kB);
}
/*--------------------------------------------
 pre run check it throw excepctions
 --------------------------------------------*/
void MDNST::pre_run_chk(AtomsMD* atoms,ForceFieldMD* ff)
{
    //check if configuration is loaded
    if(!atoms)
        throw std::string("cannot start md without initial conditions");
    
    //check if force field is loaded
    if(!ff)
        throw std::string("cannot start md without governing equations (force field)");
    
    if(std::isnan(atoms->kB))
        throw std::string("boltzmann constant should be set prior to MD");
    
    //check to see if the H_dof components are consistent with stoms->dof
    
    if(!atoms->x_dof->is_empty())
    {
        bool* dof=atoms->x_dof->begin();
        int __dof_lcl[__dim__]{DESIG(__dim__,0)};
        for(int i=0;i<atoms->natms_lcl;i++,dof+=__dim__)
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
void MDNST::pre_init()
{
    MDNVT::pre_init();
    ndof_baro=0;
    T_baro=0.0;
    Algebra::DoLT<__dim__>::func([this](int i,int j)
    {
        if(S_dof[i][j])
        {
            V_H_prefac[i][j]=1.0/(tau[i][j]*tau[i][j]*ndof_part*kB*T);
            ndof_baro++;
            T_baro+=V_H[i][j]*V_H[i][j]*tau[i][j]*tau[i][j];
            
        }
        else
        {
            S[i][j]=0.0;
            V_H[i][j]=0.0;
        }
        S_dev[i][j]=S[i][j];
        
    });
    
    T_baro*=T*ndof_part*kB/ndof_baro;
    s_hyd=Algebra::Tr(S)/static_cast<type0>(__dim__);
    Algebra::Do<__dim__>::func([this](int i)
    {
        S_dev[i][i]=S[i][i]-s_hyd;
    });
}
/*--------------------------------------------
 
 --------------------------------------------*/
void MDNST::update_V_H()
{
    Algebra::MLT_mul_MLT(B_ref,atoms->H,MLT0);
    Algebra::MUT_mul_MSY_mul_MLT(S_dev,MLT0,MLT1,MLT2);
    
    
    Algebra::Do<__dim__>::func([this](int i)
    {
        if(S_dof[i][i]) V_H[i][i]+=(kB*T_part+atoms->vol*s_hyd)*V_H_prefac[i][i]*dt2;
    });
    
    T_baro=0.0;

    Algebra::DoLT<__dim__>::func([this](const int i,const int j)
    {
        if(S_dof[i][j])
        {
            V_H[i][j]+=(vol_ref*MLT2[i][j]
            -atoms->vol*atoms->S_pe[i][j]
            +mvv[i+j*__dim__-j*(j+1)/2])*V_H_prefac[i][j]*dt2;
            
            T_baro+=V_H[i][j]*V_H[i][j]*tau[i][j]*tau[i][j];
            
        }
    });
    
    T_baro*=ndof_part*T/ndof_baro;

}
/*--------------------------------------------
 
 --------------------------------------------*/
void MDNST::init()
{
    pre_init();
    dynamic=new DynamicMD(atoms,ff,true,{},{atoms->x_d,atoms->x_dof},{});
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
void MDNST::fin()
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
void MDNST::run(int nsteps)
{
    int step=atoms->step;
    
    ff->force_calc();

    
    int nevery_xprt=xprt==NULL ? 0:xprt->nevery;
    if(nevery_xprt) xprt->write(step);
    
    ThermoDynamics thermo(6,"T",T_part,"PE",atoms->pe,
    "S[0][0]",S_part[0][0],
    "S[1][1]",S_part[1][1],
    "S[2][2]",S_part[2][2],
    "S[1][2]",S_part[2][1],
    "S[2][0]",S_part[2][0],
    "S[0][1]",S_part[1][0]);
    
    if(ntally) thermo.init();
    Algebra::DoLT<__dim__>::func([this](const int i,const int j)
    {
        S_part[i][j]=atoms->S_pe[i][j]-mvv[i+j*__dim__-j*(j+1)/2]/atoms->vol;
    });
    
    
    if(ntally) thermo.print(step);
    
    //type0* x_d;
    type0 fac,fac_x_d=1.0;
    
    // store reference configuration
    Algebra::DoLT<__dim__>::func([this](int i,int j){B_ref[i][j]=atoms->B[i][j];});
    vol_ref=atoms->vol;
    
    for(int istep=0;istep<nsteps;istep++)
    {
        // baro thermostat
        fac=thermo_baro(T_baro/T,ndof_baro);
        Algebra::DoLT<__dim__>::func([this,&fac](int i,int j){V_H[i][j]*=fac;});
        T_baro*=fac*fac;

        
        // particle thermostat
        fac_x_d*=fac=thermo_part(T_part/T,ndof_part);
        fac*=fac;
        Algebra::Do<__nvoigt__>::func([this,&fac](int i){mvv[i]*=fac;});
        T_part*=fac;
        
        update_V_H();
        
        update_x_d__x__x_d(fac_x_d);
        
        update_V_H();
        
        // particle thermostat
        fac_x_d=fac=thermo_part(T_part/T,ndof_part);
        fac*=fac;
        Algebra::Do<__nvoigt__>::func([this,&fac](int i){mvv[i]*=fac;});
        T_part*=fac;
        
        // baro thermostat
        fac=thermo_baro(T_baro/T,ndof_baro);
        Algebra::DoLT<__dim__>::func([this,&fac](int i,int j){V_H[i][j]*=fac;});
        T_baro*=fac*fac;
        
        Algebra::DoLT<__dim__>::func([this](const int i,const int j)
        {
            S_part[i][j]=atoms->S_pe[i][j]-mvv[i+j*__dim__-j*(j+1)/2]/atoms->vol;
        });
        
        
        
        if(nreset && (istep+1)%nreset==0)
        {
            Algebra::DoLT<__dim__>::func([this](int i,int j){B_ref[i][j]=atoms->B[i][j];});
            vol_ref=atoms->vol;
        }
        
        if(ntally && (istep+1)%ntally==0) thermo.print(step+istep+1);
        if(nevery_xprt && (istep+1)%nevery_xprt==0) xprt->write(step+istep+1);
    }
    
    if(ntally && nsteps%ntally) thermo.print(step+nsteps);
    if(nevery_xprt && nsteps%nevery_xprt) xprt->write(step+nsteps);
    
    update_x_d_final(fac_x_d);
    
    if(ntally) thermo.fin();
    
    atoms->step+=nsteps;
}




















/*------------------------------------------------------------------------------------------------------------------------------------
 
 ------------------------------------------------------------------------------------------------------------------------------------*/
PyObject* MDNST::__new__(PyTypeObject* type,PyObject* args,PyObject* kwds)
{
    Object* __self=reinterpret_cast<Object*>(type->tp_alloc(type,0));
    PyObject* self=reinterpret_cast<PyObject*>(__self);
    return self;
}
/*--------------------------------------------
 
 --------------------------------------------*/
int MDNST::__init__(PyObject* self,PyObject* args,PyObject* kwds)
{
    FuncAPI<symm<type0[__dim__][__dim__]>,type0,type0> f("__init__",{"S","T","dt"});
    
    f.logics<1>()[0]=VLogics("gt",0.0);
    f.logics<2>()[0]=VLogics("gt",0.0);
    
    if(f(args,kwds)==-1) return -1;
    Object* __self=reinterpret_cast<Object*>(self);
    __self->md=new MDNST(f.val<0>(),f.val<1>(),f.val<2>());
    __self->xprt=NULL;
    
    return 0;
}
/*--------------------------------------------
 
 --------------------------------------------*/
PyObject* MDNST::__alloc__(PyTypeObject* type,Py_ssize_t)
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
void MDNST::__dealloc__(PyObject* self)
{
    Object* __self=reinterpret_cast<Object*>(self);
    delete __self->md;
    __self->md=NULL;
    if(__self->xprt) Py_DECREF(__self->xprt);
    __self->xprt=NULL;
    delete __self;
}
/*--------------------------------------------*/
PyTypeObject MDNST::TypeObject ={PyObject_HEAD_INIT(NULL)};
/*--------------------------------------------*/
int MDNST::setup_tp()
{
    TypeObject.tp_name="mapp4py.md.nst";
    TypeObject.tp_doc=R"---(
    __init__(S,T,dt)
    
    :math:`N\mathbf{S}T` ensemble
    
    Molecular dynamics of isothermal-isostress ensemble
        
    Parameters
    ----------
    S : symm<double[dim][dim]>
       External stress imposed on system, here dim is the dimension of simulation
    T : double
       Temperature of the ensemble
    dt : double
       Time step for simulation
    
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
PyGetSetDef MDNST::getset[]=EmptyPyGetSetDef(7);
/*--------------------------------------------*/
void MDNST::setup_tp_getset()
{
    getset_S(getset[0]);
    getset_niters_s(getset[1]);
    getset_L_s(getset[2]);
    getset_t_relax_s(getset[3]);
    getset_tau(getset[4]);
    getset_nreset(getset[5]);
}
/*--------------------------------------------
 
 --------------------------------------------*/
void MDNST::getset_L_s(PyGetSetDef& getset)
{
    getset.name=(char*)"L_s";
    getset.doc=(char*)R"---(
    (int) length of stress NHC
    
    Number of links in NHC of stress thermostat, see :ref:`nhc-ref`
    )---";
    getset.get=[](PyObject* self,void*)->PyObject*
    {
        int L=reinterpret_cast<Object*>(self)->md->thermo_baro.L;
        return var<int>::build(L);
    };
    getset.set=[](PyObject* self,PyObject* op,void*)->int
    {
        VarAPI<int> L_s("L_s");
        L_s.logics[0]=VLogics("gt",0);
        int ichk=L_s.set(op);
        if(ichk==-1) return -1;
        if(reinterpret_cast<Object*>(self)->md->thermo_baro.L==L_s.val)
            return 0;
        
        ThermostatNHC& thermo_baro=reinterpret_cast<Object*>(self)->md->thermo_baro;
        int niters=thermo_baro.niters;
        type0 t_relax=thermo_baro.t_relax;
        type0 __dt=reinterpret_cast<Object*>(self)->md->dt2;
        thermo_baro.~ThermostatNHC();
        new (&thermo_baro) ThermostatNHC(__dt,t_relax,L_s.val,niters);
        return 0;
    };
}
/*--------------------------------------------
 
 --------------------------------------------*/
void MDNST::getset_niters_s(PyGetSetDef& getset)
{
    getset.name=(char*)"niters_s";
    getset.doc=(char*)R"---(
    (int) number of iterations in NHC for stress
    
    Number of iterations in stress thermostat per each step of evolution, see :ref:`nhc-ref`
    )---";
    getset.get=[](PyObject* self,void*)->PyObject*
    {
        int niters=reinterpret_cast<Object*>(self)->md->thermo_baro.niters;
        return var<int>::build(niters);
    };
    getset.set=[](PyObject* self,PyObject* op,void*)->int
    {
        VarAPI<int> niters_s("niters_s");
        niters_s.logics[0]=VLogics("gt",0);
        int ichk=niters_s.set(op);
        if(ichk==-1) return -1;
        if(reinterpret_cast<Object*>(self)->md->thermo_baro.niters==niters_s.val)
            return 0;
        
        ThermostatNHC& thermo_baro=reinterpret_cast<Object*>(self)->md->thermo_baro;
        int L=thermo_baro.L;
        type0 t_relax=thermo_baro.t_relax;
        type0 __dt=reinterpret_cast<Object*>(self)->md->dt2;
        thermo_baro.~ThermostatNHC();
        new (&thermo_baro) ThermostatNHC(__dt,t_relax,L,niters_s.val);
        return 0;
    };
}
/*--------------------------------------------
 
 --------------------------------------------*/
void MDNST::getset_t_relax_s(PyGetSetDef& getset)
{
    getset.name=(char*)"t_relax_s";
    getset.doc=(char*)R"---(
    (double) NHC relaxation parameter for stress
    
    NHC relaxation parameter (unit of time) for stress thermostat. Roughly equivalent to relaxation time of stress with thermal bath, see :ref:`nhc-ref`
    )---";
    getset.get=[](PyObject* self,void*)->PyObject*
    {
        type0 t_relax=reinterpret_cast<Object*>(self)->md->thermo_baro.t_relax;
        return var<type0>::build(t_relax);
    };
    getset.set=[](PyObject* self,PyObject* op,void*)->int
    {
        VarAPI<type0> t_relax_s("t_relax_s");
        t_relax_s.logics[0]=VLogics("gt",0.0);
        int ichk=t_relax_s.set(op);
        if(ichk==-1) return -1;
        if(reinterpret_cast<Object*>(self)->md->thermo_baro.t_relax==t_relax_s.val)
            return 0;
        
        ThermostatNHC& thermo_baro=reinterpret_cast<Object*>(self)->md->thermo_baro;
        int L=thermo_baro.L;
        int niters=thermo_baro.niters;
        type0 __dt=reinterpret_cast<Object*>(self)->md->dt2;
        thermo_baro.~ThermostatNHC();
        new (&thermo_baro) ThermostatNHC(__dt,t_relax_s.val,L,niters);
        return 0;
    };
}
/*--------------------------------------------
 
 --------------------------------------------*/
void MDNST::getset_S(PyGetSetDef& getset)
{
    getset.name=(char*)"S";
    getset.doc=(char*)R"---(
    (symm<double[dim][dim]>) external stress tensor
    
    External stress imposed on system, here dim is the dimension of simulation
    )---";
    getset.get=[](PyObject* self,void*)->PyObject*
    {
        return var<symm<type0[__dim__][__dim__]>>::build(reinterpret_cast<Object*>(self)->md->S);
    };
    getset.set=[](PyObject* self,PyObject* op,void*)->int
    {
        VarAPI<symm<type0[__dim__][__dim__]>> S("S");
        int ichk=S.set(op);
        if(ichk==-1) return -1;
        
        bool (&__S_dof)[__dim__][__dim__]=reinterpret_cast<Object*>(self)->md->S_dof;
        type0 (&__S)[__dim__][__dim__]=reinterpret_cast<Object*>(self)->md->S;
        type0 (&__tau)[__dim__][__dim__]=reinterpret_cast<Object*>(self)->md->tau;
        type0 __dt=reinterpret_cast<Object*>(self)->md->dt;
        Algebra::DoLT<__dim__>::func([&__S_dof,&__S,&S,&__tau,&__dt](int i,int j)
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
                if(std::isnan(__tau[i][j]))
                    __tau[i][j]=__tau[j][i]=1000.0*__dt;
            }
        });
        
        return 0;
    };
}
/*--------------------------------------------
 
 --------------------------------------------*/
void MDNST::getset_tau(PyGetSetDef& getset)
{
    getset.name=(char*)"tau";
    getset.doc=(char*)R"---(
    (symm<double[dim][dim]>) relaxation parameter for stress
    
    Relaxation parameter (unit of time) for stress. Roughly equivalent to relaxation time of internal stress with external stress, see :ref:`nst-ref`
    )---";
    getset.get=[](PyObject* self,void*)->PyObject*
    {
        return var<symm<type0[__dim__][__dim__]>>::build(reinterpret_cast<Object*>(self)->md->tau);
    };
    getset.set=[](PyObject* self,PyObject* op,void*)->int
    {
        VarAPI<symm<type0[__dim__][__dim__]>> tau("tau");
        tau.logics[0]=VLogics("gt",0.0);
        int ichk=tau.set(op);
        if(ichk==-1) return -1;
        
        type0 (&__tau)[__dim__][__dim__]=reinterpret_cast<Object*>(self)->md->tau;
        Algebra::DoLT<__dim__>::func([&__tau,&tau](int i,int j)
        {
            __tau[i][j]=__tau[j][i]=tau.val[i][j];
        });
        
        return 0;
    };
}
/*--------------------------------------------
 
 --------------------------------------------*/
void MDNST::getset_nreset(PyGetSetDef& getset)
{
    getset.name=(char*)"nreset";
    getset.doc=(char*)R"---(
    (int) configuration reset period
    
    Number of steps to reset initial configuration. If set to 0 the very inital configuration would be reference
    )---";
    getset.get=[](PyObject* self,void*)->PyObject*
    {
        int niters=reinterpret_cast<Object*>(self)->md->nreset;
        return var<int>::build(niters);
    };
    getset.set=[](PyObject* self,PyObject* op,void*)->int
    {
        VarAPI<int> nreset("nreset");
        nreset.logics[0]=VLogics("ge",0);
        int ichk=nreset.set(op);
        if(ichk==-1) return -1;
        reinterpret_cast<Object*>(self)->md->nreset=nreset.val;
        return 0;
    };
}






























/*--------------------------------------------
 
 --------------------------------------------*/
type0 MDMath::f(type0 x)
{
    if(fabs(x)<0.015625)
        return 1.0+x/2.0+x*x/6.0+x*x*x/24.0;
    return (exp(x)-1.0)/x;
}
/*--------------------------------------------
 
 --------------------------------------------*/
type0 MDMath::df(type0 x)
{
    if(fabs(x)<0.015625)
        return 0.5+x/3.0+x*x/8.0+x*x*x/30.0;
    type0 exp_x=exp(x);
    return (exp_x-(exp_x-1.0)/x)/x;
}
/*--------------------------------------------
 
 --------------------------------------------*/
type0 MDMath::ddf(type0 x)
{
    if(fabs(x)<0.015625)
        return 1.0/3.0+x/4.0+x*x/10.0+x*x*x/36.0;
    type0 exp_x=exp(x);
    return (exp_x-2.0*((exp_x-(exp_x-1.0)/x)/x))/x;
}
/*--------------------------------------------
 
 --------------------------------------------*/
type0 MDMath::dddf(type0 x)
{
    if(fabs(x)<0.015625)
        return 1.0/4.0+x/5.0+x*x/12.0+x*x*x/42.0;
    type0 exp_x=exp(x);
    return (exp_x-3.0*(exp_x-2.0*((exp_x-(exp_x-1.0)/x)/x))/x)/x;
}
/*--------------------------------------------
 
 --------------------------------------------*/
void MDMath::calc(const type0 a,const type0& dt,const type0& dt2,type0(&A)[3][3],type0(&EXP_A)[3][3],type0(&F_A0)[3][3],type0(&F_A1)[3][3])
{

    auto div_dif_exp=[](type0& x,type0& y,type0& fx,type0& fy)->type0
    {
        if(fabs(x-y)<sqrt_eps)
            return exp(0.5*(x+y));
        return (fx-fy)/(x-y);
    };
    auto div_dif_f=[](type0& x,type0& y,type0& fx,type0& fy)->type0
    {
        if(fabs(x-y)<sqrt(std::numeric_limits<type0>::epsilon()*f(0.5*(x+y))/ddf(0.5*(x+y))))
            return df(0.5*(x+y));
        return (fx-fy)/(x-y);
    };
    

    type0 x=-dt2*(A[0][0]+a),y=-dt2*(A[1][1]+a),z=-dt2*(A[2][2]+a);
    type0 __x=dt*A[0][0],__y=dt*A[1][1],__z=dt*A[2][2];
    F_A0[0][0]=f(x);
    F_A0[1][1]=f(y);
    F_A0[2][2]=f(z);
    F_A0[2][1]=F_A0[1][0]=F_A0[2][0]=0.0;
    
    F_A1[0][0]=f(__x);
    F_A1[1][1]=f(__y);
    F_A1[2][2]=f(__z);
    F_A1[2][1]=F_A1[1][0]=F_A1[2][0]=0.0;
    
    EXP_A[0][0]=exp(x);
    EXP_A[1][1]=exp(y);
    EXP_A[2][2]=exp(z);
    EXP_A[2][1]=EXP_A[1][0]=EXP_A[2][0]=0.0;
    
    
    
    type0 f_yz=0.0,f_zx=0.0,f_xy=0.0,__f_yz=0.0,__f_zx=0.0,__f_xy=0.0,g_yz=0.0,g_zx=0.0,g_xy=0.0;
    if(A[2][1]!=0.0)
    {
        f_yz=div_dif_f(y,z,F_A0[1][1],F_A0[2][2]);
        __f_yz=div_dif_f(__y,__z,F_A1[1][1],F_A1[2][2]);
        g_yz=div_dif_exp(y,z,EXP_A[1][1],EXP_A[2][2]);
    }
    if(A[2][0]!=0.0)
    {
        f_zx=div_dif_f(z,x,F_A0[2][2],F_A0[0][0]);
        __f_zx=div_dif_f(__z,__x,F_A1[2][2],F_A1[0][0]);
        g_zx=div_dif_exp(z,x,EXP_A[2][2],EXP_A[0][0]);
    }
    if(A[1][0]!=0.0)
    {
        f_xy=div_dif_f(x,y,F_A0[0][0],F_A0[1][1]);
        __f_xy=div_dif_f(__x,__y,F_A1[0][0],F_A1[1][1]);
        g_xy=div_dif_exp(x,y,EXP_A[0][0],EXP_A[1][1]);
    }
    
    
    F_A0[0][0]*=dt2;
    F_A0[1][1]*=dt2;
    F_A0[2][2]*=dt2;
    F_A0[2][1]=-A[2][1]*dt2*dt2*f_yz;
    F_A0[2][0]=-A[2][0]*dt2*dt2*f_zx;
    F_A0[1][0]=-A[1][0]*dt2*dt2*f_xy;
    
    
    F_A1[0][0]*=dt;
    F_A1[1][1]*=dt;
    F_A1[2][2]*=dt;
    F_A1[2][1]=A[2][1]*dt*dt*__f_yz;
    F_A1[2][0]=A[2][0]*dt*dt*__f_zx;
    F_A1[1][0]=A[1][0]*dt*dt*__f_xy;
    
    EXP_A[2][1]=-A[2][1]*dt2*g_yz;
    EXP_A[2][0]=-A[2][0]*dt2*g_zx;
    EXP_A[1][0]=-A[1][0]*dt2*g_xy;
    
    
    if(A[1][0]==0.0 || A[2][1]==0.0)
        return;
    auto div_dif_f3=[](type0 x,type0 y,type0 z,type0& f_yz,type0& f_zx,type0& f_xy)->type0
    {
        type0 __ddf=ddf((x+y+z)/3.0);
        type0 lim=sqrt(std::numeric_limits<type0>::epsilon()*f((x+y+z)/3.0)/__ddf);
        
        if(fabs(x-z)>lim)
            return (f_xy-f_yz)/(x-z);
        
        if(fabs(y-x)>lim)
            return (f_yz-f_zx)/(y-x);
        
        if(fabs(z-y)>lim)
            return (f_zx-f_xy)/(z-y);
        
        return __ddf;
    };
    
    auto div_dif_exp3=[](type0 x,type0 y,type0 z,type0& f_yz,type0& f_zx,type0& f_xy)->type0
    {
        if(fabs(x-z)>sqrt_eps)
            return (f_xy-f_yz)/(x-z);
        
        if(fabs(y-x)>sqrt_eps)
            return (f_yz-f_zx)/(y-x);
        
        if(fabs(z-y)>sqrt_eps)
            return (f_zx-f_xy)/(z-y);
        
        return exp((x+y+z)/3.0);
    };
    
    F_A0[2][0]+=dt2*dt2*dt2*A[2][1]*A[1][0]*div_dif_f3(x,y,z,f_yz,f_zx,f_xy);
    F_A1[2][0]+=dt*dt*dt*A[2][1]*A[1][0]*div_dif_f3(__x,__y,__z,__f_yz,__f_zx,__f_xy);
    EXP_A[2][0]+=dt2*dt2*A[2][1]*A[1][0]*div_dif_exp3(x,y,z,g_yz,g_zx,g_xy);
}
/*--------------------------------------------
 
 --------------------------------------------*/
void MDMath::calc(const type0 a,const type0& dt,const type0& dt2,type0(&A)[2][2],type0(&EXP_A)[2][2],type0(&F_A0)[2][2],type0(&F_A1)[2][2])
{

    auto div_dif_exp=[](type0& x,type0& y,type0& fx,type0& fy)->type0
    {
        if(fabs(x-y)<sqrt_eps)
            return exp(0.5*(x+y));
        return (fx-fy)/(x-y);
    };
    auto div_dif_f=[](type0& x,type0& y,type0& fx,type0& fy)->type0
    {
        if(fabs(x-y)<sqrt(std::numeric_limits<type0>::epsilon()*f(0.5*(x+y))/ddf(0.5*(x+y))))
            return df(0.5*(x+y));
        return (fx-fy)/(x-y);
    };
    

    type0 x=-dt2*(A[0][0]+a),y=-dt2*(A[1][1]+a);
    type0 __x=dt*A[0][0],__y=dt*A[1][1];
    F_A0[0][0]=f(x);
    F_A0[1][1]=f(y);
    F_A0[1][0]=0.0;
    
    F_A1[0][0]=f(__x);
    F_A1[1][1]=f(__y);
    F_A1[1][0]=0.0;
    
    EXP_A[0][0]=exp(x);
    EXP_A[1][1]=exp(y);
    EXP_A[1][0]=0.0;
    
    
    
    type0 f_xy=0.0,__f_xy=0.0,g_xy=0.0;
    if(A[1][0]!=0.0)
    {
        f_xy=div_dif_f(x,y,F_A0[0][0],F_A0[1][1]);
        __f_xy=div_dif_f(__x,__y,F_A1[0][0],F_A1[1][1]);
        g_xy=div_dif_exp(x,y,EXP_A[0][0],EXP_A[1][1]);
    }
    
    
    F_A0[0][0]*=dt2;
    F_A0[1][1]*=dt2;
    F_A0[1][0]=-A[1][0]*dt2*dt2*f_xy;
    
    
    F_A1[0][0]*=dt;
    F_A1[1][1]*=dt;
    F_A1[1][0]=A[1][0]*dt*dt*__f_xy;

    EXP_A[1][0]=-A[1][0]*dt2*g_xy;
}
/*--------------------------------------------
 
 --------------------------------------------*/
void MDMath::calc(const type0 a,const type0& dt,const type0& dt2,type0(&A)[1][1],type0(&EXP_A)[1][1],type0(&F_A0)[1][1],type0(&F_A1)[1][1])
{
    F_A0[0][0]=f(-dt2*(A[0][0]+a))*dt2;
    F_A1[0][0]=f(dt*A[0][0])*dt;
    EXP_A[0][0]=exp(-dt2*(A[0][0]+a));
}
/*--------------------------------------------
 
 --------------------------------------------*/
void MDMath::calc(type0& dt,type0(&A)[3][3],type0(&EXP_A)[3][3],type0(&F_A)[3][3])
{
    auto div_dif_f=[](type0& dt,type0& x,type0& y,type0& fx,type0& fy)->type0
    {
        if(fabs(x-y)<sqrt(std::numeric_limits<type0>::epsilon()*f(0.5*dt*(x+y))/ddf(0.5*dt*(x+y))))
            return df(0.5*dt*(x+y));
        return (fx-fy)/(x*dt-y*dt);
    };
    
    auto div_dif_exp=[](type0& dt,type0& x,type0& y,type0& fx,type0& fy)->type0
    {
        if(fabs((x-y)*dt)<sqrt_eps)
            return exp(0.5*dt*(x+y));
        return (fx-fy)/((x-y)*dt);
    };
    
    F_A[0][0]=f(dt*A[0][0]);
    F_A[1][1]=f(dt*A[1][1]);
    F_A[2][2]=f(dt*A[2][2]);
    F_A[2][1]=F_A[1][0]=F_A[2][0]=0.0;
    
    EXP_A[0][0]=exp(dt*A[0][0]);
    EXP_A[1][1]=exp(dt*A[1][1]);
    EXP_A[2][2]=exp(dt*A[2][2]);
    EXP_A[2][1]=EXP_A[1][0]=EXP_A[2][0]=0.0;
    
    
    
    type0 f_yz=0.0,f_zx=0.0,f_xy=0.0,g_yz=0.0,g_zx=0.0,g_xy=0.0;
    if(A[2][1]!=0.0)
    {
        f_yz=div_dif_f(dt,A[1][1],A[2][2],F_A[1][1],F_A[2][2]);
        g_yz=div_dif_exp(dt,A[1][1],A[2][2],EXP_A[1][1],EXP_A[2][2]);
    }
    if(A[2][0]!=0.0)
    {
        f_zx=div_dif_f(dt,A[2][2],A[0][0],F_A[2][2],F_A[0][0]);
        g_zx=div_dif_exp(dt,A[2][2],A[0][0],EXP_A[2][2],EXP_A[0][0]);
    }
    if(A[1][0]!=0.0)
    {
        f_xy=div_dif_f(dt,A[0][0],A[1][1],F_A[0][0],F_A[1][1]);
        g_xy=div_dif_exp(dt,A[0][0],A[1][1],EXP_A[0][0],EXP_A[1][1]);
    }
    
    
    F_A[0][0]*=dt;
    F_A[1][1]*=dt;
    F_A[2][2]*=dt;
    F_A[2][1]=A[2][1]*dt*dt*f_yz;
    F_A[2][0]=A[2][0]*dt*dt*f_zx;
    F_A[1][0]=A[1][0]*dt*dt*f_xy;
    
    EXP_A[2][1]=A[2][1]*dt*g_yz;
    EXP_A[2][0]=A[2][0]*dt*g_zx;
    EXP_A[1][0]=A[1][0]*dt*g_xy;
    
    
    if(A[1][0]==0.0 || A[2][1]==0.0)
        return;
    
    auto div_dif_f3=[](type0& dt,type0& x,type0& y,type0& z,type0& f_yz,type0& f_zx,type0& f_xy)->type0
    {
        type0 __ddf=ddf((x+y+z)*dt/3.0);
        type0 lim=sqrt(std::numeric_limits<type0>::epsilon()*f((x+y+z)*dt/3.0)/__ddf);
        
        if(fabs((x-z)*dt)>lim)
            return (f_xy-f_yz)/((x-z)*dt);
        
        if(fabs((y-x)*dt)>lim)
            return (f_yz-f_zx)/((y-x)*dt);
        
        if(fabs((z-y)*dt)>lim)
            return (f_zx-f_xy)/((z-y)*dt);
        
        return __ddf;
    };
    
    auto div_dif_exp3=[](type0& dt,type0& x,type0& y,type0& z,type0& f_yz,type0& f_zx,type0& f_xy)->type0
    {
        if(fabs((x-z)*dt)>sqrt_eps)
            return (f_xy-f_yz)/((x-z)*dt);
        
        if(fabs((y-x)*dt)>sqrt_eps)
            return (f_yz-f_zx)/((y-x)*dt);
        
        if(fabs((z-y)*dt)>sqrt_eps)
            return (f_zx-f_xy)/((z-y)*dt);
        
        return exp((x+y+z)*dt/3.0);
    };
    
    F_A[2][0]+=dt*dt*dt*A[2][1]*A[1][0]*div_dif_f3(dt,A[0][0],A[1][1],A[2][2],f_yz,f_zx,f_xy);
    EXP_A[2][0]+=dt*dt*A[2][1]*A[1][0]*div_dif_exp3(dt,A[0][0],A[1][1],A[2][2],g_yz,g_zx,g_xy);
}
/*--------------------------------------------
 
 --------------------------------------------*/
void MDMath::calc(type0& dt,type0(&A)[2][2],type0(&EXP_A)[2][2],type0(&F_A)[2][2])
{
    auto div_dif_f=[](type0 dt,type0& x,type0& y,type0& fx,type0& fy)->type0
    {
        if(fabs(x-y)<sqrt(std::numeric_limits<type0>::epsilon()*f(0.5*dt*(x+y))/ddf(0.5*dt*(x+y))))
            return df(0.5*dt*(x+y));
        return (fx-fy)/(x*dt-y*dt);
    };
    
    auto div_dif_exp=[](type0& dt,type0& x,type0& y,type0& fx,type0& fy)->type0
    {
        if(fabs((x-y)*dt)<sqrt_eps)
            return exp(0.5*dt*(x+y));
        return (fx-fy)/((x-y)*dt);
    };
    
    F_A[0][0]=f(dt*A[0][0]);
    F_A[1][1]=f(dt*A[1][1]);
    F_A[1][0]=0.0;
    
    EXP_A[0][0]=exp(dt*A[0][0]);
    EXP_A[1][1]=exp(dt*A[1][1]);
    EXP_A[1][0]=0.0;
    
    if(A[1][0]!=0.0)
    {
        F_A[1][0]=A[1][0]*dt*dt*div_dif_f(dt,A[0][0],A[1][1],F_A[0][0],F_A[1][1]);
        EXP_A[1][0]=A[1][0]*div_dif_exp(dt,A[0][0],A[1][1],EXP_A[0][0],EXP_A[1][1]);
    }
    F_A[0][0]*=dt;
    F_A[1][1]*=dt;
}
/*--------------------------------------------
 
 --------------------------------------------*/
void MDMath::calc(type0& dt,type0(&A)[1][1],type0(&EXP_A)[1][1],type0(&F_A)[1][1])
{
    EXP_A[0][0]=exp(A[0][0]*dt);
    F_A[0][0]=dt*f(dt*A[0][0]);
}
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    











