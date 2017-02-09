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
MDNST::MDNST(AtomsMD*& __atoms,ForceFieldMD*& __ff,
type0(&__S)[__dim__][__dim__],type0 __T,type0 __dt):
MDNVT(__atoms,__ff,__T,__dt),
thermo_baro(__dt/2.0,1000.0*__dt,3,1),
nreset(0),
MLT0{[0 ... __dim__-1][0 ... __dim__-1]=0.0},
MLT1{[0 ... __dim__-1][0 ... __dim__-1]=0.0},
MLT2{[0 ... __dim__-1][0 ... __dim__-1]=0.0},
S{[0 ... __dim__-1][0 ... __dim__-1]=0.0},
S_dev{[0 ... __dim__-1][0 ... __dim__-1]=0.0},
B_ref{[0 ... __dim__-1][0 ... __dim__-1]=0.0},
t_relax_S{[0 ... __dim__-1][0 ... __dim__-1]=0.0},
V_H{[0 ... __dim__-1][0 ... __dim__-1]=0.0},
S_dof{[0 ... __dim__-1][0 ... __dim__-1]=false}
{
    Algebra::DoLT<__dim__>::func([this,&__S,&__dt]
    (int i,int j)
    {
        S_dof[i][j]=S_dof[j][i]=true;
        S[i][j]=S[j][i]=__S[i][j];
        t_relax_S[i][j]=t_relax_S[j][i]=1000.0*__dt;
        
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
    
    int nchains=thermo_part.nchains;
    int niters=thermo_part.niters;
    type0 t_relax=thermo_part.t_relax;
    
    thermo_part.~ThermostatNHC();
    new (&thermo_part) ThermostatNHC(__dt/2.0,t_relax,nchains,niters);
    
    
    nchains=thermo_baro.nchains;
    niters=thermo_baro.niters;
    t_relax=thermo_baro.t_relax;
    
    thermo_baro.~ThermostatNHC();
    new (&thermo_baro) ThermostatNHC(__dt/2.0,t_relax,nchains,niters);
}
/*--------------------------------------------
 
 --------------------------------------------*/
void MDNST::update_x()
{
    MDMath::calc(dt,V_H,MLT0,MLT1);
    
    type0* x=atoms->x->begin();
    type0* x_d=atoms->x_d->begin();
    for(int i=0;i<atoms->natms;++i)
    {
        Algebra::V_mul_MLT(x,MLT0,v0);
        Algebra::V_mul_MLT(x_d,MLT1,v1);
        
        Algebra::Do<__dim__>::func(
        [&x,this](const int j){x[j]=v0[j]+v1[j];});
        
        x+=__dim__;
        x_d+=__dim__;
    }
    
    Algebra::MLT_mul_MLT(atoms->H,MLT0,atoms->H);
    atoms->update_H();
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
    type0* m=atoms->elements->masses;
    type0 m_i,m_inv;
    const int natms0=atoms->natms;
    for(int i=0;i<natms0;++i)
    {
        m_inv=1.0/m[*elem];
        MDMath::____NONAME0<__dim__,__dim__>::func(xi,x_d,&MLT0[0][0],m_inv,f,&MLT1[0][0]);
        MDMath::____NONAME1<__dim__,__dim__>::func(x,&V_H[0][0],x_d,v0);
        MDMath::____NONAME2<__dim__,__dim__>::func(v0,&MLT2[0][0]);
        
        Algebra::V_add<__dim__>(v0,x);
        
        f+=__dim__;
        x_d+=__dim__;
        x+=__dim__;
        ++elem;
    }
    
    Algebra::MLT_mul_MLT(MLT2,V_H,MLT2);
    Algebra::Do<__dim__>::func([this](int i){MLT2[i][i]++;});
    Algebra::MLT_mul_MLT(atoms->H,MLT2,atoms->H);
    atoms->update_H();
    
    dynamic->update(atoms->x);
    ff->reset();
    ff->force_calc_timer();
    
    
    f=ff->f->begin();
    x_d=atoms->x_d->begin();
    elem=atoms->elem->begin();
    Algebra::zero(__vec_lcl);
    const int natms1=atoms->natms;
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
    MPI_Allreduce(__vec_lcl,mvv,__nvoigt__,MPI_TYPE0,MPI_SUM,world);
    T_part=Algebra::Tr_DyadicV<__dim__>(mvv)/(ndof_part*kB);
}
/*--------------------------------------------
 
 --------------------------------------------*/
void MDNST::update_x_d(type0 fac_x_d)
{
    
    Algebra::DoLT<__dim__>::func([this](int i,int j){MLT0[i][j]=-V_H[i][j];});
    type0 tr=Algebra::Tr(V_H)/ndof_part;
    Algebra::Do<__dim__>::func([this,&tr](int i){MLT0[i][i]-=tr;});
    MDMath::calc(dt2,MLT0,MLT1,MLT2);
    if(fac_x_d!=1.0) Algebra::DoLT<__dim__>::func([this,&fac_x_d](int i,int j){MLT1[i][j]*=fac_x_d;});
    
    type0* f=ff->f->begin();
    type0* x_d=atoms->x_d->begin();
    elem_type* elem=atoms->elem->begin();
    type0* m=atoms->elements->masses;
    type0 m_i;
    Algebra::zero(__vec_lcl);
    
    for(int i=0;i<atoms->natms;++i)
    {
        Algebra::V_mul_MLT(x_d,MLT1,v0);
        Algebra::V_mul_MLT(f,MLT2,v1);
        m_i=m[*elem];
        
        Algebra::Do<__dim__>::func([&x_d,&m_i,this](const int j){x_d[j]=v0[j]+v1[j]/m_i;});
        Algebra::DyadicV<__dim__>(m_i,x_d,__vec_lcl);
        
        f+=__dim__;
        x_d+=__dim__;
        ++elem;
    }
    MPI_Allreduce(__vec_lcl,mvv,__nvoigt__,MPI_TYPE0,MPI_SUM,world);
    T_part=Algebra::Tr_DyadicV<__dim__>(mvv)/(ndof_part*kB);
}
/*--------------------------------------------
 
 --------------------------------------------*/
void MDNST::dof_consistency()
{
    if(!atoms->dof) return;
    bool* dof=atoms->dof->begin();
    int __dof_lcl[__dim__]{[0 ... __dim__-1]=0};
    for(int i=0;i<atoms->natms;i++,dof+=__dim__)
        Algebra::Do<__dim__>::func([&dof,&__dof_lcl](int i){ if(!dof[i]) __dof_lcl[i]=1;});
    
    int __dof[__dim__]{[0 ... __dim__-1]=0};
    MPI_Allreduce(__dof_lcl,__dof,__dim__,MPI_INT,MPI_MAX,world);
    for(int i=0;i<__dim__;i++)
        for(int j=i;j<__dim__;j++)
            if(S_dof[i][j] && __dof[i])
                throw "cannot impose stress component ["+std::to_string(i)+"]["+std::to_string(j)
                +"] while any of the atoms do not have degree freedom in "+std::to_string(i)
                +" direction";
    
}
/*--------------------------------------------
 
 --------------------------------------------*/
void MDNST::init()
{
    MDNVT::init();
    
    ndof_baro=0;
    T_baro=0.0;
    Algebra::DoLT<__dim__>::func([this](int i,int j)
    {
        if(S_dof[i][j])
        {
            V_H_prefac[i][j]=1.0/(t_relax_S[i][j]*t_relax_S[i][j]*ndof_part*kB*T);
            ndof_baro++;
            T_baro+=V_H[i][j]*V_H[i][j]*t_relax_S[i][j]*t_relax_S[i][j];
            
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
void MDNST::fin()
{}
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
            -atoms->vol*ff->nrgy_strss[1+i+j*__dim__-j*(j+1)/2]
            +mvv[i+j*__dim__-j*(j+1)/2])*V_H_prefac[i][j]*dt2;
            
            T_baro+=V_H[i][j]*V_H[i][j]*t_relax_S[i][j]*t_relax_S[i][j];
            
        }
    });
    
    T_baro*=ndof_part*T/ndof_baro;

}
/*--------------------------------------------
 
 --------------------------------------------*/
void MDNST::run(int nsteps)
{
    init();
    
    if(atoms->dof)
        dynamic=new DynamicMD(atoms,ff,true,{atoms->elem},{atoms->x_d,atoms->dof});
    else
        dynamic=new DynamicMD(atoms,ff,true,{atoms->elem},{atoms->x_d});
    
    dynamic->init();
    
    ff->reset();
    ff->force_calc_timer();
    ThermoDynamics thermo(6,"T",T_part,"PE",ff->nrgy_strss[0],
    "S[0][0]",S_part[0][0],
    "S[1][1]",S_part[1][1],
    "S[2][2]",S_part[2][2],
    "S[1][2]",S_part[2][1],
    "S[2][0]",S_part[2][0],
    "S[0][1]",S_part[1][0]);
    
    thermo.init();
    Algebra::DoLT<__dim__>::func([this](const int i,const int j)
    {
        S_part[i][j]=ff->nrgy_strss[1+i+j*__dim__-j*(j+1)/2]-mvv[i+j*__dim__-j*(j+1)/2]/atoms->vol;
    });
    
    
    thermo.print(0);
    
    //type0* x_d;
    type0 fac,fac_x_d=1.0;
    
    // store reference configuration
    Algebra::DoLT<__dim__>::func([this](int i,int j){B_ref[i][j]=atoms->B[i][j];});
    vol_ref=atoms->vol;
    
    for(int i=0;i<nsteps;i++)
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
            S_part[i][j]=ff->nrgy_strss[1+i+j*__dim__-j*(j+1)/2]-mvv[i+j*__dim__-j*(j+1)/2]/atoms->vol;
        });
        
        
        
        if(nreset && (i+1)%nreset==0)
        {
            Algebra::DoLT<__dim__>::func([this](int i,int j){B_ref[i][j]=atoms->B[i][j];});
            vol_ref=atoms->vol;
        }
        
        if((i+1)%ntally==0) thermo.print(i+1);
    }
    
    if(nsteps%ntally) thermo.print(nsteps);
    
    update_x_d_final(fac_x_d);
    
    thermo.fin();
    
    dynamic->fin();
    delete dynamic;
    dynamic=NULL;
    fin();
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
    FuncAPI<OP<AtomsMD>,symm<type0[__dim__][__dim__]>,type0,type0> f("__init__",{"atoms_md","S","T","dt"});
    
    f.logics<2>()[0]=VLogics("gt",0.0);
    f.logics<3>()[0]=VLogics("gt",0.0);
    
    if(f(args,kwds)==-1) return -1;
    Object* __self=reinterpret_cast<Object*>(self);
    AtomsMD::Object* atoms_md=reinterpret_cast<AtomsMD::Object*>(f.val<0>().ob);
    __self->md=new MDNST(atoms_md->atoms,atoms_md->ff,f.val<1>(),f.val<2>(),f.val<3>());
    __self->atoms_md=atoms_md;
    Py_INCREF(atoms_md);
    return 0;
}
/*--------------------------------------------
 
 --------------------------------------------*/
PyObject* MDNST::__alloc__(PyTypeObject* type,Py_ssize_t)
{
    Object* __self=new Object;
    __self->ob_type=type;
    __self->ob_refcnt=1;
    __self->md=NULL;
    __self->atoms_md=NULL;
    return reinterpret_cast<PyObject*>(__self);
}
/*--------------------------------------------
 
 --------------------------------------------*/
void MDNST::__dealloc__(PyObject* self)
{
    Object* __self=reinterpret_cast<Object*>(self);
    delete __self->md;
    __self->md=NULL;
    if(__self->atoms_md) Py_DECREF(__self->atoms_md);
    __self->atoms_md=NULL;
    delete __self;
}
/*--------------------------------------------*/
PyTypeObject MDNST::TypeObject ={PyObject_HEAD_INIT(NULL)};
/*--------------------------------------------*/
void MDNST::setup_tp()
{
    TypeObject.tp_name="md_nst";
    TypeObject.tp_doc="MD of isothermal–isobaric ensemble";
    
    TypeObject.tp_flags=Py_TPFLAGS_DEFAULT;
    TypeObject.tp_basicsize=sizeof(Object);
    
    TypeObject.tp_new=__new__;
    TypeObject.tp_init=__init__;
    TypeObject.tp_alloc=__alloc__;
    TypeObject.tp_dealloc=__dealloc__;
    setup_tp_getset();
    TypeObject.tp_getset=getset;
    
    TypeObject.tp_base=&MDNVT::TypeObject;
}
/*--------------------------------------------*/
PyGetSetDef MDNST::getset[]={[0 ... 7]={NULL,NULL,NULL,NULL,NULL}};
/*--------------------------------------------*/
void MDNST::setup_tp_getset()
{
    getset_niters_baro(getset[0]);
    getset_nchains_baro(getset[1]);
    getset_t_relax_baro(getset[2]);
    getset_S_dof(getset[3]);
    getset_S(getset[4]);
    getset_t_relax_S(getset[5]);
    getset_nreset(getset[6]);
}
/*--------------------------------------------
 
 --------------------------------------------*/
void MDNST::getset_nchains_baro(PyGetSetDef& getset)
{
    getset.name=(char*)"nchains_baro";
    getset.doc=(char*)"number of chains in Nose-Hoover chain thermostat";
    getset.get=[](PyObject* self,void*)->PyObject*
    {
        int nchains=reinterpret_cast<Object*>(self)->md->thermo_baro.nchains;
        return var<int>::build(nchains,NULL);
    };
    getset.set=[](PyObject* self,PyObject* op,void*)->int
    {
        VarAPI<int> nchains_baro("nchains_baro");
        nchains_baro.logics[0]=VLogics("gt",0);
        int ichk=nchains_baro.set(op);
        if(ichk==-1) return -1;
        if(reinterpret_cast<Object*>(self)->md->thermo_baro.nchains==nchains_baro.val)
            return 0;
        
        ThermostatNHC& thermo_baro=reinterpret_cast<Object*>(self)->md->thermo_baro;
        int niters=thermo_baro.niters;
        type0 t_relax=thermo_baro.t_relax;
        type0 __dt=reinterpret_cast<Object*>(self)->md->dt2;
        thermo_baro.~ThermostatNHC();
        new (&thermo_baro) ThermostatNHC(__dt,t_relax,nchains_baro.val,niters);
        return 0;
    };
}
/*--------------------------------------------
 
 --------------------------------------------*/
void MDNST::getset_niters_baro(PyGetSetDef& getset)
{
    getset.name=(char*)"niters_baro";
    getset.doc=(char*)"number of iterations in Nose-Hoover chain thermostat";
    getset.get=[](PyObject* self,void*)->PyObject*
    {
        int niters=reinterpret_cast<Object*>(self)->md->thermo_baro.niters;
        return var<int>::build(niters,NULL);
    };
    getset.set=[](PyObject* self,PyObject* op,void*)->int
    {
        VarAPI<int> niters_baro("niters_baro");
        niters_baro.logics[0]=VLogics("gt",0);
        int ichk=niters_baro.set(op);
        if(ichk==-1) return -1;
        if(reinterpret_cast<Object*>(self)->md->thermo_baro.niters==niters_baro.val)
            return 0;
        
        ThermostatNHC& thermo_baro=reinterpret_cast<Object*>(self)->md->thermo_baro;
        int nchains=thermo_baro.nchains;
        type0 t_relax=thermo_baro.t_relax;
        type0 __dt=reinterpret_cast<Object*>(self)->md->dt2;
        thermo_baro.~ThermostatNHC();
        new (&thermo_baro) ThermostatNHC(__dt,t_relax,nchains,niters_baro.val);
        return 0;
    };
}
/*--------------------------------------------
 
 --------------------------------------------*/
void MDNST::getset_t_relax_baro(PyGetSetDef& getset)
{
    getset.name=(char*)"t_relax_baro";
    getset.doc=(char*)"thermostat parameter, relaxation time";
    getset.get=[](PyObject* self,void*)->PyObject*
    {
        type0 t_relax=reinterpret_cast<Object*>(self)->md->thermo_baro.t_relax;
        return var<type0>::build(t_relax,NULL);
    };
    getset.set=[](PyObject* self,PyObject* op,void*)->int
    {
        VarAPI<type0> t_relax_baro("t_relax_baro");
        t_relax_baro.logics[0]=VLogics("gt",0.0);
        int ichk=t_relax_baro.set(op);
        if(ichk==-1) return -1;
        if(reinterpret_cast<Object*>(self)->md->thermo_baro.t_relax==t_relax_baro.val)
            return 0;
        
        ThermostatNHC& thermo_baro=reinterpret_cast<Object*>(self)->md->thermo_baro;
        int nchains=thermo_baro.nchains;
        int niters=thermo_baro.niters;
        type0 __dt=reinterpret_cast<Object*>(self)->md->dt2;
        thermo_baro.~ThermostatNHC();
        new (&thermo_baro) ThermostatNHC(__dt,t_relax_baro.val,nchains,niters);
        return 0;
    };
}
/*--------------------------------------------
 
 --------------------------------------------*/
void MDNST::getset_S_dof(PyGetSetDef& getset)
{
    getset.name=(char*)"S_dof";
    getset.doc=(char*)"S_dof";
    getset.get=[](PyObject* self,void*)->PyObject*
    {
        return var<symm<bool[__dim__][__dim__]>>::build(reinterpret_cast<Object*>(self)->md->S_dof,NULL);
    };
    getset.set=[](PyObject* self,PyObject* op,void*)->int
    {
        VarAPI<symm<bool[__dim__][__dim__]>> S_dof("S_dof");
        int ichk=S_dof.set(op);
        if(ichk==-1) return -1;
        
        bool (&__S_dof)[__dim__][__dim__]=reinterpret_cast<Object*>(self)->md->S_dof;
        type0 (&__S)[__dim__][__dim__]=reinterpret_cast<Object*>(self)->md->S;
        Algebra::DoLT<__dim__>::func([&__S_dof,&__S,&S_dof](int i,int j)
        {
            __S_dof[i][j]=__S_dof[j][i]=S_dof.val[i][j];
            if(!__S_dof[i][j])
                __S[i][j]=__S[j][i]=0.0;
        });
        
        return 0;
    };
}
/*--------------------------------------------
 
 --------------------------------------------*/
void MDNST::getset_S(PyGetSetDef& getset)
{
    getset.name=(char*)"S";
    getset.doc=(char*)"external stress of isothermal–isobaric ensemble";
    getset.get=[](PyObject* self,void*)->PyObject*
    {
        return var<symm<type0[__dim__][__dim__]>>::build(reinterpret_cast<Object*>(self)->md->S,NULL);
    };
    getset.set=[](PyObject* self,PyObject* op,void*)->int
    {
        VarAPI<symm<type0[__dim__][__dim__]>> S("S_dof");
        int ichk=S.set(op);
        if(ichk==-1) return -1;
        
        bool (&__S_dof)[__dim__][__dim__]=reinterpret_cast<Object*>(self)->md->S_dof;
        type0 (&__S)[__dim__][__dim__]=reinterpret_cast<Object*>(self)->md->S;
        type0 (&__t_relax_S)[__dim__][__dim__]=reinterpret_cast<Object*>(self)->md->t_relax_S;
        type0 __dt=reinterpret_cast<Object*>(self)->md->dt;
        Algebra::DoLT<__dim__>::func([&__S_dof,&__S,&S,&__t_relax_S,&__dt](int i,int j)
        {
            __S_dof[i][j]=__S_dof[j][i]=true;
            __S[i][j]=__S[j][i]=S.val[i][j];
            if(__t_relax_S[i][j]==0.0)
                __t_relax_S[i][j]=__t_relax_S[j][i]=1000.0*__dt;
        });
        
        return 0;
    };
}
/*--------------------------------------------
 
 --------------------------------------------*/
void MDNST::getset_t_relax_S(PyGetSetDef& getset)
{
    getset.name=(char*)"t_relax_S";
    getset.doc=(char*)"t_relax of external stress of isothermal–isobaric ensemble";
    getset.get=[](PyObject* self,void*)->PyObject*
    {
        return var<symm<type0[__dim__][__dim__]>>::build(reinterpret_cast<Object*>(self)->md->t_relax_S,NULL);
    };
    getset.set=[](PyObject* self,PyObject* op,void*)->int
    {
        VarAPI<symm<type0[__dim__][__dim__]>> t_relax_S("t_relax_S");
        t_relax_S.logics[0]=VLogics("gt",0.0);
        int ichk=t_relax_S.set(op);
        if(ichk==-1) return -1;
        
        type0 (&__t_relax_S)[__dim__][__dim__]=reinterpret_cast<Object*>(self)->md->t_relax_S;
        Algebra::DoLT<__dim__>::func([&__t_relax_S,&t_relax_S](int i,int j)
        {
            __t_relax_S[i][j]=__t_relax_S[j][i]=t_relax_S.val[i][j];
        });
        
        return 0;
    };
}
/*--------------------------------------------
 
 --------------------------------------------*/
void MDNST::getset_nreset(PyGetSetDef& getset)
{
    getset.name=(char*)"nreset";
    getset.doc=(char*)"number of steps to reset initial conf";
    getset.get=[](PyObject* self,void*)->PyObject*
    {
        int niters=reinterpret_cast<Object*>(self)->md->nreset;
        return var<int>::build(niters,NULL);
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
    if(A[2][1])
    {
        f_yz=div_dif_f(y,z,F_A0[1][1],F_A0[2][2]);
        __f_yz=div_dif_f(__y,__z,F_A1[1][1],F_A1[2][2]);
        g_yz=div_dif_exp(y,z,EXP_A[1][1],EXP_A[2][2]);
    }
    if(A[2][0])
    {
        f_zx=div_dif_f(z,x,F_A0[2][2],F_A0[0][0]);
        __f_zx=div_dif_f(__z,__x,F_A1[2][2],F_A1[0][0]);
        g_zx=div_dif_exp(z,x,EXP_A[2][2],EXP_A[0][0]);
    }
    if(A[1][0])
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
    if(A[1][0])
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
    if(A[2][1])
    {
        f_yz=div_dif_f(dt,A[1][1],A[2][2],F_A[1][1],F_A[2][2]);
        g_yz=div_dif_exp(dt,A[1][1],A[2][2],EXP_A[1][1],EXP_A[2][2]);
    }
    if(A[2][0])
    {
        f_zx=div_dif_f(dt,A[2][2],A[0][0],F_A[2][2],F_A[0][0]);
        g_zx=div_dif_exp(dt,A[2][2],A[0][0],EXP_A[2][2],EXP_A[0][0]);
    }
    if(A[1][0])
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
    
    if(A[1][0])
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
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    











