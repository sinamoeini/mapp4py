#include "dae_imp.h"
#include "atoms_dmd.h"
#include "ff_dmd.h"
#include "dynamic_dmd.h"
#include "gmres.h"
#include "memory.h"
#include "xmath.h"
using namespace MAPP_NS;
/*--------------------------------------------
 
 --------------------------------------------*/
DAEImplicit::DAEImplicit():
DAE(),
max_nnewton_iters(5),
max_ngmres_iters(5),
y_0(NULL),
c_0(NULL),
del_c(NULL),
F(NULL),
a(NULL),
F_ptr(NULL),
del_c_ptr(NULL),
gmres(NULL)
{
}
/*--------------------------------------------
 
 --------------------------------------------*/
DAEImplicit::~DAEImplicit()
{
    delete [] y_0;
    y_0=c_0=a=NULL;
    delete F_ptr;
    delete del_c_ptr;
}
/*--------------------------------------------
 
 --------------------------------------------*/
void DAEImplicit::init_static()
{
    DAE::init_static();
    Memory::alloc(y_0,3*ncs);
    c_0=y_0+ncs;
    a=y_0+2*ncs;
    del_c_ptr=new Vec<type0>(atoms,c_dim);
    del_c=del_c_ptr->begin();
    F_ptr=new Vec<type0>(atoms,c_dim);
    F=F_ptr->begin();
    gmres=new GMRES(atoms,max_ngmres_iters,c_dim);
}
/*--------------------------------------------
 
 --------------------------------------------*/
void DAEImplicit::fin_static()
{
    delete gmres;
    gmres=NULL;
    delete F_ptr;
    F_ptr=NULL;
    F=NULL;
    delete del_c_ptr;
    del_c_ptr=NULL;
    del_c=NULL;
    Memory::dealloc(y_0);
    y_0=c_0=a=NULL;
    DAE::fin_static();
}
/*--------------------------------------------
 
 --------------------------------------------*/
void DAEImplicit::init()
{
    DAE::init();
    
    //static related

}
/*--------------------------------------------
 
 --------------------------------------------*/
void DAEImplicit::fin()
{
    //static related 

    
    
    DAE::fin();
}
/*--------------------------------------------
 
 --------------------------------------------*/
inline type0 DAEImplicit::update_c()
{
    type0 tmp;
    type0 r,r_lcl=1.0;
    for(int i=0;i<ncs;i++)
    {
        if(c[i]>=0.0)
        {
            tmp=c[i]-r_lcl*del_c[i];
            if(tmp>1.0)
            {
                r_lcl=(c[i]-1.0)/del_c[i];
                while(c[i]-r_lcl*del_c[i]>1.0)
                    r_lcl=nextafter(r_lcl,0.0);
            }
            if(tmp<0.0)
            {
                r_lcl=c[i]/del_c[i];
                while(c[i]-r_lcl*del_c[i]<0.0)
                    r_lcl=nextafter(r_lcl,0.0);
            }
        }
        else
            del_c[i]=0.0;
    }
    
    
    MPI_Allreduce(&r_lcl,&r,1,Vec<type0>::MPI_T,MPI_MIN,atoms->world);
    volatile type0 c0;
    
    for(int i=0;i<ncs;i++)
    {
        c0=c[i]-r*del_c[i];
        --++c0;
        c[i]=c0;
        /*
        dc=fabs(c0-c_init[i]);
        
        c[i]=c0;
        if(dc>0.0)
        {
            if(c[i]==0.0)
                curr_max_dc_rel_lcl=std::numeric_limits<type0>::infinity();
            if(c[i]>0.0)
                curr_max_dc_rel_lcl=MAX(dc/c_init[i],curr_max_dc_rel_lcl);
        }*/
    }
    //MPI_Allreduce(&curr_max_dc_rel_lcl,&curr_max_dc_rel,1,Vec<type0>::MPI_T,MPI_MAX,atoms->world);
    return r;
}
/*--------------------------------------------
 solve the implicit equation
 
 nc_dofs
 err_fac
 --------------------------------------------*/
bool DAEImplicit::newton()
{
    auto F_norm=[this]()->type0
    {
        ff->c_d_calc_timer();
        type0 F_norm_lcl=0.0;
        type0 beta_inv=1.0/beta;
        for(int i=0;i<ncs;i++)
        {
            F[i]=beta_inv*c[i]+a[i]-c_d[i];
            F_norm_lcl+=F[i]*F[i];
        }
        type0 F_norm;
        MPI_Allreduce(&F_norm_lcl,&F_norm,1,Vec<type0>::MPI_T,MPI_SUM,atoms->world);
        return sqrt(F_norm)/a_tol_sqrt_nc_dofs;
    };
    
    auto Jacobain_calc=[this](Vec<type0>* x,Vec<type0>* Jx)->void
    {
        ff->J_timer(x,Jx);
        type0* __Jx=Jx->begin();
        type0* __x=x->begin();
        type0 beta_inv=1.0/beta;
        for(int i=0;i<ncs;i++)
            __Jx[i]=-__Jx[i]+beta_inv*__x[i];
    };
    
    
    
    type0 res_tol=0.005*a_tol_sqrt_nc_dofs/err_fac;
    type0 denom=0.1*a_tol_sqrt_nc_dofs/err_fac;
    type0 norm=1.0,delta=0.0,delta_prev=0.0,ratio=1.0,R=1.0;
    
    int iter=0;
    
    //memcpy(c_0,c,ncs*sizeof(type0));
    //memcpy(c,y_0,ncs*sizeof(type0));
    //dynamic->update(atoms->c);
    type0 cost=F_norm();

    bool converge=false,diverge=false;
    while(iter<max_nnewton_iters && !converge && !diverge)
    {
        gmres->solve(Jacobain_calc,F_ptr,res_tol,norm,del_c_ptr);
        
        ratio=update_c();
        dynamic->update(atoms->c);
        
        delta=fabs(ratio*norm/denom);
        
        if(iter) R=MAX(0.3*R,delta/delta_prev);
        if(iter>1 && R*delta<1.0)
        {
            cost=F_norm();
            converge=true;
            continue;
        }
        
        if(iter>1 && delta/delta_prev>2.0)
        {
            diverge=true;
            continue;
        }
        
        delta_prev=delta;
        
        cost=F_norm();
        iter++;
    }
    
    
    if(converge)
    {
        if(iter) nnonlin_acc++;
        return true;
    }
    memcpy(c,c_0,ncs*sizeof(type0));
    dynamic->update(atoms->c);
    nnonlin_rej++;
    return false;
}
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
#include <iostream>
#include "random.h"
#include "ff_eam_dmd.h"
void DAEImplicit::J_test()
{
    
    Random rand(562148245);
    constexpr int nvecs=100;
    type0 delta=1.0e-9;
    
    Vec<type0>* h =new Vec<type0>(atoms,c_dim);
    Vec<type0>* Jh =new Vec<type0>(atoms,c_dim);
    type0** dFs;
    Memory::alloc(dFs,nvecs,ncs);
    

    type0* __h=h->begin();
    for(int i=0;i<ncs;i++)
    {
        c[i]=rand.uniform();
        if(rand.uniform()>0.5)
        {
            __h[i]=rand.uniform()*((1.0-c[i])/(delta*nvecs));
        }
        else
        {
            __h[i]=-rand.uniform()*(c[i]/(delta*nvecs));
        }
    }
    
    ForceFieldEAMDMD* __ff=dynamic_cast<ForceFieldEAMDMD*>(ff);
    
    dynamic->update(atoms->c);
    __ff->prep_J();
    __ff->J(h,Jh);
    memcpy(F,c_d,sizeof(type0)*ncs);
    
    memcpy(c_0,c,ncs*sizeof(type0));
    for(int ivec=0;ivec<nvecs;ivec++)
    {
        type0 __delta=delta*ivec;
        for(int i=0;i<ncs;i++)
            c[i]=c_0[i]+__delta*__h[i];
        dynamic->update(atoms->c);
        __ff->prep_J();
        for(int i=0;i<ncs;i++)
            dFs[ivec][i]=c_d[i]-F[i];
    }
    
    type0* __Jh=Jh->begin();
    FILE* fp=fopen("/Users/sina/Desktop/data.txt","w");
    for(int i=0;i<ncs;i++)
        for(int ivec=0;ivec<nvecs;ivec++)
            fprintf(fp,"%e\t%e\t%e\n",ivec*delta,dFs[ivec][i],__Jh[i]*ivec*delta);

    
    fclose(fp);
    

    Memory::dealloc(dFs);
    delete Jh;
    delete h;
}*/
/*------------------------------------------------------------------------------------------------------------------------------------
 
 ------------------------------------------------------------------------------------------------------------------------------------*/
PyObject* DAEImplicit::__new__(PyTypeObject* type,PyObject* args,PyObject* kwds)
{
    Object* __self=reinterpret_cast<Object*>(type->tp_alloc(type,0));
    PyObject* self=reinterpret_cast<PyObject*>(__self);
    return self;
}
/*--------------------------------------------
 
 --------------------------------------------*/
int DAEImplicit::__init__(PyObject* self,PyObject* args,PyObject* kwds)
{
    FuncAPI<> f("__init__");
    
    if(f(args,kwds)==-1) return -1;
    Object* __self=reinterpret_cast<Object*>(self);
    __self->dae=new DAEImplicit();
    __self->xprt=NULL;
    return 0;
}
/*--------------------------------------------
 
 --------------------------------------------*/
PyObject* DAEImplicit::__alloc__(PyTypeObject* type,Py_ssize_t)
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
void DAEImplicit::__dealloc__(PyObject* self)
{
    Object* __self=reinterpret_cast<Object*>(self);
    delete __self->dae;
    __self->dae=NULL;
    if(__self->xprt) Py_DECREF(__self->xprt);
    __self->xprt=NULL;
    delete __self;
}
/*--------------------------------------------*/
PyTypeObject DAEImplicit::TypeObject={PyObject_HEAD_INIT(NULL)};
/*--------------------------------------------*/
int DAEImplicit::setup_tp()
{
    TypeObject.tp_name="mapp.dmd.dae_implicit";
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
PyGetSetDef DAEImplicit::getset[]={[0 ... 2]={NULL,NULL,NULL,NULL,NULL}};
/*--------------------------------------------*/
void DAEImplicit::setup_tp_getset()
{
    getset_max_ngmres_iters(getset[0]);
    getset_max_nnewton_iters(getset[1]);
}
/*--------------------------------------------*/
PyMethodDef DAEImplicit::methods[]={[0 ... 0]={NULL}};
/*--------------------------------------------*/
void DAEImplicit::setup_tp_methods()
{
}
/*--------------------------------------------
 
 --------------------------------------------*/
void DAEImplicit::getset_max_ngmres_iters(PyGetSetDef& getset)
{
    getset.name=(char*)"max_ngmres_iters";
    getset.doc=(char*)R"---(
    (int) maximim number of gmres iterations
    
    Maximum number of iterations of linear solver (GMRES)
    )---";
    getset.get=[](PyObject* self,void*)->PyObject*
    {
        return var<int>::build(reinterpret_cast<Object*>(self)->dae->max_ngmres_iters,NULL);
    };
    getset.set=[](PyObject* self,PyObject* op,void*)->int
    {
        VarAPI<int> max_ngmres_iters("max_ngmres_iters");
        max_ngmres_iters.logics[0]=VLogics("gt",0);
        int ichk=max_ngmres_iters.set(op);
        if(ichk==-1) return -1;
        reinterpret_cast<Object*>(self)->dae->max_ngmres_iters=max_ngmres_iters.val;
        return 0;
    };
}
/*--------------------------------------------
 
 --------------------------------------------*/
void DAEImplicit::getset_max_nnewton_iters(PyGetSetDef& getset)
{
    getset.name=(char*)"max_nnewton_iters";
    getset.doc=(char*)R"---(
    (int) maximim number of newton iterations
    
    Maximum number of iterations of nonlinear solver (newton)
    )---";
    getset.get=[](PyObject* self,void*)->PyObject*
    {
        return var<int>::build(reinterpret_cast<Object*>(self)->dae->max_nnewton_iters,NULL);
    };
    getset.set=[](PyObject* self,PyObject* op,void*)->int
    {
        VarAPI<int> max_nnewton_iters("max_nnewton_iters");
        max_nnewton_iters.logics[0]=VLogics("gt",0);
        int ichk=max_nnewton_iters.set(op);
        if(ichk==-1) return -1;
        reinterpret_cast<Object*>(self)->dae->max_nnewton_iters=max_nnewton_iters.val;
        return 0;
    };
}
