#include "potfit_o.h"
#include "import_cfg.h"
#include "min.h"


using namespace MAPP_NS;
const char* PotFitO::err_msgs[]=
{
    [LS_S]="",
    [LS_F_DOWNHILL]="line search failed: not downhill direction\n",
    [LS_F_GRAD0]="line search failed: gradient is zero\n",
    [LS_MIN_ALPHA]="line search failed: minimum alpha reached\n",
    [MIN_S_TOLERANCE]="minimization finished: energy tolerance reached\n",
    [MIN_F_MAX_ITER]="minimization finished: maximum iteration reached\n",
    [B_S]="",
    [B_F_MAX_ALPHA]="bracketing failed: maximum alpha reached\n",
    [B_F_DOWNHILL]="bracketing failed: not downhill direction\n"
    
};
/*--------------------------------------------
 constructor
 --------------------------------------------*/
PotFitO::PotFitO(
type0(&__A_rho_H)[nrho_H],type0(&__A_phi_FeH)[nphi_FeH],type0(&__A_phi_HH)[nphi_HH],type0(&__A_F_H)[nF_H],
std::string*&& __names_str,
std::string*& files,
int*& nprocs,
type0(*& __targets)[1+__nvoigt__],
type0(*& __coefs)[1+__nvoigt__],
PyObject** dof_funcs,
size_t __nconfigs,
MPI_Comm& __world):
atoms(NULL),
world(__world)
{
    nmin_steps=1000;
    max_ntrials=5;
    err_tol=1.0e-8;
    type0* __dx_max;
    __dx_max=dx_max+rho_H_offset;
    __dx_max[0]=100.0;
    __dx_max[1]=0.1;
    
    __dx_max=dx_max+phi_FeH_offset;
    __dx_max[0]=0.1;
    __dx_max[1]=0.1;
    __dx_max[2]=0.05;
    __dx_max[2]=0.1;
    
    __dx_max=dx_max+phi_HH_offset;
    __dx_max[0]=0.01;
    __dx_max[1]=0.01;
    __dx_max[2]=0.01;
    __dx_max[3]=0.01;
    
    __dx_max=dx_max+F_H_offset;
    __dx_max[0]=0.1;
    __dx_max[1]=1.0;
    __dx_max[2]=1.0;
    //__dx_max[3]=1.0;
    
    
    nconfigs=static_cast<int>(__nconfigs);
    names_str=__names_str;
    __names_str=NULL;
    Memory::alloc(errs,__nconfigs);
    Memory::alloc(roots,__nconfigs);
    Memory::alloc(names,__nconfigs);
    Memory::alloc(nHs,__nconfigs);
    Memory::alloc(mean_rhoHs,__nconfigs);
    MPI_Comm_rank(world,&my_rank);
    //lets figure out which configuration I am
    int tot_nproc=0;
    my_conf=0;
    my_lcl_rank=0;
    for(int i=0;i<nconfigs;i++)
    {
        names[i]=names_str[i].c_str();
        if(tot_nproc<=my_rank && my_rank<tot_nproc+nprocs[i])
        {
            my_conf=i;
            my_lcl_rank=my_rank-tot_nproc;
        }
        roots[i]=tot_nproc;
        tot_nproc+=nprocs[i];
    }
    my_world=new MPI_Comm;
    
    MPI_Comm_split(world,my_conf,my_lcl_rank,my_world);
    coef=__coefs[my_conf][0];
    target=__targets[my_conf][0];
    
    
    
    
    
    ntrial=0;

    

    
    
    /* loading the configuration */
    ImportCFGMD read(*my_world);
    atoms=read(files[my_conf].c_str());
    if(dof_funcs[my_conf]) atoms->DO(dof_funcs[my_conf]);
    
    
    ff=new ForceFieldEAMFitO(atoms,__A_rho_H,__A_phi_FeH,__A_phi_HH,__A_F_H);
    
    bool H_dof[__dim__][__dim__];
    for(int i=0;i<__dim__;i++) for(int j=0;j<__dim__;j++) H_dof[i][j]=false;
    H_dof[0][0]=H_dof[1][1]=H_dof[2][2]=true;
    
    Xorig.~VecTens();
    new (&Xorig) VecTens<type0,1>(atoms,true,__dim__);
    X0.~VecTens();
    new (&X0) VecTens<type0,1>(atoms,true,__dim__);
    
    min_ls=new LineSearchBrent();
    min=new MinCGPotFit(0.0,H_dof,false,0.4,min_ls,Xorig.vecs[0],X0.vecs[0]);
    min->atoms=atoms;
    min->ff=ff;
    
    memcpy(Xorig.vecs[0]->begin(),atoms->x->begin(),sizeof(type0)*__dim__*atoms->natms_lcl);
    memcpy(&Xorig.A[0][0],&atoms->H[0][0],sizeof(type0)*__dim__*__dim__);
    
    
}
/*--------------------------------------------
 destructor
 --------------------------------------------*/
PotFitO::~PotFitO()
{
    X0.~VecTens();
    Xorig.~VecTens();
    delete ff;
    delete min_ls;
    delete min;
    delete atoms;
    Memory::dealloc(roots);
    Memory::dealloc(errs);
    *names=NULL;
    Memory::dealloc(names);
    Memory::dealloc(names_str);
    Memory::dealloc(mean_rhoHs);
    Memory::dealloc(nHs);
}
/*--------------------------------------------
 
 --------------------------------------------*/
void PotFitO::calc_nHs()
{
    int nH_lcl=0,nH;
    elem_type* elem=atoms->elem->begin();
    const int natms_lcl=atoms->natms_lcl;
    for(int i=0;i<natms_lcl;i++) if(elem[i]==1) nH_lcl++;
    MPI_Allreduce(&nH_lcl,&nH,1,Vec<int>::MPI_T,MPI_SUM,*my_world);
    nHs[my_conf]=nH;
    for(int i=0;i<nconfigs;i++)
        MPI_Bcast(&nHs[i],1,Vec<int>::MPI_T,roots[i],world);
}
/*--------------------------------------------
 
 --------------------------------------------*/
void PotFitO::calc_mean_rho_Hs()
{
    for(int i=0;i<nconfigs;i++)
        MPI_Bcast(&mean_rhoHs[i],1,Vec<type0>::MPI_T,roots[i],world);
}
/*--------------------------------------------
 just one test run
 --------------------------------------------*/
void PotFitO::gen_test(int comp,type0 delta,int n)
{
    type0 dx[nHvars];
    Algebra::zero<nHvars>(dx);
    dx[comp]=1.0;
    
    min->init();
    min->run(0);
    ff->gradient();
    type0 pe0=atoms->pe;
    type0 deriv=0.0;
    Algebra::Do<nHvars>::func([this,&deriv,&dx](int i){deriv+=ff->dv[i]*dx[i];});
    
    for(int i=0;i<n;i++)
    {
        min->run(0);
        printf("%.16lf\t%.16lf\t%.16lf\n",i*delta,deriv*delta*i,atoms->pe-pe0);
        Algebra::Do<nHvars>::func([this,&delta,&dx](int i){ff->v[i]+=dx[i]*delta;});
    }

    min->fin();
}
/*--------------------------------------------
 just one test run
 --------------------------------------------*/
void PotFitO::store_x0()
{
    Algebra::V_eq<nHvars>(ff->v,x0);
    memcpy(X0.vecs[0]->begin(),atoms->x->begin(),sizeof(type0)*__dim__*atoms->natms_lcl);
    memcpy(&X0.A[0][0],&atoms->H[0][0],sizeof(type0)*__dim__*__dim__);
}
/*--------------------------------------------
 just one test run
 --------------------------------------------*/
void PotFitO::restore_x0()
{
    Algebra::V_eq<nHvars>(x0,ff->v);
    memcpy(atoms->x->begin(),X0.vecs[0]->begin(),sizeof(type0)*__dim__*atoms->natms_lcl);
    memcpy(&atoms->H[0][0],&X0.A[0][0],sizeof(type0)*__dim__*__dim__);
    atoms->update_H();
}
/*--------------------------------------------
 just one test run
 --------------------------------------------*/
void PotFitO::full_reset()
{
    memcpy(atoms->x->begin(),Xorig.vecs[0]->begin(),sizeof(type0)*__dim__*atoms->natms_lcl);
    memcpy(&atoms->H[0][0],&Xorig.A[0][0],sizeof(type0)*__dim__*__dim__);
    atoms->update_H();
}
/*--------------------------------------------
 just one test run
 --------------------------------------------*/
ThermoDynamics PotFitO::get_thermo()
{
    if(nconfigs==2)
        return ThermoDynamics(6,"ERR",tot_err,names[0],errs[0],names[1],errs[1]);
    else if(nconfigs==3)
        return ThermoDynamics(6,"ERR",tot_err,names[0],errs[0],names[1],errs[1],names[2],errs[2]);
    else if(nconfigs==4)
        return ThermoDynamics(6,"ERR",tot_err,names[0],errs[0],names[1],errs[1],names[2],errs[2],names[3],errs[3]);
    else if(nconfigs==5)
        return ThermoDynamics(6,"ERR",tot_err,names[0],errs[0],names[1],errs[1],names[2],errs[2],names[3],errs[3],names[4],errs[4]);
    else if(nconfigs==6)
        return ThermoDynamics(6,"ERR",tot_err,names[0],errs[0],names[1],errs[1],names[2],errs[2],names[3],errs[3],names[4],errs[4],names[5],errs[5]);
    else if(nconfigs==7)
        return ThermoDynamics(6,"ERR",tot_err,names[0],errs[0],names[1],errs[1],names[2],errs[2],names[3],errs[3],names[4],errs[4],names[5],errs[5],names[6],errs[6]);
    else if(nconfigs==8)
        return ThermoDynamics(6,"ERR",tot_err,names[0],errs[0],names[1],errs[1],names[2],errs[2],names[3],errs[3],names[4],errs[4],names[5],errs[5],names[6],errs[6],names[7],errs[7]);
    else if(nconfigs==9)
        return ThermoDynamics(6,"ERR",tot_err,names[0],errs[0],names[1],errs[1],names[2],errs[2],names[3],errs[3],names[4],errs[4],names[5],errs[5],names[6],errs[6],names[7],errs[7],names[8],errs[8]);
    else if(nconfigs==10)
        return ThermoDynamics(6,"ERR",tot_err,names[0],errs[0],names[1],errs[1],names[2],errs[2],names[3],errs[3],names[4],errs[4],names[5],errs[5],names[6],errs[6],names[7],errs[7],names[8],errs[8],names[9],errs[9]);
    return ThermoDynamics(6,"ERR",tot_err,names[0],errs[0]);
    
}
/*--------------------------------------------
 just one test run
 --------------------------------------------*/
void PotFitO::min_cg(int nsteps)
{
    type0 curr_err,prev_err,f0_f0,f_f0,f_f,f_h,ratio,alpha,tol=1.0e-12;

    curr_err=iter();
    Algebra::V_eq<nHvars>(f,h);
    f0_f0=Algebra::V_mul_V<nHvars>(f,f);

    ThermoDynamics thermo=get_thermo();

    int istep=0;
    thermo.init();
    thermo.print(0);
    int err=nsteps==0? MIN_F_MAX_ITER:LS_S;
    for(;err==LS_S;istep++)
    {
        if(f0_f0==0.0)
        {
            err=LS_F_GRAD0;
            continue;
        }

        store_x0();
        
        Algebra::V_eq<nHvars>(f,f0);
        prev_err=curr_err;
        f_h=Algebra::V_mul_V<nHvars>(f,h);
        if(f_h<0.0)
        {
            Algebra::V_eq<nHvars>(f,h);
            f_h=f0_f0;
        }


        //perform linesearch here
        //printf("%.16lf\n",curr_err);
        err=min_ls->min(this,curr_err,alpha,1);
        //check for error from linesearch here
        if(err!=LS_S)
        {
            // this was a bullshit step so we have to decrease the setup once to compensate
            // the last step was the previous one
            istep--;
            continue;
        }
        curr_err=iter();
        thermo.print(istep+1);

        if(prev_err-curr_err<tol || curr_err<err_tol) err=MIN_S_TOLERANCE;
        if(istep+1==nsteps) err=MIN_F_MAX_ITER;
        if(err) continue;

        f_f=Algebra::V_mul_V<nHvars>(f,f);
        f_f0=Algebra::V_mul_V<nHvars>(f,f0);
        ratio=(f_f-f_f0)/(f0_f0);
        
        for(int j=0;j<nHvars;j++) h[j]=ratio*h[j]+f[j];
        //advance h here
        f0_f0=f_f;
    }

    thermo.fin();
    fprintf(MAPP::mapp_out,"%s",err_msgs[err]);
}
/*--------------------------------------------
 just one test run
 --------------------------------------------*/
type0 PotFitO::iter()
{
    if(ntrial%max_ntrials==0) full_reset();
    ntrial++;
    min->init();
    min->run(nmin_steps);
    mean_rhoHs[my_conf]=ff->mean_rho_H();
    ff->gradient();
    min->fin();
    type0 err_lcl=0.0,err;
    Algebra::zero<nHvars>(f_lcl);
    if(my_lcl_rank==0)
    {
        err_lcl=coef*Algebra::pow<2>(atoms->pe-target);
        Algebra::V_eq_x_mul_V<nHvars>(-2.0*coef*(atoms->pe-target),ff->dv,f_lcl);
        errs[my_conf]=atoms->pe-target;
    }
    MPI_Allreduce(f_lcl,f,nHvars,Vec<type0>::MPI_T,MPI_SUM,world);
    MPI_Allreduce(&err_lcl,&err,1,Vec<type0>::MPI_T,MPI_SUM,world);
    for(int i=0;i<nconfigs;i++)
        MPI_Bcast(&errs[i],1,Vec<type0>::MPI_T,roots[i],world);
    tot_err=err;
    
    Algebra::Do<nHvars>::func([this](int i){if(!ff->dof[i]) f[i]=0.0;});
    
    return err;
}
/*--------------------------------------------
 find maximum h
 lets find the next sensible number
 
 x=x_0+h*alpha
 (x-x0)/alpha=sqrt(eps)/alpha
 
 --------------------------------------------*/
void PotFitO::ls_prep(type0& dfa,type0& h_norm,type0& max_a)
{
    
    h_norm=Algebra::V_mul_V<nHvars>(h,h);
    
    dfa=-Algebra::V_mul_V<nHvars>(f,h);
    
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
    
    
//    //fix this shit
    max_a=std::numeric_limits<type0>::infinity();
    max_a=MIN(max_a,find_max_alpha_A_rho_H());
    max_a=MIN(max_a,find_max_alpha_A_phi_FeH());
    max_a=MIN(max_a,find_max_alpha_A_phi_HH());
    
    
    
    max_a=MIN(max_a,find_max_alpha_A_F_H());
}
/*--------------------------------------------
 just one test run
 --------------------------------------------*/
type0 PotFitO::F(type0 alpha)
{
    Algebra::V_eq<nHvars>(x0,ff->v);
    Algebra::V_add_x_mul_V<nHvars>(alpha,h,ff->v);
    bool param_err=false;
    for(int i=0;i<nHvars;i++)
        if(std::isnan(ff->v[i]) || std::isinf(ff->v[i]))
            param_err=true;
    
    /*
    if(param_err && my_rank==0)
    {
        printf("%e\n",alpha);
        printf("%e %e %e %e %e %e %e %e %e %e %e %e %e %e\n",h[0],h[1],h[2],h[3],h[4],h[5],h[6],h[7],h[8],h[9],h[10],h[11],h[12],h[13]);
    }*/
    
    //printf("%e %e %e %e %e %e %e %e %e %e %e %e %e %e",ff->v[0],ff->v[1],ff->v[2],ff->v[3],ff->v[4],ff->v[5],ff->v[6],ff->v[7],ff->v[8],ff->v[9],ff->v[10],ff->v[11],ff->v[12],ff->v[13]);


    if(ntrial%max_ntrials==0) full_reset();
    ntrial++;
    min->init();
    min->run(nmin_steps);
    mean_rhoHs[my_conf]=ff->mean_rho_H();
    min->fin();
    type0 err_lcl=0.0,err;
    if(my_lcl_rank==0) err_lcl=coef*Algebra::pow<2>(atoms->pe-target);
    MPI_Allreduce(&err_lcl,&err,1,Vec<type0>::MPI_T,MPI_SUM,world);
    tot_err=err;
    return err;
}
/*--------------------------------------------
 just one test run
 --------------------------------------------*/
void PotFitO::F_reset()
{
    restore_x0();
    if(ntrial%max_ntrials==0) full_reset();
    ntrial++;
    min->init();
    min->run(nmin_steps);
    mean_rhoHs[my_conf]=ff->mean_rho_H();
    min->fin();
}
/*--------------------------------------------
 just one test run
 --------------------------------------------*/
type0 PotFitO::find_max_alpha_A_rho_H()
{
    type0 __max_a=std::numeric_limits<type0>::infinity();
    type0* __x0=x0+rho_H_offset;
    bool* __dof=ff->dof+rho_H_offset;
    type0* __h=h+rho_H_offset;
    type0* __dx_max=dx_max+rho_H_offset;

    __max_a=MIN(__max_a,find_max_alpha(0.0, std::numeric_limits<type0>::quiet_NaN(),__dof[0],__dx_max[0],__x0[0],__h[0]));
    __max_a=MIN(__max_a,find_max_alpha(0.0, std::numeric_limits<type0>::quiet_NaN(),__dof[1],__dx_max[1],__x0[1],__h[1]));
    
    return __max_a;
}
/*--------------------------------------------
 just one test run
 --------------------------------------------*/
type0 PotFitO::find_max_alpha_A_phi_FeH()
{
    type0 __max_a=std::numeric_limits<type0>::infinity();
    type0* __x0=x0+phi_FeH_offset;
    
    type0* __h=h+phi_FeH_offset;
    type0* __dx_max=dx_max+phi_FeH_offset;
#ifdef FeH_SPLINE
    type0 ___h,___x0,t;
    for(int i=0;i<nphi_FeH-1;i++)
    {
        ___h=___x0=0.0;
        for(int j=i+1;j<nphi_FeH;j++)
        {
            t=Algebra::pow<3>(ff->R_phi_FeH[j]-ff->R_phi_FeH[i]);
            ___h+=__h[j]*t;
            ___x0+=__x0[j]*t;
            __max_a=MIN(__max_a,find_max_alpha(std::numeric_limits<type0>::quiet_NaN(), std::numeric_limits<type0>::quiet_NaN(),true,__dx_max[i],___x0,___h));
        }
    }
    return __max_a;
#else
    bool* __dof=ff->dof+phi_FeH_offset;
    __max_a=MIN(__max_a,find_max_alpha(0.0, std::numeric_limits<type0>::quiet_NaN(),__dof[0],__dx_max[0],__x0[0],__h[0]));
    __max_a=MIN(__max_a,find_max_alpha(0.0, std::numeric_limits<type0>::quiet_NaN(),__dof[1],__dx_max[1],__x0[1],__h[1]));
    __max_a=MIN(__max_a,find_max_alpha(0.0, std::numeric_limits<type0>::quiet_NaN(),__dof[2],__dx_max[2],__x0[2],__h[2]));
    __max_a=MIN(__max_a,find_max_alpha(std::numeric_limits<type0>::quiet_NaN(), std::numeric_limits<type0>::quiet_NaN(),__dof[3],__dx_max[3],__x0[3],__h[3]));

    return __max_a;
#endif
}
/*--------------------------------------------
 just one test run
 --------------------------------------------*/
type0 PotFitO::find_max_alpha_A_phi_HH()
{
    type0 __max_a=std::numeric_limits<type0>::infinity();
    type0* __x0=x0+phi_HH_offset;
    bool* __dof=ff->dof+phi_HH_offset;
    type0* __h=h+phi_HH_offset;
    type0* __dx_max=dx_max+phi_HH_offset;

    __max_a=MIN(__max_a,find_max_alpha(std::numeric_limits<type0>::quiet_NaN(), std::numeric_limits<type0>::quiet_NaN(),__dof[0],__dx_max[0],__x0[0],__h[0]));
    __max_a=MIN(__max_a,find_max_alpha(std::numeric_limits<type0>::quiet_NaN(), std::numeric_limits<type0>::quiet_NaN(),__dof[1],__dx_max[1],__x0[1],__h[1]));
    __max_a=MIN(__max_a,find_max_alpha(0.0, std::numeric_limits<type0>::quiet_NaN(),__dof[2],__dx_max[2],__x0[2],__h[2]));
    __max_a=MIN(__max_a,find_max_alpha(0.9,1.5,__dof[3],__dx_max[3],__x0[3],__h[3]));
    
    return __max_a;
}
/*--------------------------------------------
 
 --------------------------------------------*/
type0 PotFitO::find_max_alpha_A_F_H()
{
    type0 __max_a=std::numeric_limits<type0>::infinity();
    type0* __x0=x0+F_H_offset;
    bool* __dof=ff->dof+F_H_offset;
    type0* __h=h+F_H_offset;
    type0* __dx_max=dx_max+F_H_offset;
    __max_a=MIN(__max_a,find_max_alpha(0.0, std::numeric_limits<type0>::quiet_NaN(),__dof[0],__dx_max[0],__x0[0],__h[0]));
    __max_a=MIN(__max_a,find_max_alpha(std::numeric_limits<type0>::quiet_NaN(), std::numeric_limits<type0>::quiet_NaN(),__dof[1],__dx_max[1],__x0[1],__h[1]));
    __max_a=MIN(__max_a,find_max_alpha(0.0, std::numeric_limits<type0>::quiet_NaN(),__dof[2],__dx_max[2],__x0[2],__h[2]));
    __max_a=MIN(__max_a,find_max_alpha(0.0, std::numeric_limits<type0>::quiet_NaN(),__dof[3],__dx_max[3],__x0[3],__h[3]));
    __max_a=MIN(__max_a,find_max_alpha(0.0, std::numeric_limits<type0>::quiet_NaN(),__dof[4],__dx_max[4],__x0[4],__h[4]));
    
    return __max_a;
}
/*--------------------------------------------
 
 --------------------------------------------*/
type0 PotFitO::find_max_alpha(const type0 xlo,const type0 xhi,bool dof,type0 max_dx,type0 x,type0 h)
{
    
    if(dof==false || h==0.0) return std::numeric_limits<type0>::infinity();
    type0 max_alpha=fabs(max_dx/h);
    
    if(h>0.0 && std::isnan(xhi)==false)
    {
        max_alpha=MIN(max_alpha,(xhi-x)/h);
        while(x+max_alpha*h>xhi)
        max_alpha=nextafter(max_alpha,0.0);
    }
    else if(h<0.0 && std::isnan(xlo)==false)
    {
        max_alpha=MIN(max_alpha,(xlo-x)/h);
        while(x+max_alpha*h<xlo)
        max_alpha=nextafter(max_alpha,0.0);
        
    }
    
    return max_alpha;
}
/*--------------------------------------------
 
 --------------------------------------------*/
size_t PotFitO::get_rFeH(type0*& Rs,type0*& Fs,int*& Ns)
{
    if(ntrial%max_ntrials==0) full_reset();
    ntrial++;
    min->init();
    min->run(nmin_steps);
    size_t sz=ff->get_rFeH(Rs,Fs,Ns);
    min->fin();
    
    return sz;
}
/*--------------------------------------------
 
 --------------------------------------------*/
type0* PotFitO::get_coefs()
{
    type0* coefs=NULL;
    
    Memory::alloc(coefs,nconfigs);
    coefs[my_conf]=coef;
    for(int i=0;i<nconfigs;i++)
        MPI_Bcast(&coefs[i],1,Vec<type0>::MPI_T,roots[i],world);
    
    return coefs;
}
/*--------------------------------------------
 
 --------------------------------------------*/
void PotFitO::set_coefs(type0* coefs)
{
    coef=coefs[my_conf];
}
/*------------------------------------------------------------------------------------------------------------------------------------

 ------------------------------------------------------------------------------------------------------------------------------------*/
PyObject* PotFitO::__new__(PyTypeObject* type,PyObject* args,PyObject* kwds)
{
    Object* __self=reinterpret_cast<Object*>(type->tp_alloc(type,0));
    PyObject* self=reinterpret_cast<PyObject*>(__self);
    return self;
}
/*--------------------------------------------

 --------------------------------------------*/
int PotFitO::__init__(PyObject* self,PyObject* args,PyObject* kwds)
{
    

    FuncAPI<type0[nrho_H],type0[nphi_FeH],type0[nphi_HH],type0[nF_H],
    OB<PyListObject,PyList_Type>>
    f("__init__",{
        "A_rho_H",    //0
        "A_phi_FeH",  //1
        "A_phi_HH",   //2
        "A_F_H",      //3
        "conf"        //4
    });
    if(f(args,kwds)) return -1;

    PyObject* configs=f.val<4>();
    size_t nconfigs=PyList_Size(configs);
    std::string* names;
    Memory::alloc(names,nconfigs);
    std::string* files;
    Memory::alloc(files,nconfigs);
    int* nprocs;
    Memory::alloc(nprocs,nconfigs);
    type0(*  targets)[1+__nvoigt__];
    Memory::alloc(targets,nconfigs);
    type0(*  coefs)[1+__nvoigt__];
    Memory::alloc(coefs,nconfigs);
    PyObject** funcs;
    Memory::alloc(funcs,nconfigs);


    int tot_nproc=0;
    for(size_t i=0;i<nconfigs;i++)
    {
        PyObject* config=PyList_GetItem(configs,i);
        FuncAPI<std::string,std::string,
        int,
        type0,symm<type0[__dim__][__dim__]>,
        type0,symm<type0[__dim__][__dim__]>,
        OB<PyFunctionObject,PyFunction_Type>>
        conf("configuration tuple",{"name","cfg_file","nproc","pe_target","S_target","pe_coef","S_coef","dof_func"});

        conf.val<7>()=NULL;
        conf.noptionals=1;
        if(conf(config,NULL))
        {
            Memory::dealloc(names);
            Memory::dealloc(files);
            Memory::dealloc(targets);
            Memory::dealloc(coefs);
            *funcs=NULL;
            Memory::dealloc(funcs);
            return -1;
        }
        names[i]=conf.val<0>();
        files[i]=conf.val<1>();
        nprocs[i]=conf.val<2>();
        targets[i][0]=conf.val<3>();
        Algebra::MSY_2_DyadicV(conf.val<4>(),&targets[i][1]);
        coefs[i][0]=conf.val<5>();
        Algebra::MSY_2_DyadicV(conf.val<6>(),&coefs[i][1]);
        funcs[i]=conf.val<7>();
        tot_nproc+=nprocs[i];
    }

    MPI_Comm world=MPI_COMM_WORLD;
    int comm_size;
    MPI_Comm_size(world,&comm_size);
    if(tot_nproc!=comm_size)
    {
        PyErr_Format(PyExc_TypeError,"number of processors do not match: %d vs %d",comm_size,tot_nproc);
        Memory::dealloc(names);
        Memory::dealloc(files);
        Memory::dealloc(targets);
        Memory::dealloc(coefs);
        *funcs=NULL;
        Memory::dealloc(funcs);
        return -1;
    }




    Object* __self=reinterpret_cast<Object*>(self);
    __self->potfit=new PotFitO(f.val<0>(),f.val<1>(),f.val<2>(),f.val<3>(),
                              std::move(names),files,nprocs,targets,coefs,funcs,nconfigs,world);
    Memory::dealloc(names);
    Memory::dealloc(files);
    Memory::dealloc(targets);
    Memory::dealloc(coefs);
    *funcs=NULL;
    Memory::dealloc(funcs);

    return 0;
}
/*--------------------------------------------

 --------------------------------------------*/
PyObject* PotFitO::__alloc__(PyTypeObject* type,Py_ssize_t)
{
    Object* __self=new Object;
    Py_TYPE(__self)=type;
    Py_REFCNT(__self)=1;
    __self->potfit=NULL;
    return reinterpret_cast<PyObject*>(__self);
}
/*--------------------------------------------

 --------------------------------------------*/
void PotFitO::__dealloc__(PyObject* self)
{
    Object* __self=reinterpret_cast<Object*>(self);
    delete __self->potfit;
    __self->potfit=NULL;
    delete __self;
}
/*--------------------------------------------*/
PyTypeObject PotFitO::TypeObject={PyObject_HEAD_INIT(NULL)};
/*--------------------------------------------*/
int PotFitO::setup_tp()
{
    TypeObject.tp_name="mapp.potfit_o";
    TypeObject.tp_doc="";

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
    /*
     this is a shitty hack since python does not have a slot
     for __init__, __new__, __call__, and etc. they use
     a wrapper_desriptor with a default doc here I change it
     */
    GET_WRAPPER_DOC(TypeObject,__init__)=NULL;
    return ichk;
}
/*--------------------------------------------*/
PyGetSetDef PotFitO::getset[]=EmptyPyGetSetDef(20);
/*--------------------------------------------*/
void PotFitO::setup_tp_getset()
{
    getset_A_rho_H(getset[0]);
    getset_A_phi_FeH(getset[1]);
    getset_A_phi_HH(getset[2]);
    getset_A_F_H(getset[3]);
    
    getset_dA_rho_H_max(getset[4]);
    getset_dA_phi_FeH_max(getset[5]);
    getset_dA_phi_HH_max(getset[6]);
    getset_dA_F_H_max(getset[7]);
    
    getset_A_rho_H_dof(getset[8]);
    getset_A_phi_FeH_dof(getset[9]);
    getset_A_phi_HH_dof(getset[10]);
    getset_A_F_H_dof(getset[11]);
    
    getset_mean_rho_H(getset[12]);
    getset_max_ntrials(getset[13]);
    getset_nmin_steps(getset[14]);
    getset_tol(getset[15]);
    getset_coefs(getset[16]);
    getset_errs(getset[17]);
    
    getset_RFeH(getset[18]);
}
/*--------------------------------------------*/
PyMethodDef PotFitO::methods[]=EmptyPyMethodDef(6);
/*--------------------------------------------*/
void PotFitO::setup_tp_methods()
{
    ml_test_A_rho_H(methods[0]);
    ml_test_A_phi_FeH(methods[1]);
    ml_test_A_phi_HH(methods[2]);
    ml_test_A_F_H(methods[3]);
    ml_min(methods[4]);
//    ml_reset(methods[1]);
}
/*--------------------------------------------
 
 --------------------------------------------*/
void PotFitO::getset_A_rho_H(PyGetSetDef& getset)
{
    getset.name=(char*)"A_rho_H";
    getset.doc=(char*)"";
    
    getset.get=[](PyObject* self,void*)->PyObject*
    {
        PotFitO* potfit=reinterpret_cast<Object*>(self)->potfit;
        size_t sz;
        size_t* szp=&sz;
        sz=nrho_H;
        type0* v=potfit->ff->v+rho_H_offset;
        return var<type0*>::build(v,&szp);
    };
    getset.set=[](PyObject* self,PyObject* val,void*)->int
    {
        PotFitO* potfit=reinterpret_cast<Object*>(self)->potfit;
        VarAPI<type0[nrho_H]> var("A_rho_H");
        if(var.set(val)!=0) return -1;
        memcpy(potfit->ff->v+rho_H_offset,var.val,nrho_H*sizeof(type0));
        return 0;
    };
}
/*--------------------------------------------
 
 --------------------------------------------*/
void PotFitO::getset_A_phi_FeH(PyGetSetDef& getset)
{
    getset.name=(char*)"A_phi_FeH";
    getset.doc=(char*)"";
    
    getset.get=[](PyObject* self,void*)->PyObject*
    {
        PotFitO* potfit=reinterpret_cast<Object*>(self)->potfit;
        size_t sz;
        size_t* szp=&sz;
        sz=nphi_FeH;
        type0* v=potfit->ff->v+phi_FeH_offset;
        return var<type0*>::build(v,&szp);
    };
    getset.set=[](PyObject* self,PyObject* val,void*)->int
    {
        PotFitO* potfit=reinterpret_cast<Object*>(self)->potfit;
        VarAPI<type0[nphi_FeH]> var("A_phi_FeH");
        if(var.set(val)!=0) return -1;
        memcpy(potfit->ff->v+phi_FeH_offset,var.val,nphi_FeH*sizeof(type0));
        return 0;
    };
}
/*--------------------------------------------
 
 --------------------------------------------*/
void PotFitO::getset_A_phi_HH(PyGetSetDef& getset)
{
    getset.name=(char*)"A_phi_HH";
    getset.doc=(char*)"";
    
    getset.get=[](PyObject* self,void*)->PyObject*
    {
        PotFitO* potfit=reinterpret_cast<Object*>(self)->potfit;
        size_t sz;
        size_t* szp=&sz;
        sz=nphi_HH;
        type0* v=potfit->ff->v+phi_HH_offset;
        return var<type0*>::build(v,&szp);
    };
    getset.set=[](PyObject* self,PyObject* val,void*)->int
    {
        PotFitO* potfit=reinterpret_cast<Object*>(self)->potfit;
        VarAPI<type0[nphi_HH]> var("A_phi_HH");
        if(var.set(val)!=0) return -1;
        memcpy(potfit->ff->v+phi_HH_offset,var.val,nphi_HH*sizeof(type0));
        return 0;
    };
}
/*--------------------------------------------
 
 --------------------------------------------*/
void PotFitO::getset_A_F_H(PyGetSetDef& getset)
{
    getset.name=(char*)"A_F_H";
    getset.doc=(char*)"";
    
    getset.get=[](PyObject* self,void*)->PyObject*
    {
        PotFitO* potfit=reinterpret_cast<Object*>(self)->potfit;
        size_t sz;
        size_t* szp=&sz;
        sz=nF_H;
        type0* v=potfit->ff->v+F_H_offset;
        return var<type0*>::build(v,&szp);
    };
    getset.set=[](PyObject* self,PyObject* val,void*)->int
    {
        PotFitO* potfit=reinterpret_cast<Object*>(self)->potfit;
        VarAPI<type0[nF_H]> var("A_F_H");
        if(var.set(val)!=0) return -1;
        memcpy(potfit->ff->v+F_H_offset,var.val,nF_H*sizeof(type0));
        return 0;
    };
}
/*--------------------------------------------
 
 --------------------------------------------*/
void PotFitO::getset_dA_rho_H_max(PyGetSetDef& getset)
{
    getset.name=(char*)"dA_rho_H_max";
    getset.doc=(char*)"";
    
    getset.get=[](PyObject* self,void*)->PyObject*
    {
        PotFitO* potfit=reinterpret_cast<Object*>(self)->potfit;
        size_t sz;
        size_t* szp=&sz;
        sz=nrho_H;
        type0* v=potfit->dx_max+rho_H_offset;
        return var<type0*>::build(v,&szp);
    };
    getset.set=[](PyObject* self,PyObject* val,void*)->int
    {
        PotFitO* potfit=reinterpret_cast<Object*>(self)->potfit;
        VarAPI<type0[nrho_H]> var("dA_rho_H_max");
        if(var.set(val)!=0) return -1;
        memcpy(potfit->dx_max+rho_H_offset,var.val,nrho_H*sizeof(type0));
        return 0;
    };
}
/*--------------------------------------------
 
 --------------------------------------------*/
void PotFitO::getset_dA_phi_FeH_max(PyGetSetDef& getset)
{
    getset.name=(char*)"dA_phi_FeH_max";
    getset.doc=(char*)"";
    
    getset.get=[](PyObject* self,void*)->PyObject*
    {
        PotFitO* potfit=reinterpret_cast<Object*>(self)->potfit;
        size_t sz;
        size_t* szp=&sz;
        sz=nphi_FeH;
        type0* v=potfit->dx_max+phi_FeH_offset;
        return var<type0*>::build(v,&szp);
    };
    getset.set=[](PyObject* self,PyObject* val,void*)->int
    {
        PotFitO* potfit=reinterpret_cast<Object*>(self)->potfit;
        VarAPI<type0[nphi_FeH]> var("dA_phi_FeH_max");
        if(var.set(val)!=0) return -1;
        memcpy(potfit->dx_max+phi_FeH_offset,var.val,nphi_FeH*sizeof(type0));
        return 0;
    };
}
/*--------------------------------------------
 
 --------------------------------------------*/
void PotFitO::getset_dA_phi_HH_max(PyGetSetDef& getset)
{
    getset.name=(char*)"dA_phi_HH_max";
    getset.doc=(char*)"";
    
    getset.get=[](PyObject* self,void*)->PyObject*
    {
        PotFitO* potfit=reinterpret_cast<Object*>(self)->potfit;
        size_t sz;
        size_t* szp=&sz;
        sz=nphi_HH;
        type0* v=potfit->dx_max+phi_HH_offset;
        return var<type0*>::build(v,&szp);
    };
    getset.set=[](PyObject* self,PyObject* val,void*)->int
    {
        PotFitO* potfit=reinterpret_cast<Object*>(self)->potfit;
        VarAPI<type0[nphi_HH]> var("dA_phi_HH_max");
        if(var.set(val)!=0) return -1;
        memcpy(potfit->dx_max+phi_HH_offset,var.val,nphi_HH*sizeof(type0));
        return 0;
    };
}
/*--------------------------------------------
 
 --------------------------------------------*/
void PotFitO::getset_dA_F_H_max(PyGetSetDef& getset)
{
    getset.name=(char*)"dA_F_H_max";
    getset.doc=(char*)"";
    
    getset.get=[](PyObject* self,void*)->PyObject*
    {
        PotFitO* potfit=reinterpret_cast<Object*>(self)->potfit;
        size_t sz;
        size_t* szp=&sz;
        sz=nF_H;
        type0* v=potfit->dx_max+F_H_offset;
        return var<type0*>::build(v,&szp);
    };
    getset.set=[](PyObject* self,PyObject* val,void*)->int
    {
        PotFitO* potfit=reinterpret_cast<Object*>(self)->potfit;
        VarAPI<type0[nF_H]> var("dA_F_H_max");
        if(var.set(val)!=0) return -1;
        memcpy(potfit->dx_max+F_H_offset,var.val,nF_H*sizeof(type0));
        return 0;
    };
}
/*--------------------------------------------

 --------------------------------------------*/
void PotFitO::getset_A_rho_H_dof(PyGetSetDef& getset)
{
    getset.name=(char*)"A_rho_H_dof";
    getset.doc=(char*)"";

    getset.get=[](PyObject* self,void*)->PyObject*
    {
        PotFitO* potfit=reinterpret_cast<Object*>(self)->potfit;
        size_t sz;
        size_t* szp=&sz;
        sz=nrho_H;
        bool* v=potfit->ff->dof+rho_H_offset;
        return var<bool*>::build(v,&szp);
    };
    getset.set=[](PyObject* self,PyObject* val,void*)->int
    {
        PotFitO* potfit=reinterpret_cast<Object*>(self)->potfit;
        VarAPI<bool[nrho_H]> var("A_rho_H_dof");
        if(var.set(val)!=0) return -1;
        memcpy(potfit->ff->dof+rho_H_offset,var.val,nrho_H*sizeof(bool));
        return 0;
    };
}
/*--------------------------------------------
 
 --------------------------------------------*/
void PotFitO::getset_A_phi_FeH_dof(PyGetSetDef& getset)
{
    getset.name=(char*)"A_phi_FeH_dof";
    getset.doc=(char*)"";
    
    getset.get=[](PyObject* self,void*)->PyObject*
    {
        PotFitO* potfit=reinterpret_cast<Object*>(self)->potfit;
        size_t sz;
        size_t* szp=&sz;
        sz=nphi_FeH;
        bool* v=potfit->ff->dof+phi_FeH_offset;
        return var<bool*>::build(v,&szp);
    };
    getset.set=[](PyObject* self,PyObject* val,void*)->int
    {
        PotFitO* potfit=reinterpret_cast<Object*>(self)->potfit;
        VarAPI<bool[nphi_FeH]> var("A_phi_FeH_dof");
        if(var.set(val)!=0) return -1;
        memcpy(potfit->ff->dof+phi_FeH_offset,var.val,nphi_FeH*sizeof(bool));
        return 0;
    };
}
/*--------------------------------------------
 
 --------------------------------------------*/
void PotFitO::getset_A_phi_HH_dof(PyGetSetDef& getset)
{
    getset.name=(char*)"A_phi_HH_dof";
    getset.doc=(char*)"";
    
    getset.get=[](PyObject* self,void*)->PyObject*
    {
        PotFitO* potfit=reinterpret_cast<Object*>(self)->potfit;
        size_t sz;
        size_t* szp=&sz;
        sz=nphi_HH;
        bool* v=potfit->ff->dof+phi_HH_offset;
        return var<bool*>::build(v,&szp);
    };
    getset.set=[](PyObject* self,PyObject* val,void*)->int
    {
        PotFitO* potfit=reinterpret_cast<Object*>(self)->potfit;
        VarAPI<bool[nphi_HH]> var("A_phi_HH_dof");
        if(var.set(val)!=0) return -1;
        memcpy(potfit->ff->dof+phi_HH_offset,var.val,nphi_HH*sizeof(bool));
        return 0;
    };
}
/*--------------------------------------------
 
 --------------------------------------------*/
void PotFitO::getset_A_F_H_dof(PyGetSetDef& getset)
{
    getset.name=(char*)"A_F_H_dof";
    getset.doc=(char*)"";
    
    getset.get=[](PyObject* self,void*)->PyObject*
    {
        PotFitO* potfit=reinterpret_cast<Object*>(self)->potfit;
        size_t sz;
        size_t* szp=&sz;
        sz=nF_H;
        bool* v=potfit->ff->dof+F_H_offset;
        return var<bool*>::build(v,&szp);
    };
    getset.set=[](PyObject* self,PyObject* val,void*)->int
    {
        PotFitO* potfit=reinterpret_cast<Object*>(self)->potfit;
        VarAPI<bool[nF_H]> var("A_F_H_dof");
        if(var.set(val)!=0) return -1;
        memcpy(potfit->ff->dof+F_H_offset,var.val,nF_H*sizeof(bool));
        return 0;
    };
}
/*--------------------------------------------
 
 --------------------------------------------*/
void PotFitO::getset_mean_rho_H(PyGetSetDef& getset)
{
    getset.name=(char*)"mean_rho_H";
    getset.doc=(char*)"";
    
    getset.get=[](PyObject* self,void*)->PyObject*
    {
        
        PotFitO* potfit=reinterpret_cast<Object*>(self)->potfit;
        size_t sz;
        size_t* szp=&sz;
        sz=potfit->nconfigs;
        potfit->calc_mean_rho_Hs();
        type0* v=potfit->mean_rhoHs;
        return var<type0*>::build(v,&szp);
    };
    getset.set=[](PyObject*,PyObject*,void*)->int
    {
        PyErr_SetString(PyExc_TypeError,"readonly attribute");
        return -1;
    };
}
/*--------------------------------------------
 
 --------------------------------------------*/
void PotFitO::getset_max_ntrials(PyGetSetDef& getset)
{
    getset.name=(char*)"max_ntrials";
    getset.doc=(char*)"";
    
    getset.get=[](PyObject* self,void*)->PyObject*
    {
        return var<int>::build(reinterpret_cast<Object*>(self)->potfit->max_ntrials);
    };
    getset.set=[](PyObject* self,PyObject* val,void*)->int
    {
        VarAPI<int> var("max_ntrials");
        if(var.set(val)!=0) return -1;
        reinterpret_cast<Object*>(self)->potfit->max_ntrials=var.val;
        return 0;
    };
}
/*--------------------------------------------
 
 --------------------------------------------*/
void PotFitO::getset_nmin_steps(PyGetSetDef& getset)
{
    getset.name=(char*)"nmin_steps";
    getset.doc=(char*)"";
    
    getset.get=[](PyObject* self,void*)->PyObject*
    {
        return var<int>::build(reinterpret_cast<Object*>(self)->potfit->nmin_steps);
    };
    getset.set=[](PyObject* self,PyObject* val,void*)->int
    {
        VarAPI<int> var("nmin_steps");
        if(var.set(val)!=0) return -1;
        reinterpret_cast<Object*>(self)->potfit->nmin_steps=var.val;
        return 0;
    };
}
/*--------------------------------------------
 
 --------------------------------------------*/
void PotFitO::getset_tol(PyGetSetDef& getset)
{
    getset.name=(char*)"tol";
    getset.doc=(char*)"";
    
    getset.get=[](PyObject* self,void*)->PyObject*
    {
        return var<type0>::build(reinterpret_cast<Object*>(self)->potfit->err_tol);
    };
    getset.set=[](PyObject* self,PyObject* val,void*)->int
    {
        VarAPI<type0> var("tol");
        if(var.set(val)!=0) return -1;
        reinterpret_cast<Object*>(self)->potfit->err_tol=var.val;
        return 0;
    };
}
/*--------------------------------------------
 
 --------------------------------------------*/
void PotFitO::getset_errs(PyGetSetDef& getset)
{
    getset.name=(char*)"errs";
    getset.doc=(char*)"";
    
    getset.get=[](PyObject* self,void*)->PyObject*
    {
        PotFitO* potfit=reinterpret_cast<Object*>(self)->potfit;
        size_t sz;
        size_t* szp=&sz;
        sz=potfit->nconfigs;
        PyObject* op=var<type0*>::build(potfit->errs,&szp);
        return op;
    };
    getset.set=[](PyObject* self,PyObject* val,void*)->int
    {
        PyErr_SetString(PyExc_TypeError,"readonly attribute");
        return -1;
    };
}
/*--------------------------------------------
 
 --------------------------------------------*/
void PotFitO::getset_coefs(PyGetSetDef& getset)
{
    getset.name=(char*)"coefs";
    getset.doc=(char*)"";
    
    getset.get=[](PyObject* self,void*)->PyObject*
    {
        PotFitO* potfit=reinterpret_cast<Object*>(self)->potfit;
        size_t sz;
        size_t* szp=&sz;
        sz=potfit->nconfigs;
        type0* v=potfit->get_coefs();
        PyObject* op=var<type0*>::build(v,&szp);
        Memory::dealloc(v);
        return op;
    };
    getset.set=[](PyObject* self,PyObject* val,void*)->int
    {
        VarAPI<type0*> var("coefs");
        if(var.set(val)!=0) return -1;
        PotFitO* potfit=reinterpret_cast<Object*>(self)->potfit;
        if(var.__var__.size!=static_cast<size_t>(potfit->nconfigs))
        {
            PyErr_SetString(PyExc_TypeError,"size mismatch");
            return -1;
        }
        potfit->set_coefs(var.val);
        return 0;
    };
}
/*--------------------------------------------
 
 --------------------------------------------*/
void PotFitO::getset_RFeH(PyGetSetDef& getset)
{
    getset.name=(char*)"RFeH";
    getset.doc=(char*)"";
 
    getset.get=[](PyObject* self,void*)->PyObject*
    {
 
        PotFitO* potfit=reinterpret_cast<Object*>(self)->potfit;
        if(potfit->nconfigs!=1)
        {
            PyErr_SetString(PyExc_TypeError,"cannot calculate RFeH for multiple configurations");
            return NULL;
        }
        
        type0* Rs=NULL;
        type0* Fs=NULL;
        int* Ns=NULL;
        size_t sz=potfit->get_rFeH(Rs,Fs,Ns);
        size_t* szp=&sz;

        PyObject* py_obj=PyList_New(3);
        PyList_SET_ITEM(py_obj,0,var<type0*>::build(Rs,&szp));
        PyList_SET_ITEM(py_obj,1,var<type0*>::build(Fs,&szp));
        PyList_SET_ITEM(py_obj,2,var<int*>::build(Ns,&szp));
 
        Memory::dealloc(Rs);
        Memory::dealloc(Fs);
        Memory::dealloc(Ns);
 
        return py_obj;
 
 
    };
    getset.set=[](PyObject*,PyObject*,void*)->int
    {
        PyErr_SetString(PyExc_TypeError,"readonly attribute");
        return -1;
    };
}
/*--------------------------------------------

 --------------------------------------------*/
void PotFitO::ml_min(PyMethodDef& tp_methods)
{
    tp_methods.ml_flags=METH_VARARGS | METH_KEYWORDS;
    tp_methods.ml_name="min";
    tp_methods.ml_meth=(PyCFunction)(PyCFunctionWithKeywords)(
    [](PyObject* self,PyObject* args,PyObject* kwds)->PyObject*
    {

        FuncAPI<int>f("min",{"nsteps"});
        if(f(args,kwds)) return NULL;

        reinterpret_cast<Object*>(self)->potfit->min_cg(f.val<0>());

        Py_RETURN_NONE;
    });


    tp_methods.ml_doc=(char*)R"---(
    quick function for calculateing phonon freq
    use with caution

    )---";
}
///*--------------------------------------------
//
// --------------------------------------------*/
//void PotFitO::ml_reset(PyMethodDef& tp_methods)
//{
//    tp_methods.ml_flags=METH_VARARGS | METH_KEYWORDS;
//    tp_methods.ml_name="reset";
//    tp_methods.ml_meth=(PyCFunction)(PyCFunctionWithKeywords)(
//    [](PyObject* self,PyObject* args,PyObject* kwds)->PyObject*
//    {
//
//        FuncAPI<>f("reset");
//        if(f(args,kwds)) return NULL;
//
//        reinterpret_cast<Object*>(self)->potfit->full_reset();
//
//        Py_RETURN_NONE;
//    });
//
//
//    tp_methods.ml_doc="";
//}
//
//
/*--------------------------------------------

 --------------------------------------------*/
void PotFitO::ml_test_A_rho_H(PyMethodDef& tp_methods)
{
    tp_methods.ml_flags=METH_VARARGS | METH_KEYWORDS;
    tp_methods.ml_name="test_A_rho_H";
    tp_methods.ml_meth=(PyCFunction)(PyCFunctionWithKeywords)(
    [](PyObject* self,PyObject* args,PyObject* kwds)->PyObject*
    {
        FuncAPI<int,type0,int>f("test_A_rho_H",{"comp","delta","n"});
        if(f(args,kwds)) return NULL;
        reinterpret_cast<Object*>(self)->potfit->gen_test(rho_H_offset+f.val<0>(),f.val<1>(),f.val<2>());
        Py_RETURN_NONE;
    });
    tp_methods.ml_doc="";
}
/*--------------------------------------------

 --------------------------------------------*/
void PotFitO::ml_test_A_phi_FeH(PyMethodDef& tp_methods)
{
    tp_methods.ml_flags=METH_VARARGS | METH_KEYWORDS;
    tp_methods.ml_name="test_A_phi_FeH";
    tp_methods.ml_meth=(PyCFunction)(PyCFunctionWithKeywords)(
    [](PyObject* self,PyObject* args,PyObject* kwds)->PyObject*
    {
        FuncAPI<int,type0,int>f("test_A_phi_FeH",{"comp","delta","n"});
        if(f(args,kwds)) return NULL;
        reinterpret_cast<Object*>(self)->potfit->gen_test(phi_HH_offset+f.val<0>(),f.val<1>(),f.val<2>());
        Py_RETURN_NONE;
    });
    tp_methods.ml_doc="";
}
/*--------------------------------------------

 --------------------------------------------*/
void PotFitO::ml_test_A_phi_HH(PyMethodDef& tp_methods)
{
    tp_methods.ml_flags=METH_VARARGS | METH_KEYWORDS;
    tp_methods.ml_name="test_A_phi_HH";
    tp_methods.ml_meth=(PyCFunction)(PyCFunctionWithKeywords)(
    [](PyObject* self,PyObject* args,PyObject* kwds)->PyObject*
    {
        FuncAPI<int,type0,int>f("test_A_phi_HH",{"comp","delta","n"});
        if(f(args,kwds)) return NULL;
        reinterpret_cast<Object*>(self)->potfit->gen_test(phi_FeH_offset+f.val<0>(),f.val<1>(),f.val<2>());
        Py_RETURN_NONE;
    });
    tp_methods.ml_doc="";
}
/*--------------------------------------------

 --------------------------------------------*/
void PotFitO::ml_test_A_F_H(PyMethodDef& tp_methods)
{
    tp_methods.ml_flags=METH_VARARGS | METH_KEYWORDS;
    tp_methods.ml_name="test_A_F_H";
    tp_methods.ml_meth=(PyCFunction)(PyCFunctionWithKeywords)(
    [](PyObject* self,PyObject* args,PyObject* kwds)->PyObject*
    {
        FuncAPI<int,type0,int>f("test_A_F_H",{"comp","delta","n"});
        if(f(args,kwds)) return NULL;
        reinterpret_cast<Object*>(self)->potfit->gen_test(F_H_offset+f.val<0>(),f.val<1>(),f.val<2>());
        Py_RETURN_NONE;
    });
    tp_methods.ml_doc="";
}
