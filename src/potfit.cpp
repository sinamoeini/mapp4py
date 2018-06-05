#include "potfit.h"
#include "import_cfg.h"
#include "min.h"
#include "hydride_const.h"


using namespace MAPP_NS;
const char* PotFit::err_msgs[]=
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
PotFit::PotFit(
type0(***&& __rho_AR)[2],
size_t**&& __rho_sz,
type0(***&& __phi_AR)[2],
size_t**&& __phi_sz,
type0(*&& __F_A)[3],
               
size_t*&& __phi_const_list,size_t __phi_const_list_sz,
size_t*&& __phi_nconst_list,size_t __phi_nconst_list_sz,
type0*&& __phi_const0,type0**&& __A,type0*&& __rho0s,type0**&& __B,type0**&& __L,type0**&& __D0,type0***&& __C,
               
               
               
std::string*&& __names_str,
std::string*& files,
int*& nprocs,
type0(*& __targets)[1+__nvoigt__],
type0(*& __coefs)[1+__nvoigt__],
PyObject** dof_funcs,
size_t __nconfigs,
MPI_Comm& __world):
atoms(NULL),
rho_sz(__rho_sz),
phi_sz(__phi_sz),
world(__world)
{
    max_ntrials=5;
    max_AF=0.5;
    max_rhobF=0.5;
    max_alphaF=0.1;
    max_Arho=0.5;
    max_Aphi=0.1;
    
    
    
    
    hyd_const=new HydrideConst(this,__phi_AR[1][0],__phi_sz[1][0],__rho_AR[1][0],__rho_sz[1][0],
                               std::move(__phi_const_list),__phi_const_list_sz,
                               std::move(__phi_nconst_list),__phi_nconst_list_sz,
                               std::move(__phi_const0),std::move(__A),std::move(__rho0s),std::move(__B),std::move(__L),std::move(__D0),std::move(__C));
    
    
    nconfigs=static_cast<int>(__nconfigs);
    names_str=__names_str;
    __names_str=NULL;
    Memory::alloc(errs,__nconfigs);
    Memory::alloc(roots,__nconfigs);
    Memory::alloc(names,__nconfigs);
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
    
    
    
    
    
    
    
    
    ntrial=0;
    __rho_sz=NULL;
    __phi_sz=NULL;
    
    Algebra::V_eq<1+__nvoigt__>(__targets[my_conf],target);
    Algebra::V_eq<1+__nvoigt__>(__coefs[my_conf],coef);
    

    
    
    alloc(__rho_AR,__phi_AR,__F_A);
    
    
    /* loading the configuration */
    ImportCFGMD read(*my_world);
    atoms=read(files[my_conf].c_str());
    if(dof_funcs[my_conf]) atoms->DO(dof_funcs[my_conf]);
    
    
    ff=new ForceFieldEAMFit(atoms,rho_A,rho_R,rho_sz,phi_A,phi_R,phi_sz,F_A);
    
    bool H_dof[__dim__][__dim__];
    for(int i=0;i<__dim__;i++) for(int j=0;j<__dim__;j++) H_dof[i][j]=false;
    H_dof[0][0]=H_dof[1][1]=H_dof[2][2]=true;
    
    
    Xorig.~VecTens();
    new (&Xorig) VecTens<type0,1>(atoms,true,__dim__);
    X0.~VecTens();
    new (&X0) VecTens<type0,1>(atoms,true,__dim__);
    
    min_ls=new LineSearchBrent();
    min=new MinCGFit(0.0,H_dof,false,0.4,min_ls,Xorig.vecs[0],X0.vecs[0]);
    min->atoms=atoms;
    min->ff=ff;
    
    memcpy(Xorig.vecs[0]->begin(),atoms->x->begin(),sizeof(type0)*__dim__*atoms->natms_lcl);
    memcpy(&Xorig.A[0][0],&atoms->H[0][0],sizeof(type0)*__dim__*__dim__);
    
    hyd_const->prep();
    hyd_const->update();
    
    
     
}
/*--------------------------------------------
 destructor
 --------------------------------------------*/
PotFit::~PotFit()
{
    X0.~VecTens();
    Xorig.~VecTens();
    delete ff;
    delete min_ls;
    delete min;
    delete atoms;
    dealloc();
    Memory::dealloc(rho_sz);
    Memory::dealloc(phi_sz);
    Memory::dealloc(roots);
    Memory::dealloc(errs);
    *names=NULL;
    Memory::dealloc(names);
    Memory::dealloc(names_str);
}
/*--------------------------------------------
 
 --------------------------------------------*/
void PotFit::alloc(
type0(***& __rho_AR)[2],
type0(***& __phi_AR)[2],
type0(*& __F_A)[3])
{
    size_t nelems=2;
    size_t phi_tot_sz=0;
    size_t rho_tot_sz=0;
    
    for(size_t ielem=0;ielem<nelems;ielem++)
        for(size_t jelem=0;jelem<nelems;jelem++)
        {
            
            rho_tot_sz+=rho_sz[ielem][jelem];
            if(jelem<=ielem)
                phi_tot_sz+=phi_sz[ielem][jelem];
        }
    size_t tot_sz=rho_tot_sz+phi_tot_sz;
    nvars=static_cast<int>(rho_tot_sz+phi_tot_sz+3*2);
    
    type0* __Rs=NULL;
    Memory::alloc(__Rs,tot_sz);
    Memory::alloc(rho_R,nelems,nelems);
    Memory::alloc(phi_R,nelems,nelems);
    
    type0* __As=NULL;
    Memory::alloc(__As,nvars);
    Memory::alloc(rho_A,nelems,nelems);
    Memory::alloc(phi_A,nelems,nelems);
    Memory::alloc(F_A,nelems);
    
    type0* __A0s=NULL;
    Memory::alloc(__A0s,nvars);
    Memory::alloc(rho_A0,nelems,nelems);
    Memory::alloc(phi_A0,nelems,nelems);
    Memory::alloc(F_A0,nelems);
    
    type0* __As_h=NULL;
    Memory::alloc(__As_h,nvars);
    Memory::alloc(rho_A_h,nelems,nelems);
    Memory::alloc(phi_A_h,nelems,nelems);
    Memory::alloc(F_A_h,nelems);
    
    
    bool* __dofs=NULL;
    Memory::alloc(__dofs,nvars);
    for(int i=0;i<nvars;i++) __dofs[i]=true;
    Memory::alloc(rho_A_dof,nelems,nelems);
    Memory::alloc(phi_A_dof,nelems,nelems);
    Memory::alloc(F_A_dof,nelems);
    
    type0* __dAs_tot=NULL;
    Memory::alloc(__dAs_tot,tot_sz+6);
    Memory::alloc(rho_A_f,nelems,nelems);
    Memory::alloc(phi_A_f,nelems,nelems);
    Memory::alloc(F_A_f,nelems);
    
    
    type0 (* __dAs)[1+__nvoigt__]=NULL;
    Memory::alloc(__dAs,nvars);
    Memory::alloc(drho_A,nelems,nelems);
    Memory::alloc(dphi_A,nelems,nelems);
    Memory::alloc(dF_A,nelems);
    
    
    
    

    
    size_t __i=0;
    for(size_t ielem=0;ielem<nelems;ielem++)
        for(size_t jelem=0;jelem<nelems;jelem++)
        {
            rho_R[ielem][jelem]=__Rs+__i;
            rho_A[ielem][jelem]=__As+__i;
            rho_A0[ielem][jelem]=__A0s+__i;
            rho_A_h[ielem][jelem]=__As_h+__i;
            rho_A_f[ielem][jelem]=__dAs_tot+__i;
            rho_A_dof[ielem][jelem]=__dofs+__i;
            drho_A[ielem][jelem]=__dAs+__i;
            for(size_t i=0;i<rho_sz[ielem][jelem];i++)
            {
                __Rs[__i]=__rho_AR[ielem][jelem][i][1];
                __As[__i]=__rho_AR[ielem][jelem][i][0];
                __i++;
            }
            sort_AR_ij(rho_A[ielem][jelem],rho_R[ielem][jelem],rho_A_dof[ielem][jelem],rho_sz[ielem][jelem]);
        }
    
    for(size_t ielem=0;ielem<nelems;ielem++)
        for(size_t jelem=0;jelem<ielem+1;jelem++)
        {
            phi_R[jelem][ielem]=phi_R[ielem][jelem]=__Rs+__i;
            phi_A[jelem][ielem]=phi_A[ielem][jelem]=__As+__i;
            phi_A0[ielem][jelem]=__A0s+__i;
            phi_A_h[ielem][jelem]=__As_h+__i;
            phi_A_f[jelem][ielem]=phi_A_f[ielem][jelem]=__dAs_tot+__i;
            phi_A_dof[jelem][ielem]=phi_A_dof[ielem][jelem]=__dofs+__i;
            dphi_A[jelem][ielem]=dphi_A[ielem][jelem]=__dAs+__i;
            
            for(size_t i=0;i<phi_sz[ielem][jelem];i++)
            {
                __Rs[__i]=__phi_AR[ielem][jelem][i][1];
                __As[__i]=__phi_AR[ielem][jelem][i][0];
                __i++;
            }
            sort_AR_ij(phi_A[ielem][jelem],phi_R[ielem][jelem],phi_A_dof[ielem][jelem],phi_sz[ielem][jelem]);
        }
    
    for(size_t ielem=0;ielem<nelems;ielem++)
    {
        F_A[ielem]=__As+__i;
        F_A0[ielem]=__A0s+__i;
        F_A_h[ielem]=__As_h+__i;
        F_A_f[ielem]=__dAs_tot+__i;
        F_A_dof[ielem]=__dofs+__i;
        dF_A[ielem]=__dAs+__i;
        for(size_t i=0;i<3;i++)
        {
            __As[__i]=__F_A[ielem][i];
            __i++;
        }
    }
    
    
    
    Memory::dealloc(__rho_AR);
    __rho_AR=NULL;
    Memory::dealloc(__phi_AR);
    __phi_AR=NULL;
    Memory::dealloc(__F_A);
    __F_A=NULL;
    
    
    x=**rho_A;
    x0=**rho_A0;
    h=**rho_A_h;
    f_lcl=**rho_A_f;
    Memory::alloc(f,nvars);
    dof=**rho_A_dof;
    dA=**drho_A;
    Memory::alloc(f0,nvars);
    Memory::alloc(h0,nvars);
}
/*--------------------------------------------
 
 --------------------------------------------*/
void PotFit::dealloc()
{
    **rho_R=NULL;
    **rho_A=NULL;
    **rho_A0=NULL;
    **rho_A_h=NULL;
    **rho_A_dof=NULL;
    **rho_A_f=NULL;
    **drho_A=NULL;
    
    **phi_R=NULL;
    **phi_A=NULL;
    **phi_A0=NULL;
    **phi_A_h=NULL;
    **phi_A_dof=NULL;
    **phi_A_f=NULL;
    **dphi_A=NULL;
    
    
    *F_A=NULL;
    *F_A0=NULL;
    *F_A_h=NULL;
    *F_A_dof=NULL;
    *F_A_f=NULL;
    *dF_A=NULL;
    
    
    
    Memory::dealloc(rho_R);
    Memory::dealloc(rho_A);
    Memory::dealloc(rho_A0);
    Memory::dealloc(rho_A_h);
    Memory::dealloc(rho_A_f);
    Memory::dealloc(rho_A_dof);
    Memory::dealloc(drho_A);

    Memory::dealloc(phi_R);
    Memory::dealloc(phi_A);
    Memory::dealloc(phi_A0);
    Memory::dealloc(phi_A_f);
    Memory::dealloc(phi_A_h);
    Memory::dealloc(phi_A_dof);
    Memory::dealloc(dphi_A);
    
    Memory::dealloc(F_A);
    Memory::dealloc(F_A0);
    Memory::dealloc(F_A_h);
    Memory::dealloc(F_A_f);
    Memory::dealloc(F_A_dof);
    Memory::dealloc(dF_A);
    
    
    Memory::dealloc(h0);
    Memory::dealloc(f0);
    Memory::dealloc(x);
    Memory::dealloc(x0);
    Memory::dealloc(f_lcl);
    Memory::dealloc(f);
    Memory::dealloc(h);
    Memory::dealloc(dof);
    Memory::dealloc(dA);
}
/*--------------------------------------------
 destructor
 --------------------------------------------*/
void PotFit::sort_AR_ij(type0*& A,type0*& R,bool*& dof,size_t sz)
{
    if(sz==0) return;
    size_t* key=NULL;
    Memory::alloc(key,sz);
    for(size_t i=0;i<sz;i++) key[i]=i;
    
    XMath::quicksort(key,key+sz,
    [&R](size_t* rank_i,size_t* rank_j){return (R[*rank_i]>R[*rank_j]);},
    [&A,&R,&dof](size_t* rank_i,size_t* rank_j)
    {
        std::swap(R[*rank_i],R[*rank_j]);
        std::swap(A[*rank_i],A[*rank_j]);
        std::swap(dof[*rank_i],dof[*rank_j]);
    });
    
    Memory::dealloc(key);
}
/*--------------------------------------------
 just one test run
 --------------------------------------------*/
void PotFit::store_x0()
{
    memcpy(x0,x,sizeof(type0)*nvars);
    memcpy(X0.vecs[0]->begin(),atoms->x->begin(),sizeof(type0)*__dim__*atoms->natms_lcl);
    memcpy(&X0.A[0][0],&atoms->H[0][0],sizeof(type0)*__dim__*__dim__);
}
/*--------------------------------------------
 just one test run
 --------------------------------------------*/
void PotFit::restore_x0()
{
    memcpy(x,x0,sizeof(type0)*nvars);
    memcpy(atoms->x->begin(),X0.vecs[0]->begin(),sizeof(type0)*__dim__*atoms->natms_lcl);
    memcpy(&atoms->H[0][0],&X0.A[0][0],sizeof(type0)*__dim__*__dim__);
    atoms->update_H();
}
/*--------------------------------------------
 just one test run
 --------------------------------------------*/
void PotFit::full_reset()
{
    memcpy(atoms->x->begin(),Xorig.vecs[0]->begin(),sizeof(type0)*__dim__*atoms->natms_lcl);
    memcpy(&atoms->H[0][0],&Xorig.A[0][0],sizeof(type0)*__dim__*__dim__);
    atoms->update_H();
}
/*--------------------------------------------
 just one test run
 --------------------------------------------*/
ThermoDynamics PotFit::get_thermo()
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
    return ThermoDynamics(6,"ERR",tot_err,names[0],errs[0]);
    
}
/*--------------------------------------------
 just one test run
 --------------------------------------------*/
void PotFit::min_cg(int nsteps)
{
    type0 curr_err,prev_err,f0_f0,f_f0,f_f,ratio,alpha,tol=1.0e-12;
    
    
    
    curr_err=iter();
    memcpy(h,f,sizeof(type0)*nvars);
    f0_f0=dot(f,f);
    
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
        memcpy(f0,f,sizeof(type0)*nvars);
        prev_err=curr_err;
        f_h=dot(f,h);
        if(f_h<0.0)
        {
            memcpy(h,f,sizeof(type0)*nvars);
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
        
        if(prev_err-curr_err<tol) err=MIN_S_TOLERANCE;
        if(istep+1==nsteps) err=MIN_F_MAX_ITER;
        if(err) continue;
        
        f_f=dot(f,f);
        f_f0=dot(f,f0);
        ratio=(f_f-f_f0)/(f0_f0);
        for(int j=0;j<nvars;j++) h[j]=ratio*h[j]+f[j];
        //advance h here 
        f0_f0=f_f;
    }

    thermo.fin();
    fprintf(MAPP::mapp_out,"%s",err_msgs[err]);
}
/*--------------------------------------------
 just one test run
 --------------------------------------------*/
type0 PotFit::iter()
{
    if(ntrial%max_ntrials==0) full_reset();
    ntrial++;
    min->init();
    min->run(1000);
    ff->calc_deriv(rho_A_dof,drho_A,phi_A_dof,dphi_A,F_A_dof,dF_A);
    min->fin();
    curr_diff[0]=atoms->pe;
    Algebra::MSY_2_DyadicV(atoms->S_pe,&curr_diff[1]);
    type0 err,err_lcl;
    err_lcl=val(curr_diff,target,coef);
    calc_f(curr_diff);
    hyd_const->prep();
    hyd_const->adj_deriv();
    Algebra::V_sub<1+__nvoigt__>(target,curr_diff);
    
    
    if(my_lcl_rank!=0)
    {
        for(int i=0;i<nvars;i++) f_lcl[i]=0;
        err_lcl=0.0;
    }
    else
        errs[my_conf]=curr_diff[0];
    
    for(int i=0;i<nconfigs;i++)
        MPI_Bcast(&errs[i],1,Vec<type0>::MPI_T,roots[i],world);
    
    MPI_Allreduce(f_lcl,f,nvars,Vec<type0>::MPI_T,MPI_SUM,world);
    MPI_Allreduce(&err_lcl,&err,1,Vec<type0>::MPI_T,MPI_SUM,world);
    tot_err=err;
    return err;
}
/*--------------------------------------------
 just one test run
 --------------------------------------------*/
type0 PotFit::F(type0 alpha)
{
    for(int i=0;i<nvars;i++)
    {
        
        x[i]=x0[i]+alpha*h[i];
    }
    hyd_const->prep();
    hyd_const->update();
    
    if(ntrial%max_ntrials==0) full_reset();
    ntrial++;
    min->init();
    min->run(1000);
    min->fin();
    curr_diff[0]=atoms->pe;
    Algebra::MSY_2_DyadicV(atoms->S_pe,&curr_diff[1]);
    type0 err,err_lcl=val(curr_diff,target,coef);
    Algebra::V_sub<1+__nvoigt__>(target,curr_diff);
    if(my_lcl_rank!=0) err_lcl=0.0;
    MPI_Allreduce(&err_lcl,&err,1,Vec<type0>::MPI_T,MPI_SUM,world);
    tot_err=err;
    return err;
}
/*--------------------------------------------
 just one test run
 --------------------------------------------*/
void PotFit::F_reset()
{
    restore_x0();
    if(ntrial%max_ntrials==0) full_reset();
    ntrial++;
    min->init();
    min->run(1000);
    min->fin();
}
/*--------------------------------------------
 find maximum h
 lets find the next sensible number
 
 x=x_0+h*alpha
 (x-x0)/alpha=sqrt(eps)/alpha

 --------------------------------------------*/
void PotFit::ls_prep(type0& dfa,type0& h_norm,type0& max_a)
{
    
    h_norm=dot(h,h);
    
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
    
    
    //fix this shit
    max_a=std::numeric_limits<type0>::infinity();
    max_a=MIN(max_a,find_max_alpha_rho());
    max_a=MIN(max_a,find_max_alpha_phi());
    max_a=MIN(max_a,find_max_alpha_F());
}
/*--------------------------------------------
 just one test run
 --------------------------------------------*/
type0 PotFit::find_max_alpha_F()
{
    type0 __max_a=std::numeric_limits<type0>::infinity();
    if(F_A_dof[1][0])
    {
        
        if(F_A_h[1][0]<0.0)
        {
            type0 max_a_h=-F_A[1][0]/F_A_h[1][0];
            while(F_A[1][0]+max_a_h*F_A_h[1][0]<0.0)
                max_a_h=nextafter(max_a_h,0.0);
            __max_a=MIN(__max_a,max_a_h);
        }
        
        if(F_A_h[1][0]!=0.0)
        __max_a=MIN(__max_a,max_AF/fabs(F_A_h[1][0]));
        
    }
    
    if(F_A_dof[1][1])
    {
        if(F_A_h[1][1]<0.0)
        {
            type0 max_a_h=-F_A[1][1]/F_A_h[1][1];
            while(F_A[1][1]+max_a_h*F_A_h[1][1]<0.0)
                max_a_h=nextafter(max_a_h,0.0);
            
            __max_a=MIN(__max_a,max_a_h);
        }
        
        if(F_A_h[1][1]!=0.0)
            __max_a=MIN(__max_a,max_rhobF/fabs(F_A_h[1][1]));
        
        
        
    }
    
    
    if(F_A_dof[1][2])
    {
        if(F_A_h[1][2]>0.0)
            __max_a=MIN(__max_a,(0.5-F_A[1][2])/F_A_h[1][2]);
        if(F_A_h[1][2]<0.0)
            __max_a=MIN(__max_a,-F_A[1][2]/F_A_h[1][2]);
    
        if(F_A_h[1][2]!=0.0)
            __max_a=MIN(__max_a,max_alphaF/fabs(F_A_h[1][2]));
    }
    
    
    
    return __max_a;
}
/*--------------------------------------------
 just one test run
 --------------------------------------------*/
type0 PotFit::find_max_alpha_rho()
{
    size_t nelems=2;
    type0 __max_a=std::numeric_limits<type0>::infinity();
    bool rho_dof=false;
    type0 max_h=0.0;
    for(size_t ielem=0;ielem<nelems;ielem++)
        for(size_t jelem=0;jelem<nelems;jelem++)
        {
            spline2(rho_A[ielem][jelem],rho_A_h[ielem][jelem],rho_R[ielem][jelem],rho_sz[ielem][jelem],__max_a);
            for(size_t i=0;i<rho_sz[ielem][jelem];i++)
                if(rho_A_dof[ielem][jelem][i])
                {
                    rho_dof=true;
                    max_h=MAX(max_h,fabs(rho_A_h[ielem][jelem][i]));
                    
                }
        }
    
    if(rho_dof && max_h!=0.0)
        __max_a=MIN(__max_a,max_Arho/max_h);
    return __max_a;
}
/*--------------------------------------------
 just one test run
 --------------------------------------------*/
type0 PotFit::find_max_alpha_phi()
{
    size_t nelems=2;
    bool phi_dof=false;
    type0 max_h=0.0;
    for(size_t ielem=0;ielem<nelems;ielem++)
        for(size_t jelem=0;jelem<ielem+1;jelem++)
            for(size_t i=0;i<phi_sz[ielem][jelem];i++)
                if(phi_A_dof[ielem][jelem][i])
                    {
                        phi_dof=true;
                        max_h=MAX(max_h,fabs(phi_A_h[ielem][jelem][i]));
                    }
    if(phi_dof && max_h!=0.0)
         return max_Aphi/max_h;
    
    return std::numeric_limits<type0>::infinity();
}
/*--------------------------------------------
 just one test run
 --------------------------------------------*/
size_t PotFit::get_rFeH(type0*& Rs,int*& Ns)
{
    if(ntrial%max_ntrials==0) full_reset();
    ntrial++;
    min->init();
    min->run(1000);
    size_t sz=ff->get_rFeH(Rs,Ns);
    min->fin();
    
    return sz;
}
/*------------------------------------------------------------------------------------------------------------------------------------
 
 ------------------------------------------------------------------------------------------------------------------------------------*/
PyObject* PotFit::__new__(PyTypeObject* type,PyObject* args,PyObject* kwds)
{
    Object* __self=reinterpret_cast<Object*>(type->tp_alloc(type,0));
    PyObject* self=reinterpret_cast<PyObject*>(__self);
    return self;
}
/*--------------------------------------------
 
 --------------------------------------------*/
int PotFit::__init__(PyObject* self,PyObject* args,PyObject* kwds)
{
    FuncAPI<> func("__init__");
    
    FuncAPI<type0(***)[2],symm<type0(***)[2]>,type0(*)[3],
    size_t*,size_t*,type0*,type0**,type0*,type0**,type0**,type0**,type0***,
    OB<PyListObject,PyList_Type>>
    f("__init__",{
        "rho_AR", //0
        "phi_AR", //1
        "F_A", //2
        "phi_const_list", //3
        "phi_nconst_list", //4
        "phi0", //5
        "A", //6
        "rho0", //7
        "B", //8
        "L", //9
        "D0", //10
        "C", //11
        "conf" //12
    });
    if(f(args,kwds)) return -1;
    
    
    
    
    
    size_t nelems=2;
    size_t** rho_AR_sz=NULL;
    Memory::alloc(rho_AR_sz,nelems,nelems);
    size_t** phi_AR_sz=NULL;
    Memory::alloc(phi_AR_sz,nelems,nelems);
    for(size_t ielem=0;ielem<nelems;ielem++)
        for(size_t jelem=0;jelem<nelems;jelem++)
        {
            rho_AR_sz[ielem][jelem]=f.v<0>()[ielem][jelem].size;
            phi_AR_sz[ielem][jelem]=f.v<1>()[ielem][jelem].size;
        }
    
    
    
    
    
    
    
    PyObject* configs=f.val<12>();
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
    __self->potfit=new PotFit(f.mov<0>(),std::move(rho_AR_sz),
                              f.mov<1>(),std::move(phi_AR_sz),
                              f.mov<2>(),
                              f.mov<3>(),f.v<3>().size,
                              f.mov<4>(),f.v<4>().size,
                              f.mov<5>(),f.mov<6>(),f.mov<7>(),f.mov<8>(),f.mov<9>(),f.mov<10>(),f.mov<11>(),
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
PyObject* PotFit::__alloc__(PyTypeObject* type,Py_ssize_t)
{
    Object* __self=new Object;
    __self->ob_type=type;
    __self->ob_refcnt=1;
    __self->potfit=NULL;
    return reinterpret_cast<PyObject*>(__self);
}
/*--------------------------------------------
 
 --------------------------------------------*/
void PotFit::__dealloc__(PyObject* self)
{
    Object* __self=reinterpret_cast<Object*>(self);
    delete __self->potfit;
    __self->potfit=NULL;
    delete __self;
}
/*--------------------------------------------*/
PyTypeObject PotFit::TypeObject={PyObject_HEAD_INIT(NULL)};
/*--------------------------------------------*/
int PotFit::setup_tp()
{
    TypeObject.tp_name="mapp.potfit";
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
PyGetSetDef PotFit::getset[]=EmptyPyGetSetDef(15);
/*--------------------------------------------*/
void PotFit::setup_tp_getset()
{
    getset_rho_A(getset[0]);
    getset_rho_R(getset[1]);
    getset_rho_dof(getset[2]);
    getset_phi_A(getset[3]);
    getset_phi_R(getset[4]);
    getset_phi_dof(getset[5]);
    getset_F_A(getset[6]);
    getset_F_dof(getset[7]);
    getset_max_ntrials(getset[8]);
    getset_max_rhobF(getset[9]);
    getset_max_AF(getset[10]);
    getset_max_alphaF(getset[11]);
    getset_max_Arho(getset[12]);
    getset_max_Aphi(getset[13]);
    
}
/*--------------------------------------------*/
PyMethodDef PotFit::methods[]=EmptyPyMethodDef(3);
/*--------------------------------------------*/
void PotFit::setup_tp_methods()
{
    ml_min(methods[0]);
    ml_reset(methods[1]);
}
/*--------------------------------------------
 
 --------------------------------------------*/
void PotFit::getset_rho_A(PyGetSetDef& getset)
{
    getset.name=(char*)"rho_A";
    getset.doc=(char*)"";
    
    getset.get=[](PyObject* self,void*)->PyObject*
    {
        
        PotFit* potfit=reinterpret_cast<Object*>(self)->potfit;
        size_t nelems=2;
        size_t sz;
        size_t* szp=&sz;
        
        
        PyObject* py_obj=PyList_New(nelems);
        for(size_t i=0;i<nelems;i++)
        {
            PyObject* __py_obj=PyList_New(nelems);
            for(size_t j=0;j<nelems;j++)
            {
                sz=potfit->rho_sz[i][j];
                PyList_SET_ITEM(__py_obj,j,var<type0*>::build(potfit->rho_A[i][j],&szp));
            }
            PyList_SET_ITEM(py_obj,i,__py_obj);
        }

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
void PotFit::getset_rho_R(PyGetSetDef& getset)
{
    getset.name=(char*)"rho_R";
    getset.doc=(char*)"";
    
    getset.get=[](PyObject* self,void*)->PyObject*
    {
        
        PotFit* potfit=reinterpret_cast<Object*>(self)->potfit;
        size_t nelems=2;
        size_t sz;
        size_t* szp=&sz;
        
        
        PyObject* py_obj=PyList_New(nelems);
        for(size_t i=0;i<nelems;i++)
        {
            PyObject* __py_obj=PyList_New(nelems);
            for(size_t j=0;j<nelems;j++)
            {
                sz=potfit->rho_sz[i][j];
                PyList_SET_ITEM(__py_obj,j,var<type0*>::build(potfit->rho_R[i][j],&szp));
            }
            PyList_SET_ITEM(py_obj,i,__py_obj);
        }

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
void PotFit::getset_rho_dof(PyGetSetDef& getset)
{
    getset.name=(char*)"rho_dof";
    getset.doc=(char*)"";
    
    getset.get=[](PyObject* self,void*)->PyObject*
    {
        
        PotFit* potfit=reinterpret_cast<Object*>(self)->potfit;
        size_t nelems=2;
        size_t sz;
        size_t* szp=&sz;
        
        
        PyObject* py_obj=PyList_New(nelems);
        for(size_t i=0;i<nelems;i++)
        {
            PyObject* __py_obj=PyList_New(nelems);
            for(size_t j=0;j<nelems;j++)
            {
                sz=potfit->rho_sz[i][j];
                PyList_SET_ITEM(__py_obj,j,var<bool*>::build(potfit->rho_A_dof[i][j],&szp));
            }
            PyList_SET_ITEM(py_obj,i,__py_obj);
        }

        return py_obj;
    };
    getset.set=[](PyObject* self,PyObject* val,void*)->int
    {
        PotFit* potfit=reinterpret_cast<Object*>(self)->potfit;
        VarAPI<bool***> rho_dof("rho_dof");
        if(rho_dof.set(val)!=0) return -1;
        size_t nelems=2;
        for(size_t ielem=0;ielem<nelems;ielem++)
            for(size_t jelem=0;jelem<nelems;jelem++)
            {
                if(rho_dof.__var__[ielem][jelem].size!=potfit->rho_sz[ielem][jelem])
                {
                    PyErr_SetString(PyExc_TypeError,"size mismatch");
                    return -1;
                }
            }
        
        for(size_t ielem=0;ielem<nelems;ielem++)
            for(size_t jelem=0;jelem<nelems;jelem++)
                for(size_t i=0;i<potfit->rho_sz[ielem][jelem];i++)
                    potfit->rho_A_dof[ielem][jelem][i]=rho_dof.val[ielem][jelem][i];
        
        return 0;
    };
}
/*--------------------------------------------
 
 --------------------------------------------*/
void PotFit::getset_phi_A(PyGetSetDef& getset)
{
    getset.name=(char*)"phi_A";
    getset.doc=(char*)"";
    
    getset.get=[](PyObject* self,void*)->PyObject*
    {
        
        PotFit* potfit=reinterpret_cast<Object*>(self)->potfit;
        size_t nelems=2;
        size_t sz;
        size_t* szp=&sz;
        
        
        PyObject* py_obj=PyList_New(nelems);
        for(size_t i=0;i<nelems;i++)
        {
            PyObject* __py_obj=PyList_New(nelems);
            for(size_t j=0;j<nelems;j++)
            {
                sz=potfit->phi_sz[i][j];
                PyList_SET_ITEM(__py_obj,j,var<type0*>::build(potfit->phi_A[i][j],&szp));
            }
            PyList_SET_ITEM(py_obj,i,__py_obj);
        }

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
void PotFit::getset_phi_R(PyGetSetDef& getset)
{
    getset.name=(char*)"phi_R";
    getset.doc=(char*)"";
    
    getset.get=[](PyObject* self,void*)->PyObject*
    {
        
        PotFit* potfit=reinterpret_cast<Object*>(self)->potfit;
        size_t nelems=2;
        size_t sz;
        size_t* szp=&sz;
        
        
        PyObject* py_obj=PyList_New(nelems);
        for(size_t i=0;i<nelems;i++)
        {
            PyObject* __py_obj=PyList_New(nelems);
            for(size_t j=0;j<nelems;j++)
            {
                sz=potfit->phi_sz[i][j];
                PyList_SET_ITEM(__py_obj,j,var<type0*>::build(potfit->phi_R[i][j],&szp));
            }
            PyList_SET_ITEM(py_obj,i,__py_obj);
        }

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
void PotFit::getset_phi_dof(PyGetSetDef& getset)
{
    getset.name=(char*)"phi_dof";
    getset.doc=(char*)"";
    
    getset.get=[](PyObject* self,void*)->PyObject*
    {
        
        PotFit* potfit=reinterpret_cast<Object*>(self)->potfit;
        size_t nelems=2;
        size_t sz;
        size_t* szp=&sz;
        
        
        PyObject* py_obj=PyList_New(nelems);
        for(size_t i=0;i<nelems;i++)
        {
            PyObject* __py_obj=PyList_New(nelems);
            for(size_t j=0;j<nelems;j++)
            {
                sz=potfit->phi_sz[i][j];
                PyList_SET_ITEM(__py_obj,j,var<bool*>::build(potfit->phi_A_dof[i][j],&szp));
            }
            PyList_SET_ITEM(py_obj,i,__py_obj);
        }

        return py_obj;
    };
    getset.set=[](PyObject* self,PyObject* val,void*)->int
    {
        PotFit* potfit=reinterpret_cast<Object*>(self)->potfit;
        VarAPI<symm<bool***>> phi_dof("phi_dof");
        if(phi_dof.set(val)!=0) return -1;
        size_t nelems=2;
        for(size_t ielem=0;ielem<nelems;ielem++)
            for(size_t jelem=0;jelem<nelems;jelem++)
            {
                if(phi_dof.__var__[ielem][jelem].size!=potfit->phi_sz[ielem][jelem])
                {
                    PyErr_SetString(PyExc_TypeError,"size mismatch");
                    return -1;
                }
            }
        
        for(size_t ielem=0;ielem<nelems;ielem++)
            for(size_t jelem=0;jelem<nelems;jelem++)
                for(size_t i=0;i<potfit->phi_sz[ielem][jelem];i++)
                    potfit->phi_A_dof[ielem][jelem][i]=phi_dof.val[ielem][jelem][i];
        potfit->hyd_const->readj_dof();
        return 0;
    };
}
/*--------------------------------------------
 
 --------------------------------------------*/
void PotFit::getset_F_A(PyGetSetDef& getset)
{
    getset.name=(char*)"F_A";
    getset.doc=(char*)"";
    
    getset.get=[](PyObject* self,void*)->PyObject*
    {
        
        PotFit* potfit=reinterpret_cast<Object*>(self)->potfit;
        size_t i=2;
        size_t j=3;
        size_t*szp[2]={&i,&j};
        return var<type0**>::build(potfit->F_A,szp);
        
    };
    getset.set=[](PyObject*,PyObject*,void*)->int
    {
        PyErr_SetString(PyExc_TypeError,"readonly attribute");
        return -1;
    };
}
/*--------------------------------------------
 
 --------------------------------------------*/
void PotFit::getset_F_dof(PyGetSetDef& getset)
{
    getset.name=(char*)"F_dof";
    getset.doc=(char*)"";
    
    getset.get=[](PyObject* self,void*)->PyObject*
    {
        
        PotFit* potfit=reinterpret_cast<Object*>(self)->potfit;
        size_t nelems=2;
        size_t sz[2]={nelems,3};
        size_t* szp=sz;
        
        return var<bool**>::build(potfit->F_A_dof,&szp);
        
    };
    getset.set=[](PyObject* self,PyObject* val,void*)->int
    {
        PotFit* potfit=reinterpret_cast<Object*>(self)->potfit;

        VarAPI<bool**> F_dof("F_dof");
        if(F_dof.set(val)!=0) return -1;
                size_t nelems=2;
        for(size_t ielem=0;ielem<nelems;ielem++)
            if(F_dof.__var__[ielem].size!=3)
            {
                PyErr_SetString(PyExc_TypeError,"size mismatch");
                return -1;
            }

        for(size_t ielem=0;ielem<nelems;ielem++)
            for(size_t i=0;i<3;i++)
                potfit->F_A_dof[ielem][i]=F_dof.__var__[ielem][i];
        
        return 0;
    };
}
/*--------------------------------------------
 
 --------------------------------------------*/
/*
void PotFit::getset_RFeH(PyGetSetDef& getset)
{
    getset.name=(char*)"RFeH";
    getset.doc=(char*)"";
    
    getset.get=[](PyObject* self,void*)->PyObject*
    {
        
        PotFit* potfit=reinterpret_cast<Object*>(self)->potfit;
        type0* Rs=NULL;
        int* Ns=NULL;
        size_t sz=potfit->get_rFeH(Rs,Ns);
        size_t* szp=&sz;

        PyObject* py_obj=PyList_New(2);
        PyList_SET_ITEM(py_obj,0,var<type0*>::build(Rs,&szp));
        PyList_SET_ITEM(py_obj,1,var<int*>::build(Ns,&szp));
        
        Memory::dealloc(Rs);
        Memory::dealloc(Ns);
        
        return py_obj;
        
        
    };
    getset.set=[](PyObject*,PyObject*,void*)->int
    {
        PyErr_SetString(PyExc_TypeError,"readonly attribute");
        return -1;
    };
}
 */
/*--------------------------------------------
 
 --------------------------------------------*/
void PotFit::getset_max_ntrials(PyGetSetDef& getset)
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
void PotFit::getset_max_rhobF(PyGetSetDef& getset)
{
    getset.name=(char*)"max_rhobF";
    getset.doc=(char*)"";
    
    getset.get=[](PyObject* self,void*)->PyObject*
    {
        return var<type0>::build(reinterpret_cast<Object*>(self)->potfit->max_rhobF);
    };
    getset.set=[](PyObject* self,PyObject* val,void*)->int
    {
        VarAPI<type0> var("max_rhobF");
        if(var.set(val)!=0) return -1;
        reinterpret_cast<Object*>(self)->potfit->max_rhobF=var.val;
        return 0;
    };
}
/*--------------------------------------------
 
 --------------------------------------------*/
void PotFit::getset_max_AF(PyGetSetDef& getset)
{
    getset.name=(char*)"max_AF";
    getset.doc=(char*)"";
    
    getset.get=[](PyObject* self,void*)->PyObject*
    {
        return var<type0>::build(reinterpret_cast<Object*>(self)->potfit->max_AF);
    };
    getset.set=[](PyObject* self,PyObject* val,void*)->int
    {
        VarAPI<type0> var("max_AF");
        if(var.set(val)!=0) return -1;
        reinterpret_cast<Object*>(self)->potfit->max_AF=var.val;
        return 0;
    };
}
/*--------------------------------------------
 
 --------------------------------------------*/
void PotFit::getset_max_alphaF(PyGetSetDef& getset)
{
    getset.name=(char*)"max_alphaF";
    getset.doc=(char*)"";
    
    getset.get=[](PyObject* self,void*)->PyObject*
    {
        return var<type0>::build(reinterpret_cast<Object*>(self)->potfit->max_alphaF);
    };
    getset.set=[](PyObject* self,PyObject* val,void*)->int
    {
        VarAPI<type0> var("max_alphaF");
        if(var.set(val)!=0) return -1;
        reinterpret_cast<Object*>(self)->potfit->max_alphaF=var.val;
        return 0;
    };
}
/*--------------------------------------------
 
 --------------------------------------------*/
void PotFit::getset_max_Arho(PyGetSetDef& getset)
{
    getset.name=(char*)"max_Arho";
    getset.doc=(char*)"";
    
    getset.get=[](PyObject* self,void*)->PyObject*
    {
        return var<type0>::build(reinterpret_cast<Object*>(self)->potfit->max_Arho);
    };
    getset.set=[](PyObject* self,PyObject* val,void*)->int
    {
        VarAPI<type0> var("max_Arho");
        if(var.set(val)!=0) return -1;
        reinterpret_cast<Object*>(self)->potfit->max_Arho=var.val;
        return 0;
    };
}
/*--------------------------------------------
 
 --------------------------------------------*/
void PotFit::getset_max_Aphi(PyGetSetDef& getset)
{
    getset.name=(char*)"max_Aphi";
    getset.doc=(char*)"";
    
    getset.get=[](PyObject* self,void*)->PyObject*
    {
        return var<type0>::build(reinterpret_cast<Object*>(self)->potfit->max_Aphi);
    };
    getset.set=[](PyObject* self,PyObject* val,void*)->int
    {
        VarAPI<type0> var("max_Aphi");
        if(var.set(val)!=0) return -1;
        reinterpret_cast<Object*>(self)->potfit->max_Aphi=var.val;
        return 0;
    };
}
/*--------------------------------------------
 
 --------------------------------------------*/
void PotFit::ml_min(PyMethodDef& tp_methods)
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
/*--------------------------------------------
 
 --------------------------------------------*/
void PotFit::ml_reset(PyMethodDef& tp_methods)
{
    tp_methods.ml_flags=METH_VARARGS | METH_KEYWORDS;
    tp_methods.ml_name="reset";
    tp_methods.ml_meth=(PyCFunction)(PyCFunctionWithKeywords)(
    [](PyObject* self,PyObject* args,PyObject* kwds)->PyObject*
    {
        
        FuncAPI<>f("reset");
        if(f(args,kwds)) return NULL;
        
        reinterpret_cast<Object*>(self)->potfit->full_reset();
        
        Py_RETURN_NONE;
    });
    
    
    tp_methods.ml_doc="";
}



