#ifndef __MAPP__potfit__
#define __MAPP__potfit__
#include "ff_eam_potfit.h"
#include "atoms_md.h"
#include "min_cg.h"
#include "thermo_dynamics.h"
#include "import_cfg.h"
#include "random.h"
namespace MAPP_NS
{
    template<class FF,size_t NELEMS>
    class PotFit
    {
    private:
        type0* f;
        type0* h;
        type0* f0;
        size_t nvars;
        int ntally;
        type0 err_tol;
        int nmin_steps;
        static const char* err_msgs[];
        type0* get_en_coefs();
        type0* get_f_coefs();
        type0* get_S_coefs();
        bool* get_H_dofs();
        type0* get_mean_rho(elem_type);
        Random* random;
    protected:
    public:
        ForceFieldEAMPotFit<NELEMS>* ff;
        AtomsMD* atoms;
        MinCG* min;
        LineSearchBrent* min_ls;
        
        int nconfigs;
        type0* errs;
        int* roots;
        std::string* names_str;
        const char** names;
        
        MPI_Comm world;
        MPI_Comm* my_world;
        int my_rank,my_conf,my_lcl_rank;
        
        type0 f_coef;
        type0 coef[1+__nvoigt__];
        type0 target[1+__nvoigt__];
        type0 err;
        
        
        
        
        VecTens<type0,1> Xorig;
        
        PotFit(PyObject*,
        std::string*&&,std::string*&,int*&,
        type0(*&)[1+__nvoigt__],
        MinCG**,
        PyObject**,size_t,MPI_Comm&);
        ~PotFit();
        
        void init();
        void fin();
        void assign_min(MinCG*);
        void store_x0();
        void restore_x0();
        void full_reset();
        void swtch2dflt();
        type0 val();
        type0 deriv();
        type0 iter();
        
        void test_deriv(ptrdiff_t,type0,int);
        void min_cg(int);
        void min_struc();
        void new_conf();
        void min_sa(int,int,int,type0,type0);
        void ls_prep(type0&,type0&,type0&);
        type0 F(type0);
        void F_reset();
        
        
        
        ThermoDynamics get_thermo(){RET_POTFIT_THRMO};
        
        typedef struct
        {
            PyObject_HEAD
            class PotFit<FF,NELEMS>* potfit;
            MinCG::Object** mins;
        }Object;
        
        static PyTypeObject TypeObject;
        static int setup_tp();
        
        static PyObject* __new__(PyTypeObject*,PyObject*, PyObject*);
        static int __init__(PyObject*, PyObject*,PyObject*);
        static PyObject* __alloc__(PyTypeObject*,Py_ssize_t);
        static void __dealloc__(PyObject*);
        
        static PyGetSetDef getset[];
        static void setup_tp_getset();
        
        static void getset_A_phi(PyGetSetDef& getset);
        static void getset_dA_phi_max(PyGetSetDef& getset);
        static void getset_A_phi_dof(PyGetSetDef& getset);
        static void getset_dA_phi(PyGetSetDef& getset);
        
        static void getset_A_rho(PyGetSetDef& getset);
        static void getset_dA_rho_max(PyGetSetDef& getset);
        static void getset_A_rho_dof(PyGetSetDef& getset);
        static void getset_dA_rho(PyGetSetDef& getset);
        
        static void getset_A_F(PyGetSetDef& getset);
        static void getset_dA_F_max(PyGetSetDef& getset);
        static void getset_A_F_dof(PyGetSetDef& getset);
        static void getset_dA_F(PyGetSetDef& getset);
        
        static void getset_ntally(PyGetSetDef&);
        static void getset_nmin_steps(PyGetSetDef&);
        static void getset_tol(PyGetSetDef&);
        static void getset_en_coefs(PyGetSetDef&);
        static void getset_f_coefs(PyGetSetDef&);
        static void getset_S_coefs(PyGetSetDef&);
        static void getset_mins(PyGetSetDef&);
        static void getset_errs(PyGetSetDef&);
        static void getset_err(PyGetSetDef&);
        
        
        static PyMethodDef methods[];
        static void setup_tp_methods();
        static void ml_val(PyMethodDef&);
        static void ml_deriv(PyMethodDef&);
        static void ml_min_struc(PyMethodDef&);
        static void ml_reset(PyMethodDef&);
        static void ml_swtch2dflt(PyMethodDef&);
        static void ml_min_cg(PyMethodDef&);
        static void ml_min_sa(PyMethodDef&);
        static void ml_mean_rho(PyMethodDef&);
        static void ml_test_A_phi(PyMethodDef&);
        static void ml_test_A_rho(PyMethodDef&);
        static void ml_test_A_F(PyMethodDef&);
    };
    
    
}
/*--------------------------------------------

 --------------------------------------------*/
template<class FF,size_t NELEMS>
const char* PotFit<FF,NELEMS>::err_msgs[]=
{
    //[LS_S]=
    "",
    //[LS_F_DOWNHILL]=
    "line search failed: not downhill direction\n",
    //[LS_F_GRAD0]=
    "line search failed: gradient is zero\n",
    //[LS_MIN_ALPHA]=
    "line search failed: minimum alpha reached\n",
    //[MIN_S_TOLERANCE]=
    "minimization finished: energy tolerance reached\n",
    //[MIN_F_MAX_ITER]=
    "minimization finished: maximum iteration reached\n",
    //[B_S]=
    "",
    //[B_F_MAX_ALPHA]=
    "bracketing failed: maximum alpha reached\n",
    //[B_F_DOWNHILL]=
    "bracketing failed: not downhill direction\n"
    
};
/*--------------------------------------------
 constructor
 --------------------------------------------*/
template<class FF,size_t NELEMS>
PotFit<FF,NELEMS>::PotFit(PyObject* ff_args,
std::string*&& __names_str,std::string*& files,int*& nprocs,
type0(*& __targets)[1+__nvoigt__],
MinCG** __mins,
PyObject** dof_funcs,size_t __nconfigs,
MPI_Comm& __world):
atoms(NULL),
world(__world)
{
    //defaults
    ntally=1;
    nmin_steps=0;
    err_tol=1.0e-8;
    
    nconfigs=static_cast<int>(__nconfigs);
    names_str=__names_str;
    __names_str=NULL;
    Memory::alloc(errs,__nconfigs);
    Memory::alloc(roots,__nconfigs);
    Memory::alloc(names,__nconfigs);
    MPI_Comm_rank(world,&my_rank);
    
    //lets figure out which configuration I am
    // based on my rank in the external world
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
    
    // create my local world
    my_world=new MPI_Comm;
    MPI_Comm_split(world,my_conf,my_lcl_rank,my_world);
    
    // get my coefficient and target
    // for now we just focus on energy
    Algebra::Do<1+__nvoigt__>::func([this](int i){coef[i]=1.0;});
    f_coef=0.0;
    Algebra::zero<1+__nvoigt__>(coef);
    coef[0]=1.0;
    Algebra::V_eq<1+__nvoigt__>(__targets[my_conf],target);
    
    //now lets load my conf
    ImportCFGMD read(*my_world);
    atoms=read(files[my_conf].c_str());
    if(dof_funcs[my_conf]) atoms->DO(dof_funcs[my_conf]);
    ff=FF::get_new_ff(atoms,ff_args,NULL);
    
    bool H_dof[__dim__][__dim__];
    for(int i=0;i<__dim__;i++) for(int j=0;j<__dim__;j++) H_dof[i][j]=false;
    H_dof[0][0]=H_dof[1][1]=H_dof[2][2]=true;
    
    Xorig.~VecTens();
    new (&Xorig) VecTens<type0,1>(atoms,true,__dim__);
    
    min_ls=new LineSearchBrent();
    
    
    min=__mins[my_conf];
    min->atoms=atoms;
    min->ff=ff;
    
    
    memcpy(Xorig.vecs[0]->begin(),atoms->x->begin(),sizeof(type0)*__dim__*atoms->natms_lcl);
    memcpy(&Xorig.A[0][0],&atoms->H[0][0],sizeof(type0)*__dim__*__dim__);

    f=ff->fvs;
    f0=ff->f0vs;
    h=ff->hvs;
    nvars=ff->nvs;
    init();
}
/*--------------------------------------------
 destructor
 --------------------------------------------*/
template<class FF,size_t NELEMS>
PotFit<FF,NELEMS>::~PotFit()
{
    min->atoms=NULL;
    min->ff=NULL;
    fin();
    Xorig.~VecTens();
    
    delete min_ls;
    
    
    delete ff;
    delete atoms;
    MPI_Comm_free(my_world);
    delete my_world;
    
    
    Memory::dealloc(roots);
    Memory::dealloc(errs);
    *names=NULL;
    Memory::dealloc(names);
    Memory::dealloc(names_str);
}
/*--------------------------------------------
 init
 --------------------------------------------*/
template<class FF,size_t NELEMS>
void PotFit<FF,NELEMS>::init()
{
    min->init(Xorig.vecs[0]);
}
/*--------------------------------------------
 init
 --------------------------------------------*/
template<class FF,size_t NELEMS>
void PotFit<FF,NELEMS>::fin()
{
    min->fin();
}
/*--------------------------------------------
 init
 --------------------------------------------*/
template<class FF,size_t NELEMS>
void PotFit<FF,NELEMS>::assign_min(MinCG* __min)
{
    if(min)
    {
        fin();
        min->atoms=NULL;
        min->ff=NULL;
    }
    min=__min;
    min->atoms=atoms;
    min->ff=ff;
    init();
}
/*--------------------------------------------
 new configuration
 --------------------------------------------*/
template<class FF,size_t NELEMS>
void PotFit<FF,NELEMS>::new_conf()
{
    bool valid=false;
    while(!valid)
    {
        for(size_t i=0;i<ff->nuniq_phi_ptrs;i++)
            ff->uniq_phi_ptr[i]->random_neigh(random);
        for(size_t i=0;i<ff->nuniq_F_ptrs;i++)
            ff->uniq_F_ptr[i]->random_neigh(random);
        for(size_t i=0;i<ff->nuniq_rho_ptrs;i++)
            ff->uniq_rho_ptr[i]->random_neigh(random);
        
        for(size_t i=0;i<nvars;i++)
            if(!ff->dofs[i]) h[i]=0.0;
        
        for(size_t i=0;i<nvars;i++)
            ff->vs[i]=ff->v0s[i]+h[i];
        
        valid=true;
        for(size_t i=0;i<ff->nuniq_phi_ptrs && valid;i++)
            valid=ff->uniq_phi_ptr[i]->validate();
        for(size_t i=0;i<ff->nuniq_F_ptrs && valid;i++)
            valid=ff->uniq_F_ptr[i]->validate();
        for(size_t i=0;i<ff->nuniq_rho_ptrs && valid;i++)
            valid=ff->uniq_rho_ptr[i]->validate();
    }
    
}
/*--------------------------------------------
 minimize error function simulated anealling
 --------------------------------------------*/
template<class FF,size_t NELEMS>
void PotFit<FF,NELEMS>::min_sa(int seed,int nsteps,int ntrials,type0 temp,type0 gamma)
{

    
    random=new Random (seed);
    
    //type0 old_skin=atoms->comm.skin;
    //atoms->comm.skin=0.0001;
    type0 * lwst_vars=NULL;
    Memory::alloc(lwst_vars,nvars);
    type0 lwst_val;
    
    
    
    type0 dc,c_n,c_o,dpe;
    lwst_val=err=val();
    memcpy(lwst_vars,ff->vs,nvars*sizeof(type0));
    c_o=err;
    dpe=atoms->pe-target[0];
    
    ThermoDynamics thermo=get_thermo();
    thermo.init();
    if(ntally)
    {
        errs[my_conf]=dpe;
        for(int i=0;i<nconfigs;i++)
            MPI_Bcast(&errs[i],1,Vec<type0>::MPI_T,roots[i],world);
        thermo.print(0);
    }
    for(int istep=0;istep<nsteps;istep++)
    {
        
        for(int j=0;j<ntrials;j++)
        {
            memcpy(ff->v0s,ff->vs,nvars*sizeof(type0));
            // do the randomization here
            new_conf();
            c_n=val();
            dc=c_n-c_o;
            //printf("%lf\n",dc);
            if(random->uniform()>=exp(-dc/temp))
                memcpy(ff->vs,ff->v0s,nvars*sizeof(type0));
            else
            {
                c_o=c_n;
                err=c_o;
                dpe=atoms->pe-target[0];
                if(c_n<lwst_val)
                {
                    lwst_val=c_n;
                    memcpy(lwst_vars,ff->vs,nvars*sizeof(type0));
                }
            }
        }
        
        
        if(ntally && (istep+1)%ntally==0)
        {
            errs[my_conf]=dpe;
            for(int i=0;i<nconfigs;i++)
                MPI_Bcast(&errs[i],1,Vec<type0>::MPI_T,roots[i],world);
            thermo.print(istep+1);
        }
        temp*=gamma;
    }
    
    if(ntally && nsteps%ntally)
    {
        errs[my_conf]=atoms->pe-target[0];
        for(int i=0;i<nconfigs;i++)
            MPI_Bcast(&errs[i],1,Vec<type0>::MPI_T,roots[i],world);
        thermo.print(nsteps);
    }
    thermo.fin();
    
    delete random;
    random=NULL;
    //atoms->comm.skin=old_skin;
    err=lwst_val;
    memcpy(ff->vs,lwst_vars,nvars*sizeof(type0));
    Memory::dealloc(lwst_vars);
}
/*--------------------------------------------
 minimize error function
 --------------------------------------------*/
template<class FF,size_t NELEMS>
void PotFit<FF,NELEMS>::min_struc()
{
    int __ntally=min->ntally;
    min->ntally=0;
    min->run(nmin_steps);
    min->ntally=__ntally;
    err=val();
    errs[my_conf]=atoms->pe-target[0];
    for(int i=0;i<nconfigs;i++)
        MPI_Bcast(&errs[i],1,Vec<type0>::MPI_T,roots[i],world);
    
    if(ntally)
    {
        ThermoDynamics thermo=get_thermo();
        thermo.init();
        thermo.print(0);
        thermo.fin();
    }
    
}
/*--------------------------------------------
 minimize error function
 --------------------------------------------*/
template<class FF,size_t NELEMS>
void PotFit<FF,NELEMS>::min_cg(int nsteps)
{
    ThermoDynamics thermo=get_thermo();
    
    type0 curr_err,prev_err,f0_f0,f_f0,f_f,f_h,ratio,alpha,tol=1.0e-12;
    err=iter();
    curr_err=err;
    memcpy(h,f,nvars*sizeof(type0));
    f0_f0=0.0;
    for(size_t i=0;i<nvars;i++) f0_f0+=f[i]*f[i];
    
    int istep=0;
    
    if(ntally)
    {
        thermo.init();
        errs[my_conf]=atoms->pe-target[0];
        for(int i=0;i<nconfigs;i++)
            MPI_Bcast(&errs[i],1,Vec<type0>::MPI_T,roots[i],world);
        thermo.print(0);
    }
    int ERR=nsteps==0? MIN_F_MAX_ITER:LS_S;
    for(;ERR==LS_S;istep++)
    {
        if(f0_f0==0.0)
        {
            ERR=LS_F_GRAD0;
            continue;
        }
        
        store_x0();
        
        memcpy(f0,f,nvars*sizeof(type0));
        prev_err=curr_err;
        f_h=0.0;
        for(size_t i=0;i<nvars;i++) f_h+=f[i]*h[i];
        if(f_h<0.0)
        {
            memcpy(h,f,nvars*sizeof(type0));
            f_h=f0_f0;
        }
        
        
        //perform linesearch here
        //printf("%.16lf\n",curr_err);
        ERR=min_ls->min(this,curr_err,alpha,1);
        //check for error from linesearch here
        if(ERR!=LS_S)
        {
            // this was a bullshit step so we have to decrease the setup once to compensate
            // the last step was the previous one
            istep--;
            continue;
        }
        err=iter();
        curr_err=err;
        if(ntally && (istep+1)%ntally==0)
        {
            errs[my_conf]=atoms->pe-target[0];
            for(int i=0;i<nconfigs;i++)
                MPI_Bcast(&errs[i],1,Vec<type0>::MPI_T,roots[i],world);
            thermo.print(istep+1);
        }
        
        if(prev_err-curr_err<tol || curr_err<err_tol) ERR=MIN_S_TOLERANCE;
        if(istep+1==nsteps) ERR=MIN_F_MAX_ITER;
        if(ERR) continue;
        
        f_f=0.0;
        for(size_t i=0;i<nvars;i++) f_f+=f[i]*f[i];
        f_f0=0.0;
        for(size_t i=0;i<nvars;i++) f_f0+=f[i]*f0[i];
        ratio=(f_f-f_f0)/(f0_f0);
        for(size_t j=0;j<nvars;j++) h[j]=ratio*h[j]+f[j];
        //advance h here
        f0_f0=f_f;
    }
    if(ntally && istep%ntally)
    {
        errs[my_conf]=atoms->pe-target[0];
        for(int i=0;i<nconfigs;i++)
            MPI_Bcast(&errs[i],1,Vec<type0>::MPI_T,roots[i],world);
        thermo.print(istep);
        
        
    }
    if(ntally) thermo.fin();
    
}
/*--------------------------------------------
 find maximum h
 lets find the next sensible number
 
 x=x_0+h*alpha
 (x-x0)/alpha=sqrt(eps)/alpha
  --------------------------------------------*/
template<class FF,size_t NELEMS>
void PotFit<FF,NELEMS>::ls_prep(type0& dfa,type0& h_norm,type0& max_a)
{
    
    //fix this shit
    max_a=std::numeric_limits<type0>::infinity();
    for(size_t i=0;i<ff->nuniq_phi_ptrs;i++)
    ff->uniq_phi_ptr[i]->find_max_alpha(max_a);
    for(size_t i=0;i<ff->nuniq_F_ptrs;i++)
    ff->uniq_F_ptr[i]->find_max_alpha(max_a);
    for(size_t i=0;i<ff->nuniq_rho_ptrs;i++)
    ff->uniq_rho_ptr[i]->find_max_alpha(max_a);
    
    h_norm=0.0;
    for(size_t i=0;i<nvars;i++) h_norm+=h[i]*h[i];
    
    dfa=0.0;
    for(size_t i=0;i<nvars;i++) dfa-=f[i]*h[i];
    
    
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
    
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<class FF,size_t NELEMS>
type0 PotFit<FF,NELEMS>::F(type0 alpha)
{
    for(size_t i=0;i<nvars;i++)
        ff->vs[i]=ff->v0s[i]+alpha*h[i];
    err=val();
    return err;
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<class FF,size_t NELEMS>
void PotFit<FF,NELEMS>::test_deriv(ptrdiff_t offset,type0 delta,int nsteps)
{
    int __nmin_steps=nmin_steps;
    nmin_steps=0;
    
    full_reset();
    store_x0();
    err=iter();
    type0 err0=err;
    type0 deriv=-f[offset];
    for(size_t i=0;i<nvars;i++)
        h[i]=0.0;
    h[offset]=1.0;
    type0 DeltaX=0.0,DeltaERR,DERRDX;
    ThermoDynamics thermo(6,"DeltaX",DeltaX,"DeltaERR",DeltaERR,"DERR*DX",DERRDX);
    thermo.init();
    for(int i=0;i<nsteps;i++)
    {
        DeltaERR=F(DeltaX)-err0;
        DERRDX=DeltaX*deriv;
        thermo.print(i);
        DeltaX+=delta;
    }
    thermo.fin();
    restore_x0();
    
    nmin_steps=__nmin_steps;
    
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<class FF,size_t NELEMS>
type0 PotFit<FF,NELEMS>::val()
{
    type0 err_lcl,err_lcl_lcl;
    type0* curr_val=ff->derivative();
    type0* f_vec=ff->f->begin();
    err_lcl_lcl=0.0;
    for(int i=0;i<atoms->natms_lcl*__dim__;i++) err_lcl_lcl+=f_vec[i]*f_vec[i];
    MPI_Allreduce(&err_lcl_lcl,&err_lcl,1,Vec<type0>::MPI_T,MPI_SUM,*my_world);
    err_lcl*=f_coef/static_cast<type0>(atoms->natms*__dim__);
    Algebra::Do<1+__nvoigt__>::func([this,&err_lcl,&curr_val](int i)
    {err_lcl+=coef[i]*(curr_val[i]-target[i])*(curr_val[i]-target[i]);});
    if(my_lcl_rank!=0) err_lcl=0.0;
    type0 err=0.0;
    MPI_Allreduce(&err_lcl,&err,1,Vec<type0>::MPI_T,MPI_SUM,world);
    return 0.5*err;
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<class FF,size_t NELEMS>
type0 PotFit<FF,NELEMS>::deriv()
{
    type0 err=val();
    ff->force_gradient(target,coef,f_coef);
    MPI_Allreduce(ff->dvs_lcl,ff->dvs,
    static_cast<int>(nvars),Vec<type0>::MPI_T,MPI_SUM,world);
    return err;
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<class FF,size_t NELEMS>
type0 PotFit<FF,NELEMS>::iter()
{
    type0 err=deriv();
    for(size_t i=0;i<nvars;i++)
        f[i]= (ff->dofs[i]) ? -ff->dvs[i]:0.0;
    return err;
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<class FF,size_t NELEMS>
void PotFit<FF,NELEMS>::F_reset()
{
    memcpy(ff->vs,ff->v0s,nvars*sizeof(type0));
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<class FF,size_t NELEMS>
void PotFit<FF,NELEMS>::store_x0()
{
    memcpy(ff->v0s,ff->vs,nvars*sizeof(type0));
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<class FF,size_t NELEMS>
void PotFit<FF,NELEMS>::restore_x0()
{
    memcpy(ff->vs,ff->v0s,nvars*sizeof(type0));
}
/*--------------------------------------------
 just one test run
 --------------------------------------------*/
template<class FF,size_t NELEMS>
void PotFit<FF,NELEMS>::full_reset()
{
    memcpy(atoms->x->begin(),Xorig.vecs[0]->begin(),sizeof(type0)*__dim__*atoms->natms_lcl);
    memcpy(&atoms->H[0][0],&Xorig.A[0][0],sizeof(type0)*__dim__*__dim__);
    atoms->update_H();
    min->dynamic->update<true>();
}
/*--------------------------------------------
 just one test run
 --------------------------------------------*/
template<class FF,size_t NELEMS>
void PotFit<FF,NELEMS>::swtch2dflt()
{
    memcpy(Xorig.vecs[0]->begin(),atoms->x->begin(),sizeof(type0)*__dim__*atoms->natms_lcl);
    memcpy(&Xorig.A[0][0],&atoms->H[0][0],sizeof(type0)*__dim__*__dim__);
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<class FF,size_t NELEMS>
type0* PotFit<FF,NELEMS>::get_en_coefs()
{
    type0* __en_coefs=NULL;
    
    Memory::alloc(__en_coefs,nconfigs);
    __en_coefs[my_conf]=coef[0];
    for(int i=0;i<nconfigs;i++)
        MPI_Bcast(&__en_coefs[i],1,Vec<type0>::MPI_T,roots[i],world);
    
    return __en_coefs;
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<class FF,size_t NELEMS>
type0* PotFit<FF,NELEMS>::get_f_coefs()
{
    type0* __f_coefs=NULL;
    
    Memory::alloc(__f_coefs,nconfigs);
    __f_coefs[my_conf]=f_coef;
    for(int i=0;i<nconfigs;i++)
        MPI_Bcast(&__f_coefs[i],1,Vec<type0>::MPI_T,roots[i],world);
    
    return __f_coefs;
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<class FF,size_t NELEMS>
type0* PotFit<FF,NELEMS>::get_S_coefs()
{
    type0(* __S_coefs)[__dim__][__dim__]=NULL;
    
    Memory::alloc(__S_coefs,nconfigs);

    Algebra::DyadicV_2_MSY(coef+1,__S_coefs[my_conf]);
    for(int i=0;i<nconfigs;i++)
        MPI_Bcast(&__S_coefs[i][0][0],__dim__*__dim__,Vec<type0>::MPI_T,roots[i],world);
    
    return &__S_coefs[0][0][0];
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<class FF,size_t NELEMS>
bool* PotFit<FF,NELEMS>::get_H_dofs()
{
    bool(* __H_dofs)[__dim__][__dim__]=NULL;
    Memory::alloc(__H_dofs,nconfigs);
    memcpy(&__H_dofs[my_conf][0][0],&(min->H_dof[0][0]),__dim__*__dim__*sizeof(bool));
    //min->get_H_dof(__H_dofs[my_conf]);
    for(int i=0;i<nconfigs;i++)
        MPI_Bcast(&__H_dofs[i][0][0],__dim__*__dim__,Vec<bool>::MPI_T,roots[i],world);
    
    return &__H_dofs[0][0][0];
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<class FF,size_t NELEMS>
type0* PotFit<FF,NELEMS>::get_mean_rho(elem_type ielem)
{
    type0* __mean_rhos=NULL;
    
    Memory::alloc(__mean_rhos,nconfigs);
    __mean_rhos[my_conf]=ff->mean_rho(ielem);
    for(int i=0;i<nconfigs;i++)
        MPI_Bcast(&__mean_rhos[i],1,Vec<type0>::MPI_T,roots[i],world);
    
    return __mean_rhos;
}
/*------------------------------------------------------------------------------------------------------------------------------------
 
 ------------------------------------------------------------------------------------------------------------------------------------*/
template<class FF,size_t NELEMS>
PyObject* PotFit<FF,NELEMS>::__new__(PyTypeObject* type,PyObject* args,PyObject* kwds)
{
    Object* __self=reinterpret_cast<Object*>(type->tp_alloc(type,0));
    PyObject* self=reinterpret_cast<PyObject*>(__self);
    return self;
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<class FF,size_t NELEMS>
int PotFit<FF,NELEMS>::__init__(PyObject* self,PyObject* args,PyObject* kwds)
{
    
    
    FuncAPI<OB<PyListObject,PyList_Type>,OB<PyListObject,PyList_Type>>
    f("__init__",{
        "ff_param",    //0
        "conf"        //1
    });
    if(f(args,kwds)) return -1;
    
    PyObject* configs=f.val<1>();
    size_t nconfigs=PyList_Size(configs);
    std::string* names;
    Memory::alloc(names,nconfigs);
    std::string* files;
    Memory::alloc(files,nconfigs);
    int* nprocs;
    Memory::alloc(nprocs,nconfigs);
    type0(*  targets)[1+__nvoigt__];
    Memory::alloc(targets,nconfigs);
    PyObject** funcs;
    Memory::alloc(funcs,nconfigs);
    MinCG::Object** min_objs;
    Memory::alloc(min_objs,nconfigs);
    
    
    int tot_nproc=0;
    for(size_t i=0;i<nconfigs;i++)
    {
        PyObject* config=PyList_GetItem(configs,i);
        FuncAPI<std::string,std::string,
        int,
        type0,symm<type0[__dim__][__dim__]>,
        OP<MinCG>,
        OB<PyFunctionObject,PyFunction_Type>>
        conf("configuration tuple",{"name","cfg_file","nproc","pe_target","S_target","min_obj","dof_func"});
        
        conf.val<6>()=NULL;
        conf.noptionals=1;
        if(conf(config,NULL))
        {
            Memory::dealloc(names);
            Memory::dealloc(files);
            Memory::dealloc(targets);
            *min_objs=NULL;
            Memory::dealloc(min_objs);
            *funcs=NULL;
            Memory::dealloc(funcs);
            return -1;
        }
        names[i]=conf.val<0>();
        files[i]=conf.val<1>();
        nprocs[i]=conf.val<2>();
        targets[i][0]=conf.val<3>();
        Algebra::MSY_2_DyadicV(conf.val<4>(),&targets[i][1]);
        min_objs[i]=reinterpret_cast<MinCG::Object*>(conf.val<5>().ob);
        funcs[i]=conf.val<6>();
        tot_nproc+=nprocs[i];
    }
    
    MPI_Comm world=MPI_COMM_WORLD;
    int comm_size;
    MPI_Comm_size(world,&comm_size);
    if(tot_nproc!=comm_size)
    {
        PyErr_Format(PyExc_TypeError,"number of processors do not match: available: %d vs requested: %d",comm_size,tot_nproc);
        Memory::dealloc(names);
        Memory::dealloc(files);
        Memory::dealloc(targets);
        *min_objs=NULL;
        Memory::dealloc(min_objs);
        *funcs=NULL;
        Memory::dealloc(funcs);
        return -1;
    }
    
    
    
    MinCG** mins;
    Memory::alloc(mins,nconfigs);
    for(size_t i=0;i<nconfigs;i++) mins[i]=min_objs[i]->min;
    Object* __self=reinterpret_cast<Object*>(self);
    __self->potfit=new PotFit<FF,NELEMS>(f.val<0>(),std::move(names),
    files,nprocs,targets,mins,funcs,nconfigs,world);
    *mins=NULL;
    Memory::dealloc(mins);
    
    __self->mins=min_objs;
    for(size_t i=0;i<nconfigs;i++)
        Py_INCREF(min_objs[i]);
    
    min_objs=NULL;
    
    
    Memory::dealloc(names);
    Memory::dealloc(files);
    Memory::dealloc(targets);
    *funcs=NULL;
    Memory::dealloc(funcs);
    return 0;
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<class FF,size_t NELEMS>
PyObject* PotFit<FF,NELEMS>::__alloc__(PyTypeObject* type,Py_ssize_t)
{
    Object* __self=new Object;
    Py_TYPE(__self)=type;
    Py_REFCNT(__self)=1;
    __self->potfit=NULL;
    __self->mins=NULL;
    return reinterpret_cast<PyObject*>(__self);
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<class FF,size_t NELEMS>
void PotFit<FF,NELEMS>::__dealloc__(PyObject* self)
{
    Object* __self=reinterpret_cast<Object*>(self);
    if(__self->potfit)
    {
        int __nconfigs=__self->potfit->nconfigs;
        for(int i=0;i<__nconfigs;i++)
        Py_DECREF(__self->mins[i]);
    }
    delete __self->potfit;
    __self->potfit=NULL;
    *(__self->mins)=NULL;
    Memory::dealloc(__self->mins);
    __self->mins=NULL;
    delete __self;
}
/*--------------------------------------------*/
template<class FF,size_t NELEMS>
PyTypeObject PotFit<FF,NELEMS>::TypeObject={PyObject_HEAD_INIT(NULL)};
/*--------------------------------------------*/
template<class FF,size_t NELEMS>
int PotFit<FF,NELEMS>::setup_tp()
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
template<class FF,size_t NELEMS>
PyGetSetDef PotFit<FF,NELEMS>::getset[]=EmptyPyGetSetDef(22);
/*--------------------------------------------*/
template<class FF,size_t NELEMS>
void PotFit<FF,NELEMS>::setup_tp_getset()
{
    getset_A_phi(getset[0]);
    getset_dA_phi_max(getset[1]);
    getset_A_phi_dof(getset[2]);
    getset_dA_phi(getset[3]);
    
    getset_A_rho(getset[4]);
    getset_dA_rho_max(getset[5]);
    getset_A_rho_dof(getset[6]);
    getset_dA_rho(getset[7]);
    
    getset_A_F(getset[8]);
    getset_dA_F_max(getset[9]);
    getset_A_F_dof(getset[10]);
    getset_dA_F(getset[11]);
    
    
    getset_ntally(getset[12]);
    getset_nmin_steps(getset[13]);
    getset_tol(getset[14]);
    getset_en_coefs(getset[15]);
    getset_f_coefs(getset[16]);
    getset_S_coefs(getset[17]);
    getset_mins(getset[18]);
    getset_errs(getset[19]);
    getset_err(getset[20]);
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<class FF,size_t NELEMS>
void PotFit<FF,NELEMS>::getset_A_phi(PyGetSetDef& getset)
{
    getset.name=(char*)"A_phi";
    getset.doc=(char*)"";
    getset.get=[](PyObject* self,void*)->PyObject*
    {
        PotFit<FF,NELEMS>* potfit=reinterpret_cast<Object*>(self)->potfit;
        PyObject* py_obj=PyList_New(NELEMS);
        for(size_t i=0;i<NELEMS;i++)
            PyList_SET_ITEM(py_obj,i,PyList_New(NELEMS));
        for(size_t i=0;i<NELEMS;i++)
            for(size_t j=0;j<NELEMS;j++)
                PyList_SET_ITEM(PyList_GET_ITEM(py_obj,i),j,potfit->ff->phi_ptr[i][j]->get_A());
        
        return py_obj;
    };
    getset.set=[](PyObject* self,PyObject* val,void*)->int
    {
        PotFit<FF,NELEMS>* potfit=reinterpret_cast<Object*>(self)->potfit;
        VarAPI<OB<PyListObject,PyList_Type>> var("A_phi");
        if(var.set(val)!=0) return -1;
        PyObject* py_ob=var.val;
        size_t nelems=PyList_Size(py_ob);
        if(nelems!=NELEMS)
        {
            PyErr_Format(PyExc_TypeError,"size mismatch in A_phi");
            return -1;
        }
        
        for(size_t i=0;i<nelems;i++)
            if(PyList_Size(PyList_GET_ITEM(py_ob,i))!=NELEMS)
            {
                PyErr_Format(PyExc_TypeError,"size mismatch in A_phi");
                return -1;
            }
        
        for(size_t i=0;i<nelems;i++)
            for(size_t j=0;j<nelems;j++)
                if(potfit->ff->phi_ptr[i][j]->set_A(PyList_GET_ITEM(PyList_GET_ITEM(py_ob,i),j))!=0)
                   return -1;
            
        
        return 0;
    };
    
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<class FF,size_t NELEMS>
void PotFit<FF,NELEMS>::getset_dA_phi(PyGetSetDef& getset)
{
    getset.name=(char*)"dA_phi";
    getset.doc=(char*)"";
    getset.get=[](PyObject* self,void*)->PyObject*
    {
        PotFit<FF,NELEMS>* potfit=reinterpret_cast<Object*>(self)->potfit;
        PyObject* py_obj=PyList_New(NELEMS);
        for(size_t i=0;i<NELEMS;i++)
            PyList_SET_ITEM(py_obj,i,PyList_New(NELEMS));
        for(size_t i=0;i<NELEMS;i++)
            for(size_t j=0;j<NELEMS;j++)
                PyList_SET_ITEM(PyList_GET_ITEM(py_obj,i),j,potfit->ff->phi_ptr[i][j]->get_dA());
        
        return py_obj;
    };
    getset.set=[](PyObject* self,PyObject* val,void*)->int
    {
        PyErr_SetString(PyExc_TypeError,"readonly attribute");
        return -1;
    };
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<class FF,size_t NELEMS>
void PotFit<FF,NELEMS>::getset_dA_phi_max(PyGetSetDef& getset)
{
    getset.name=(char*)"dA_phi_max";
    getset.doc=(char*)"";
    getset.get=[](PyObject* self,void*)->PyObject*
    {
        PotFit<FF,NELEMS>* potfit=reinterpret_cast<Object*>(self)->potfit;
        PyObject* py_obj=PyList_New(NELEMS);
        for(size_t i=0;i<NELEMS;i++)
            PyList_SET_ITEM(py_obj,i,PyList_New(NELEMS));
        for(size_t i=0;i<NELEMS;i++)
            for(size_t j=0;j<NELEMS;j++)
                PyList_SET_ITEM(PyList_GET_ITEM(py_obj,i),j,potfit->ff->phi_ptr[i][j]->get_dA_max());
        
        return py_obj;
    };
    getset.set=[](PyObject* self,PyObject* val,void*)->int
    {
        PotFit<FF,NELEMS>* potfit=reinterpret_cast<Object*>(self)->potfit;
        VarAPI<OB<PyListObject,PyList_Type>> var("dA_phi_max");
        if(var.set(val)!=0) return -1;
        PyObject* py_ob=var.val;
        size_t nelems=PyList_Size(py_ob);
        if(nelems!=NELEMS)
        {
            PyErr_Format(PyExc_TypeError,"size mismatch in dA_phi_max");
            return -1;
        }
        
        for(size_t i=0;i<nelems;i++)
            if(PyList_Size(PyList_GET_ITEM(py_ob,i))!=NELEMS)
            {
                PyErr_Format(PyExc_TypeError,"size mismatch in dA_phi_max");
                return -1;
            }
        
        for(size_t i=0;i<nelems;i++)
            for(size_t j=0;j<nelems;j++)
                if(potfit->ff->phi_ptr[i][j]->set_dA_max(PyList_GET_ITEM(PyList_GET_ITEM(py_ob,i),j))!=0)
                   return -1;
        
        
        return 0;
    };
    
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<class FF,size_t NELEMS>
void PotFit<FF,NELEMS>::getset_A_phi_dof(PyGetSetDef& getset)
{
    getset.name=(char*)"A_phi_dof";
    getset.doc=(char*)"";
    getset.get=[](PyObject* self,void*)->PyObject*
    {
        PotFit<FF,NELEMS>* potfit=reinterpret_cast<Object*>(self)->potfit;
        PyObject* py_obj=PyList_New(NELEMS);
        for(size_t i=0;i<NELEMS;i++)
            PyList_SET_ITEM(py_obj,i,PyList_New(NELEMS));
        for(size_t i=0;i<NELEMS;i++)
            for(size_t j=0;j<NELEMS;j++)
                PyList_SET_ITEM(PyList_GET_ITEM(py_obj,i),j,potfit->ff->phi_ptr[i][j]->get_A_dof());
        
        return py_obj;
    };
    getset.set=[](PyObject* self,PyObject* val,void*)->int
    {
        PotFit<FF,NELEMS>* potfit=reinterpret_cast<Object*>(self)->potfit;
        VarAPI<OB<PyListObject,PyList_Type>> var("A_phi_dof");
        if(var.set(val)!=0) return -1;
        PyObject* py_ob=var.val;
        size_t nelems=PyList_Size(py_ob);
        if(nelems!=NELEMS)
        {
            PyErr_Format(PyExc_TypeError,"size mismatch in A_phi_dof");
            return -1;
        }
        
        for(size_t i=0;i<nelems;i++)
            if(PyList_Size(PyList_GET_ITEM(py_ob,i))!=NELEMS)
            {
                PyErr_Format(PyExc_TypeError,"size mismatch in A_phi_dof");
                return -1;
            }
        
        for(size_t i=0;i<nelems;i++)
            for(size_t j=0;j<nelems;j++)
                if(potfit->ff->phi_ptr[i][j]->set_A_dof(PyList_GET_ITEM(PyList_GET_ITEM(py_ob,i),j))!=0)
                   return -1;
        
        
        return 0;
    };
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<class FF,size_t NELEMS>
void PotFit<FF,NELEMS>::getset_A_rho(PyGetSetDef& getset)
{
    getset.name=(char*)"A_rho";
    getset.doc=(char*)"";
    getset.get=[](PyObject* self,void*)->PyObject*
    {
        PotFit<FF,NELEMS>* potfit=reinterpret_cast<Object*>(self)->potfit;
        PyObject* py_obj=PyList_New(NELEMS);
        for(size_t i=0;i<NELEMS;i++)
            PyList_SET_ITEM(py_obj,i,PyList_New(NELEMS));
        for(size_t i=0;i<NELEMS;i++)
            for(size_t j=0;j<NELEMS;j++)
                PyList_SET_ITEM(PyList_GET_ITEM(py_obj,i),j,potfit->ff->rho_ptr[i][j]->get_A());
        
        return py_obj;
    };
    getset.set=[](PyObject* self,PyObject* val,void*)->int
    {
        PotFit<FF,NELEMS>* potfit=reinterpret_cast<Object*>(self)->potfit;
        VarAPI<OB<PyListObject,PyList_Type>> var("A_rho");
        if(var.set(val)!=0) return -1;
        PyObject* py_ob=var.val;
        size_t nelems=PyList_Size(py_ob);
        if(nelems!=NELEMS)
        {
            PyErr_Format(PyExc_TypeError,"size mismatch in A_rho");
            return -1;
        }
        
        for(size_t i=0;i<nelems;i++)
            if(PyList_Size(PyList_GET_ITEM(py_ob,i))!=NELEMS)
            {
                PyErr_Format(PyExc_TypeError,"size mismatch in A_rho");
                return -1;
            }
        
        for(size_t i=0;i<nelems;i++)
            for(size_t j=0;j<nelems;j++)
                if(potfit->ff->rho_ptr[i][j]->set_A(PyList_GET_ITEM(PyList_GET_ITEM(py_ob,i),j))!=0)
                   return -1;
        
        
        return 0;
    };
    
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<class FF,size_t NELEMS>
void PotFit<FF,NELEMS>::getset_dA_rho(PyGetSetDef& getset)
{
    getset.name=(char*)"dA_rho";
    getset.doc=(char*)"";
    getset.get=[](PyObject* self,void*)->PyObject*
    {
        PotFit<FF,NELEMS>* potfit=reinterpret_cast<Object*>(self)->potfit;
        PyObject* py_obj=PyList_New(NELEMS);
        for(size_t i=0;i<NELEMS;i++)
            PyList_SET_ITEM(py_obj,i,PyList_New(NELEMS));
        for(size_t i=0;i<NELEMS;i++)
            for(size_t j=0;j<NELEMS;j++)
                PyList_SET_ITEM(PyList_GET_ITEM(py_obj,i),j,potfit->ff->rho_ptr[i][j]->get_dA());
        
        return py_obj;
    };
    getset.set=[](PyObject* self,PyObject* val,void*)->int
    {
        PyErr_SetString(PyExc_TypeError,"readonly attribute");
        return -1;
    };
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<class FF,size_t NELEMS>
void PotFit<FF,NELEMS>::getset_dA_rho_max(PyGetSetDef& getset)
{
    getset.name=(char*)"dA_rho_max";
    getset.doc=(char*)"";
    getset.get=[](PyObject* self,void*)->PyObject*
    {
        PotFit<FF,NELEMS>* potfit=reinterpret_cast<Object*>(self)->potfit;
        PyObject* py_obj=PyList_New(NELEMS);
        for(size_t i=0;i<NELEMS;i++)
            PyList_SET_ITEM(py_obj,i,PyList_New(NELEMS));
        for(size_t i=0;i<NELEMS;i++)
            for(size_t j=0;j<NELEMS;j++)
                PyList_SET_ITEM(PyList_GET_ITEM(py_obj,i),j,potfit->ff->rho_ptr[i][j]->get_dA_max());
        
        return py_obj;
    };
    getset.set=[](PyObject* self,PyObject* val,void*)->int
    {
        PotFit<FF,NELEMS>* potfit=reinterpret_cast<Object*>(self)->potfit;
        VarAPI<OB<PyListObject,PyList_Type>> var("dA_rho_max");
        if(var.set(val)!=0) return -1;
        PyObject* py_ob=var.val;
        size_t nelems=PyList_Size(py_ob);
        if(nelems!=NELEMS)
        {
            PyErr_Format(PyExc_TypeError,"size mismatch in dA_rho_max");
            return -1;
        }
        
        for(size_t i=0;i<nelems;i++)
            if(PyList_Size(PyList_GET_ITEM(py_ob,i))!=NELEMS)
            {
                PyErr_Format(PyExc_TypeError,"size mismatch in dA_rho_max");
                return -1;
            }
        
        for(size_t i=0;i<nelems;i++)
            for(size_t j=0;j<nelems;j++)
                if(potfit->ff->rho_ptr[i][j]->set_dA_max(PyList_GET_ITEM(PyList_GET_ITEM(py_ob,i),j))!=0)
                   return -1;
        
        
        return 0;
    };
    
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<class FF,size_t NELEMS>
void PotFit<FF,NELEMS>::getset_A_rho_dof(PyGetSetDef& getset)
{
    getset.name=(char*)"A_rho_dof";
    getset.doc=(char*)"";
    getset.get=[](PyObject* self,void*)->PyObject*
    {
        PotFit<FF,NELEMS>* potfit=reinterpret_cast<Object*>(self)->potfit;
        PyObject* py_obj=PyList_New(NELEMS);
        for(size_t i=0;i<NELEMS;i++)
            PyList_SET_ITEM(py_obj,i,PyList_New(NELEMS));
        for(size_t i=0;i<NELEMS;i++)
            for(size_t j=0;j<NELEMS;j++)
                PyList_SET_ITEM(PyList_GET_ITEM(py_obj,i),j,potfit->ff->rho_ptr[i][j]->get_A_dof());
        
        return py_obj;
    };
    getset.set=[](PyObject* self,PyObject* val,void*)->int
    {
        PotFit<FF,NELEMS>* potfit=reinterpret_cast<Object*>(self)->potfit;
        VarAPI<OB<PyListObject,PyList_Type>> var("A_rho_dof");
        if(var.set(val)!=0) return -1;
        PyObject* py_ob=var.val;
        size_t nelems=PyList_Size(py_ob);
        if(nelems!=NELEMS)
        {
            PyErr_Format(PyExc_TypeError,"size mismatch in A_rho_dof");
            return -1;
        }
        
        for(size_t i=0;i<nelems;i++)
            if(PyList_Size(PyList_GET_ITEM(py_ob,i))!=NELEMS)
            {
                PyErr_Format(PyExc_TypeError,"size mismatch in A_rho_dof");
                return -1;
            }
        
        for(size_t i=0;i<nelems;i++)
            for(size_t j=0;j<nelems;j++)
                if(potfit->ff->rho_ptr[i][j]->set_A_dof(PyList_GET_ITEM(PyList_GET_ITEM(py_ob,i),j))!=0)
                   return -1;
        
        
        return 0;
    };
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<class FF,size_t NELEMS>
void PotFit<FF,NELEMS>::getset_A_F(PyGetSetDef& getset)
{
    getset.name=(char*)"A_F";
    getset.doc=(char*)"";
    getset.get=[](PyObject* self,void*)->PyObject*
    {
        PotFit<FF,NELEMS>* potfit=reinterpret_cast<Object*>(self)->potfit;
        PyObject* py_obj=PyList_New(NELEMS);
        for(size_t i=0;i<NELEMS;i++)
            PyList_SET_ITEM(py_obj,i,potfit->ff->F_ptr[i]->get_A());
        return py_obj;
    };
    getset.set=[](PyObject* self,PyObject* val,void*)->int
    {
        PotFit<FF,NELEMS>* potfit=reinterpret_cast<Object*>(self)->potfit;
        VarAPI<OB<PyListObject,PyList_Type>> var("A_F");
        if(var.set(val)!=0) return -1;
        PyObject* py_ob=var.val;
        size_t nelems=PyList_Size(py_ob);
        if(nelems!=NELEMS)
        {
            PyErr_Format(PyExc_TypeError,"size mismatch in A_F");
            return -1;
        }
        
        for(size_t i=0;i<nelems;i++)
            if(potfit->ff->F_ptr[i]->set_A(PyList_GET_ITEM(py_ob,i))!=0)
                   return -1;
        return 0;
    };
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<class FF,size_t NELEMS>
void PotFit<FF,NELEMS>::getset_dA_F(PyGetSetDef& getset)
{
    getset.name=(char*)"dA_F";
    getset.doc=(char*)"";
    getset.get=[](PyObject* self,void*)->PyObject*
    {
        PotFit<FF,NELEMS>* potfit=reinterpret_cast<Object*>(self)->potfit;
        PyObject* py_obj=PyList_New(NELEMS);
        for(size_t i=0;i<NELEMS;i++)
            PyList_SET_ITEM(py_obj,i,potfit->ff->F_ptr[i]->get_dA());
        return py_obj;
    };
    getset.set=[](PyObject* self,PyObject* val,void*)->int
    {
        PyErr_SetString(PyExc_TypeError,"readonly attribute");
        return -1;
    };
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<class FF,size_t NELEMS>
void PotFit<FF,NELEMS>::getset_dA_F_max(PyGetSetDef& getset)
{
    getset.name=(char*)"dA_F_max";
    getset.doc=(char*)"";
    getset.get=[](PyObject* self,void*)->PyObject*
    {
        PotFit<FF,NELEMS>* potfit=reinterpret_cast<Object*>(self)->potfit;
        PyObject* py_obj=PyList_New(NELEMS);
        for(size_t i=0;i<NELEMS;i++)
            PyList_SET_ITEM(py_obj,i,potfit->ff->F_ptr[i]->get_dA_max());
        return py_obj;
    };
    getset.set=[](PyObject* self,PyObject* val,void*)->int
    {
        PotFit<FF,NELEMS>* potfit=reinterpret_cast<Object*>(self)->potfit;
        VarAPI<OB<PyListObject,PyList_Type>> var("dA_F_max");
        if(var.set(val)!=0) return -1;
        PyObject* py_ob=var.val;
        size_t nelems=PyList_Size(py_ob);
        if(nelems!=NELEMS)
        {
            PyErr_Format(PyExc_TypeError,"size mismatch in dA_F_max");
            return -1;
        }
        
        for(size_t i=0;i<nelems;i++)
            if(potfit->ff->F_ptr[i]->set_dA_max(PyList_GET_ITEM(py_ob,i))!=0)
                return -1;
        return 0;
    };
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<class FF,size_t NELEMS>
void PotFit<FF,NELEMS>::getset_A_F_dof(PyGetSetDef& getset)
{
    getset.name=(char*)"A_F_dof";
    getset.doc=(char*)"";
    getset.get=[](PyObject* self,void*)->PyObject*
    {
        PotFit<FF,NELEMS>* potfit=reinterpret_cast<Object*>(self)->potfit;
        PyObject* py_obj=PyList_New(NELEMS);
        for(size_t i=0;i<NELEMS;i++)
            PyList_SET_ITEM(py_obj,i,potfit->ff->F_ptr[i]->get_A_dof());
        return py_obj;
    };
    getset.set=[](PyObject* self,PyObject* val,void*)->int
    {
        PotFit<FF,NELEMS>* potfit=reinterpret_cast<Object*>(self)->potfit;
        VarAPI<OB<PyListObject,PyList_Type>> var("A_F_dof");
        if(var.set(val)!=0) return -1;
        PyObject* py_ob=var.val;
        size_t nelems=PyList_Size(py_ob);
        if(nelems!=NELEMS)
        {
            PyErr_Format(PyExc_TypeError,"size mismatch in A_F_dof");
            return -1;
        }
        
        for(size_t i=0;i<nelems;i++)
            if(potfit->ff->F_ptr[i]->set_A_dof(PyList_GET_ITEM(py_ob,i))!=0)
                return -1;
        return 0;
    };
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<class FF,size_t NELEMS>
void PotFit<FF,NELEMS>::getset_ntally(PyGetSetDef& getset)
{
    getset.name=(char*)"ntally";
    getset.doc=(char*)"";
    
    getset.get=[](PyObject* self,void*)->PyObject*
    {
        return var<int>::build(reinterpret_cast<Object*>(self)->potfit->ntally);
    };
    getset.set=[](PyObject* self,PyObject* val,void*)->int
    {
        VarAPI<int> var("ntally");
        var.logics[0]=VLogics("ge",0);
        int ichk=var.set(val);
        if(ichk==-1) return -1;
        reinterpret_cast<Object*>(self)->potfit->ntally=var.val;
        return 0;
    };
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<class FF,size_t NELEMS>
void PotFit<FF,NELEMS>::getset_nmin_steps(PyGetSetDef& getset)
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
template<class FF,size_t NELEMS>
void PotFit<FF,NELEMS>::getset_tol(PyGetSetDef& getset)
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
template<class FF,size_t NELEMS>
void PotFit<FF,NELEMS>::getset_errs(PyGetSetDef& getset)
{
    getset.name=(char*)"errs";
    getset.doc=(char*)"";
    
    getset.get=[](PyObject* self,void*)->PyObject*
    {
        PotFit<FF,NELEMS>* potfit=reinterpret_cast<Object*>(self)->potfit;
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
template<class FF,size_t NELEMS>
void PotFit<FF,NELEMS>::getset_err(PyGetSetDef& getset)
{
    getset.name=(char*)"err";
    getset.doc=(char*)"";
    
    getset.get=[](PyObject* self,void*)->PyObject*
    {
        return var<type0>::build(reinterpret_cast<Object*>(self)->potfit->err);
    };
    getset.set=[](PyObject* self,PyObject* val,void*)->int
    {
        PyErr_SetString(PyExc_TypeError,"readonly attribute");
        return -1;
    };
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<class FF,size_t NELEMS>
void PotFit<FF,NELEMS>::getset_en_coefs(PyGetSetDef& getset)
{
    getset.name=(char*)"en_coefs";
    getset.doc=(char*)"";
    
    getset.get=[](PyObject* self,void*)->PyObject*
    {
        PotFit<FF,NELEMS>* potfit=reinterpret_cast<Object*>(self)->potfit;
        size_t sz;
        size_t* szp=&sz;
        sz=potfit->nconfigs;
        type0* v=potfit->get_en_coefs();
        PyObject* op=var<type0*>::build(v,&szp);
        Memory::dealloc(v);
        return op;
    };
    getset.set=[](PyObject* self,PyObject* val,void*)->int
    {
        VarAPI<type0*> var("en_coefs");
        if(var.set(val)!=0) return -1;
        PotFit<FF,NELEMS>* potfit=reinterpret_cast<Object*>(self)->potfit;
        if(var.__var__.size!=static_cast<size_t>(potfit->nconfigs))
        {
            PyErr_SetString(PyExc_TypeError,"size mismatch");
            return -1;
        }
        
        int __my_conf=potfit->my_conf;
        potfit->coef[0]=var.val[__my_conf];
        return 0;
    };
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<class FF,size_t NELEMS>
void PotFit<FF,NELEMS>::getset_f_coefs(PyGetSetDef& getset)
{
    getset.name=(char*)"f_coefs";
    getset.doc=(char*)"";
    
    getset.get=[](PyObject* self,void*)->PyObject*
    {
        PotFit<FF,NELEMS>* potfit=reinterpret_cast<Object*>(self)->potfit;
        size_t sz;
        size_t* szp=&sz;
        sz=potfit->nconfigs;
        type0* v=potfit->get_f_coefs();
        PyObject* op=var<type0*>::build(v,&szp);
        Memory::dealloc(v);
        return op;
    };
    getset.set=[](PyObject* self,PyObject* val,void*)->int
    {
        VarAPI<type0*> var("f_coefs");
        if(var.set(val)!=0) return -1;
        PotFit<FF,NELEMS>* potfit=reinterpret_cast<Object*>(self)->potfit;
        if(var.__var__.size!=static_cast<size_t>(potfit->nconfigs))
        {
            PyErr_SetString(PyExc_TypeError,"size mismatch");
            return -1;
        }
        
        int __my_conf=potfit->my_conf;
        potfit->f_coef=var.val[__my_conf];
        return 0;
    };
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<class FF,size_t NELEMS>
void PotFit<FF,NELEMS>::getset_S_coefs(PyGetSetDef& getset)
{
    getset.name=(char*)"S_coefs";
    getset.doc=(char*)"";
    
    getset.get=[](PyObject* self,void*)->PyObject*
    {
        PotFit<FF,NELEMS>* potfit=reinterpret_cast<Object*>(self)->potfit;
        size_t sz;
        size_t* szp=&sz;
        sz=potfit->nconfigs;
        type0(* v)[__dim__][__dim__]=reinterpret_cast<type0(*)[__dim__][__dim__]>(potfit->get_S_coefs());
        
        PyObject* op=var<type0(*)[__dim__][__dim__]>::build(v,&szp);
        Memory::dealloc(v);
        return op;
    };
    getset.set=[](PyObject* self,PyObject* val,void*)->int
    {
        VarAPI<symm<type0[__dim__][__dim__]>*> var("S_coefs");
        if(var.set(val)!=0) return -1;
        PotFit<FF,NELEMS>* potfit=reinterpret_cast<Object*>(self)->potfit;
        if(var.__var__.size!=static_cast<size_t>(potfit->nconfigs))
        {
            PyErr_SetString(PyExc_TypeError,"size mismatch");
            return -1;
        }
        int __my_conf=potfit->my_conf;
        Algebra::MSY_2_DyadicV(var.val[__my_conf],potfit->coef+1);
        
        return 0;
    };
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<class FF,size_t NELEMS>
void PotFit<FF,NELEMS>::getset_mins(PyGetSetDef& getset)
{
    getset.name=(char*)"mins";
    getset.doc=(char*)"";
    
    getset.get=[](PyObject* self,void*)->PyObject*
    {
        
        size_t __nconfigs=static_cast<size_t>(reinterpret_cast<Object*>(self)->potfit->nconfigs);
        MinCG::Object** mins=reinterpret_cast<Object*>(self)->mins;
        
        PyObject* py_obj=PyList_New(__nconfigs);
        for(size_t i=0;i<__nconfigs;i++)
            PyList_SET_ITEM(py_obj,i,reinterpret_cast<PyObject*>(mins[i]));
    
        return py_obj;
    };
    getset.set=[](PyObject* self,PyObject* val,void*)->int
    {
        VarAPI<OP<MinCG>*> var("mins");
        if(var.set(val)!=0) return -1;
        Object* __self=reinterpret_cast<Object*>(self);
        PotFit<FF,NELEMS>* potfit=__self->potfit;
        int __nconfigs=potfit->nconfigs;
        int __my_conf=potfit->my_conf;
        if(var.__var__.size!=static_cast<size_t>(__nconfigs))
        {
            PyErr_SetString(PyExc_TypeError,"size mismatch");
            return -1;
        }
        OP<MinCG>* __mins=var.val;
        
        for(int i=0;i<__nconfigs;i++)
        {
            if(i==__my_conf)
                __self->potfit->assign_min(reinterpret_cast<MinCG::Object*>(__mins[i].ob)->min);
            
            Py_DECREF(__self->mins[i]);
           __self->mins[i]=reinterpret_cast<MinCG::Object*>(__mins[i].ob);
           Py_INCREF(__self->mins[i]);
        }

        
        return 0;
    };
}
/*--------------------------------------------*/
template<class FF,size_t NELEMS>
PyMethodDef PotFit<FF,NELEMS>::methods[]=EmptyPyMethodDef(12);
/*--------------------------------------------*/
template<class FF,size_t NELEMS>
void PotFit<FF,NELEMS>::setup_tp_methods()
{
    ml_reset(methods[0]);
    ml_swtch2dflt(methods[1]);
    ml_min_struc(methods[2]);
    ml_min_cg(methods[3]);
    ml_min_sa(methods[4]);
    ml_mean_rho(methods[5]);
    ml_test_A_phi(methods[6]);
    ml_test_A_rho(methods[7]);
    ml_test_A_F(methods[8]);
    ml_val(methods[9]);
    ml_deriv(methods[10]);
}
/*--------------------------------------------

 --------------------------------------------*/
template<class FF,size_t NELEMS>
void PotFit<FF,NELEMS>::ml_val(PyMethodDef& tp_methods)
{
    tp_methods.ml_flags=METH_VARARGS | METH_KEYWORDS;
    tp_methods.ml_name="val";
    tp_methods.ml_meth=(PyCFunction)(PyCFunctionWithKeywords)(
    [](PyObject* self,PyObject* args,PyObject* kwds)->PyObject*
    {
        FuncAPI<>f("val");
        if(f(args,kwds)) return NULL;
        type0 ans=reinterpret_cast<Object*>(self)->potfit->val();
        return var<type0>::build(ans);
    });
    tp_methods.ml_doc=(char*)"";
}
/*--------------------------------------------

 --------------------------------------------*/
template<class FF,size_t NELEMS>
void PotFit<FF,NELEMS>::ml_deriv(PyMethodDef& tp_methods)
{
    tp_methods.ml_flags=METH_VARARGS | METH_KEYWORDS;
    tp_methods.ml_name="deriv";
    tp_methods.ml_meth=(PyCFunction)(PyCFunctionWithKeywords)(
    [](PyObject* self,PyObject* args,PyObject* kwds)->PyObject*
    {
        FuncAPI<>f("deriv");
        if(f(args,kwds)) return NULL;
        type0 ans=reinterpret_cast<Object*>(self)->potfit->deriv();
        return var<type0>::build(ans);
    });
    tp_methods.ml_doc=(char*)"";
}
/*--------------------------------------------

 --------------------------------------------*/
template<class FF,size_t NELEMS>
void PotFit<FF,NELEMS>::ml_reset(PyMethodDef& tp_methods)
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
    tp_methods.ml_doc=(char*)"";
}
/*--------------------------------------------

 --------------------------------------------*/
template<class FF,size_t NELEMS>
void PotFit<FF,NELEMS>::ml_swtch2dflt(PyMethodDef& tp_methods)
{
    tp_methods.ml_flags=METH_VARARGS | METH_KEYWORDS;
    tp_methods.ml_name="swtch2dflt";
    tp_methods.ml_meth=(PyCFunction)(PyCFunctionWithKeywords)(
    [](PyObject* self,PyObject* args,PyObject* kwds)->PyObject*
    {
        FuncAPI<>f("swtch2dflt");
        if(f(args,kwds)) return NULL;
        reinterpret_cast<Object*>(self)->potfit->swtch2dflt();
        Py_RETURN_NONE;
    });
    tp_methods.ml_doc=(char*)"";
}
/*--------------------------------------------

 --------------------------------------------*/
template<class FF,size_t NELEMS>
void PotFit<FF,NELEMS>::ml_min_struc(PyMethodDef& tp_methods)
{
    tp_methods.ml_flags=METH_VARARGS | METH_KEYWORDS;
    tp_methods.ml_name="min_struc";
    tp_methods.ml_meth=(PyCFunction)(PyCFunctionWithKeywords)(
    [](PyObject* self,PyObject* args,PyObject* kwds)->PyObject*
    {
        FuncAPI<>f("min_struc");
        if(f(args,kwds)) return NULL;
        reinterpret_cast<Object*>(self)->potfit->min_struc();
        Py_RETURN_NONE;
    });
    tp_methods.ml_doc=(char*)"";
}
/*--------------------------------------------

 --------------------------------------------*/
template<class FF,size_t NELEMS>
void PotFit<FF,NELEMS>::ml_min_cg(PyMethodDef& tp_methods)
{
    tp_methods.ml_flags=METH_VARARGS | METH_KEYWORDS;
    tp_methods.ml_name="min_cg";
    tp_methods.ml_meth=(PyCFunction)(PyCFunctionWithKeywords)(
    [](PyObject* self,PyObject* args,PyObject* kwds)->PyObject*
    {

        FuncAPI<int>f("min_cg",{"nsteps"});
        if(f(args,kwds)) return NULL;

        reinterpret_cast<Object*>(self)->potfit->min_cg(f.val<0>());

        Py_RETURN_NONE;
    });


    tp_methods.ml_doc=(char*)"";
}
/*--------------------------------------------

 --------------------------------------------*/
template<class FF,size_t NELEMS>
void PotFit<FF,NELEMS>::ml_min_sa(PyMethodDef& tp_methods)
{
    tp_methods.ml_flags=METH_VARARGS | METH_KEYWORDS;
    tp_methods.ml_name="min_sa";
    tp_methods.ml_meth=(PyCFunction)(PyCFunctionWithKeywords)(
    [](PyObject* self,PyObject* args,PyObject* kwds)->PyObject*
    {

        FuncAPI<int,int,int,type0,type0>f("min_sa",{"seed","nsteps","ntrials","temp","gamma"});
        if(f(args,kwds)) return NULL;

        reinterpret_cast<Object*>(self)->potfit->min_sa(f.val<0>(),f.val<1>(),
        f.val<2>(),f.val<3>(),f.val<4>());

        Py_RETURN_NONE;
    });


    tp_methods.ml_doc=(char*)"";
}
/*--------------------------------------------

 --------------------------------------------*/
template<class FF,size_t NELEMS>
void PotFit<FF,NELEMS>::ml_mean_rho(PyMethodDef& tp_methods)
{
    tp_methods.ml_flags=METH_VARARGS | METH_KEYWORDS;
    tp_methods.ml_name="mean_rho";
    tp_methods.ml_meth=(PyCFunction)(PyCFunctionWithKeywords)(
    [](PyObject* self,PyObject* args,PyObject* kwds)->PyObject*
    {

        FuncAPI<elem_type>f("mean_rho",{"ielem"});
        if(f(args,kwds)) return NULL;
        PotFit<FF,NELEMS>* potfit=reinterpret_cast<Object*>(self)->potfit;
        
        size_t sz;
        size_t* szp=&sz;
        sz=potfit->nconfigs;
        
        potfit->ff->derivative();
        type0* v=potfit->get_mean_rho(f.val<0>());
        
        PyObject* op=var<type0*>::build(v,&szp);
        Memory::dealloc(v);
        return op;
    });


    tp_methods.ml_doc=(char*)"";
}
/*--------------------------------------------

 --------------------------------------------*/
template<class FF,size_t NELEMS>
void PotFit<FF,NELEMS>::ml_test_A_phi(PyMethodDef& tp_methods)
{
    tp_methods.ml_flags=METH_VARARGS | METH_KEYWORDS;
    tp_methods.ml_name="test_A_phi";
    tp_methods.ml_meth=(PyCFunction)(PyCFunctionWithKeywords)(
    [](PyObject* self,PyObject* args,PyObject* kwds)->PyObject*
    {

        
        FuncAPI<int,int,int,type0,int>f("test_A_phi",{"i","j","k","delta","nsteps"});
        if(f(args,kwds)) return NULL;
        
        PotFit<FF,NELEMS>* potfit=reinterpret_cast<Object*>(self)->potfit;
        ptrdiff_t offset=potfit->ff->phi_ptr[f.val<0>()][f.val<1>()]->vars-potfit->ff->vs;
        offset+=f.val<2>();
        potfit->test_deriv(offset,f.val<3>(),f.val<4>());
        
        
        Py_RETURN_NONE;
    });


    tp_methods.ml_doc=(char*)"";
}
/*--------------------------------------------

 --------------------------------------------*/
template<class FF,size_t NELEMS>
void PotFit<FF,NELEMS>::ml_test_A_rho(PyMethodDef& tp_methods)
{
    tp_methods.ml_flags=METH_VARARGS | METH_KEYWORDS;
    tp_methods.ml_name="test_A_rho";
    tp_methods.ml_meth=(PyCFunction)(PyCFunctionWithKeywords)(
    [](PyObject* self,PyObject* args,PyObject* kwds)->PyObject*
    {

        
        FuncAPI<int,int,int,type0,int>f("test_A_rho",{"i","j","k","delta","nsteps"});
        if(f(args,kwds)) return NULL;
        
        PotFit<FF,NELEMS>* potfit=reinterpret_cast<Object*>(self)->potfit;
        ptrdiff_t offset=potfit->ff->rho_ptr[f.val<0>()][f.val<1>()]->vars-potfit->ff->vs;
        offset+=f.val<2>();
        potfit->test_deriv(offset,f.val<3>(),f.val<4>());
        
        
        Py_RETURN_NONE;
    });


    tp_methods.ml_doc=(char*)"";
}
/*--------------------------------------------

 --------------------------------------------*/
template<class FF,size_t NELEMS>
void PotFit<FF,NELEMS>::ml_test_A_F(PyMethodDef& tp_methods)
{
    tp_methods.ml_flags=METH_VARARGS | METH_KEYWORDS;
    tp_methods.ml_name="test_A_F";
    tp_methods.ml_meth=(PyCFunction)(PyCFunctionWithKeywords)(
    [](PyObject* self,PyObject* args,PyObject* kwds)->PyObject*
    {

        
        FuncAPI<int,int,type0,int>f("test_A_F",{"i","j","delta","nsteps"});
        if(f(args,kwds)) return NULL;
        
        PotFit<FF,NELEMS>* potfit=reinterpret_cast<Object*>(self)->potfit;
        ptrdiff_t offset=potfit->ff->F_ptr[f.val<0>()]->vars-potfit->ff->vs;
        offset+=f.val<1>();
        potfit->test_deriv(offset,f.val<2>(),f.val<3>());
        
        
        Py_RETURN_NONE;
    });


    tp_methods.ml_doc=(char*)"";
}
#endif
