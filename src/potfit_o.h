#ifndef __MAPP__potfit_o__
#define __MAPP__potfit_o__
#include "atoms_md.h"
#include "min_cg_fit.h"
#include "ff_eam_fit_o.h"
#include "thermo_dynamics.h"
namespace MAPP_NS
{
    class PotFitO
    {
    private:
        int ntrial;
        type0 err_tol;
        int nmin_steps;
        static const char* err_msgs[];
        
    protected:
    public:
        PotFitO(
               type0(&)[nrho_H],type0(&)[nphi_FeH],type0(&)[nphi_HH],type0(&)[nF_H],
               std::string*&&,std::string*&,int*&,type0(*&)[1+__nvoigt__],type0(*&)[1+__nvoigt__],PyObject**,size_t,MPI_Comm&);
        
        ~PotFitO();
        
        ThermoDynamics get_thermo();
        type0 F(type0);
        void ls_prep(type0&,type0&,type0&);
        void F_reset();
        size_t get_rFeH(type0*&,int*&);
        
        
        int nconfigs;
        type0* errs;
        int* roots;
        std::string* names_str;
        const char** names;
        MPI_Comm world;
        MPI_Comm* my_world;
        int my_rank,my_conf,my_lcl_rank;
        AtomsMD* atoms;
        MinCGFit* min;
        LineSearchBrent* min_ls;
        ForceFieldEAMFitO* ff;
        int* nHs;
        type0* mean_rhoHs;
        
        VecTens<type0,1> X0;
        VecTens<type0,1> Xorig;
        
        
        type0 find_max_alpha_A_rho_H();
        type0 find_max_alpha_A_phi_FeH();
        type0 find_max_alpha_A_phi_HH();
        type0 find_max_alpha_A_F_H();
        
       
        
        type0 tot_err;
        void calc_nHs();
        void calc_mean_rho_Hs();
        type0 iter();
        void min_cg(int);
        void store_x0();
        void restore_x0();
        void full_reset();
        void gen_test(int,type0,int);
        
        int max_ntrials;
        type0 coef;
        type0 target;
        type0 f_lcl[nHvars];
        type0 f[nHvars];
        type0 f0[nHvars];
        type0 h[nHvars];
        type0 x0[nHvars];
        type0 dx_max[nHvars];
        
        type0* get_coefs();
        void set_coefs(type0*);
        
        
        typedef struct
        {
            PyObject_HEAD
            class PotFitO* potfit;
        }Object;
        
        
        
        static PyTypeObject TypeObject;
        static PyObject* __new__(PyTypeObject*,PyObject*, PyObject*);
        static int __init__(PyObject*, PyObject*,PyObject*);
        static PyObject* __alloc__(PyTypeObject*,Py_ssize_t);
        static void __dealloc__(PyObject*);
        
        
        static PyGetSetDef getset[];
        static void setup_tp_getset();

        
        
        static void getset_A_rho_H(PyGetSetDef&);
        static void getset_A_phi_FeH(PyGetSetDef&);
        static void getset_A_phi_HH(PyGetSetDef&);
        static void getset_A_F_H(PyGetSetDef&);
        
        static void getset_dA_rho_H_max(PyGetSetDef&);
        static void getset_dA_phi_FeH_max(PyGetSetDef&);
        static void getset_dA_phi_HH_max(PyGetSetDef&);
        static void getset_dA_F_H_max(PyGetSetDef&);
        
        static void getset_A_rho_H_dof(PyGetSetDef&);
        static void getset_A_phi_FeH_dof(PyGetSetDef&);
        static void getset_A_phi_HH_dof(PyGetSetDef&);
        static void getset_A_F_H_dof(PyGetSetDef&);
        static void getset_mean_rho_H(PyGetSetDef&);
        static void getset_max_ntrials(PyGetSetDef&);
        static void getset_nmin_steps(PyGetSetDef&);
        static void getset_tol(PyGetSetDef&);
        static void getset_coefs(PyGetSetDef&);
        
        
        static void getset_RFeH(PyGetSetDef&);
        
        
        
        
        static PyMethodDef methods[];
        static void setup_tp_methods();
        static void ml_reset(PyMethodDef&);
        static void ml_min(PyMethodDef&);
        static void ml_test_A_rho_H(PyMethodDef&);
        static void ml_test_A_F_H(PyMethodDef&);
        static void ml_test_A_phi_FeH(PyMethodDef&);
        static void ml_test_A_phi_HH(PyMethodDef&);
        
        
        static int setup_tp();
        
        
        
        type0 find_max_alpha(const type0,const type0,bool,type0,type0,type0);
        
        

        


    };
    
    
}


#endif
