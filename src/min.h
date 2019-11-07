#ifndef __MAPP__min__
#define __MAPP__min__
#include "api.h"
#include "ls.h"
#include "global.h"
#include "atoms_styles.h"
#include "export_styles.h"
#include "xmath.h"
#include "thermo_dynamics.h"
#include "MAPP.h"
namespace MAPP_NS
{
    namespace MinHelper
    {
        template<bool...Bs0>
        class CondB
        {
        public:
            template<class F,class LS>
            static void run(F& f,int n,LS* ls)
            {
                return f.template __run<Bs0...>(ls,n);
            }
            
            template<class F,class LS, class... Bs>
            static void run(F& f,int n,LS* ls,bool b0,Bs...bs)
            {
                if(b0==true)
                    return CondB<Bs0...,true>::run(f,n,ls,bs...);
                else
                    return CondB<Bs0...,false>::run(f,n,ls,bs...);
            }
            
            
            template<class F>
            static void init(F& f)
            {
                return f.template __init<Bs0...>();
            }
            
            template<class F,class... Vs>
            static void init_var(F& f,Vs*&... vs)
            {
                return f.template __init<Bs0...>(vs...);
            }
            
            template<class F, class... Bs>
            static void init(F& f,bool b0,Bs...bs)
            {
                if(b0==true)
                    return CondB<Bs0...,true>::init(f,bs...);
                else
                    return CondB<Bs0...,false>::init(f,bs...);
            }
            
            template<class F, class... Bs,class... Vs>
            static void init_var(F& f,bool b0,Bs...bs,Vs*&... vs)
            {
                if(b0==true)
                    return CondB<Bs0...,true>::init_var(f,bs...,vs...);
                else
                    return CondB<Bs0...,false>::init_var(f,bs...,vs...);
            }
            
            template<class F>
            static void fin(F& f)
            {
                return f.template __fin<Bs0...>();
            }
            
            template<class F, class... Bs>
            static void fin(F& f,bool b0,Bs...bs)
            {
                if(b0==true)
                    return CondB<Bs0...,true>::fin(f,bs...);
                else
                    return CondB<Bs0...,false>::fin(f,bs...);
            }
            
        };
        
        
        template<class LS0,class... LSs>
        class CondLS
        {
        public:
            template<class F,class LS00, class... Bs>
            static void run(F& f,int n,LS00* ls,Bs...bs)
            {
                if(dynamic_cast<LS0*>(ls)!=NULL)
                    return CondB<>::run(f,n,dynamic_cast<LS0*>(ls),bs...);
                return CondLS<LSs...>::run(f,n,ls,bs...);
                
            }
        };
        
        template<class LS0>
        class CondLS<LS0>
        {
        public:
            template<class F,class LS00, class... Bs>
            static void run(F& f,int n,LS00* ls,Bs...bs)
            {
                return CondB<>::run(f,n,dynamic_cast<LS0*>(ls),bs...);
            }
        };
        
        template<template<bool ...> class F,int N,bool ...BS>
        class MaxSizeAlign
        {
        public:
            static constexpr int MaxS=(MaxSizeAlign<F,N-1,BS...,true>::MaxS > MaxSizeAlign<F,N-1,BS...,false>::MaxS) ? MaxSizeAlign<F,N-1,BS...,true>::MaxS : MaxSizeAlign<F,N-1,BS...,false>::MaxS;
            static constexpr int MaxA=(MaxSizeAlign<F,N-1,BS...,true>::MaxA > MaxSizeAlign<F,N-1,BS...,false>::MaxA) ? MaxSizeAlign<F,N-1,BS...,true>::MaxA : MaxSizeAlign<F,N-1,BS...,false>::MaxA;
        };
        
        template<template<bool ...> class F,bool ...BS>
        class MaxSizeAlign<F,0,BS...>
        {
        public:
            static constexpr int MaxS=sizeof(F<BS...>);
            static constexpr int MaxA=alignof(F<BS...>);
        };
        
        inline void prep_x_d(const int natms_lcl,
        type0* x,type0 (*x_A)[__dim__],
        type0* h,type0 (*h_A)[__dim__],
        type0* x_d,type0 (*x_d_A)[__dim__])
        {

            Algebra::MLT_mul_MLT(x_A,h_A,x_d_A);
            memcpy(x_d,h,natms_lcl*__dim__*sizeof(type0));
            for(int iatm=0;iatm<natms_lcl;iatm++,x+=__dim__,x_d+=__dim__)
                Algebra::V_mul_MLT_add_in(x,h_A,x_d);
        }
        inline void prep_x_d_affine(const int natms_lcl,
        type0* x,type0 (*x_A)[__dim__],
        type0 (*h_A)[__dim__],
        type0* x_d,type0 (*x_d_A)[__dim__])
        {

            Algebra::MLT_mul_MLT(x_A,h_A,x_d_A);
            for(int iatm=0;iatm<natms_lcl;iatm++,x+=__dim__,x_d+=__dim__)
                Algebra::V_mul_MLT(x,h_A,x_d);
        }

        
        template<bool>
        class ThermoHandler
        {
        public:
            template<class F>
            static void init(F&,int&)
            {}
            template<class F>
            static void print(F&,int&,int&)
            {}
            template<class F>
            static void fin(F&,int&,int&,const char*&)
            {}
        };
        
        template<>
        class ThermoHandler<true>
        {
        public:
            template<class F>
            static void init(F& f,int& step)
            {
                f.create_thermo();
                f.thermo->init();
                f.thermo->print(step);
            }
            template<class F>
            static void print(F& f,int& step,int& istep)
            {
                if((istep+1)%f.ntally==0)
                    f.thermo->print(step+istep+1);
            }
            
            template<class F>
            static void fin(F& f,int& step,int& istep,const char*& __err)
            {
                if(istep%f.ntally)
                    f.thermo->print(step+istep);
                f.thermo->fin();
                MAPP::print_stdout("%s",__err);
            }
            
        };
        
        template<bool>
        class ExportHandler
        {
        public:
            template<class F>
            static void init(F&,int&)
            {}
            template<class F>
            static void print(F&,int&,int&)
            {}
            template<class F>
            static void fin(F&,int&,int&)
            {}
            
        };
        
        template<>
        class ExportHandler<true>
        {
        public:
            template<class F>
            static void init(F& f,int& step)
            {
                f.xprt->write(step);
            }
            template<class F>
            static void print(F& f,int& step,int& istep)
            {
                if((istep+1)%f.xprt->nevery==0) f.xprt->write(step+istep+1);
            }
            
            template<class F>
            static void fin(F& f,int& step,int& istep)
            {
                if(istep%f.xprt->nevery) f.xprt->write(step+istep);
            }
            
        };
        
    }
    class Min
    {
    template<bool>
    friend class MinHelper::ThermoHandler;
    template<bool>
    friend class MinHelper::ExportHandler;
    private:
    protected:
        static const char* err_msgs[];
        type0 max_dx;
        bool chng_box;
        type0 f_h;
        bool affine;
        type0 e_tol;
        bool H_dof[__dim__][__dim__];
        int ntally;
        class LineSearch* ls;
        typedef std::aligned_storage<sizeof(ThermoDynamics),alignof(ThermoDynamics)>::type MemThermoDynamics;
        MemThermoDynamics thermo_mem;
        ThermoDynamics* thermo;
        
        
        void pre_run_chk(Atoms*,ForceField*);
    public:
        Min();
        Min(type0,bool(&)[__dim__][__dim__],bool,type0,class LineSearch*);
        virtual ~Min();
        
        typedef struct
        {
            PyObject_HEAD
            Min* min;
            LineSearch::Object* ls;
            Export::Object* xprt;
        }Object;
        
        static PyTypeObject TypeObject;
        static PyObject* __new__(PyTypeObject*,PyObject*, PyObject*);
        static int __init__(PyObject*, PyObject*,PyObject*);
        static PyObject* __alloc__(PyTypeObject*,Py_ssize_t);
        static void __dealloc__(PyObject*);
        
        static PyMethodDef methods[];
        static void setup_tp_methods();
        
        static PyGetSetDef getset[];
        static void setup_tp_getset();
        static void getset_e_tol(PyGetSetDef&);
        static void getset_affine(PyGetSetDef&);
        static void getset_H_dof(PyGetSetDef&);
        static void getset_max_dx(PyGetSetDef&);
        static void getset_ntally(PyGetSetDef&);
        static void getset_ls(PyGetSetDef&);
        
        static int setup_tp();
    };
}
#endif
