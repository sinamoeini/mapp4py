#ifndef __MAPP__md_nst__
#define __MAPP__md_nst__
#include "md_nvt.h"



namespace MAPP_NS
{
    
    class MDNST:public MDNVT
    {
    private:
        type0 MLT0[__dim__][__dim__];
        type0 MLT1[__dim__][__dim__];
        type0 MLT2[__dim__][__dim__];
        type0 S_dev[__dim__][__dim__];
        type0 s_hyd;
        type0 B_ref[__dim__][__dim__];
        
        
        type0 V_H_prefac[__dim__][__dim__];
        
    
        type0 v0[__dim__];
        type0 v1[__dim__];
        
        
        type0 V_H[__dim__][__dim__];
        type0 vol_ref;
    protected:
        
        
        
        type0 T_baro;
        int ndof_baro;
        ThermostatNHC thermo_baro;
        
        
        
                
        bool S_dof[__dim__][__dim__];
        type0 S[__dim__][__dim__];
        type0 t_relax_S[__dim__][__dim__];
        int nreset;
        
        void update_x();
        void update_x_d(type0=1.0);
        void update_x_d__x__x_d(type0);
        void update_V_H();
        void change_dt(type0);
        void pre_run_chk(AtomsMD*,ForceFieldMD*);
        void pre_init();
    public:
        MDNST(type0(&)[__dim__][__dim__],type0,type0);
        virtual ~MDNST();
        
        void init();
        void run(int);
        void fin();
        
        
        
        
        typedef struct
        {
            PyObject_HEAD
            MDNST* md;
            ExportMD::Object* xprt;
        }Object;
    
        
        static PyTypeObject TypeObject;
        static PyObject* __new__(PyTypeObject*,PyObject*, PyObject*);
        static int __init__(PyObject*, PyObject*,PyObject*);
        static PyObject* __alloc__(PyTypeObject*,Py_ssize_t);
        static void __dealloc__(PyObject*);
        
        static PyGetSetDef getset[];
        static void setup_tp_getset();
        static void getset_niters_baro(PyGetSetDef&);
        static void getset_nlinks_baro(PyGetSetDef&);
        static void getset_t_relax_baro(PyGetSetDef&);
        static void getset_S_dof(PyGetSetDef&);
        static void getset_S(PyGetSetDef&);
        static void getset_t_relax_S(PyGetSetDef&);
        static void getset_nreset(PyGetSetDef&);
        
        static int setup_tp();
    };
}

/*--------------------------------------------
 
 --------------------------------------------*/
#include "xmath.h"
namespace MAPP_NS
{
    namespace MDMath
    {
        
        static const type0 sqrt_eps=sqrt(std::numeric_limits<type0>::epsilon());
        template<const int dim>
        void calc(type0&,type0(&)[dim][dim],type0(&)[dim][dim],type0(&)[dim][dim]){};
        void calc(const type0,const type0&,const type0&,type0(&)[3][3],type0(&)[3][3],type0(&)[3][3],type0(&)[3][3]);
        void calc(const type0,const type0&,const type0&,type0(&)[2][2],type0(&)[2][2],type0(&)[2][2],type0(&)[2][2]);
        void calc(const type0,const type0&,const type0&,type0(&)[1][1],type0(&)[1][1],type0(&)[1][1],type0(&)[1][1]);
        void calc(type0&,type0(&)[3][3],type0(&)[3][3],type0(&)[3][3]);
        void calc(type0&,type0(&)[2][2],type0(&)[2][2],type0(&)[2][2]);
        void calc(type0&,type0(&)[1][1],type0(&)[1][1],type0(&)[1][1]);
        
        type0 f(type0);
        type0 df(type0);
        type0 ddf(type0);
        type0 dddf(type0);
        
        template<const int i,const int dim>
        class ____NONAME0
        {
        public:
            template<class T>
            static inline void func(const T& xi,T* x_d,T* MLT_x_d,const T& m_inv,T* f,T* MLT_f)
            {
                *x_d=xi*Algebra::__V_strd_mul_V<i,dim>::func(MLT_x_d,x_d)+m_inv*Algebra::__V_strd_mul_V<i,dim>::func(MLT_f,f);
                ____NONAME0<i-1,dim>::func(xi,x_d+1,MLT_x_d+dim+1,m_inv,f+1,MLT_f+dim+1);
            }
            
            template<class T>
            static inline void func(T* x_d,T* MLT_x_d,const T& m_inv,T* f,T* MLT_f)
            {
                *x_d=Algebra::__V_strd_mul_V<i,dim>::func(MLT_x_d,x_d)+m_inv*Algebra::__V_strd_mul_V<i,dim>::func(MLT_f,f);
                ____NONAME0<i-1,dim>::func(x_d+1,MLT_x_d+dim+1,m_inv,f+1,MLT_f+dim+1);
            }
        };
        
        template<const int dim>
        class ____NONAME0<1,dim>
        {
        public:
            template<class T>
            static inline void func(const T& xi,T* x_d,T* MLT_x_d,const T& m_inv,T* f,T* MLT_f)
            {
                *x_d=xi**MLT_x_d**x_d+m_inv**MLT_f**f;
            }
            template<class T>
            static inline void func(T* x_d,T* MLT_x_d,const T& m_inv,T* f,T* MLT_f)
            {
                *x_d=*MLT_x_d**x_d+m_inv**MLT_f**f;
            }
        };
        
        
        template<const int i,const int dim>
        class ____NONAME1
        {
        public:
            template<class T>
            static inline void func(T* x,T* MLT_x,T* x_d,T* dx)
            {
                *dx=Algebra::__V_strd_mul_V<i,dim>::func(MLT_x,x)+*x_d;
                ____NONAME1<i-1,dim>::func(x+1,MLT_x+dim+1,x_d+1,dx+1);
            }
        };
        
        template<const int dim>
        class ____NONAME1<1,dim>
        {
        public:
            template<class T>
            static inline void func(T* x,T* MLT_x,T* x_d,T* dx)
            {
                *dx=*MLT_x**x+*x_d;
            }
        };
        
        template<const int i,const int dim>
        class ____NONAME2
        {
        public:
            template<class T>
            static inline void func(T* dx,T* MLT)
            {
                *dx=Algebra::__V_strd_mul_V<i,dim>::func(MLT,dx);
                ____NONAME2<i-1,dim>::func(dx+1,MLT+dim+1);
            }
        };
        
        template<const int dim>
        class ____NONAME2<1,dim>
        {
        public:
            template<class T>
            static inline void func(T* dx,T* MLT)
            {
                *dx*=*MLT;
            }
        };
    }
}


#endif













