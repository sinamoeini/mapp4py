#ifndef __MAPP__dmd_bdf__
#define __MAPP__dmd_bdf__
#include "dmd_implicit.h"
namespace MAPP_NS
{
    template<typename> class Vec;
    class DMDBDF:public DMDImplicit
    {
    private:
    protected:
        void nonlin_fail();
    public:
        
        static constexpr int max_q=5;
        int q;
        int dq;
        
        int nconst_q;
        int nconst_dt;
        
        
        type0 eta[3];
        
        type0 lo_err_fac[2];
        type0 hi_err_fac[2];
        
        
        type0 t[max_q+1];
        type0 l[max_q+1];
        
        type0 A[max_q+1][max_q+1];
        type0 A_bar[max_q+1][max_q+1];
        
        type0* z;
        type0* dy;
        

        
        
        
        
        
        
        
        
        
        
        DMDBDF();
        ~DMDBDF();
        
        void init_static();
        void fin_static();
        void run_static(type0);
        bool integrate();
        void integrate_fail();
        bool interpolate();
        void interpolate_fail();
        void update_z();
        void prep_for_next();
        void reset();

        
        
        
        
        
        
        
        typedef struct
        {
            PyObject_HEAD
            DMDBDF* dmd;
        }Object;
        
        static PyTypeObject TypeObject;
        static PyObject* __new__(PyTypeObject*,PyObject*, PyObject*);
        static int __init__(PyObject*, PyObject*,PyObject*);
        static PyObject* __alloc__(PyTypeObject*,Py_ssize_t);
        static void __dealloc__(PyObject*);
        
        
        static PyMethodDef methods[];
        static void setup_tp_methods();
        static void ml_run(PyMethodDef&);
        
        
        static PyGetSetDef getset[];
        static void setup_tp_getset();
        
        static void setup_tp();

    };
}
using namespace MAPP_NS;
#include "xmath.h"
/*--------------------------------------------
 
 --------------------------------------------*/
namespace MAPP_NS
{
    namespace DMDBDFMath
    {

        template<const int i,const int dim>
        class A_z_y_l
        {
        public:
            template<class T>
            static inline void func0(T* RESTRICT A,T* RESTRICT z,T dy,T* RESTRICT l)
            {
                *z=Algebra::V_mul_V<i>(A,z)+dy**l;
                A_z_y_l<i-1,dim>::func0(A+DMDBDF::max_q+2,z+1,dy,l+1);
            }
        };
        
        template<const int dim>
        class A_z_y_l<1,dim>
        {
            public:
            template<class T>
            static inline void func0(T* RESTRICT A,T* RESTRICT z,T dy,T* RESTRICT l)
            {
                *z=*A**z+dy**l;
            }
        };
        
        
        
        

        
        
        template<const int dim>
        void update_z(const int n,type0* RESTRICT A,type0* RESTRICT z,type0* RESTRICT dy,type0* RESTRICT l)
        {
            for(int i=0;i<n;i++)
            {
                A_z_y_l<dim+1,dim+1>::func0(A,z,dy[i],l);
                z+=DMDBDF::max_q+1;
            }
            
        }
        
        template<const int dim>
        void update_z_inc(const int n,type0* RESTRICT A,type0* RESTRICT z,type0* RESTRICT dy,type0* RESTRICT l)
        {
            for(int i=0;i<n;i++)
            {
                z[dim+1]=l[dim+1]*dy[i];
                A_z_y_l<dim+1,dim+1>::func0(A,z,dy[i],l);
                z+=DMDBDF::max_q+1;
            }
        }
        
        
        
        template<const int dim>
        bool interpolate(const int n,const type0 beta,const type0* RESTRICT A0,const type0* RESTRICT A1,const type0* RESTRICT z,type0* RESTRICT y_0,type0* RESTRICT a)
        {
            volatile type0 y0;
            for(int i=0;i<n;i++)
            {
                y0=Algebra::V_mul_V<dim+1>(A0,z);
                --++y0;
                if(y0<0.0 || y0>1.0) return false;
                
                y_0[i]=y0;
                a[i]=Algebra::V_mul_V<dim>(A1,z+1)*beta-y0;
                z+=DMDBDF::max_q+1;
            }
            
            return true;
        }
        
        
        
        
        
        
        
        
        
        
        
        
        
        inline void l_calc(int q,type0& dt,type0 (&t)[DMDBDF::max_q+1],type0 (&l)[DMDBDF::max_q+1])
        {
            l[0]=1.0;
            for(int i=1;i<q+1;i++)
                l[i]=0.0;
            
            type0 k0;
            
            for(int i=0;i<q;i++)
            {
                k0=1.0/(dt-t[i]);
                for(int j=i+1;j>0;j--)
                    l[j]+=k0*l[j-1];
            }
        }
        
        
        // f_flc
        inline type0 beta_calc(int& q,type0& dt,type0 (&t)[DMDBDF::max_q+1])
        {
            type0 iq=0.0;
            type0 beta_inv=0.0;
            for(int i=0;i<q;i++,iq++)
                beta_inv+=1.0/(iq+1.0);
            beta_inv/=dt;
            return 1.0/beta_inv;
        }
        
        
        
        // f_flc
        inline void prep_A_bar_l(int& q,type0& dt,type0 (&t)[DMDBDF::max_q+1],int dq,type0 (&l)[DMDBDF::max_q+1],type0 (&A_bar)[DMDBDF::max_q+1][DMDBDF::max_q+1])
        {
            if(dq==1)
            {
                type0 k2=0.0,iq=0.0;
                for(int i=0;i<q;i++,iq++)
                    k2+=1.0/(dt*(iq+1.0))-1.0/(dt-t[i]);
                
                
                l_calc(q,dt,t,l);
                l[q+1]=0.0;
                for(int i=q+1;i>0;i--)
                    l[i]+=k2*l[i-1];
            }
            else if(dq==-1)
            {
                l_calc(q-2,dt,t,l);
                type0 p_q_2=1.0;
                for(int i=0;i<q-2;i++)
                    p_q_2*=dt-t[i];
                for(int i=2;i<q;i++)
                    A_bar[i][q]-=p_q_2*l[i-2];
                
                
                type0 k1=0.0,iq=0.0;
                for(int i=0;i<q-2;i++,iq++)
                    k1+=1.0/(dt*(iq+1.0))-1.0/(dt-t[i]);
                k1+=1.0/(dt*(iq+1.0))+1.0/(dt*(iq+2.0));
                
                
                l[q-1]=0.0;
                for(int i=q-1;i>0;i--)
                    l[i]+=k1*l[i-1];
            }
            else
            {
                l_calc(q-1,dt,t,l);
                
                type0 k0=0.0,iq=0.0;
                for(int i=0;i<q-1;i++,iq++)
                    k0+=1.0/(dt*(iq+1.0))-1.0/(dt-t[i]);
                k0+=1.0/(dt*(iq+1.0));
                
                l[q]=0.0;
                for(int i=q;i>0;i--)
                    l[i]+=k0*l[i-1];
            }
        }
        
        
        
        
        
        
        
        // f_flc
        inline type0  err_fac_calc(int& q,type0& dt,type0 (&t)[DMDBDF::max_q+1],type0 (&lo_err_fac)[2],type0 (&hi_err_fac)[2])
        {
            type0 u_q=0.0,s_bar_q=0.0;
            type0 iq=0.0;
            
            for(int i=0;i<q;i++,iq++)
            {
                s_bar_q+=1.0/(iq+1.0);
                u_q+=1.0/(dt-t[i]);
            }
            
            
            u_q*=dt;
            
            u_q++;
            u_q-=s_bar_q;
            
            type0 err_fac=fabs(u_q/(s_bar_q*(1.0+iq*u_q)));
            
            if(q<DMDBDF::max_q)
            {
                type0 r=1.0;
                for(int i=1;i<q;i++)
                    r*=(t[i]-dt)/t[i];
                type0 u_p_q=0.0;
                for(int i=0;i<q;i++)
                    u_p_q+=1.0/t[i+1];
                
                u_p_q*=t[1];
                u_p_q++;
                u_p_q-=s_bar_q;
                
                type0 tmp=(1.0-t[q]/dt)*(u_q-1.0/(iq+1.0)+dt/(dt-t[q]));
                tmp/=(1.0+iq*u_q)*(s_bar_q+1.0/(iq+1.0))*(iq+2.0);
                hi_err_fac[0]=fabs(tmp);
                hi_err_fac[1]=-dt*dt*r*(1.0+iq*u_q)/(t[1]*t[q]*(1.0+iq*u_p_q));
            }
            
            if(q>1)
            {
                type0 p_q=1.0;
                for(int i=0;i<q;i++,iq++)
                    p_q*=dt-t[i];
                
                type0 tmp=u_q-dt/(dt-t[q-1])+1.0/iq;
                tmp/=s_bar_q-1.0/iq;
                tmp/=1.0-t[q-1]/dt;
                
                lo_err_fac[0]=fabs(tmp);
                lo_err_fac[1]=p_q;
            }
            return err_fac;
        }
        

        
        
    }
}

#endif
