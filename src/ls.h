#ifndef __MAPP__ls__
#define __MAPP__ls__
#include <Python.h>
#include "global.h"
#include <limits>
#include <stdio.h>
namespace MAPP_NS
{
    enum
    {
        LS_S,
        LS_F_DOWNHILL,
        LS_F_GRAD0,
        LS_MIN_ALPHA,
        
        MIN_S_TOLERANCE,
        MIN_F_MAX_ITER,
        
        B_S,
        B_F_MAX_ALPHA,
        B_F_DOWNHILL
    };
}
#include <cmath>
using namespace MAPP_NS;
/*------------------------------------------------------------------------------------------------------------------------------------
 
 ------------------------------------------------------------------------------------------------------------------------------------*/
namespace MAPP_NS
{
    class LineSearch
    {
    private:
    protected:
    
        
        const type0 epsilon_3_4;
        const type0 sqrt_epsilon;
        const type0 golden;
        type0 prev_val;
    public:
        template<class Func>
        int bracket(Func*,type0,type0,type0&,type0&,type0&,type0&,type0&,type0&);
    
        LineSearch();
        virtual ~LineSearch(){};
        //template<class Func>int min(Func*,type0&,type0&,int){return 0;}
        
        struct Object;
        
        static PyTypeObject TypeObject;
        static PyObject* __new__(PyTypeObject*,PyObject*, PyObject*);
        static int __init__(PyObject*, PyObject*,PyObject*);
        static PyObject* __alloc__(PyTypeObject*,Py_ssize_t);
        static void __dealloc__(PyObject*);
        
        static PyMethodDef methods[];
        static void setup_tp_methods();
        
        static PyGetSetDef getset[];
        static void setup_tp_getset();
        
        static int setup_tp();
    };
    
    struct LineSearch::Object
    {
        PyObject_HEAD
        LineSearch ls;
    };
}
/*--------------------------------------------
 
 --------------------------------------------*/
namespace MAPP_NS
{
    class LineSearchGoldenSection:public LineSearch
    {
    private:
    protected:
    public:
        type0 tol;
        int max_iter;
        bool brack;
    
        LineSearchGoldenSection();
        ~LineSearchGoldenSection(){};
        template<class Func>
        int min(Func*,type0&,type0&,int);
        
        struct Object;
        
        static PyTypeObject TypeObject;
        static PyObject* __new__(PyTypeObject*,PyObject*, PyObject*);
        static int __init__(PyObject*, PyObject*,PyObject*);
        static PyObject* __alloc__(PyTypeObject*,Py_ssize_t);
        static void __dealloc__(PyObject*);
        
        static PyMethodDef methods[];
        static void setup_tp_methods();
        
        static PyGetSetDef getset[];
        static void setup_tp_getset();
        static void getset_bracket(PyGetSetDef&);
        static void getset_max_iter(PyGetSetDef&);
        static void getset_tol(PyGetSetDef&);
        static int setup_tp();
    };
    
    struct LineSearchGoldenSection::Object
    {
        PyObject_HEAD
        LineSearchGoldenSection ls;
    };
}
/*--------------------------------------------
 
 --------------------------------------------*/
namespace MAPP_NS
{
    class LineSearchBrent:public LineSearch
    {
    private:
    protected:
    public:
        type0 tol,zeps;
        int max_iter;
        bool brack;
    
        LineSearchBrent();
        ~LineSearchBrent(){};
        template<class Func>
        int min(Func*,type0&,type0&,int);
        
        struct Object;
        
        static PyTypeObject TypeObject;
        static PyObject* __new__(PyTypeObject*,PyObject*, PyObject*);
        static int __init__(PyObject*, PyObject*,PyObject*);
        static PyObject* __alloc__(PyTypeObject*,Py_ssize_t);
        static void __dealloc__(PyObject*);
        
        static PyMethodDef methods[];
        static void setup_tp_methods();
        
        static PyGetSetDef getset[];
        static void setup_tp_getset();
        static void getset_bracket(PyGetSetDef&);
        static void getset_max_iter(PyGetSetDef&);
        static void getset_tol(PyGetSetDef&);
        static void getset_zeps(PyGetSetDef&);
        static int setup_tp();
    };
    
    struct LineSearchBrent::Object
    {
        PyObject_HEAD
        LineSearchBrent ls;
    };
}
/*--------------------------------------------
 
 --------------------------------------------*/
namespace MAPP_NS
{
    class LineSearchBackTrack:public LineSearch
    {
    private:
    protected:
    public:
        type0 c,rho,min_alpha;
    
        LineSearchBackTrack();
        ~LineSearchBackTrack(){};
        template<class Func>
        int min(Func*,type0&,type0&,int);
        
        struct Object;
        struct Object1;
        
        static PyTypeObject TypeObject;
        static PyObject* __new__(PyTypeObject*,PyObject*, PyObject*);
        static int __init__(PyObject*, PyObject*,PyObject*);
        static PyObject* __alloc__(PyTypeObject*,Py_ssize_t);
        static void __dealloc__(PyObject*);
        
        static PyMethodDef methods[];
        static void setup_tp_methods();
        
        static PyGetSetDef getset[];
        static void setup_tp_getset();
        static void getset_rho(PyGetSetDef&);
        static void getset_c(PyGetSetDef&);
        static void getset_min_alpha(PyGetSetDef&);
        static int setup_tp();
    };
    
    struct LineSearchBackTrack::Object
    {
        PyObject_HEAD
        LineSearchBackTrack ls;
    };
    
    struct LineSearchBackTrack::Object1
    {
        PyObject_HEAD
        LineSearchBackTrack* ls1;
    };
}
/*--------------------------------------------
 bracketing routine
 --------------------------------------------*/
template<class Func>
int LineSearch::bracket(Func* func,type0 dfa,type0 max_a,type0& a,
type0& b,type0& c,type0& fa,type0& fb,type0& fc)
{
    type0 u,fu,r,q,ulim;
    
    q=1.0e-14/(max_a*fabs(dfa));
    if(0.01<q && q<1.0)
        b=q*max_a;
    else if(q>=1.0)
        b=max_a;
    else
        b=0.01*max_a;
    fb=func->F(b);
    
    if(fb>=fa)
    {
        c=b;
        fc=fb;
        int iter=20;
        type0 r=(fa-fc)/(c*dfa);
        b=c*0.5/(1.0+r);
        while(fb>=fa && iter && b>std::numeric_limits<type0>::epsilon())
        {
            fb=func->F(b);
            if(fb<fa)
                continue;
            c=b;
            fc=fb;
            iter--;
            r=(fa-fc)/(c*dfa);
            b=c*0.5/(1.0+r);
        }
        
        if(fb>=fa)
        {
            // last ditch effort
            b=sqrt_epsilon<max_a ? sqrt_epsilon:max_a;
            fb=func->F(b);
            if(fb>=fa)
                return B_F_DOWNHILL;
        }
        
        return B_S;
    }

    fc=fb;
    while(fb>=fc)
    {
        
        c=b+golden*(b-a);
        if(c>=max_a)
        {
            c=max_a;
            fc=func->F(c);
            return B_S;
        }
        
        fc=func->F(c);
        if(fc>fb)
            continue;
        
        ulim=MIN(b+(golden+2.0)*(b-a),max_a);
        ulim=max_a;
        
        r=(b-a)*(fb-fc);
        q=(b-c)*(fb-fa);
        
        u=b-((b-c)*q-(b-a)*r)/(2.0*(q-r));
        if(b<u && u<c)
        {
            fu=func->F(u);
            if(fu<fc)
            {
                a=b;
                b=u;
                fa=fb;
                fb=fu;
                return B_S;
            }
            else if(fu>fb)
            {
                c=u;
                fc=fu;
                return B_S;
            }
            
            a=b;
            b=c;
            
            fa=fb;
            fb=fc;
        }
        else if(u>c)
        {
            
            u=MIN(u,ulim);
            fu=func->F(u);
            a=b;
            b=c;
            c=u;
            
            fa=fb;
            fb=fc;
            fc=fu;
            
            if(fu>fc)
            {
                return B_S;
            }
        }
        else
        {
            a=b;
            b=c;
            
            fa=fb;
            fb=fc;
        }
    }
    
    return B_S;
}
/*--------------------------------------------
 minimize line
 --------------------------------------------*/
template<class Func>
int LineSearchGoldenSection::min(Func* func,type0& nrgy
,type0& alpha,int init_flag)
{
    int calc=0;
    type0 a,b,c;
    type0 fa,fb,fc,h_norm;
    type0 max_a;
    type0 x0,x1,x2,x3,f1,f2;
    type0 dfa;
    constexpr type0 gold=0.61803399;
    constexpr type0 cgold=0.38196601;
    
    func->ls_prep(dfa,h_norm,max_a);
    
    if(dfa==0.0) return LS_F_GRAD0;
    if(max_a==0.0) return LS_MIN_ALPHA;
    if(dfa>=0.0) return LS_F_DOWNHILL;
    
    
    
    a=0.0;
    fa=nrgy;
    
    
    if(brack)
    {
        int chk_bracket=bracket(func,dfa,max_a,a,b,c,fa,fb,fc);
        if(chk_bracket!=B_S)
        {
            func->F_reset();
            return chk_bracket;
        }
        x0=a;
        x3=c;
        
        
        if(c-b>b-a)
        {
            x1=b;
            x2=b+cgold*(c-b);
            
            f1=fb;
            f2=func->F(x2);
            calc=2;
        }
        else
        {
            x1=b+cgold*(a-b);
            x2=b;
            
            f1=func->F(x1);
            f2=fb;
            calc=1;
        }
        
        
        int iter=max_iter;
        while(x3-x0>tol*(x1+x2) && x3-x0>std::numeric_limits<type0>::epsilon() && iter)
        {
            if(f2<f1)
            {
                x0=x1;
                x1=x2;
                x2=gold*x2+cgold*x3;
                
                f1=f2;
                f2=func->F(x2);
                calc=2;
            }
            else
            {
                x3=x2;
                x2=x1;
                x1=gold*x1+cgold*x0;
                f2=f1;
                f1=func->F(x1);
                calc=1;
            }
            iter--;
        }
    }
    else
    {
        type0 delta,f0,f3;
        f1=f2=0.0;
        x0=x1=x2=0.0;
        f0=fa;
        x3=max_a;
        f3=func->F(x3);
        
        int iter=max_iter;
        
        bool set_left=false;
        bool set_right=false;
        while(x3-x0>tol*(x1+x2) && x3-x0>std::numeric_limits<type0>::epsilon() && iter)
        {
            delta=cgold*(x3-x0);
            if(!set_left)
            {
                x1=x0+delta;
                f1=func->F(x1);
                calc=1;
                set_left=1;
            }
            
            if(!set_right)
            {
                x2=x3-delta;
                f2=func->F(x2);
                calc=2;
                set_right=1;
            }
            
            
            if((f0<f2 && f0<f3) || (f1<f2 && f1<f3))
            {
                set_left=0;
                x3=x2;
                x2=x1;
                f3=f2;
                f2=f1;
            }
            else
            {
                set_right=0;
                x0=x1;
                x1=x2;
                f0=f1;
                f1=f2;
            }
            
            iter--;
        }
        
        if(!set_left)
        {
            if(f0<f1)
            {
                if(x0==0.0)
                {
                    func->F_reset();
                    alpha=0.0;
                    return LS_MIN_ALPHA;
                }
                
                x1=x0;
                f1=func->F(x1);
                calc=1;
            }
        }
        if(!set_right)
        {
            if(f3<f2)
            {
                
                x2=x3;
                f2=func->F(x2);
                calc=2;
            }
        }
    }
    
    if(f1<f2)
    {
        if(calc==2)
            f1=func->F(x1);
        nrgy=f1;
        alpha=x1;
    }
    else
    {
        if(calc==1)
            f2=func->F(x2);
        nrgy=f2;
        alpha=x2;
    }

    prev_val=-dfa*alpha;
    return LS_S;
}
/*--------------------------------------------
 minimize line
 --------------------------------------------*/
template<class Func>
int LineSearchBrent::min(Func* func,type0& nrgy
,type0& alpha,int init_flag)
{
    type0 a,b;
    type0 fa;
    type0 max_a,h_norm;
    type0 cgold,gold;
    type0 dfa;

    type0 d=0.0,etemp,fu,fv,fw,fx;
    type0 r,q,p,tol1,tol2,u,v,w,x,xm;
    type0 e=0.0;
    u=fu=0.0;
    
    gold=0.61803399;
    cgold=0.38196601;
    
    func->ls_prep(dfa,h_norm,max_a);
    
    if(dfa==0.0) return LS_F_GRAD0;
    if(max_a==0.0) return LS_MIN_ALPHA;
    if(dfa>=0.0) return LS_F_DOWNHILL;
    
    a=0.0;
    fa=nrgy;
    
    
    type0 fb;
    
    if(brack)
    {
        int chk_bracket=bracket(func,dfa,max_a,a,x,b,fa,fx,fb);
        
        if(chk_bracket!=B_S)
        {
            func->F_reset();
            return chk_bracket;
        }
    }
    else
    {
        a=0.0;
        b=max_a;
        x=cgold*(b-a);
        fx=func->F(x);
    }
    
    w=v=x;
    fw=fv=fx;

    for(int iter=0;iter<max_iter;iter++)
    {
        xm=0.5*(a+b);
        tol1=tol*fabs(x)+zeps;
        tol2=2.0*tol1;
        
        if(fabs(x-xm)<=(tol2-0.5*(b-a)))
        {
            if(u!=x)
                func->F(x);
            nrgy=fx;
            alpha=x;
            prev_val=alpha;
            return LS_S;
        }
        
        if(fabs(e)>tol1)
        {
            r=(x-w)*(fx-fv);
            q=(x-v)*(fx-fw);
            p=(x-v)*q-(x-w)*r;
            q=2.0*(q-r);
            if(q>0.0) p=-p;
            q=fabs(q);
            etemp=e;
            e=d;
            if(fabs(p)>=fabs(0.5*q*etemp)
               || p<=q*(a-x)
               || p>=q*(b-x))
            {
                if(x>=xm)
                    e=(a-x);
                else
                    e=b-x;
                d=cgold*e;
            }
            else
            {
                d=p/q;
                u=x+d;
                if(u-a<tol2 || b-u<tol2)
                {
                    if(xm-x>=0.0)
                        d=fabs(tol1);
                    else
                        d=-fabs(tol1);
                }
            }
        }
        else
        {
            if(x>=xm)
                e=(a-x);
            else
                e=b-x;
            
            d=cgold*e;
        }
        
        if(fabs(d)>=tol1)
        {
            u=x+d;
        }
        else
        {
            if(d>=0.0)
            {
                u=x+fabs(tol1);
            }
            else
            {
                u=x-fabs(tol1);
            }
        }
        
        fu=func->F(u);
        
        if(fu<=fx)
        {
            if(u>=x)
                a=x;
            else
                b=x;
            
            v=w;
            w=x;
            x=u;
            
            fv=fw;
            fw=fx;
            fx=fu;
        }
        else
        {
            if(u<x)
                a=u;
            else
                b=u;
            
            if(fu<=fw
            || w==x)
            {
                v=w;
                w=u;
                fv=fw;
                fw=fu;
            }
            else if(fu<=fv
            || v==x
            || v==w)
            {
                v=u;
                fv=fu;
            }
        }
    }

    //make sure that the result is less than initial value
    if(fa<fx)
    {
        func->F_reset();
        alpha=0.0;
        return LS_MIN_ALPHA;
    }
    
    if(u!=x)
        func->F(x);
    
    nrgy=fx;
    alpha=x;
    
    prev_val=alpha;
    return LS_S;
}
/*--------------------------------------------
 minimize line
 --------------------------------------------*/
template<class Func>
int LineSearchBackTrack::min(Func* func,type0& nrgy
,type0& alpha,int init_flag)
{
    type0 max_a,h_norm;
    type0 dfa,init_energy;
    type0 current_energy,ideal_energy;
    
    init_energy=nrgy;
    
    func->ls_prep(dfa,h_norm,max_a);

    if(dfa==0.0) return LS_F_GRAD0;
    if(max_a==0.0) return LS_MIN_ALPHA;
    if(dfa>=0.0) return LS_F_DOWNHILL;
    
    if(init_flag==0)
    {
        max_a=MIN(max_a,1.0);
    }
    else if(init_flag==1)
    {
        if(prev_val>0.0 && dfa<0.0)
            max_a=MIN(max_a,-prev_val/dfa);
        prev_val=-dfa;
    }
    else if(init_flag==2)
    {
        if(prev_val>0.0 && dfa<0.0)
            max_a=MIN(max_a,MIN(1.0,2.02*(prev_val+nrgy)/dfa));
        prev_val=-nrgy;
    }
    
    if(max_a<=min_alpha)
        return LS_MIN_ALPHA;
    
    
    alpha=max_a;
    while(1)
    {
        ideal_energy=nrgy+alpha*c*dfa;
        current_energy=func->F(alpha);
        if(current_energy<=ideal_energy)
        {
            if(init_flag==1)
                prev_val*=alpha;
            nrgy=current_energy;
            return LS_S;
        }
        alpha*=rho;
        
        if(alpha<=min_alpha)
        {
            nrgy=init_energy;
            func->F_reset();
            return LS_MIN_ALPHA;
        }
    }
}
#endif 
