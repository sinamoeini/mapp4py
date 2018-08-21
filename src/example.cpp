#include "api.h"
#include "example.h"
#include <structmember.h>
#include "global.h"
using namespace MAPP_NS;

/*--------------------------------------------*/
PyMethodDef ExamplePython::tp_methods[]=EmptyPyMethodDef(2);
/*--------------------------------------------*/
void ExamplePython::setup_tp_methods()
{
    tp_methods[0].ml_name="func";
    tp_methods[0].ml_flags=METH_VARARGS | METH_KEYWORDS;
    tp_methods[0].ml_meth=(PyCFunction)(PyCFunctionWithKeywords)(
    [](PyObject* self,PyObject* args,PyObject* kwds)->PyObject*
    {

        Py_RETURN_NONE;
    });
    
    //tp_methods[0].ml_doc="this is a test function that I created";
    tp_methods[0].ml_doc=(char *)
    R"---(
func(eps)

A one-line summary that does not use variable names or the function name.

Parameters
----------
eps : array_like
    hi
)---";
}
/*--------------------------------------------*/
PyMemberDef ExamplePython::tp_members[]=EmptyPyMemberDef(2);
/*--------------------------------------------*/
void ExamplePython::setup_tp_members()
{
    
}
/*--------------------------------------------*/
PyGetSetDef ExamplePython::tp_getset[]=EmptyPyGetSetDef(2);
/*--------------------------------------------*/
void ExamplePython::setup_tp_getset()
{/*
    tp_getset[0].name=(char*)"mass";
    tp_getset[0].get=get_mass;
    tp_getset[0].set=set_mass;
    tp_getset[0].doc=(char*)"   defines the masses the simulation";*/
}
/*--------------------------------------------
 
 --------------------------------------------*/
int ExamplePython::set_mass(PyObject* self,PyObject* op,void*)
{
    return 0;
}
/*--------------------------------------------
 
 --------------------------------------------*/
PyObject* ExamplePython::get_mass(PyObject* self,void* closure)
{
    return NULL;
}
/*--------------------------------------------
 
 --------------------------------------------*/
PyObject* ExamplePython::tp_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
    Object *self;
    self = (Object *)type->tp_alloc(type,0);
    return (PyObject*)self;
}
/*--------------------------------------------
 
 --------------------------------------------*/
int ExamplePython::tp_init(PyObject *self_, PyObject *args, PyObject *kwds)
{
//    Object* self=(Object*) self_;
//    delete self->dt;
//    self->dt=new VarAPI(self->__dt,"dt");
//    delete self->mass;
//    delete [] self->__mass;
//    self->__mass=NULL;
//    self->mass=new VarAPI(self->__mass,"mass");
    return 0;
}
/*--------------------------------------------
 
 --------------------------------------------*/
void ExamplePython::tp_dealloc(PyObject* self_)
{
    Object* self=(Object*)self_;
    Py_TYPE(self)->tp_free((PyObject*)self);
}
/*--------------------------------------------
 
 --------------------------------------------*/
PyObject* ExamplePython::test_func(PyObject* self,PyObject* args,PyObject* kwds)
{
    Py_RETURN_NONE;
}
/*--------------------------------------------*/
PyTypeObject ExamplePython::TypeObject ={PyObject_HEAD_INIT(NULL)};
/*--------------------------------------------*/
int ExamplePython::setup_tp()
{
    TypeObject.tp_name="xmpl.obj";
    TypeObject.tp_doc="I will add doc here";
    
    TypeObject.tp_flags=Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE;
    TypeObject.tp_basicsize=sizeof(Object);
    
    TypeObject.tp_new=tp_new;
    TypeObject.tp_init=tp_init;
    TypeObject.tp_dealloc=tp_dealloc;
    
    setup_tp_members();
    TypeObject.tp_members=tp_members;
    setup_tp_getset();
    TypeObject.tp_getset=tp_getset;
    setup_tp_methods();
    TypeObject.tp_methods=tp_methods;
    
    int ichk=PyType_Ready(&TypeObject);
    if(ichk<0) return ichk;
    Py_INCREF(&TypeObject);
    return ichk;
}

/*--------------------------------------------*/
#ifdef IS_PY3K
PyModuleDef ExamplePython::module=EmptyModule;
#endif
MOD_INIT(xmpl,MAPP_NS::ExamplePython::init_module())
/*--------------------------------------------*/
PyObject* ExamplePython::init_module(void)
{
#ifdef IS_PY3K
    ExamplePython::module.m_name="xmpl";
    ExamplePython::module.m_doc="MIT Atomistic Parallel Package";
    ExamplePython::module.m_methods=NULL;
    PyObject* module_ob=PyModule_Create(&ExamplePython::module);
#else
    PyObject* module_ob=Py_InitModule3("xmpl",NULL,"MIT Atomistic Parallel Package");
#endif
    ExamplePython::setup_tp();
    if(PyType_Ready(&ExamplePython::TypeObject)<0) return NULL;
    Py_INCREF(&ExamplePython::TypeObject);
    PyModule_AddObject(module_ob,"obj",reinterpret_cast<PyObject*>(&ExamplePython::TypeObject));
    return module_ob;
}
/*--------------------------------------------
 
 --------------------------------------------*/
#include <iostream>
#include <frameobject.h>
#include <pyerrors.h>
#include "xmath.h"
#include "atoms_styles.h"

void ExamplePython::ml_test(PyMethodDef& tp_methods)
{
    tp_methods.ml_flags=METH_VARARGS | METH_KEYWORDS;
    tp_methods.ml_name="test";
    tp_methods.ml_doc="run simulation for n steps";
    
    tp_methods.ml_meth=(PyCFunction)(PyCFunctionWithKeywords)(
    [](PyObject* self,PyObject* args,PyObject* kwds)->PyObject*
    {
        FuncAPI<OP<AtomsDMD>>f("test",{"atoms"});
        if(f(args,kwds)) return NULL;
        AtomsDMD* atoms=reinterpret_cast<AtomsDMD::Object*>(f.val<0>().ob)->atoms;
        int nelems=static_cast<int>(atoms->elements.nelems);
        type0* alpha_ave_lcl=new type0[nelems];
        type0* alpha_ave=new type0[nelems];
        type0* n_lcl=new type0[nelems];
        type0* n=new type0[nelems];
        for(size_t i=0;i<nelems;i++) alpha_ave_lcl[i]=n_lcl[i]=0.0;
        type0* c=atoms->c->begin();
        type0* alpha=atoms->alpha->begin();
        elem_type* elem=atoms->elem->begin();
        int natms_lcl=atoms->natms_lcl;
        int c_dim=atoms->c_dim;
        for(int i=0;i<natms_lcl;i++)
        {
            for(int j=0;j<c_dim;j++)
            {
                if(c[j]>=0.0)
                {
                    alpha_ave_lcl[elem[j]]+=alpha[j];
                    n_lcl[elem[j]]++;
                }
                    
            }
            c+=c_dim;
            alpha+=c_dim;
        }
        
        MPI_Allreduce(alpha_ave_lcl,alpha_ave,nelems,Vec<type0>::MPI_T,MPI_SUM,atoms->world);
        MPI_Allreduce(n_lcl,n,nelems,Vec<type0>::MPI_T,MPI_SUM,atoms->world);

        for(int i=0;i<nelems;i++)
            alpha_ave[i]/=n[i];
    
        size_t* __nelems=&(atoms->elements.nelems);
        PyObject* op=var<type0*>::build(alpha_ave,&__nelems);
        delete [] n;
        delete [] n_lcl;
        delete [] alpha_ave;
        delete [] alpha_ave_lcl;
        
        return op;
        
        /*
        FuncAPI<int> f("run",{"N"});
        f.logics<0>()[0]=VLogics("gt",0);
        if(f(args,kwds)) return NULL;
        
        int n=f.val<0>();
        
        type0* xi=new type0[n];
        type0* wi=new type0[n];
        
        XMath::quadrature_hg(n,xi,wi);
        
        
        for(int i=0;i<n;i++)
            printf("%d\t\t%e\t\t%e\n",i,xi[i],wi[i]);
        
        delete [] xi;
        delete [] wi;
        
        */
        
        
        //PyFunctionObject* fp=(PyFunctionObject*) PyTuple_GetItem(args,0);
        
        
        
        /*

        PyObject* temp=PyTuple_GetItem(args,0);
        
        if (!PyCallable_Check(temp))
        {
            PyErr_SetString(PyExc_TypeError, "parameter must be callable");
            
            return NULL;
        }
        Py_XINCREF(temp);
        
        PyObject* x=PyString_FromString("x");
        PyObject* __x=PyString_FromString("xooo");
        PyObject* y=PyString_FromString("y");
        PyObject* __y=PyFloat_FromDouble(2.0);
        PyObject* z=PyString_FromString("z");
        PyObject* __z=PyFloat_FromDouble(0.6);
        
        PyObject* dict=PyDict_New();
        PyDict_SetItem(dict,x,__x);
        PyDict_SetItem(dict,y,__y);
        PyDict_SetItem(dict,z,__z);
        
        
        
        PyCodeObject* co=(PyCodeObject *)PyFunction_GET_CODE(temp);
        PyObject* co_varnames=co->co_varnames;
        size_t sz=PyTuple_Size(co_varnames);
        
        for(size_t i=0;i<sz;i++)
        {
            char* str=PyString_AsString(PyTuple_GetItem(co_varnames,i));
            printf("%s\n",str);
        }
        
        
        char* str=PyString_AsString(co->co_name);
        printf("%s\n",str);
        
    
        
        PyObject* ptype;
        PyObject* pvalue;
        PyObject* ptraceback;
        PyErr_Fetch(&ptype,&pvalue,&ptraceback);
        
        if ( ptype != NULL )
        {
            PyObject* pRepr = PyObject_Repr( ptype ) ;
            std::cout << "- EXC type: " << PyString_AsString(pRepr) << std::endl ;
            Py_DECREF( pRepr ) ;
            Py_DECREF( ptype ) ;
        }
        if ( pvalue != NULL )
        {
            PyObject* pRepr = PyObject_Repr( pvalue ) ;
            std::cout << "- EXC value: " << PyString_AsString(pRepr) << std::endl ;
            Py_DECREF( pRepr ) ;
            Py_DECREF(pvalue) ;
        }
        if ( ptraceback != NULL )
        {
            PyObject* pRepr = PyObject_Repr( pvalue ) ;
            std::cout << "- EXC traceback: " << PyString_AsString(pRepr) << std::endl ;
            Py_DECREF( pRepr ) ;
            Py_DECREF( ptraceback ) ;
        }
        
         str=PyString_AsString(pvalue);
         printf("%s\n",str);
        
        
        
        PyErr_Clear();
        */
        //Py_RETURN_NONE;
        
        
        /*
         keywords: 
            x         [__dim__]
            x_d       [__dim__]
            dof       [__dim__]
            elem
            id        <readonly>

         
         keywords:
            x         [__dim__]
            dof       [__dim__]
            alpha     [];
            dof_alpha [];
            c         [];
            c_dof     [];
            elem      [];
            id        <readonly>;
         
         */
        
    });
}

/*--------------------------------------------
 
 --------------------------------------------*/
#include "dynamic_md.h"
#include "ff_styles.h"
void ExamplePython::ml_phonon(PyMethodDef& tp_methods)
{
    tp_methods.ml_flags=METH_VARARGS | METH_KEYWORDS;
    tp_methods.ml_name="phonon";
    tp_methods.ml_meth=(PyCFunction)(PyCFunctionWithKeywords)(
    [](PyObject* self,PyObject* args,PyObject* kwds)->PyObject*
    {
        /*

        FuncAPI<OP<AtomsMD>,type0,int> f("phonon",{"atoms","max_disp","N"});
        f.logics<1>()[0]=VLogics("gt",0.0);
        f.logics<2>()[0]=VLogics("gt",0);
        if(f(args,kwds)) return NULL;
        
    
        type0 disp=-f.v<1>();
        int n=f.v<2>();
        type0 delta=2.0*f.v<1>()/n;
        
        
        AtomsMD* atoms=reinterpret_cast<AtomsMD::Object*>(f.val<0>().ob)->atoms;
        ForceFieldMD* ff=reinterpret_cast<AtomsMD::Object*>(f.val<0>().ob)->ff;
        
        
        DynamicMD* dynamic=new DynamicMD(atoms,ff,false,{},{},{});
        dynamic->init();
        int rank=atoms->comm_rank;
        type0 x0=0.0;
        if(!rank)
            x0=atoms->x->begin()[0];
    
        type0 sum=0.0;
        type0 no=0.0;
        if(!rank)
            printf("----------------------------------------------------\n");
        
        for(int i=0;i<n+1;i++)
        {
            if(!rank)
                atoms->x->begin()[0]=x0+disp;
            dynamic->update(atoms->x);
            
            ff->derivative_timer();
            if(!rank)
            {
                if(fabs(disp)>1.0e-8)
                {
                    printf("%.12lf\t%.12lf\n",disp,-ff->f->begin()[0]);
                    sum+=ff->f->begin()[0]/disp;
                    no++;
                }
                disp+=delta;
            }
            
        }
        
        if(!rank)
            printf("----------------------------------------------------\n");
        
        if(!rank)
            atoms->x->begin()[0]=x0;
        dynamic->update(atoms->x);
        
        dynamic->fin();
        delete dynamic;
        
        sum/=no;
        MPI_Bcast(&sum,1,Vec<type0>::MPI_T,0,atoms->world);
        
        PyObject* op=var<type0>::build(sum);
        
        return op;*/
        
        
        FuncAPI<OP<AtomsMD>,type0> f("phonon",{"atoms","max_disp"});
        f.logics<1>()[0]=VLogics("gt",0.0);
        if(f(args,kwds)) return NULL;
        
        
        type0 disp=f.v<1>();

        
        
        AtomsMD* atoms=reinterpret_cast<AtomsMD::Object*>(f.val<0>().ob)->atoms;
        ForceFieldMD* ff=reinterpret_cast<AtomsMD::Object*>(f.val<0>().ob)->ff;
        
        
        DynamicMD* dynamic=new DynamicMD(atoms,ff,false,{},{},{});
        dynamic->init();
        
        
        Vec<type0> v0(atoms,__dim__);
        Vec<type0> v1(atoms,__dim__);
        Vec<type0> v2(atoms,__dim__);
        type0 x0;
        x0=atoms->x->begin()[0];
        atoms->x->begin()[0]+=disp;
        dynamic->update(atoms->x);
        ff->derivative();
        memcpy(v0.begin(),ff->f->begin(),atoms->natms_lcl*__dim__*sizeof(type0));
        atoms->x->begin()[0]=x0;
        
        
        x0=atoms->x->begin()[1];
        atoms->x->begin()[1]+=disp;
        dynamic->update(atoms->x);
        ff->derivative();
        memcpy(v1.begin(),ff->f->begin(),atoms->natms_lcl*__dim__*sizeof(type0));
        atoms->x->begin()[1]=x0;
        
        
        x0=atoms->x->begin()[2];
        atoms->x->begin()[2]+=disp;
        dynamic->update(atoms->x);
        ff->derivative();
        memcpy(v2.begin(),ff->f->begin(),atoms->natms_lcl*__dim__*sizeof(type0));
        atoms->x->begin()[2]=x0;
        
        
        
        
        atoms->x2s_lcl();
        
        type0 B[__dim__][__dim__]={{-0.5,0.5,0.5},{0.5,-0.5,0.5},{0.5,0.5,-0.5}};
        type0 s[__dim__];
        
        
        printf("----------------------------------------------------\n");
        
        type0* __s=atoms->x->begin();
        type0* __v0=v0.begin();
        type0* __v1=v1.begin();
        type0* __v2=v2.begin();
        type0 min[__dim__]={0,0,0};
        type0 max[__dim__]={0,0,0};
        
        for(int i=0;i<atoms->natms_lcl;i++)
        {
            for(int j=0;j<__dim__;j++)
            {
                __v0[j]/=-disp;
                __v1[j]/=-disp;
                __v2[j]/=-disp;
                
                s[j]=0.0;
                for(int k=0;k<__dim__;k++)
                    s[j]+=(__s[k]-0.5)*B[k][j]*8.0;
                
                min[j]=MIN(min[j],s[j]);
                max[j]=MAX(max[j],s[j]);
            }
            
            printf("Exp[(%.0lf*k[[1]]+%.0lf*k[[2]]+%.0lf*k[[3]])*2*Pi*I]*{{%0.7lf,%0.7lf,%0.7lf},{%0.7lf,%0.7lf,%0.7lf},{%0.7lf,%0.7lf,%0.7lf}}",s[0],s[1],s[2],
                   __v0[0],__v0[1],__v0[2],
                   __v1[0],__v1[1],__v1[2],
                   __v2[0],__v2[1],__v2[2]);


            __v0+=__dim__;
            __v1+=__dim__;
            __v2+=__dim__;
            __s+=__dim__;
            
            if(i!=atoms->natms_lcl-1)
                printf("+\n");
            else
                printf(";\n");
            
        }
        
        printf("----------------------------------------------------\n");
        printf("min %.0lf %.0lf %.0lf\n",min[0],min[1],min[2]);
        printf("max %.0lf %.0lf %.0lf\n",max[0],max[1],max[2]);
        printf("----------------------------------------------------\n");
        dynamic->update(atoms->x);
        
        dynamic->fin();
        delete dynamic;
        Py_RETURN_NONE;
    });


    tp_methods.ml_doc=(char*)R"---(
    quick function for calculateing phonon freq
    use with caution

    )---";
}
/*--------------------------------------------
 
 --------------------------------------------*/
void ExamplePython::ml_phonon_1d(PyMethodDef& tp_methods)
{
    tp_methods.ml_flags=METH_VARARGS | METH_KEYWORDS;
    tp_methods.ml_name="phonon_1d";
    tp_methods.ml_meth=(PyCFunction)(PyCFunctionWithKeywords)(
    [](PyObject* self,PyObject* args,PyObject* kwds)->PyObject*
    {
        FuncAPI<OP<AtomsMD>,type0> f("phonon",{"atoms","max_disp"});
        f.logics<1>()[0]=VLogics("ge",0.0);
        if(f(args,kwds)) return NULL;
        
        
        type0 disp=f.v<1>();

        
        
        AtomsMD* atoms=reinterpret_cast<AtomsMD::Object*>(f.val<0>().ob)->atoms;
        ForceFieldMD* ff=reinterpret_cast<AtomsMD::Object*>(f.val<0>().ob)->ff;
        
        
        DynamicMD* dynamic=new DynamicMD(atoms,ff,false,{},{},{});
        dynamic->init();
        
        
        type0 x0;
        int N=5;
        
        x0=atoms->x->begin()[N*3];
        atoms->x->begin()[N*3]+=disp;
        dynamic->update(atoms->x);
        ff->derivative();
        atoms->x->begin()[1]=x0;
        
        for(int i=0;i<N+1;i++)
            printf("%lf\n",-ff->f->begin()[i*3]/disp);

        dynamic->update(atoms->x);
        
        dynamic->fin();
        delete dynamic;
        Py_RETURN_NONE;
    });


    tp_methods.ml_doc=(char*)R"---(
    quick function for calculateing phonon freq
    use with caution

    )---";
}
/*--------------------------------------------
 
 --------------------------------------------*/
void ExamplePython::ml_phonon_1dd(PyMethodDef& tp_methods)
{
    tp_methods.ml_flags=METH_VARARGS | METH_KEYWORDS;
    tp_methods.ml_name="phonon_1dd";
    tp_methods.ml_meth=(PyCFunction)(PyCFunctionWithKeywords)(
    [](PyObject* self,PyObject* args,PyObject* kwds)->PyObject*
    {
        FuncAPI<OP<AtomsMD>,type0> f("phonon",{"atoms","max_disp"});
        f.logics<1>()[0]=VLogics("ge",0.0);
        if(f(args,kwds)) return NULL;
        
        
        type0 disp=f.v<1>();
        
        
        
        AtomsMD* atoms=reinterpret_cast<AtomsMD::Object*>(f.val<0>().ob)->atoms;
        ForceFieldMD* ff=reinterpret_cast<AtomsMD::Object*>(f.val<0>().ob)->ff;
        
        
        DynamicMD* dynamic=new DynamicMD(atoms,ff,false,{},{},{});
        dynamic->init();
        
        
        type0 x0;

        
        for(int i=0;i<5;i++)
        {
            x0=atoms->x->begin()[i*3];
            atoms->x->begin()[i*3]+=disp;
            dynamic->update(atoms->x);
            ff->derivative();
            atoms->x->begin()[i*3]=x0;
            
            for(int j=0;j<10;j++)
            {
                printf("J[[%d,%d]]=%.10lf;\n",i+1,j+1,-ff->f->begin()[j*3]/disp);
                printf("J[[%d,%d]]=%.10lf;\n",11-(i+1),11-(j+1),-ff->f->begin()[j*3]/disp);
            }
        }
        
        
        
        dynamic->update(atoms->x);
        
        dynamic->fin();
        delete dynamic;
        Py_RETURN_NONE;
    });
    
    
    tp_methods.ml_doc=(char*)R"---(
    quick function for calculateing phonon freq
    use with caution
    
    )---";
}
/*--------------------------------------------
 
 --------------------------------------------*/
#include "dynamic_dmd.h"
void ExamplePython::ml_alpha(PyMethodDef& tp_methods)
{
    tp_methods.ml_flags=METH_VARARGS | METH_KEYWORDS;
    tp_methods.ml_name="alpha";
    tp_methods.ml_meth=(PyCFunction)(PyCFunctionWithKeywords)(
    [](PyObject* self,PyObject* args,PyObject* kwds)->PyObject*
    {
        //Object* __self=reinterpret_cast<Object*>(self);
        FuncAPI<OP<AtomsDMD>,type0,type0,int> f("alpha",{"atoms","dalpha0","dalpha1","N"});
        f.logics<1>()[0]=VLogics("lt",0.0);
        f.logics<2>()[0]=VLogics("gt",0.0);
        f.logics<3>()[0]=VLogics("gt",0);
        if(f(args,kwds)) return NULL;
        
    
        type0 disp=f.v<1>();
        int n=f.v<3>();
        type0 delta=(f.v<2>()-disp)/n;
        
        
        AtomsDMD* atoms=reinterpret_cast<AtomsDMD::Object*>(f.val<0>().ob)->atoms;
        ForceFieldDMD* ff=reinterpret_cast<AtomsDMD::Object*>(f.val<0>().ob)->ff;
        
        
        Vec<type0> alpha0(atoms,atoms->c_dim);
        DynamicDMD* dynamic=new DynamicDMD(atoms,ff,false,{},{&alpha0},{});
        dynamic->init();
        int rank=atoms->comm_rank;

        memcpy(alpha0.begin(),atoms->alpha->begin(),atoms->natms_lcl*atoms->c_dim*sizeof(type0));
#ifdef NEW_UPDATE
#else
        vec* uvecs[2];
        uvecs[0]=atoms->x;
        uvecs[1]=atoms->alpha;
#endif
        
        
        if(!rank)
            printf("----------------------------------------------------\n");
        type0 en;
        for(int i=0;i<n+1;i++)
        {
            for(int j=0;j<atoms->natms_lcl;j++) atoms->alpha->begin()[j]=alpha0.begin()[j]+disp;
#ifdef NEW_UPDATE
            dynamic->update<true,true>();
#else
            dynamic->update(uvecs,2);
#endif
            en=ff->value();
            if(!rank)
            {
                
                printf("%.12lf\t%.12lf\n",disp,en);
                    
            
                disp+=delta;
            }
            
        }
        
        if(!rank)
            printf("----------------------------------------------------\n");
        
        memcpy(atoms->alpha->begin(),alpha0.begin(),atoms->natms_lcl*atoms->c_dim*sizeof(type0));
#ifdef NEW_UPDATE
        dynamic->update<true,true>();
#else
        dynamic->update(uvecs,2);
#endif
        
        dynamic->fin();
        delete dynamic;
        
        Py_RETURN_NONE;
    });


    tp_methods.ml_doc=(char*)R"---(
    quick function for calculateing phonon freq
    use with caution

    )---";
}

/*--------------------------------------------
 
 --------------------------------------------*/
void ExamplePython::ml_delta_c(PyMethodDef& tp_methods)
{
    tp_methods.ml_flags=METH_VARARGS | METH_KEYWORDS;
    tp_methods.ml_name="delta_c";
    tp_methods.ml_meth=(PyCFunction)(PyCFunctionWithKeywords)(
    [](PyObject* self,PyObject* args,PyObject* kwds)->PyObject*
    {
        //Object* __self=reinterpret_cast<Object*>(self);
        FuncAPI<OP<AtomsDMD>,type0> f("delta_c",{"atoms","dc"});
        if(f(args,kwds)) return NULL;
        
    
        AtomsDMD* atoms=reinterpret_cast<AtomsDMD::Object*>(f.val<0>().ob)->atoms;
        type0 disp=f.v<1>();
        type0* c=atoms->c->begin();
        c[0]+=disp;
        
        Py_RETURN_NONE;
    });


    tp_methods.ml_doc=(char*)R"---(
    quick function for calculateing phonon freq
    use with caution

    )---";
}
/*--------------------------------------------
 
 --------------------------------------------*/
void ExamplePython::ml_mv_c(PyMethodDef& tp_methods)
{
    tp_methods.ml_flags=METH_VARARGS | METH_KEYWORDS;
    tp_methods.ml_name="mv_c";
    tp_methods.ml_meth=(PyCFunction)(PyCFunctionWithKeywords)(
    [](PyObject* self,PyObject* args,PyObject* kwds)->PyObject*
    {
        //Object* __self=reinterpret_cast<Object*>(self);
        FuncAPI<OP<AtomsDMD>,type0> f("mv_c",{"atoms","dc"});
        if(f(args,kwds)) return NULL;
        
    
        AtomsDMD* atoms=reinterpret_cast<AtomsDMD::Object*>(f.val<0>().ob)->atoms;
        type0 x=f.v<1>();
        type0 __x=1.0-x;
        type0* c=atoms->c->begin();
        int natms_lcl=atoms->natms_lcl;
        for(int i=0;i<natms_lcl;i++,c+=2)
        {
            c[0]=x;
            c[1]=__x;
        }
        
        Py_RETURN_NONE;
    });


    tp_methods.ml_doc=(char*)R"---(
    quick function for calculateing phonon freq
    use with caution

    )---";
}
/*--------------------------------------------
 
 --------------------------------------------*/
void ExamplePython::ml_mv_c2(PyMethodDef& tp_methods)
{
    tp_methods.ml_flags=METH_VARARGS | METH_KEYWORDS;
    tp_methods.ml_name="mv_c2";
    tp_methods.ml_meth=(PyCFunction)(PyCFunctionWithKeywords)(
    [](PyObject* self,PyObject* args,PyObject* kwds)->PyObject*
    {
        //Object* __self=reinterpret_cast<Object*>(self);
        FuncAPI<OP<AtomsDMD>,type0> f("mv_c2",{"atoms","dc"});
        if(f(args,kwds)) return NULL;
        
    
        AtomsDMD* atoms=reinterpret_cast<AtomsDMD::Object*>(f.val<0>().ob)->atoms;
        type0 x=f.v<1>();
        
        
        if(x>0.25)
        {
            
            type0 __ni=(1.0-x)*4.0/3.0;
            type0 __al=1.0-__ni;
            
            type0* c=atoms->c->begin();
            int natms_lcl=atoms->natms_lcl;
            c[0]=0.0;
            c[1]=1.0;
            c+=2;
            
            for(int i=1;i<natms_lcl;i++,c+=2)
            {
                c[0]=__ni;
                c[1]=__al;
            }
        }
        else
        {
            type0* c=atoms->c->begin();
            int natms_lcl=atoms->natms_lcl;
            c[0]=1.0-4.0*x;
            c[1]=4.0*x;
            c+=2;
            
            for(int i=1;i<natms_lcl;i++,c+=2)
            {
                c[0]=1.0;
                c[1]=0.0;
            }
        }
        
        Py_RETURN_NONE;
    });


    tp_methods.ml_doc=(char*)R"---(
    quick function for calculateing phonon freq
    use with caution

    )---";
}
/*--------------------------------------------
 
 --------------------------------------------*/
void ExamplePython::ml_mv_c3(PyMethodDef& tp_methods)
{
    tp_methods.ml_flags=METH_VARARGS | METH_KEYWORDS;
    tp_methods.ml_name="mv_c3";
    tp_methods.ml_meth=(PyCFunction)(PyCFunctionWithKeywords)(
    [](PyObject* self,PyObject* args,PyObject* kwds)->PyObject*
    {
        //Object* __self=reinterpret_cast<Object*>(self);
        FuncAPI<OP<AtomsDMD>,type0,type0> f("mv_c3",{"atoms","x","a"});
        if(f(args,kwds)) return NULL;
        
    
        AtomsDMD* atoms=reinterpret_cast<AtomsDMD::Object*>(f.val<0>().ob)->atoms;
        type0 x=f.v<1>();
        type0 a=f.v<2>();
        type0 b=(4.0*x-a)/3.0;
        
        //printf(">>>>>>>>>>>>>>>>>>>> %e %e\n",a,b);
        type0* c=atoms->c->begin();
        int natms_lcl=atoms->natms_lcl;
        c[0]=a;
        c[1]=1.0-a;
        c+=2;
        
        for(int i=1;i<natms_lcl;i++,c+=2)
        {
            c[0]=b;
            c[1]=1.0-b;
        }
        
        
        Py_RETURN_NONE;
    });


    tp_methods.ml_doc=(char*)R"---(
    quick function for calculateing phonon freq
    use with caution

    )---";
}
/*--------------------------------------------
 
 --------------------------------------------*/
void ExamplePython::ml_prt(PyMethodDef& tp_methods)
{
    tp_methods.ml_flags=METH_VARARGS | METH_KEYWORDS;
    tp_methods.ml_name="prt";
    tp_methods.ml_meth=(PyCFunction)(PyCFunctionWithKeywords)(
    [](PyObject* self,PyObject* args,PyObject* kwds)->PyObject*
    {
        //Object* __self=reinterpret_cast<Object*>(self);
        FuncAPI<OP<AtomsDMD>> f("prt",{"atoms"});
        if(f(args,kwds)) return NULL;
        
        
        AtomsDMD* atoms=reinterpret_cast<AtomsDMD::Object*>(f.val<0>().ob)->atoms;
        type0* c=atoms->c->begin();
        type0* alpha=atoms->alpha->begin();
        printf("%.30lf\t%.18lf\n",c[0],alpha[0]);
        
        Py_RETURN_NONE;
    });
    
    
    tp_methods.ml_doc=(char*)R"---(
    quick function for calculateing phonon freq
    use with caution
    
    )---";
}



