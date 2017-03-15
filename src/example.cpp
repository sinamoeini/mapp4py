#include "example.h"
#include <structmember.h>
using namespace MAPP_NS;

/*--------------------------------------------*/
PyMethodDef ExamplePython::tp_methods[]={[0 ... 1]={NULL,NULL,0,NULL}};
/*--------------------------------------------*/
void ExamplePython::setup_tp_methods()
{
    tp_methods[0].ml_name="func";
    tp_methods[0].ml_flags=METH_VARARGS | METH_KEYWORDS;
    tp_methods[0].ml_meth=(PyCFunction)(PyCFunctionWithKeywords)
    [](PyObject* self,PyObject* args,PyObject* kwds)->PyObject*
    {

        Py_RETURN_NONE;
    };
    
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
PyMemberDef ExamplePython::tp_members[]={[0 ... 1]={NULL}};
/*--------------------------------------------*/
void ExamplePython::setup_tp_members()
{
    
}
/*--------------------------------------------*/
PyGetSetDef ExamplePython::tp_getset[]={[0 ... 1]={NULL,NULL,NULL,NULL}};
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
    self->ob_type->tp_free((PyObject*)self);
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


PyMODINIT_FUNC initxmpl(void)
{
    PyObject* module=Py_InitModule3("xmpl",NULL,"MIT Atomistic Parallel Package");
    ExamplePython::setup_tp();
    if(PyType_Ready(&ExamplePython::TypeObject)<0) return;
    Py_INCREF(&ExamplePython::TypeObject);
    PyModule_AddObject(module,"obj",reinterpret_cast<PyObject*>(&ExamplePython::TypeObject));
}
/*--------------------------------------------
 
 --------------------------------------------*/
#include <iostream>
#include <frameobject.h>
#include <pyerrors.h>
void ExamplePython::ml_test(PyMethodDef& tp_methods)
{
    tp_methods.ml_flags=METH_VARARGS | METH_KEYWORDS;
    tp_methods.ml_name="test";
    tp_methods.ml_doc="run simulation for n steps";
    
    tp_methods.ml_meth=(PyCFunction)(PyCFunctionWithKeywords)
    [](PyObject* self,PyObject* args,PyObject* kwds)->PyObject*
    {
        
        PyObject* op;
        if(PyArg_ParseTuple(args, "O:set_callback", &op))
        {
         
            if (!PyFunction_Check(op))
            {
                PyErr_SetString(PyExc_TypeError, "parameter must be callable");
                return NULL;
            }
        }
        else
            return NULL;
        
        
        
        
        
        
        
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
        Py_RETURN_NONE;
        
        
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
        
    };
}



