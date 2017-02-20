#include "example.h"
#include <iostream>
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
        FuncAPI<double[3],double> f("func",{"x","u"});
        f.logics<1>()[0]=VLogics("ge",0.0);
        
        if(f(args,kwds)==-1)
            return NULL;
        Py_RETURN_NONE;
    };
    
    //tp_methods[0].ml_doc="this is a test function that I created";
    tp_methods[0].ml_doc=
R"----(This is test function
This is test function
This is test function
This is test function
This is test function
This is test function
This is test function
    )----";
}
/*--------------------------------------------*/
PyMemberDef ExamplePython::tp_members[]={[0 ... 1]={NULL}};
/*--------------------------------------------*/
void ExamplePython::setup_tp_members()
{
    tp_members[0].name=(char*)"mass";
    tp_members[0].type=T_OBJECT_EX;
    tp_members[0].offset=offsetof(ExamplePython::Object,mass);
    tp_members[0].flags=0;
    tp_members[0].doc=(char*)"defines the masses the simulation";
}
/*--------------------------------------------*/
PyGetSetDef ExamplePython::tp_getset[]={[0 ... 1]={NULL,NULL,NULL,NULL}};
/*--------------------------------------------*/
void ExamplePython::setup_tp_getset()
{
    tp_getset[0].name=(char*)"mass";
    tp_getset[0].get=get_mass;
    tp_getset[0].set=set_mass;
    tp_getset[0].doc=(char*)"defines the masses the simulation";
}
/*--------------------------------------------
 
 --------------------------------------------*/
int ExamplePython::set_mass(PyObject* self,PyObject* op,void*)
{

    VarAPI<size_t> ntypes("ntypes");
    ntypes.val=2;
    
    VarAPI<symm<double[2][2]>> mass("mass",ntypes.val,ntypes.val);
    mass.logics[0]=VLogics("gt",0.0)*VLogics("lt",1.0);

    int i=mass.set(op);
    

    return i;
}
/*--------------------------------------------
 
 --------------------------------------------*/
PyObject* ExamplePython::get_mass(PyObject* self,void* closure)
{
    printf("get mass\n");
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
//    delete self->dt;
    self->ob_type->tp_free((PyObject*)self);
}
/*--------------------------------------------
 
 --------------------------------------------*/
PyObject* ExamplePython::test_func(PyObject* self,PyObject* args,PyObject* kwds)
{
    FuncAPI<double[3],double> f("f",{"x","u"});
    f.logics<1>()[0]=VLogics("ge",0.0);
    
    if(f(args,kwds)==-1)
        return NULL;
    Py_RETURN_NONE;
}
/*--------------------------------------------*/
PyTypeObject ExamplePython::TypeObject ={PyObject_HEAD_INIT(NULL)};
/*--------------------------------------------*/
void ExamplePython::setup_tp()
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
}

