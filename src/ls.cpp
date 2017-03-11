#include "min_cg.h"
#include "min_cg_dmd.h"
#include "ls.h"
/*--------------------------------------------
 constructor
 --------------------------------------------*/
LineSearch::LineSearch():
golden(0.5+0.5*sqrt(5.0)),
epsilon_3_4(pow(std::numeric_limits<type0>::epsilon(),0.75)),
sqrt_epsilon(sqrt(std::numeric_limits<type0>::epsilon())),
prev_val(-1.0)
{    
}
/*--------------------------------------------
 destructor
 --------------------------------------------*/
//LineSearch::~LineSearch(){}
/*------------------------------------------------------------------------------------------------------------------------------------
 
 ------------------------------------------------------------------------------------------------------------------------------------*/
PyObject* LineSearch::__new__(PyTypeObject* type,PyObject* args,PyObject* kwds)
{
    Object* __self=reinterpret_cast<Object*>(type->tp_alloc(type,0));
    PyObject* self=reinterpret_cast<PyObject*>(__self);
    return self;
}
/*--------------------------------------------
 
 --------------------------------------------*/
int LineSearch::__init__(PyObject* self,PyObject* args,PyObject* kwds)
{
    FuncAPI<> f("__init__");
    
    if(f(args,kwds)==-1) return -1;
    return 0;
}
/*--------------------------------------------
 
 --------------------------------------------*/
PyObject* LineSearch::__alloc__(PyTypeObject* type,Py_ssize_t)
{
    Object* __self=new Object;
    __self->ob_type=type;
    __self->ob_refcnt=1;
    return reinterpret_cast<PyObject*>(__self);
}
/*--------------------------------------------
 
 --------------------------------------------*/
void LineSearch::__dealloc__(PyObject* self)
{
    Object* __self=reinterpret_cast<Object*>(self);
    delete __self;
}
/*--------------------------------------------*/
PyTypeObject LineSearch::TypeObject={PyObject_HEAD_INIT(NULL)};
/*--------------------------------------------*/
void LineSearch::setup_tp()
{
    TypeObject.tp_name="ls";
    TypeObject.tp_doc="line search";
    
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
}
/*--------------------------------------------*/
PyGetSetDef LineSearch::getset[]={[0 ... 0]={NULL,NULL,NULL,NULL,NULL}};
/*--------------------------------------------*/
void LineSearch::setup_tp_getset()
{
}
/*--------------------------------------------*/
PyMethodDef LineSearch::methods[]={[0 ... 0]={NULL}};
/*--------------------------------------------*/
void LineSearch::setup_tp_methods()
{
}
/*--------------------------------------------
 constructor
 --------------------------------------------*/
LineSearchGoldenSection::LineSearchGoldenSection():
LineSearch()
{
    tol=sqrt(2.0*std::numeric_limits<type0>::epsilon());
    max_iter=5;
    brack=true;
}
/*--------------------------------------------
 destructor
 --------------------------------------------*/
//LineSearchGoldenSection::~LineSearchGoldenSection(){}
/*------------------------------------------------------------------------------------------------------------------------------------
 
 ------------------------------------------------------------------------------------------------------------------------------------*/
PyObject* LineSearchGoldenSection::__new__(PyTypeObject* type,PyObject* args,PyObject* kwds)
{
    Object* __self=reinterpret_cast<Object*>(type->tp_alloc(type,0));
    PyObject* self=reinterpret_cast<PyObject*>(__self);
    return self;
}
/*--------------------------------------------
 
 --------------------------------------------*/
int LineSearchGoldenSection::__init__(PyObject* self,PyObject* args,PyObject* kwds)
{
    FuncAPI<> f("__init__");
    
    if(f(args,kwds)==-1) return -1;
    return 0;
}
/*--------------------------------------------
 
 --------------------------------------------*/
PyObject* LineSearchGoldenSection::__alloc__(PyTypeObject* type,Py_ssize_t)
{
    Object* __self=new Object;
    __self->ob_type=type;
    __self->ob_refcnt=1;
    return reinterpret_cast<PyObject*>(__self);
}
/*--------------------------------------------
 
 --------------------------------------------*/
void LineSearchGoldenSection::__dealloc__(PyObject* self)
{
    Object* __self=reinterpret_cast<Object*>(self);
    delete __self;
}
/*--------------------------------------------*/
PyTypeObject LineSearchGoldenSection::TypeObject={PyObject_HEAD_INIT(NULL)};
/*--------------------------------------------*/
void LineSearchGoldenSection::setup_tp()
{
    TypeObject.tp_name="ls_golden";
    TypeObject.tp_doc="line search golden section";
    
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
    
    TypeObject.tp_base=&LineSearch::TypeObject;
}
/*--------------------------------------------*/
PyGetSetDef LineSearchGoldenSection::getset[]={[0 ... 3]={NULL,NULL,NULL,NULL,NULL}};
/*--------------------------------------------*/
void LineSearchGoldenSection::setup_tp_getset()
{
    getset_bracket(getset[0]);
    getset_max_iter(getset[1]);
    getset_tol(getset[2]);
}
/*--------------------------------------------*/
PyMethodDef LineSearchGoldenSection::methods[]={[0 ... 0]={NULL}};
/*--------------------------------------------*/
void LineSearchGoldenSection::setup_tp_methods()
{
}
/*--------------------------------------------
 
 --------------------------------------------*/
void LineSearchGoldenSection::getset_bracket(PyGetSetDef& getset)
{
    getset.name=(char*)"bracket";
    getset.doc=(char*)"set to true if perform bracketing prior to line search";
    getset.get=[](PyObject* self,void*)->PyObject*
    {
        return var<bool>::build(reinterpret_cast<Object*>(self)->ls.brack,NULL);
    };
    getset.set=[](PyObject* self,PyObject* op,void*)->int
    {
        VarAPI<bool> bracket("bracket");
        int ichk=bracket.set(op);
        if(ichk==-1) return -1;
        reinterpret_cast<Object*>(self)->ls.brack=bracket.val;
        return 0;
    };
}
/*--------------------------------------------
 
 --------------------------------------------*/
void LineSearchGoldenSection::getset_max_iter(PyGetSetDef& getset)
{
    getset.name=(char*)"max_iter";
    getset.doc=(char*)"maximum number of iterations";
    getset.get=[](PyObject* self,void*)->PyObject*
    {
        return var<int>::build(reinterpret_cast<Object*>(self)->ls.max_iter,NULL);
    };
    getset.set=[](PyObject* self,PyObject* op,void*)->int
    {
        VarAPI<int> max_iter("max_iter");
        max_iter.logics[0]=VLogics("gt",0);
        int ichk=max_iter.set(op);
        if(ichk==-1) return -1;
        reinterpret_cast<Object*>(self)->ls.max_iter=max_iter.val;
        return 0;
    };
}
/*--------------------------------------------
 
 --------------------------------------------*/
void LineSearchGoldenSection::getset_tol(PyGetSetDef& getset)
{
    getset.name=(char*)"tol";
    getset.doc=(char*)"tolerance to stop minimization";
    getset.get=[](PyObject* self,void*)->PyObject*
    {
        return var<type0>::build(reinterpret_cast<Object*>(self)->ls.tol,NULL);
    };
    getset.set=[](PyObject* self,PyObject* op,void*)->int
    {
        VarAPI<type0> tol("tol");
        tol.logics[0]=VLogics("gt",0.0);
        int ichk=tol.set(op);
        if(ichk==-1) return -1;
        reinterpret_cast<Object*>(self)->ls.tol=tol.val;
        return 0;
    };
}
/*--------------------------------------------
 constructor
 --------------------------------------------*/
LineSearchBrent::LineSearchBrent():LineSearch()
{
    max_iter=5;
    tol=sqrt(2.0*std::numeric_limits<type0>::epsilon());
    zeps=std::numeric_limits<type0>::epsilon()*1.0e-3;
    brack=true;
}
/*--------------------------------------------
 destructor
 --------------------------------------------*/
//LineSearchBrent::~LineSearchBrent(){}
/*------------------------------------------------------------------------------------------------------------------------------------
 
 ------------------------------------------------------------------------------------------------------------------------------------*/
PyObject* LineSearchBrent::__new__(PyTypeObject* type,PyObject* args,PyObject* kwds)
{
    Object* __self=reinterpret_cast<Object*>(type->tp_alloc(type,0));
    PyObject* self=reinterpret_cast<PyObject*>(__self);
    return self;
}
/*--------------------------------------------
 
 --------------------------------------------*/
int LineSearchBrent::__init__(PyObject* self,PyObject* args,PyObject* kwds)
{
    FuncAPI<> f("__init__");
    
    if(f(args,kwds)==-1) return -1;
    return 0;
}
/*--------------------------------------------
 
 --------------------------------------------*/
PyObject* LineSearchBrent::__alloc__(PyTypeObject* type,Py_ssize_t)
{
    Object* __self=new Object;
    __self->ob_type=type;
    __self->ob_refcnt=1;
    return reinterpret_cast<PyObject*>(__self);
}
/*--------------------------------------------
 
 --------------------------------------------*/
void LineSearchBrent::__dealloc__(PyObject* self)
{
    Object* __self=reinterpret_cast<Object*>(self);
    delete __self;
}
/*--------------------------------------------*/
PyTypeObject LineSearchBrent::TypeObject={PyObject_HEAD_INIT(NULL)};
/*--------------------------------------------*/
void LineSearchBrent::setup_tp()
{
    TypeObject.tp_name="ls_brent";
    TypeObject.tp_doc="line search brent";
    
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
    
    TypeObject.tp_base=&LineSearch::TypeObject;
}
/*--------------------------------------------*/
PyGetSetDef LineSearchBrent::getset[]={[0 ... 4]={NULL,NULL,NULL,NULL,NULL}};
/*--------------------------------------------*/
void LineSearchBrent::setup_tp_getset()
{
    getset_bracket(getset[0]);
    getset_max_iter(getset[1]);
    getset_tol(getset[2]);
    getset_zeps(getset[3]);
}
/*--------------------------------------------*/
PyMethodDef LineSearchBrent::methods[]={[0 ... 0]={NULL}};
/*--------------------------------------------*/
void LineSearchBrent::setup_tp_methods()
{
}
/*--------------------------------------------
 
 --------------------------------------------*/
void LineSearchBrent::getset_bracket(PyGetSetDef& getset)
{
    getset.name=(char*)"bracket";
    getset.doc=(char*)"set to true if perform bracketing prior to line search";
    getset.get=[](PyObject* self,void*)->PyObject*
    {
        return var<bool>::build(reinterpret_cast<Object*>(self)->ls.brack,NULL);
    };
    getset.set=[](PyObject* self,PyObject* op,void*)->int
    {
        VarAPI<bool> bracket("bracket");
        int ichk=bracket.set(op);
        if(ichk==-1) return -1;
        reinterpret_cast<Object*>(self)->ls.brack=bracket.val;
        return 0;
    };
}
/*--------------------------------------------
 
 --------------------------------------------*/
void LineSearchBrent::getset_max_iter(PyGetSetDef& getset)
{
    getset.name=(char*)"max_iter";
    getset.doc=(char*)"maximum number of iterations";
    getset.get=[](PyObject* self,void*)->PyObject*
    {
        return var<int>::build(reinterpret_cast<Object*>(self)->ls.max_iter,NULL);
    };
    getset.set=[](PyObject* self,PyObject* op,void*)->int
    {
        VarAPI<int> max_iter("max_iter");
        max_iter.logics[0]=VLogics("gt",0);
        int ichk=max_iter.set(op);
        if(ichk==-1) return -1;
        reinterpret_cast<Object*>(self)->ls.max_iter=max_iter.val;
        return 0;
    };
}
/*--------------------------------------------
 
 --------------------------------------------*/
void LineSearchBrent::getset_tol(PyGetSetDef& getset)
{
    getset.name=(char*)"tol";
    getset.doc=(char*)"tolerance to stop minimization";
    getset.get=[](PyObject* self,void*)->PyObject*
    {
        return var<type0>::build(reinterpret_cast<Object*>(self)->ls.tol,NULL);
    };
    getset.set=[](PyObject* self,PyObject* op,void*)->int
    {
        VarAPI<type0> tol("tol");
        tol.logics[0]=VLogics("gt",0.0);
        int ichk=tol.set(op);
        if(ichk==-1) return -1;
        reinterpret_cast<Object*>(self)->ls.tol=tol.val;
        return 0;
    };
}
/*--------------------------------------------
 
 --------------------------------------------*/
void LineSearchBrent::getset_zeps(PyGetSetDef& getset)
{
    getset.name=(char*)"zeps";
    getset.doc=(char*)
    R"---(small numnber to protect against trying to achieve fractional accuracy for a minimum that happens to be exactly 0.0)---";
    getset.get=[](PyObject* self,void*)->PyObject*
    {
        return var<type0>::build(reinterpret_cast<Object*>(self)->ls.zeps,NULL);
    };
    getset.set=[](PyObject* self,PyObject* op,void*)->int
    {
        VarAPI<type0> zeps("zeps");
        zeps.logics[0]=VLogics("gt",0.0);
        int ichk=zeps.set(op);
        if(ichk==-1) return -1;
        reinterpret_cast<Object*>(self)->ls.zeps=zeps.val;
        return 0;
    };
}
/*--------------------------------------------
 constructor
 --------------------------------------------*/
LineSearchBackTrack::LineSearchBackTrack():LineSearch()
{
    min_alpha=0.0;
    c=0.4;
    rho=0.5;
}
/*--------------------------------------------
 destructor
 --------------------------------------------*/
//LineSearchBackTrack::~LineSearchBackTrack(){}
/*------------------------------------------------------------------------------------------------------------------------------------
 
 ------------------------------------------------------------------------------------------------------------------------------------*/
PyObject* LineSearchBackTrack::__new__(PyTypeObject* type,PyObject* args,PyObject* kwds)
{
    Object* __self=reinterpret_cast<Object*>(type->tp_alloc(type,0));
    PyObject* self=reinterpret_cast<PyObject*>(__self);
    return self;
}
/*--------------------------------------------
 
 --------------------------------------------*/
int LineSearchBackTrack::__init__(PyObject* self,PyObject* args,PyObject* kwds)
{
    FuncAPI<> f("__init__");
    
    if(f(args,kwds)==-1) return -1;
    return 0;
}
/*--------------------------------------------
 
 --------------------------------------------*/
PyObject* LineSearchBackTrack::__alloc__(PyTypeObject* type,Py_ssize_t)
{
    Object* __self=new Object;
    __self->ob_type=type;
    __self->ob_refcnt=1;
    return reinterpret_cast<PyObject*>(__self);
}
/*--------------------------------------------
 
 --------------------------------------------*/
void LineSearchBackTrack::__dealloc__(PyObject* self)
{
    Object* __self=reinterpret_cast<Object*>(self);
    delete __self;
}
/*--------------------------------------------*/
PyTypeObject LineSearchBackTrack::TypeObject={PyObject_HEAD_INIT(NULL)};
/*--------------------------------------------*/
void LineSearchBackTrack::setup_tp()
{
    TypeObject.tp_name="ls_bt";
    TypeObject.tp_doc="line search backtrack";
    
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
    
    TypeObject.tp_base=&LineSearch::TypeObject;
}
/*--------------------------------------------*/
PyGetSetDef LineSearchBackTrack::getset[]={[0 ... 3]={NULL,NULL,NULL,NULL,NULL}};
/*--------------------------------------------*/
void LineSearchBackTrack::setup_tp_getset()
{
    getset_c(getset[0]);
    getset_rho(getset[1]);
    getset_min_alpha(getset[2]);
}
/*--------------------------------------------*/
PyMethodDef LineSearchBackTrack::methods[]={[0 ... 0]={NULL}};
/*--------------------------------------------*/
void LineSearchBackTrack::setup_tp_methods()
{
}
/*--------------------------------------------
 
 --------------------------------------------*/
void LineSearchBackTrack::getset_c(PyGetSetDef& getset)
{
    getset.name=(char*)"c";
    getset.doc=(char*)"parameter to determine sufficent decrease";
    getset.get=[](PyObject* self,void*)->PyObject*
    {
        return var<type0>::build(reinterpret_cast<Object*>(self)->ls.c,NULL);
    };
    getset.set=[](PyObject* self,PyObject* op,void*)->int
    {
        VarAPI<type0> c("c");
        c.logics[0]=VLogics("gt",0.0)*VLogics("lt",1.0);
        int ichk=c.set(op);
        if(ichk==-1) return -1;
        reinterpret_cast<Object*>(self)->ls.c=c.val;
        return 0;
    };
}
/*--------------------------------------------
 
 --------------------------------------------*/
void LineSearchBackTrack::getset_rho(PyGetSetDef& getset)
{
    getset.name=(char*)"rho";
    getset.doc=(char*)"step reduction in every iterations";
    getset.get=[](PyObject* self,void*)->PyObject*
    {
        return var<type0>::build(reinterpret_cast<Object*>(self)->ls.rho,NULL);
    };
    getset.set=[](PyObject* self,PyObject* op,void*)->int
    {
        VarAPI<type0> rho("rho");
        rho.logics[0]=VLogics("gt",0.0)*VLogics("lt",1.0);
        int ichk=rho.set(op);
        if(ichk==-1) return -1;
        reinterpret_cast<Object*>(self)->ls.rho=rho.val;
        return 0;
    };
}
/*--------------------------------------------
 
 --------------------------------------------*/
void LineSearchBackTrack::getset_min_alpha(PyGetSetDef& getset)
{
    getset.name=(char*)"min_alpha";
    getset.doc=(char*)"minimum alpha";
    getset.get=[](PyObject* self,void*)->PyObject*
    {
        return var<type0>::build(reinterpret_cast<Object*>(self)->ls.min_alpha,NULL);
    };
    getset.set=[](PyObject* self,PyObject* op,void*)->int
    {
        VarAPI<type0> min_alpha("min_alpha");
        min_alpha.logics[0]=VLogics("ge",0.0);
        int ichk=min_alpha.set(op);
        if(ichk==-1) return -1;
        reinterpret_cast<Object*>(self)->ls.min_alpha=min_alpha.val;
        return 0;
    };
}
