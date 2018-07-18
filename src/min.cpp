#include "min.h"
#include "ls.h"
#include "ff_styles.h"
#include "xmath.h"
using namespace MAPP_NS;
const char* Min::err_msgs[]=
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
Min::Min():
H_dof{DESIG2(__dim__,__dim__,false)},
ls(NULL),
chng_box(false),
e_tol(sqrt(std::numeric_limits<type0>::epsilon())),
affine(false),
max_dx(1.0),
ntally(1000)
{
}
/*--------------------------------------------
 constructor
 --------------------------------------------*/
Min::Min(type0 __e_tol,bool(&__H_dof)[__dim__][__dim__],bool __affine,type0 __max_dx,LineSearch* __ls):
H_dof{DESIG2(__dim__,__dim__,false)},
ls(__ls),
chng_box(false),
e_tol(__e_tol),
affine(__affine),
max_dx(__max_dx),
ntally(1000)
{
    for(int i=0;i<__dim__;i++)
        for(int j=0;j<__dim__;j++)
        {
            H_dof[i][j]=__H_dof[i][j];
            if(H_dof[i][j]) chng_box=true;
        }
}
/*--------------------------------------------
 destructor
 --------------------------------------------*/
Min::~Min()
{
}
/*--------------------------------------------
 pre run check it throw excepctions
 --------------------------------------------*/
void Min::pre_run_chk(Atoms* atoms,ForceField* ff)
{
    //check if configuration is loaded
    if(!atoms)
        throw std::string("cannot start minimization without initial conditions");
    
    //check if force field is loaded
    if(!ff)
        throw std::string("cannot start minimization without governing equations (force field)");
    
    //check to see if the H_dof components are consistent with stoms->dof
    if(chng_box && !atoms->dof->is_empty())
    {
        bool* dof=atoms->dof->begin();
        int __dof_lcl[__dim__]{DESIG(__dim__,0)};
        for(int i=0;i<atoms->natms_lcl;i++,dof+=__dim__)
            Algebra::Do<__dim__>::func([&dof,&__dof_lcl](int i){ if(!dof[i]) __dof_lcl[i]=1;});
        
        int __dof[__dim__]{DESIG(__dim__,0)};
        MPI_Allreduce(__dof_lcl,__dof,__dim__,MPI_INT,MPI_MAX,atoms->world);
        std::string err_msg=std::string();
        for(int i=0;i<__dim__;i++)
            for(int j=i;j<__dim__;j++)
                if(H_dof[i][j] && __dof[i])
                {
                    if(!err_msg.empty()) err_msg+="\n";
                    err_msg+="cannot impose stress component ["+Print::to_string(i)+"]["+Print::to_string(j)
                    +"] while any of the atoms do not have degree freedom in "+Print::to_string(i)
                    +" direction";
                }
        
        if(!err_msg.empty()) throw err_msg;
    }
}
/*------------------------------------------------------------------------------------------------------------------------------------
 
 ------------------------------------------------------------------------------------------------------------------------------------*/
PyObject* Min::__new__(PyTypeObject* type,PyObject* args,PyObject* kwds)
{
    Object* __self=reinterpret_cast<Object*>(type->tp_alloc(type,0));
    PyObject* self=reinterpret_cast<PyObject*>(__self);
    return self;
}
/*--------------------------------------------
 
 --------------------------------------------*/
int Min::__init__(PyObject* self,PyObject* args,PyObject* kwds)
{
    FuncAPI<> f("__init__");
    
    if(f(args,kwds)==-1) return -1;
    Object* __self=reinterpret_cast<Object*>(self);
    __self->min=new Min();
    __self->ls=NULL;
    __self->xprt=NULL;
    return 0;
}
/*--------------------------------------------
 
 --------------------------------------------*/
PyObject* Min::__alloc__(PyTypeObject* type,Py_ssize_t)
{
    Object* __self=new Object;
    Py_TYPE(__self)=type;
    Py_REFCNT(__self)=1;
    __self->min=NULL;
    __self->ls=NULL;
    __self->xprt=NULL;
    return reinterpret_cast<PyObject*>(__self);
}
/*--------------------------------------------
 
 --------------------------------------------*/
void Min::__dealloc__(PyObject* self)
{
    Object* __self=reinterpret_cast<Object*>(self);
    delete __self->min;
    __self->min=NULL;
    if(__self->ls) Py_DECREF(__self->ls);
    __self->ls=NULL;
    if(__self->xprt) Py_DECREF(__self->xprt);
    __self->xprt=NULL;
    delete __self;
}
/*--------------------------------------------*/
PyTypeObject Min::TypeObject={PyObject_HEAD_INIT(NULL)};
/*--------------------------------------------*/
int Min::setup_tp()
{
    TypeObject.tp_name="min_cg";
    TypeObject.tp_doc="conjugate gradient minimization";
    
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
    return ichk;
}
/*--------------------------------------------*/
PyGetSetDef Min::getset[]=EmptyPyGetSetDef(6);
/*--------------------------------------------*/
void Min::setup_tp_getset()
{
    getset_max_dx(getset[0]);
    getset_e_tol(getset[1]);
    getset_affine(getset[2]);
    getset_H_dof(getset[3]);
    getset_ntally(getset[4]);
}
/*--------------------------------------------*/
PyMethodDef Min::methods[]=EmptyPyMethodDef(1);
/*--------------------------------------------*/
void Min::setup_tp_methods()
{
}
/*--------------------------------------------
 
 --------------------------------------------*/
void Min::getset_e_tol(PyGetSetDef& getset)
{
    getset.name=(char*)"e_tol";
    getset.doc=(char*)R"---(
    (double) tolerance
    
    Energy tolerance criterion for stopping minimization
    )---";
    getset.get=[](PyObject* self,void*)->PyObject*
    {
        return var<type0>::build(reinterpret_cast<Object*>(self)->min->e_tol);
    };
    getset.set=[](PyObject* self,PyObject* op,void*)->int
    {
        VarAPI<type0> e_tol("e_tol");
        e_tol.logics[0]=VLogics("ge",0.0);
        int ichk=e_tol.set(op);
        if(ichk==-1) return -1;
        reinterpret_cast<Object*>(self)->min->e_tol=e_tol.val;
        return 0;
    };
}
/*--------------------------------------------
 
 --------------------------------------------*/
void Min::getset_H_dof(PyGetSetDef& getset)
{
    getset.name=(char*)"H_dof";
    getset.doc=(char*)R"---(
    (symm<bool[dim][dim]>) unitcell DOFs
    
    Unitcell degrees of freedom during minimization, here dim is the dimension of simulation
    )---";
    getset.get=[](PyObject* self,void*)->PyObject*
    {
        return var<symm<bool[__dim__][__dim__]>>::build(reinterpret_cast<Object*>(self)->min->H_dof);
    };
    getset.set=[](PyObject* self,PyObject* op,void*)->int
    {
        VarAPI<symm<bool[__dim__][__dim__]>> H_dof("H_dof");
        int ichk=H_dof.set(op);
        if(ichk==-1) return -1;
        
        bool (&__H_dof)[__dim__][__dim__]=reinterpret_cast<Object*>(self)->min->H_dof;
        bool chng_box=false;
        Algebra::DoLT<__dim__>::func([&H_dof,&__H_dof,&chng_box](int i,int j)
        {
            __H_dof[i][j]=__H_dof[j][i]=H_dof.val[i][j];
            if(__H_dof[i][j]) chng_box=true;
        });
        reinterpret_cast<Object*>(self)->min->chng_box=chng_box;
        return 0;
    };
}
/*--------------------------------------------
 
 --------------------------------------------*/
void Min::getset_affine(PyGetSetDef& getset)
{
    getset.name=(char*)"affine";
    getset.doc=(char*)R"---(
    (bool) affine transformation switch
    
    If set to True atomic displacements would be affine
    )---";
    getset.get=[](PyObject* self,void*)->PyObject*
    {
        return var<bool>::build(reinterpret_cast<Object*>(self)->min->affine);
    };
    getset.set=[](PyObject* self,PyObject* op,void*)->int
    {
        VarAPI<bool> affine("affine");
        int ichk=affine.set(op);
        if(ichk==-1) return -1;
        reinterpret_cast<Object*>(self)->min->affine=affine.val;
        return 0;
    };
}
/*--------------------------------------------
 
 --------------------------------------------*/
void Min::getset_max_dx(PyGetSetDef& getset)
{
    getset.name=(char*)"max_dx";
    getset.doc=(char*)R"---(
    (double) mximum displacemant
    
    Maximum displacement of any atom in one step of minimization
    )---";
    getset.get=[](PyObject* self,void*)->PyObject*
    {
        return var<type0>::build(reinterpret_cast<Object*>(self)->min->max_dx);
    };
    getset.set=[](PyObject* self,PyObject* op,void*)->int
    {
        VarAPI<type0> max_dx("max_dx");
        max_dx.logics[0]=VLogics("gt",0.0);
        int ichk=max_dx.set(op);
        if(ichk==-1) return -1;
        reinterpret_cast<Object*>(self)->min->max_dx=max_dx.val;
        return 0;
    };
}
/*--------------------------------------------
 
 --------------------------------------------*/
void Min::getset_ls(PyGetSetDef& getset)
{
    getset.name=(char*)"ls";
    getset.doc=(char*)R"---(
    (mapp.ls) line search algorithm
    
    Line search method to find the energy minimum in one dimension
    )---";
    getset.get=[](PyObject* self,void*)->PyObject*
    {
        LineSearch::Object* ls=reinterpret_cast<Object*>(self)->ls;
        Py_INCREF(ls);
        return reinterpret_cast<PyObject*>(ls);
    };
    getset.set=[](PyObject* self,PyObject* op,void*)->int
    {
        VarAPI<OP<LineSearch>> ls("ls");
        int ichk=ls.set(op);
        if(ichk==-1) return -1;
        LineSearch* __ls=&reinterpret_cast<LineSearch::Object*>(ls.val.ob)->ls;
        Py_DECREF(reinterpret_cast<Object*>(self)->ls);
        Py_INCREF(ls.val.ob);
        reinterpret_cast<Object*>(self)->ls=reinterpret_cast<LineSearch::Object*>(ls.val.ob);
        reinterpret_cast<Object*>(self)->min->ls=__ls;
        return 0;
    };
}
/*--------------------------------------------
 
 --------------------------------------------*/
void Min::getset_ntally(PyGetSetDef& getset)
{
    getset.name=(char*)"ntally";
    getset.doc=(char*)R"---(
    (int) thermodynamic tallying period
    
    Number of steps to be taken from one thermodynamics output to the next.
    )---";
    getset.get=[](PyObject* self,void*)->PyObject*
    {
        return var<int>::build(reinterpret_cast<Object*>(self)->min->ntally);
    };
    getset.set=[](PyObject* self,PyObject* op,void*)->int
    {
        VarAPI<int> ntally("ntally");
        ntally.logics[0]=VLogics("ge",0);
        int ichk=ntally.set(op);
        if(ichk==-1) return -1;
        reinterpret_cast<Object*>(self)->min->ntally=ntally.val;
        return 0;
    };
}



