#include "atoms_dmd.h"
#include "xmath.h"
/*--------------------------------------------
 
 --------------------------------------------*/
AtomsDMD::AtomsDMD(MPI_Comm& world,int __c_dim,int __N):
c_dim(__c_dim),
N(__N),
Atoms(world),
xi(new type0[__N]),
wi(new type0[__N]),
temp(NAN),
S_fe{DESIG2(__dim__,__dim__,NAN)},
fe(NAN),
s(NAN)
{
    XMath::quadrature_hg(N,xi,wi);
    elem=new Vec<elem_type>(this,c_dim,"elem");
    alpha=new DMDVec<type0>(this,0.0,"alpha");
    c=new DMDVec<type0>(this,-1.0,"c");
    dof_alpha=new DMDVec<bool>(this,true,"dof_alpha");
    dof_c=new DMDVec<bool>(this,true,"c_dof");

    dof_alpha->empty(true);
    dof_c->empty(true);
}
/*--------------------------------------------
 
 --------------------------------------------*/
AtomsDMD::~AtomsDMD()
{
    delete dof_alpha;
    delete dof_c;
    delete c;
    delete alpha;
    delete elem;
    
    delete [] wi;
    delete [] xi;
}
/*--------------------------------------------
 
 --------------------------------------------*/
#include "elements.h"
AtomsDMD& AtomsDMD::operator=(const Atoms& r)
{
    elements=r.elements;
    comm=r.comm;
    natms_lcl=r.natms_lcl;
    natms_ph=r.natms_ph;
    natms=r.natms;
    step=r.step;
    
    max_cut=r.max_cut;
    kB=r.kB;
    hP=r.hP;
    
    vol=r.vol;
    memcpy(depth_inv,r.depth_inv,__dim__*sizeof(type0));
    memcpy(__h,r.__h,__nvoigt__*sizeof(type0));
    memcpy(__b,r.__b,__nvoigt__*sizeof(type0));
    memcpy(&H[0][0],&r.H[0][0],__dim__*__dim__*sizeof(type0));
    memcpy(&B[0][0],&r.B[0][0],__dim__*__dim__*sizeof(type0));
    
    for(int i=0;i<nvecs;i++)
        if(!vecs[i]->is_empty())
            vecs[i]->resize(natms_lcl);
    
    
    memcpy(x->begin(),r.x->begin(),natms_lcl*__dim__*sizeof(type0));
    memcpy(id->begin(),r.id->begin(),natms_lcl*sizeof(unsigned int));
    return* this;
}
/*--------------------------------------------
 
 --------------------------------------------*/
type0 AtomsDMD::vac_msd()
{
    type0* __c=c->begin();
    type0* __x=x->begin();
    type0 cv_vals_lcl[__dim__+2];
    Algebra::zero<__dim__+2>(cv_vals_lcl);
    type0 __xi[__dim__];
    type0 cv;
    for(int i=0;i<natms_lcl;i++,__c+=c_dim,__x+=__dim__)
    {
        cv=1.0;
        for(int j=0;j<c_dim;j++)
            if(__c[j]>0.0)
                cv-=__c[j];
        
        Algebra::V_eq<__dim__>(__x,__xi);
        Algebra::X2S<__dim__>(B[0],__xi);
        Algebra::S2X<__dim__>(H[0],__xi);
        cv_vals_lcl[__dim__]+=cv*Algebra::V_mul_V<__dim__>(__xi,__xi);
        Algebra::V_add_x_mul_V<__dim__>(cv,__xi,cv_vals_lcl);
        cv_vals_lcl[__dim__+1]+=cv;
    }
    
    type0 cv_vals[__dim__+2];
    Algebra::zero<__dim__+2>(cv_vals);
    
    MPI_Allreduce(cv_vals_lcl,cv_vals,__dim__+2,Vec<type0>::MPI_T,MPI_SUM,world);
    
    
    Algebra::Do<__dim__+1>::func([&cv_vals](int i)
    {
        cv_vals[i]/=cv_vals[__dim__+1];
    });

    return cv_vals[__dim__]-Algebra::V_mul_V<__dim__>(cv_vals,cv_vals);
}
/*------------------------------------------------------------------------------------------------------------------------------------
 
 ------------------------------------------------------------------------------------------------------------------------------------*/
#include "ff_styles.h"
#include "import_styles.h"
PyObject* AtomsDMD::__new__(PyTypeObject* type,PyObject* args,PyObject* kwds)
{
    Object* __self=reinterpret_cast<Object*>(type->tp_alloc(type,0));
    PyObject* self=reinterpret_cast<PyObject*>(__self);
    return self;
}
/*--------------------------------------------
 
 --------------------------------------------*/
int AtomsDMD::__init__(PyObject* self,PyObject* args,PyObject* kwds)
{
    FuncAPI<> func("__init__");
    if(func(args,kwds)==-1) return -1;
    Object* __self=reinterpret_cast<Object*>(self);
    __self->atoms=NULL;
    __self->ff=NULL;
    return 0;
}
/*--------------------------------------------
 
 --------------------------------------------*/
PyObject* AtomsDMD::__alloc__(PyTypeObject* type,Py_ssize_t)
{
    Object* __self=new Object;
    __self->ob_type=type;
    __self->ob_refcnt=1;
    __self->atoms=NULL;
    __self->ff=NULL;
    return reinterpret_cast<PyObject*>(__self);
}
/*--------------------------------------------
 
 --------------------------------------------*/
void AtomsDMD::__dealloc__(PyObject* self)
{
    Object* __self=reinterpret_cast<Object*>(self);
    delete __self->ff;
    delete __self->atoms;
    delete __self;
}
/*--------------------------------------------*/
PyTypeObject AtomsDMD::TypeObject ={PyObject_HEAD_INIT(NULL)};
/*--------------------------------------------*/
int AtomsDMD::setup_tp()
{
    TypeObject.tp_name="mapp.dmd.atoms";
    TypeObject.tp_doc="container class";
    TypeObject.tp_flags=Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE;
    TypeObject.tp_basicsize=sizeof(Object);
    
    TypeObject.tp_new=__new__;
    TypeObject.tp_init=__init__;
    TypeObject.tp_alloc=__alloc__;
    TypeObject.tp_dealloc=__dealloc__;
    
    setup_tp_getset();
    TypeObject.tp_getset=getset;
    setup_tp_methods();
    TypeObject.tp_methods=methods;
    
    int ichk=PyType_Ready(&TypeObject);
    if(ichk<0) return ichk;
    Py_INCREF(&TypeObject);
    return ichk;
}
/*--------------------------------------------*/

#ifdef SC_DMD
PyGetSetDef AtomsDMD::getset[]=EmptyPyGetSetDef(19);
#else
PyGetSetDef AtomsDMD::getset[]=EmptyPyGetSetDef(17);
#endif

/*--------------------------------------------*/
void AtomsDMD::setup_tp_getset()
{
    getset_step(getset[0]);
    getset_hP(getset[1]);
    getset_kB(getset[2]);
    getset_H(getset[3]);
    getset_B(getset[4]);
    getset_vol(getset[5]);
    getset_elems(getset[6]);
    getset_skin(getset[7]);
    getset_comm_rank(getset[8]);
    getset_comm_size(getset[9]);
    getset_comm_coords(getset[10]);
    getset_comm_dims(getset[11]);
    getset_temp(getset[12]);
    getset_fe(getset[13]);
    getset_S_fe(getset[14]);
    getset_s(getset[15]);
#ifdef SC_DMD
    getset_BB(getset[16]);
    getset_delta(getset[17]);
#endif
}
/*--------------------------------------------*/

#ifdef SC_DMD
PyMethodDef AtomsDMD::methods[]=EmptyPyMethodDef(15);
#else
PyMethodDef AtomsDMD::methods[]=EmptyPyMethodDef(6);
#endif
/*--------------------------------------------*/
void AtomsDMD::setup_tp_methods()
{
    ml_strain(methods[0]);
    ImportCFGDMD::ml_import(methods[1]);
    ForceFieldEAMDMD::ml_new(methods[2],methods[3],methods[4]);
#ifdef SC_DMD
    ForceFieldEAMDMDSC::ml_new(methods[5],methods[6],methods[7]);
    ForceFieldEAMDMDSCC::ml_new(methods[8],methods[9],methods[10]);
    ForceFieldEAMDMDCLUSTER::ml_new(methods[11],methods[12],methods[13]);
#endif
}
/*--------------------------------------------
 
 --------------------------------------------*/
void AtomsDMD::getset_temp(PyGetSetDef& getset)
{
    getset.name=(char*)"temp";
    getset.doc=(char*)R"---(
    (double) temperature
    
    Temperature of the system
    )---";
    getset.get=[](PyObject* self,void*)->PyObject*
    {
        return var<type0>::build(reinterpret_cast<Object*>(self)->atoms->temp);
    };
    getset.set=[](PyObject* self,PyObject* val,void*)->int
    {
        Object* __self=reinterpret_cast<Object*>(self);
        if(!__self->atoms)
        {
            PyErr_Format(PyExc_TypeError,"cannot set 'temp' prior to loading system configuration");
            return -1;
        }
        
        VarAPI<type0> temp("temp");
        temp.logics[0]=VLogics("gt",0.0);
        if(temp.set(val)==-1) return -1;
        reinterpret_cast<Object*>(self)->atoms->temp=temp.val;
        return 0;
    };
}
/*--------------------------------------------
 
 --------------------------------------------*/
void AtomsDMD::getset_S_fe(PyGetSetDef& getset)
{
    getset.name=(char*)"S_fe";
    getset.doc=(char*)R"---(
    (double) stress
    
    Virial stress from free energy
    )---";
    getset.get=[](PyObject* self,void*)->PyObject*
    {
        return var<type0[__dim__][__dim__]>::build(reinterpret_cast<Object*>(self)->atoms->S_fe);
    };
    getset.set=[](PyObject* self,PyObject* val,void*)->int
    {
        PyErr_SetString(PyExc_TypeError,"readonly attribute");
        return -1;
    };
}
/*--------------------------------------------
 
 --------------------------------------------*/
void AtomsDMD::getset_fe(PyGetSetDef& getset)
{
    getset.name=(char*)"fe";
    getset.doc=(char*)R"---(
    (double) free energy
    
    Free energy
    )---";
    getset.get=[](PyObject* self,void*)->PyObject*
    {
        return var<type0>::build(reinterpret_cast<Object*>(self)->atoms->fe);
    };
    getset.set=[](PyObject* self,PyObject* val,void*)->int
    {
        PyErr_SetString(PyExc_TypeError,"readonly attribute");
        return -1;
    };
}
/*--------------------------------------------
 
 --------------------------------------------*/
void AtomsDMD::getset_s(PyGetSetDef& getset)
{
    getset.name=(char*)"s";
    getset.doc=(char*)R"---(
    (double) free energy
    
    Entropy
    )---";
    getset.get=[](PyObject* self,void*)->PyObject*
    {
        return var<type0>::build(reinterpret_cast<Object*>(self)->atoms->s);
    };
    getset.set=[](PyObject* self,PyObject* val,void*)->int
    {
        PyErr_SetString(PyExc_TypeError,"readonly attribute");
        return -1;
    };
}
#ifdef SC_DMD
/*--------------------------------------------
 
 --------------------------------------------*/
void AtomsDMD::getset_BB(PyGetSetDef& getset)
{
    getset.name=(char*)"BB";
    getset.doc=(char*)R"---(
    (double) free energy
    
    BB
    )---";
    getset.get=[](PyObject* self,void*)->PyObject*
    {
        return var<type0>::build(reinterpret_cast<Object*>(self)->atoms->BB);
    };
    getset.set=[](PyObject* self,PyObject* val,void*)->int
    {
        PyErr_SetString(PyExc_TypeError,"readonly attribute");
        return -1;
    };
}
/*--------------------------------------------
 
 --------------------------------------------*/
void AtomsDMD::getset_delta(PyGetSetDef& getset)
{
    getset.name=(char*)"delta";
    getset.doc=(char*)R"---(
    (double) temperature
    
    Temperature of the system
    )---";
    getset.get=[](PyObject* self,void*)->PyObject*
    {
        return var<type0>::build(reinterpret_cast<Object*>(self)->atoms->delta);
    };
    getset.set=[](PyObject* self,PyObject* val,void*)->int
    {
        Object* __self=reinterpret_cast<Object*>(self);
        if(!__self->atoms)
        {
            PyErr_Format(PyExc_TypeError,"cannot set 'delta' prior to loading system configuration");
            return -1;
        }
        
        VarAPI<type0> temp("delta");
        if(temp.set(val)==-1) return -1;
        reinterpret_cast<Object*>(self)->atoms->delta=temp.val;
        return 0;
    };
}
#endif


