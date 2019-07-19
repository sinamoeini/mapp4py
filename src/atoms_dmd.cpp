#include "atoms_dmd.h"
#include "xmath.h"
/*--------------------------------------------
 
 --------------------------------------------*/
AtomsDMD::AtomsDMD(MPI_Comm& world,int __c_dim,int __N):
Atoms(world),
S_fe{DESIG2(__dim__,__dim__,NAN)},
fe(NAN),
gp(NAN),
pe(NAN),
s(NAN),
c_dim(__c_dim),
N(__N),
xi(new type0[__N]),
wi(new type0[__N]),
temp(NAN)
{
    XMath::quadrature_hg(N,xi,wi);
    elem=new Vec<elem_type>(this,c_dim,"elem");
    alpha=new DMDVec<type0>(this,0.0,"alpha");
    c=new DMDVec<type0>(this,-1.0,"c");
    alpha_dof=new DMDVec<bool>(this,true,"alpha_dof");
    c_dof=new DMDVec<bool>(this,true,"c_dof");

    alpha_dof->empty(true);
    c_dof->empty(true);
}
/*--------------------------------------------
 
 --------------------------------------------*/
AtomsDMD::~AtomsDMD()
{
    delete alpha_dof;
    delete c_dof;
    delete c;
    delete alpha;
    delete elem;
    
    delete [] wi;
    delete [] xi;
}
/*--------------------------------------------
 use this function after changing H
 --------------------------------------------*/
void AtomsDMD::update_max_alpha()
{
    int n=natms_lcl*c_dim;
    type0 max_alpha_lcl=0.0;
    type0* alpha_vec=alpha->begin();
    for(int i=0;i<n;i++)
        max_alpha_lcl=MAX(max_alpha_lcl,alpha_vec[i]);
    MPI_Allreduce(&max_alpha_lcl,&max_alpha,1,Vec<type0>::MPI_T,MPI_MAX,world);
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
    memcpy(id->begin(),r.id->begin(),natms_lcl*sizeof(id_type));
    return* this;
}
/*--------------------------------------------
 
 --------------------------------------------*/
type0* AtomsDMD::ave_comp()
{
    size_t nelems=elements.nelems;
    type0* sum=NULL;
    type0* sum_lcl=NULL;
    type0* c_sum=NULL;
    type0* c_sum_lcl=NULL;
    if(nelems)
    {
        sum=new type0[nelems];
        sum_lcl=new type0[nelems];
        c_sum=new type0[nelems];
        c_sum_lcl=new type0[nelems];
    }
    for(size_t i=0;i<nelems;i++) sum_lcl[i]=c_sum_lcl[i]=0.0;
    
    type0* __c=c->begin();
    elem_type* __elem=elem->begin();
    for(int i=0;i<natms_lcl;i++,__c+=c_dim,__elem+=c_dim)
    {
        for(int j=0;j<c_dim;j++)
        {
            if(__c[j]<0.0) continue;
            c_sum_lcl[__elem[j]]+=__c[j];
            ++sum_lcl[__elem[j]];
        }
    }
    
    MPI_Allreduce(c_sum_lcl,c_sum,static_cast<int>(nelems),Vec<type0>::MPI_T,MPI_SUM,world);
    MPI_Allreduce(sum_lcl,sum,static_cast<int>(nelems),Vec<type0>::MPI_T,MPI_SUM,world);
    
    for(size_t i=0;i<nelems;i++) if(sum[i]!=0.0) c_sum[i]/=sum[i];
    
    delete [] c_sum_lcl;
    delete [] sum_lcl;
    delete [] sum;
    return c_sum;
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
/*--------------------------------------------
 
 --------------------------------------------*/
void AtomsDMD::DO(PyObject* op)
{
    class Func_id
    {
    public:
        void operator()(id_type* old_val,id_type* new_val)
        {
            if(*old_val!=*new_val)
                throw std::string("id of atoms cannot be changed");
        };
    };
    Func_id func_id;
    VecPy<id_type,Func_id> id_vec_py(id,func_id);
    
    class Func_x
    {
    public:
        void operator()(type0*,type0*)
        {};
    };
    Func_x func_x;
    VecPy<type0,Func_x> x_vec_py(x,func_x);
    
    class Func_dof
    {
    public:
        void operator()(bool*,bool*)
        {
        };
    };
    Func_dof func_dof;
    VecPy<bool,Func_dof> dof_vec_py(x_dof,func_dof);
    
    
    class Func_alpha
    {
    public:
        int c_dim;
        Func_alpha(int __c_dim):c_dim(__c_dim){}
        void operator()(type0* old_val,type0* new_val)
        {
            for(int i=0;i<c_dim;i++)
                if(new_val[i]<=0.0 && new_val[i]!=old_val[i])
                    throw std::string("all values of alpha corresponding to present elements (i.e. c!=-1.0) should be greater than 0.0");
        };
    };
    Func_alpha func_alpha(c_dim);
    DMDVecPy<type0,Func_alpha> alpha_vec_py(alpha,func_alpha);
    
    
    
    class Func_c
    {
    public:
        int c_dim;
        Func_c(int __c_dim):c_dim(__c_dim){}
        void operator()(type0* old_val,type0* new_val)
        {
            
            volatile type0 c0;
            for(int i=0;i<c_dim;i++)
            {
                c0=new_val[i];
                --++c0;
                new_val[i]=c0;
            }
            type0 sum_c=0.0;
            for(int i=0;i<c_dim;i++)
            {
                if(new_val[i]!=-1.0 && (new_val[i]<0.0 || new_val[i]>1.0))
                    throw std::string("all values of c should be greater or equal to 0.0 and less than or equal to 1.0, or equal to -1.0");
                if(new_val[i]!=-1.0) sum_c+=new_val[i];
                
            }
            if(sum_c>1.0)
                throw std::string("sum of present elements of one site should not exceed 1.0");
        };
    };
    Func_c func_c(c_dim);
    DMDVecPy<type0,Func_c> c_vec_py(c,func_c);
    
    
    
    DMDVecPy<bool,Func_dof> c_dof_vec_py(c_dof,func_dof);
    DMDVecPy<bool,Func_dof> alpha_dof_vec_py(alpha_dof,func_dof);
    try
    {
        VecPyFunc::Do(this,op,id_vec_py,x_vec_py,alpha_vec_py,c_vec_py,dof_vec_py,c_dof_vec_py,alpha_dof_vec_py);
    }
    catch(std::string& err_msg)
    {
        throw err_msg;
    }
    
    if(x_vec_py.inc) this->reset_domain();
    
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
    Py_TYPE(__self)=type;
    Py_REFCNT(__self)=1;
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


PyGetSetDef AtomsDMD::getset[]=EmptyPyGetSetDef(23);

/*--------------------------------------------*/
void AtomsDMD::setup_tp_getset()
{
    getset_natms(getset[0]);
    getset_step(getset[1]);
    getset_hP(getset[2]);
    getset_kB(getset[3]);
    getset_H(getset[4]);
    getset_B(getset[5]);
    getset_vol(getset[6]);
    getset_elems(getset[7]);
    getset_skin(getset[8]);
    getset_comm_rank(getset[9]);
    getset_comm_size(getset[10]);
    getset_comm_coords(getset[11]);
    getset_comm_dims(getset[12]);
    getset_temp(getset[13]);
    getset_fe(getset[14]);
    getset_pe(getset[15]);
    getset_gp(getset[16]);
    getset_S_fe(getset[17]);
    getset_s(getset[18]);
    getset_ave_mu(getset[19]);
    getset_ext_mu(getset[20]);
    getset_ave_comp(getset[21]);

}
/*--------------------------------------------*/
PyMethodDef AtomsDMD::methods[]=EmptyPyMethodDef(9);
/*--------------------------------------------*/
void AtomsDMD::setup_tp_methods()
{
    ml_do(methods[0]);
    ml_strain(methods[1]);
    ml_mul(methods[2]);
    ml_autogrid(methods[3]);
    ImportCFGDMD::ml_import(methods[4]);
    ForceFieldEAMDMD::ml_new(methods[5],methods[6],methods[7]);
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
void AtomsDMD::getset_pe(PyGetSetDef& getset)
{
    getset.name=(char*)"pe";
    getset.doc=(char*)R"---(
    (double) free energy
    
    Free energy
    )---";
    getset.get=[](PyObject* self,void*)->PyObject*
    {
        return var<type0>::build(reinterpret_cast<Object*>(self)->atoms->pe);
    };
    getset.set=[](PyObject* self,PyObject* val,void*)->int
    {
        PyErr_SetString(PyExc_TypeError,"readonly attribute");
        return -1;
    };
}
/*--------------------------------------------
 
 --------------------------------------------*/
void AtomsDMD::getset_gp(PyGetSetDef& getset)
{
    getset.name=(char*)"gp";
    getset.doc=(char*)R"---(
    (double) grand potential
    
    Grand potential
    )---";
    getset.get=[](PyObject* self,void*)->PyObject*
    {
        return var<type0>::build(reinterpret_cast<Object*>(self)->atoms->gp);
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
/*--------------------------------------------
 
 --------------------------------------------*/
void AtomsDMD::getset_ave_mu(PyGetSetDef& getset)
{
    getset.name=(char*)"ave_mu";
    getset.doc=(char*)R"---(
    (double[nelems]) average mu per element
    
    Average mu per element
    )---";
    getset.get=[](PyObject* __self,void*)->PyObject*
    {
        ForceFieldDMD* __ff=reinterpret_cast<Object*>(__self)->ff;
        AtomsDMD* __atoms=reinterpret_cast<Object*>(__self)->atoms;
        //check if force field is loaded
        if(!__ff)
        {
            PyErr_SetString(PyExc_TypeError,"cannot calculate ave_mu without governing equations (force field)");
            return NULL;
        }
        
        if(std::isnan(__atoms->kB))
        {
            PyErr_SetString(PyExc_TypeError,"boltzmann constant should be set prior to calculation of ave_mu");
            return NULL;
        }
        
        if(std::isnan(__atoms->hP))
        {
            PyErr_SetString(PyExc_TypeError,"planck constant should be set prior to calculation of ave_mu");
            return NULL;
        }
        
        __ff->calc_thermo();
        size_t* szp=&(__atoms->elements.nelems);
        return var<type0*>::build(__ff->ave_mu,&szp);
    };
    getset.set=[](PyObject* self,PyObject* val,void*)->int
    {
        PyErr_SetString(PyExc_TypeError,"readonly attribute");
        return -1;
    };
}
/*--------------------------------------------
 
 --------------------------------------------*/
void AtomsDMD::getset_ave_comp(PyGetSetDef& getset)
{
    getset.name=(char*)"ave_comp";
    getset.doc=(char*)R"---(
    (double[nelems]) average mu per element
    
    Average mu per element
    )---";
    getset.get=[](PyObject* __self,void*)->PyObject*
    {
        
        AtomsDMD* __atoms=reinterpret_cast<Object*>(__self)->atoms;
        type0* ave_comp=__atoms->ave_comp();
        size_t* szp=&(__atoms->elements.nelems);
        PyObject* op=var<type0*>::build(ave_comp,&szp);
        delete [] ave_comp;
        return op;
    };
    getset.set=[](PyObject* self,PyObject* val,void*)->int
    {
        PyErr_SetString(PyExc_TypeError,"readonly attribute");
        return -1;
    };
}
/*--------------------------------------------
 
 --------------------------------------------*/
void AtomsDMD::getset_ext_mu(PyGetSetDef& getset)
{
    getset.name=(char*)"ext_mu";
    getset.doc=(char*)R"---(
    (double[nelems]) external mu per element
    
    External mu per element
    )---";
    getset.get=[](PyObject* __self,void*)->PyObject*
    {
        ForceFieldDMD* __ff=reinterpret_cast<Object*>(__self)->ff;
        
        //check if force field is loaded
        if(!__ff)
        {
            PyErr_SetString(PyExc_TypeError,"cannot get ext_mu without governing equations (force field)");
            return NULL;
        }
        
        AtomsDMD* __atoms=reinterpret_cast<Object*>(__self)->atoms;
        size_t* szp=&(__atoms->elements.nelems);
        return var<type0*>::build(__ff->ext_mu,&szp);
    };
    getset.set=[](PyObject* __self,PyObject* val,void*)->int
    {
        ForceFieldDMD* __ff=reinterpret_cast<Object*>(__self)->ff;
        if(!__ff)
        {
            PyErr_SetString(PyExc_TypeError,"cannot set ext_mu without governing equations (force field)");
            return -1;
        }
        VarAPI<type0*> var("ext_mu");
        if(var.set(val)==-1) return -1;
        
        AtomsDMD* __atoms=reinterpret_cast<Object*>(__self)->atoms;
        
        if(var.__var__.size!=__atoms->elements.nelems)
        {
            PyErr_Format(PyExc_TypeError,"size mismatch");
            return -1;
        }
        
        memcpy(__ff->ext_mu,var.val,sizeof(type0)*__atoms->elements.nelems);
        return 0;
    };
}

