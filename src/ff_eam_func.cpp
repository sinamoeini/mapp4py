#include "ff_eam_func.h"
#include "neighbor_md.h"
#include "atoms_md.h"
#include "elements.h"
#include "memory.h"
#include "dynamic_md.h"
using namespace MAPP_NS;
/*--------------------------------------------
 This is for my personal and fitting of iron
 and hydrogen
 --------------------------------------------*/
/*--------------------------------------------
 constructor
 Fe: 0 H: 1
 list of inputs
    2x2
    for electron density:
        size_t nrho_ij: number of spilines
        type0** ARrho_ij: {a,r0}x n
    2
    for phi_00 phi_01
        size_t nphi_ij: number of spilines
        type0** ARphi_ij: {a,r0}x n
    2
    for embedded functions
 
 
 --------------------------------------------*/
ForceFieldEAMFunc::ForceFieldEAMFunc(AtomsMD* __atoms):
ForceFieldMD(__atoms)
{
}
/*--------------------------------------------
 destructor
 --------------------------------------------*/
ForceFieldEAMFunc::~ForceFieldEAMFunc()
{
}
/*--------------------------------------------
 force and energy calculation
 --------------------------------------------*/
void ForceFieldEAMFunc::force_calc()
{
}
/*--------------------------------------------
 only energy calculation this is useful for
 minimization/linesearch methods that do not
 use derivatives of energy
 --------------------------------------------*/
void ForceFieldEAMFunc::energy_calc()
{
    
}
/*--------------------------------------------
 initiate before a run
 --------------------------------------------*/
void ForceFieldEAMFunc::init()
{
    pre_init();
    
}
/*--------------------------------------------
 after a run
 --------------------------------------------*/
void ForceFieldEAMFunc::fin()
{
    
    post_fin();
}
/*--------------------------------------------
 init xchng
 --------------------------------------------*/
void ForceFieldEAMFunc::init_xchng()
{
}
/*--------------------------------------------
 fin xchng
 --------------------------------------------*/
void ForceFieldEAMFunc::fin_xchng()
{
}
/*--------------------------------------------
 pre xchng energy
 --------------------------------------------*/
void ForceFieldEAMFunc::pre_xchng_energy(GCMC*)
{
}
/*--------------------------------------------
 xchng energy
 --------------------------------------------*/
type0 ForceFieldEAMFunc::xchng_energy(GCMC*)
{
    return 0.0;
}
/*--------------------------------------------
 post xchng energy
 --------------------------------------------*/
void ForceFieldEAMFunc::post_xchng_energy(GCMC*)
{
}
/*--------------------------------------------
 python constructor
 --------------------------------------------*/


#include <dlfcn.h>
void ForceFieldEAMFunc::ml_new(PyMethodDef& tp_methods)
{
    tp_methods.ml_flags=METH_VARARGS | METH_KEYWORDS;
    tp_methods.ml_name="ff_eam_func";
    
    tp_methods.ml_meth=(PyCFunction)(PyCFunctionWithKeywords)(
    [](PyObject* self,PyObject* args,PyObject* kwds)->PyObject*
    {
        
        FuncAPI<std::string> f("ff_eam_func",{"so_file"});
        if(f(args,kwds)) return NULL;
        
        
        
        void* so=dlopen(f.val<0>().c_str(),RTLD_NOW);
        int nelem=*reinterpret_cast<int*>( dlsym(so,"nelem"));
        const char** elems=reinterpret_cast<const char **>( dlsym(so,"elems"));
        
        double(*fff)(double)=reinterpret_cast<double(*)(double)>(dlsym(so,"phi"));
        
        for(int i=0;i<nelem;i++)
            printf("%s\n",elems[i]);
        
        printf("%lf\n",fff(12.0));
        
        

        Py_RETURN_NONE;
    });
    
    tp_methods.ml_doc=(char*)R"---(
    ff_fs(A,t1,t2,k1,k2,k3,r_c_phi,r_c_rho,elems=None)
   
    Finnis-Sinclair EAM
    
    Assigns Finnis-Sinclair EAM force field to system. For explanation of the parameter see the Notes section.
    
    Parameters
    ----------
    A : double[nelems]
        :math:`A`
    t1 : double[nelems][nelems]
        :math:`t_1`
    t2 : double[nelems][nelems]
        :math:`t_2`
    k1 : symmetric double[nelems][nelems]
        :math:`k_1`
    k2 : symmetric double[nelems][nelems]
        :math:`k_2`
    k3 : symmetric double[nelems][nelems]
        :math:`k_3`
    r_c_phi : symmetric double[nelems][nelems]
        :math:`r_{c,\phi}`
    r_c_rho : symmetric double[nelems][nelems]
        :math:`r_{c,\rho}`
    elems : string[nelems]
        mapping elements
    
    Returns
    -------
    None
   
    Notes
    -----
    This is the analytical form of Finnis-Sinclair Embedded Atom Method (EAM) potential
    
    .. math::
        U=\sum_{i}\left( -A_{\alpha}\sqrt{\sum_{j\neq i} \rho_{\beta\alpha}(r_{ij})}  + \frac{1}{2}\sum_{j\neq i} \phi_{\beta\alpha}(r_{ij}) \right),
    
    where
    
    .. math::
        \rho_{\beta\alpha}(r)=
        \left\{\begin{array}{ll}
        t^{\alpha\beta}_1(r-r^{\alpha\beta}_{c,\rho})^2+t^{\alpha\beta}_2(r-r^{\alpha\beta}_{c,\rho})^3, \quad & r<r^{\alpha\beta}_{c,\rho}\\
        0 & r>r^{\alpha\beta}_{c,\rho}\
        \end{array}\right.

            
    and
        
    .. math::
        \phi_{\beta\alpha}(r)=
        \left\{\begin{array}{ll}
        (r-r^{\alpha\beta}_{c,\phi})^2(k^{\alpha\beta}_1+k^{\alpha\beta}_2 r+k^{\alpha\beta}_3 r^2), \quad & r<r^{\alpha\beta}_{c,\phi}\\
        0 & r>r^{\alpha\beta}_{c,\phi}\
        \end{array}\right.

    
    Examples
    --------
    Iron Carbon mixture
    ::
     
        >>> from mapp import md
        >>> sim=md.cfg("configs/Cementite.cfg")
        >>> sim.ff_fs(A=[1.8289905,2.9588787],
                      t1=[[1.0,10.024001],[10.482408,0.0]],
                      t2=[[0.504238,1.638980],[3.782595,-7.329211]],
                      k1=[[1.237115],[8.972488,22.061824]],
                      k2=[[-0.35921],[-4.086410,-17.468518]],
                      k3=[[-0.038560],[1.483233,4.812639]],
                      r_c_phi=[[3.40],[2.468801,2.875598]],
                      r_c_rho=[[3.569745],[2.545937,2.892070]],
                      elems=['Fe','C'])

    )---";
}
