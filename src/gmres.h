#ifndef __MAPP__gmres__
#define __MAPP__gmres__
#include "atoms.h"
/*--------------------------------------------
 
 --------------------------------------------*/
namespace MAPP_NS
{
    class __GMRES
    {
    private:
        const int m;
        const int n;
        type0* Q;
        type0** A_hat;
        type0* Ax_hat;
        type0(* cos_sin)[2];
        
        
        type0* x_hat;
        

        type0 calc(type0*,type0*);
        type0 calc(int,type0*,type0*);
        
        type0 solve_y(int,type0*);
        
        
        MPI_Comm& world;
    protected:
    public:
        __GMRES(Atoms*,int,int);
        ~__GMRES();
        template<class KERNEL>
        bool solve(KERNEL&,Vec<type0>*,type0,type0&,Vec<type0>*);
        
    };
}
using namespace MAPP_NS;
/*--------------------------------------------
 
 --------------------------------------------*/
template<class KERNEL>
bool __GMRES::solve(KERNEL& A,Vec<type0>* Ax,type0 tol,type0& norm,Vec<type0>* x)
{
    type0* __Ax=Ax->begin();
    type0* __x=x->begin();
    calc(__Ax,__x);
    
    for(int i=0;i<m;i++)
    {
        A(x,Ax);
        if(calc(i,__Ax,__x)<tol)
        {
            norm=solve_y(i+1,__x);
            return true;
        }
    }
    
    norm=solve_y(m,__x);
    return false;
}
#endif
