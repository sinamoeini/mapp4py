#ifndef __MAPP__comm__
#define __MAPP__comm__
#include "api.h"
#include "mpi_compat.h"
#include "global.h"
/*----------------------------------------------------------------------------------------------------------------
 _____   _____       ___  ___       ___  ___   _   _   __   _   _   _____       ___   _____   _   _____   __   _  
/  ___| /  _  \     /   |/   |     /   |/   | | | | | |  \ | | | | /  ___|     /   | |_   _| | | /  _  \ |  \ | | 
| |     | | | |    / /|   /| |    / /|   /| | | | | | |   \| | | | | |        / /| |   | |   | | | | | | |   \| | 
| |     | | | |   / / |__/ | |   / / |__/ | | | | | | | |\   | | | | |       / / | |   | |   | | | | | | | |\   | 
| |___  | |_| |  / /       | |  / /       | | | |_| | | | \  | | | | |___   / /  | |   | |   | | | |_| | | | \  | 
\_____| \_____/ /_/        |_| /_/        |_| \_____/ |_|  \_| |_| \_____| /_/   |_|   |_|   |_| \_____/ |_|  \_|
 
 ----------------------------------------------------------------------------------------------------------------*/
namespace MAPP_NS
{
    class Communication
    {
    private:
        template<const int dim>
        void seq_analysis(const int,const int*,const int,int(&)[dim],int(&)[dim],int(&)[dim][2]);
    protected:
    public:
        static void network_analysis(MPI_Comm&,int&,int&,int&,int*&,int**&,bool&);
        static int get_rank();
        static int get_size();
        static int get_rank(MPI_Comm&);
        static int get_size(MPI_Comm&);
        
        
        const int rank;
        const int size;
        type0 skin;
        unsigned long xchng_id;
        
        MPI_Comm world;
        
        
        
        int coords[__dim__];
        int dims[__dim__];
        int neigh[__dim__][2];
        type0 s_lo[__dim__];
        type0 s_hi[__dim__];
        
        void grid(int(&)[__dim__],type0(&)[__dim__][__dim__]);
        void grid(type0(&)[__dim__][__dim__]);
        void grid(int(&)[__dim__]);
        Communication(MPI_Comm,int(&)[__dim__],type0(&)[__dim__][__dim__],type0=0.5);
        Communication(MPI_Comm,type0(&)[__dim__][__dim__],type0=0.5);
        Communication(MPI_Comm,type0=0.5);
        Communication(const Communication&);
        Communication(Communication&&);
        Communication& operator=(const Communication&);
    };
}
using namespace MAPP_NS;
/*--------------------------------------------
 
 --------------------------------------------*/
template<const int dim>
void Communication::seq_analysis(const int my_id,const int* seq,const int seq_sz,int(&dims)[dim],int(&coord)[dim],int(&neigh)[dim][2])
{
    int base[dim];
    base[0]=1;
    for(int i=1;i<dim;i++)
        base[i]=base[i-1]*dims[i-1];
    
    int my_rank=0;
    while(seq[my_rank]!=my_id) my_rank++;
    for(int i=dim-1;i>-1;i--)
    {
        coord[i]=my_rank/base[i];
        my_rank-=coord[i]*base[i];
    }
        
    for(int i=0;i<dim;i++)
    {
        neigh[i][0]=neigh[i][1]=0;
        for(int j=0;j<dim;j++)
        {
            if(j==i)
            {
                neigh[i][0]+=(coord[j]==0? dims[j]-1:coord[j]-1)*base[j];
                neigh[i][1]+=(coord[j]==dims[j]-1? 0:coord[j]+1)*base[j];
            }
            else
            {
                neigh[i][0]+=coord[j]*base[j];
                neigh[i][1]+=coord[j]*base[j];
            }
        }
        
        neigh[i][0]=seq[neigh[i][0]];
        neigh[i][1]=seq[neigh[i][1]];
    }
    
}
/*------------------------------------------------------------------------------------------------------------------------------------
 
 ------------------------------------------------------------------------------------------------------------------------------------*/
namespace MAPP_NS
{
    
    class MAPP_MPI
    {
    public:
        typedef struct
        {
            PyObject_HEAD
            MPI_Comm world;
            int rank;
            int size;
            
        }Object;
        static PyTypeObject TypeObject;
        static PyObject* __new__(PyTypeObject*,PyObject*, PyObject*);
        static int __init__(PyObject*, PyObject*,PyObject*);
        static PyObject* __alloc__(PyTypeObject*,Py_ssize_t);
        static void __dealloc__(PyObject*);
        
        static PyMemberDef tp_members[];
        static void setup_tp_members();
        static int setup_tp();
    };
}
#endif
