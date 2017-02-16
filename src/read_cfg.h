#ifndef __MAPP__read_cfg__
#define __MAPP__read_cfg__
#include "read.h"
#include "global.h"
#include <stdio.h>

typedef struct PyMethodDef PyMethodDef;
namespace MAPP_NS
{
    class ReadCFG : public Read
    {
    private:
    protected:

        
        type0 R;
        int entry_count;
        bool is_ext;
        bool vel_xst;
        
        class vec* vec_list[4];
        int nvecs;
        
        void read_header(class Atoms*,FileReader&,char*&,size_t&);
        void read_body_ext(class Atoms*,FileReader&,char*&,size_t&);
        void read_body_std(class Atoms*,FileReader&,char*&,size_t&);
    public:
        ReadCFG(MPI_Comm&);
        ~ReadCFG();
    };
}

namespace MAPP_NS
{
    class ReadCFGMD : public ReadCFG
    {
    private:
    protected:
    public:
        ReadCFGMD(MPI_Comm&);
        ~ReadCFGMD();
        class AtomsMD* operator()(const char*);
        static void ml_cfg(PyMethodDef&);
    };
}

namespace MAPP_NS
{
    class ReadCFGDMD : public ReadCFG
    {
    private:
    protected:
    public:
        ReadCFGDMD(MPI_Comm&);
        ~ReadCFGDMD();
        class AtomsDMD* operator()(int,const char*);
        static void ml_cfg(PyMethodDef&);
    };
}

#endif


