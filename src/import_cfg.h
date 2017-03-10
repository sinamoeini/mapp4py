#ifndef __MAPP__import_cfg__
#define __MAPP__import_cfg__
#include "import.h"
#include "global.h"
#include <stdio.h>

typedef struct PyMethodDef PyMethodDef;
namespace MAPP_NS
{
    class ImportCFG : public Import
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
        ImportCFG(MPI_Comm&);
        ~ImportCFG();
    };
}

namespace MAPP_NS
{
    class ImportCFGMD : public ImportCFG
    {
    private:
    protected:
    public:
        ImportCFGMD(MPI_Comm&);
        ~ImportCFGMD();
        class AtomsMD* operator()(const char*);
        static void ml_import(PyMethodDef&);
    };
}

namespace MAPP_NS
{
    class ImportCFGDMD : public ImportCFG
    {
    private:
    protected:
    public:
        ImportCFGDMD(MPI_Comm&);
        ~ImportCFGDMD();
        class AtomsDMD* operator()(int,const char*);
        static void ml_import(PyMethodDef&);
    };
}

#endif


