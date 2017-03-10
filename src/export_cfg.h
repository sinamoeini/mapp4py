#ifndef __MAPP__export_cfg__
#define __MAPP__export_cfg__
#include "export.h"

namespace MAPP_NS
{
    class ExportCFGMD:public Export
    {
    private:
        
        std::string pattern;
        bool sort;
        bool x_d_inc;
    protected:
    public:
        ExportCFGMD(const std::string&,std::string*,size_t,bool);
        ~ExportCFGMD();
        void write_header(FILE*);
        void write_body_sort(FILE*);
        void write_body(FILE*);
        void write(int);
        void init();
        void fin();
        class AtomsMD* atoms;
    };
}
namespace MAPP_NS
{
    class ExportCFGDMD:public Export
    {
    private:
        
        std::string pattern;
        bool sort;
    protected:
    public:
        ExportCFGDMD(const std::string&,std::string*,size_t,bool);
        ~ExportCFGDMD();
        void write_header(FILE*);
        void write_body_sort(FILE*);
        void write_body(FILE*);
        void write(int);
        void init();
        void fin();
        class AtomsDMD* atoms;
    };
}

#endif 
