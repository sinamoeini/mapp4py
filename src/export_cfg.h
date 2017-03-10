#ifndef __MAPP__export_cfg__
#define __MAPP__export_cfg__
#include "export.h"

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
        class AtomsDMD* atoms;
    };
}

#endif 
