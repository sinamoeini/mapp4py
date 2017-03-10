#ifndef __MAPP__write_cfg_dmd__
#define __MAPP__write_cfg_dmd__
#include "write.h"

namespace MAPP_NS
{
    class WriteCFGDMD:public Write
    {
    private:
        
        std::string pattern;
        bool sort;
    protected:
    public:
        WriteCFGDMD(const std::string&,std::string*,size_t,bool);
        ~WriteCFGDMD();
        void write_header(FILE*);
        void write_body_sort(FILE*);
        void write_body(FILE*);
        void write(int);
        class AtomsDMD* atoms;
    };
}

#endif 
