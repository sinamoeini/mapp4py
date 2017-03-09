#ifndef __MAPP__write_cfg_dmd__
#define __MAPP__write_cfg_dmd__
#include <iostream>

namespace MAPP_NS
{
    class WriteCFGDMD
    {
    private:
        class AtomsDMD* atoms;
        const char* pattern;
        bool sort;
    protected:
    public:
        WriteCFGDMD(class AtomsDMD*,const char*,bool);
        ~WriteCFGDMD();
        void write_header(FILE*);
        void write_body_sort(FILE*);
        void write_body(FILE*);
        void write(int);
    };
}

#endif 
