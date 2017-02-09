#ifndef __MAPP__rand_engine__
#define __MAPP__rand_engine__
#include "global.h"
namespace MAPP_NS
{
    class Random
    {
    private:
        type0 reserved;
        bool gauss_chk;
        int seed;
    protected:
    public:
        Random(int);
        ~Random();
        type0 uniform();
        type0 gaussian();
        type0 gaussian(type0,type0);
    };
}
#endif 
