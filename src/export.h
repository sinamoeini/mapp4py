#ifndef __MAPP__export__
#define __MAPP__export__
#include <iostream>
namespace MAPP_NS
{
    class Export
    {
    private:
        
    protected:
        char** vec_names;
        class vec** vecs;
        int nvecs;
        int ndef_vecs;
        int nusr_vecs;
        int ndims;
    public:
        Export(std::initializer_list<const char*>,std::string*,size_t);
        void add_to_default(const char*);
        ~Export();
        class Atoms* atoms;
        void init();
        void fin();
        static void gather(class vec**,int);
        static void release(class vec**,int);
    };
}

#endif 
