#ifndef __MAPP__write__
#define __MAPP__write__
#include <iostream>
namespace MAPP_NS
{
    class Write
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
        Write(class Atoms*);
        Write(std::initializer_list<const char*>,std::string*,size_t);
        void add_to_default(const char*);
        ~Write();
        class Atoms* atoms;
        void init();
        void fin();
        static void gather(class vec**,int);
        static void release(class vec**,int);
    };
}

#endif 
