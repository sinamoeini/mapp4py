#ifndef __MAPP__export__
#define __MAPP__export__
#include "api.h"
#include <iostream>
namespace MAPP_NS
{
    class Export
    {
    private:
        
    protected:
        std::string* vec_names;
        class vec** vecs;
        int nvecs;
        int ndef_vecs;
        int nusr_vecs;
        int ndims;
        int nevery;
    public:
        Export(int,std::initializer_list<const char*>,std::string*,size_t);
        void add_to_default(const char*);
        ~Export();
        class Atoms* atoms;
        static void gather(class vec**,int);
        static void release(class vec**,int);
        void find_vecs();
        
        typedef struct
        {
            PyObject_HEAD
            Export* xprt;
        }Object;
        
        static void getset_deafult_vecs(PyGetSetDef&);
        static void getset_extra_vecs(PyGetSetDef&);
        static void getset_nevery(PyGetSetDef&);
    };
}

#endif 
