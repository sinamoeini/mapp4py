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
    public:
        int nevery;
        Export(std::initializer_list<const char*>,int,std::string*,size_t);
        void add_to_default(const char*);
        virtual ~Export();
        static void gather(class vec**,int);
        static void release(class vec**,int);
        void find_vecs(class Atoms*);
        
        static bool open(const int,MPI_Comm&,const char*,const char*,FILE*&);
        static void close(const int,MPI_Comm&,FILE*&);
        
        
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



namespace MAPP_NS
{
    class ExportMD: public Export
    {
    private:
        
    protected:
    public:
        class AtomsMD* atoms;
        
        ExportMD(std::initializer_list<const char*>,int,std::string*,size_t);
        virtual ~ExportMD();
        virtual void init(){};
        virtual void write(int){};
        virtual void write(const char*){};
        virtual void fin(){};
        
        typedef struct
        {
            PyObject_HEAD
            ExportMD* xprt;
        }Object;
        
        static PyTypeObject TypeObject;
        static PyObject* __new__(PyTypeObject*,PyObject*, PyObject*);
        static int __init__(PyObject*, PyObject*,PyObject*);
        static PyObject* __alloc__(PyTypeObject*,Py_ssize_t);
        static void __dealloc__(PyObject*);
        static PyObject* __call__(PyObject*,PyObject*,PyObject*);
        
        static PyMethodDef methods[];
        static void setup_tp_methods();
        
        static PyGetSetDef getset[];
        static void setup_tp_getset();
        
        static int setup_tp();
    };
}
namespace MAPP_NS
{
    class ExportDMD: public Export
    {
    private:
        
    protected:
    public:
        class AtomsDMD* atoms;
        
        ExportDMD(std::initializer_list<const char*>,int,std::string*,size_t);
        virtual ~ExportDMD();
        virtual void init(){};
        virtual void write(int){};
        virtual void write(const char*){};
        virtual void fin(){};
        
        typedef struct
        {
            PyObject_HEAD
            ExportDMD* xprt;
        }Object;
        
        static PyTypeObject TypeObject;
        static PyObject* __new__(PyTypeObject*,PyObject*, PyObject*);
        static int __init__(PyObject*, PyObject*,PyObject*);
        static PyObject* __alloc__(PyTypeObject*,Py_ssize_t);
        static void __dealloc__(PyObject*);
        static PyObject* __call__(PyObject*,PyObject*,PyObject*);
        
        static PyMethodDef methods[];
        static void setup_tp_methods();
        
        static PyGetSetDef getset[];
        static void setup_tp_getset();
        
        static int setup_tp();
    };
}

#endif 
