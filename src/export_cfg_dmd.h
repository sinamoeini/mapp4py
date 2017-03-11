#ifndef __MAPP__export_cfg_dmd__
#define __MAPP__export_cfg_dmd__
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
        ExportCFGDMD(const std::string&,int,std::string*,size_t,bool);
        ~ExportCFGDMD();
        void write_header(FILE*);
        void write_body_sort(FILE*);
        void write_body(FILE*);
        void write(int);
        void init();
        void fin();
        class AtomsDMD* atoms;
        
        typedef struct
        {
            PyObject_HEAD
            ExportCFGDMD* xprt;
        }Object;
        
        static PyTypeObject TypeObject;
        static PyObject* __new__(PyTypeObject*,PyObject*, PyObject*);
        static int __init__(PyObject*, PyObject*,PyObject*);
        static PyObject* __alloc__(PyTypeObject*,Py_ssize_t);
        static void __dealloc__(PyObject*);
        
        static PyMethodDef methods[];
        static void setup_tp_methods();
        
        static PyGetSetDef getset[];
        static void setup_tp_getset();
        
        static void setup_tp();
    };
}

#endif 
