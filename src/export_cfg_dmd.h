#ifndef __MAPP__export_cfg_dmd__
#define __MAPP__export_cfg_dmd__
#include "export.h"

namespace MAPP_NS
{
    class ExportCFGDMD:public ExportDMD
    {
    private:
        int max_pos(type0*& c,const int& nelems)
        {
            type0 max_c=*c;
            int ielem=0;
            for(int i=1;i<nelems;i++)
                if(max_c<c[i])
                {
                    max_c=c[i];
                    ielem=i;
                }
            
            return ielem;
        }
        
        std::string pattern;
        bool sort;
        
        void write_header(FILE*);
        void write_body_sort(FILE*);
        void write_body(FILE*);
    protected:
    public:
        ExportCFGDMD(const std::string&,int,std::string*,size_t,bool);
        ~ExportCFGDMD();

        void init();
        void write(int);
        void write(const char*);
        void fin();
        
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
        
        static int setup_tp();
    };
}
#endif 
