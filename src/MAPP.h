#ifndef __MAPP__MAPP__
#define __MAPP__MAPP__
#include <stdio.h>
struct _object;
typedef _object PyObject;
typedef struct PyMethodDef PyMethodDef;
typedef struct PyModuleDef PyModuleDef;
namespace MAPP_NS
{
    class MAPP
    {
    private:
        static const char* devnull_path;
        static int glbl_rank;
    protected:
    public:
        static PyModuleDef module;
        static PyObject* init_module(void);
        
        static FILE* __stdout__;
        static FILE* __stderr__;
        static FILE* __devnull__;
        
        static FILE* mapp_out;
        static FILE* mapp_err;
        
        static PyObject* pause_out(PyObject* =NULL);
        static PyObject* resume_out(PyObject* =NULL);
        
        static PyMethodDef methods[];
        static void setup_methods();
        
        class MD
        {
        private:
        protected:
        public:
            static PyModuleDef module;
            static PyObject* init_module(void);
            static PyMethodDef methods[];
            static void setup_methods();
        };
        
        class DMD
        {
        private:
        protected:
        public:
            static PyModuleDef module;
            static PyObject* init_module(void);
            static PyMethodDef methods[];
            static void setup_methods();
        };
    };
}
#endif 
