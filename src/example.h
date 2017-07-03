#ifndef __MAPP__example__
#define __MAPP__example__
namespace MAPP_NS
{
    
    class ExamplePython
    {
    public:
        typedef struct
        {
            PyObject_HEAD
            PyObject* mass;
            
        }Object;
        static PyTypeObject TypeObject;
        
        static int set_mass(PyObject*,PyObject*,void*);
        static PyObject* get_mass(PyObject* self,void*);
        static PyObject* tp_new(PyTypeObject*,PyObject*, PyObject*);
        static int tp_init(PyObject*, PyObject*,PyObject*);
        static void tp_dealloc(PyObject*);
        
        static PyObject* test_func(PyObject*,PyObject*,PyObject*);
        
        
        static PyMemberDef tp_members[];
        static void setup_tp_members();
        
        static PyGetSetDef tp_getset[];
        static void setup_tp_getset();
        
        static PyMethodDef tp_methods[];
        static void setup_tp_methods();
        static int setup_tp();

        static void ml_test(PyMethodDef&);
        static void ml_phonon(PyMethodDef&);
        static void ml_phonon_1d(PyMethodDef&);
        static void ml_phonon_1dd(PyMethodDef&);
        static void ml_alpha(PyMethodDef&);
        static void ml_delta_c(PyMethodDef&);
        static void ml_mv_c(PyMethodDef&);
        static void ml_mv_c2(PyMethodDef&);
        static void ml_prt(PyMethodDef&);
    };
    
    
}

#endif
