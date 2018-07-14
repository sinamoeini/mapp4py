#if PY_MAJOR_VERSION >= 3
#define IS_PY3K
#endif

#ifndef Py_TYPE
#define Py_TYPE(ob) (((PyObject*)(ob))->ob_type)
#endif



#ifdef IS_PY3K
#define PyLong_FromLong PyLong_FromLong
#define PyLong_FromSize_t PyLong_FromSize_t
#define GET_FILE_FROM_PY(file_name) file_name
#define MOD_INIT(name,func) PyMODINIT_FUNC PyInit_##name(void){PyObject* __module__=func; return __module__;}

#define __PyString_AsString(v) PyUnicode_AsUTF8AndSize(v,NULL)
#define __PyString_FromString(v) PyUnicode_FromString(v)

#else

#define PyLong_FromLong PyInt_FromLong
#define PyLong_FromSize_t PyInt_FromSize_t
#define GET_FILE_FROM_PY(file_name) reinterpret_cast<PyFileObject*>(PySys_GetObject((char*)#file_name))->f_fp
#define MOD_INIT(name,func) PyMODINIT_FUNC init##name(void){func;}

#define __PyString_AsString(v) PyString_AsString(v)
#define __PyString_FromString(v) PyString_FromString(v)
#endif

