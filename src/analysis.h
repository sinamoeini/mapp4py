#ifndef __MAPP__analysis__
#define __MAPP__analysis__
#include "global.h"
typedef struct PyMethodDef PyMethodDef;
namespace MAPP_NS
{
    class BCTPolarity
    {
    private:
    protected:
    public:
        static void ml_bct_polarity(PyMethodDef&);
        static type0* polarity(class AtomsMD*,type0**&&,int);
        static int calc_polarity(type0*,int,int*,int,type0*);
        static type0 test(type0 (&)[3][3],type0(&)[3]);
        static void QR(type0(&)[3][3],type0(&)[3][3],type0(&)[3][3]);
    };
}
#endif
