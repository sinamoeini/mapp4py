#ifndef __MAPP__potfit_funcs__
#define __MAPP__potfit_funcs__
#include "global.h"
#include "xmath.h"
#include "api.h"
namespace MAPP_NS
{
    class PotFitAux
    {
    private:
    protected:
    public:
        static type0 find_max_alpha(const type0,const type0,bool,type0,type0,type0);
        static type0 find_max_alpha(const type0,const type0,type0,type0,type0);
    };
    class PotFitPairFunc
    {
    private:
    protected:
    public:
        PotFitPairFunc(const char*);
        virtual ~PotFitPairFunc();
        const char* name;
        type0* vars;
        type0* dvars_lcl;
        type0 rc;
        size_t nvars;
        virtual type0 F(type0)=0;
        virtual type0 dF(type0)=0;
        virtual type0 ddF(type0)=0;
        virtual void DF(type0,type0,type0*)=0;
        virtual void DdF(type0,type0,type0*)=0;
        virtual void find_max_alpha(type0&,type0*,type0*)=0;
        virtual int set_init(PyObject*,type0*&,size_t&)=0;
        virtual int set(PyObject*);
        virtual PyObject* get();
        void DF(type0 coef,type0 r){DF(coef,r,dvars_lcl);};
        void DdF(type0 coef,type0 r){DdF(coef,r,dvars_lcl);};
    };
    
    class PotFitEmbFunc
    {
    private:
    protected:
    public:
        PotFitEmbFunc(const char*);
        virtual ~PotFitEmbFunc();
        const char* name;
        type0* vars;
        type0* dvars_lcl;
        size_t nvars;
        virtual type0 F(type0)=0;
        virtual type0 dF(type0)=0;
        virtual type0 ddF(type0)=0;
        virtual void DF(type0,type0*)=0;
        virtual void DdF(type0,type0*)=0;
        virtual void find_max_alpha(type0&,type0*,type0*)=0;
        virtual int set_init(PyObject*,type0*&,size_t&)=0;
        virtual int set(PyObject*);
        virtual PyObject* get();
        void DF(type0 rho){DF(rho,dvars_lcl);};
        void DdF(type0 rho){DdF(rho,dvars_lcl);};
    };
    
    
    class PotFitRhoSpl: public PotFitPairFunc
    {
    private:
        type0* R;
        int quadratic(type0*,type0*);
        int cubic(type0*,type0*);
        int quartic(type0*,type0*);
        

        void slpine_max_alpha(type0 (&)[4],type0 (&)[4],const type0&,const type0&,type0&);
        void slpine_max_alpha(type0*,type0 (&)[4],type0 (&)[4],const type0&,const type0&,type0&);

        type0 F(type0*,type0);
        type0 F(type0,type0*,type0);
    protected:
    public:
        static type0 sort(type0*,type0*,size_t);
        PotFitRhoSpl(const char*);
        ~PotFitRhoSpl();
        type0 F(type0);
        type0 dF(type0);
        type0 ddF(type0);
        void DF(type0,type0,type0*);
        void DdF(type0,type0,type0*);
        void find_max_alpha(type0&,type0*,type0*);
        int set_init(PyObject*,type0*&,size_t&);
        int set(PyObject*);
    };
    
    
    class PotFitRhoAng: public PotFitPairFunc
    {
    private:
    protected:
    public:
        PotFitRhoAng(const char*);
        ~PotFitRhoAng();
        type0 F(type0);
        type0 dF(type0);
        type0 ddF(type0);
        void DF(type0,type0,type0*);
        void DdF(type0,type0,type0*);
        void find_max_alpha(type0&,type0*,type0*);
        int set_init(PyObject*,type0*&,size_t&);
        int set(PyObject*);
    };

    
    class PotFitPhiSpl: public PotFitPairFunc
    {
    private:
        type0* R;
        type0 F(type0*,type0);
    protected:
    public:
        
        PotFitPhiSpl(const char*);
        ~PotFitPhiSpl();
        type0 F(type0);
        type0 dF(type0);
        type0 ddF(type0);
        void DF(type0,type0,type0*);
        void DdF(type0,type0,type0*);
        void find_max_alpha(type0&,type0*,type0*);
        int set_init(PyObject*,type0*&,size_t&);
        int set(PyObject*);
    };
    
    
    class PotFitEmbAck:public PotFitEmbFunc
    {
    private:
    protected:
    public:
        PotFitEmbAck(const char*);
        ~PotFitEmbAck();
        type0 F(type0);
        type0 dF(type0);
        type0 ddF(type0);
        void DF(type0,type0*);
        void DdF(type0,type0*);
        void find_max_alpha(type0&,type0*,type0*);
        int set_init(PyObject*,type0*&,size_t&);
        int set(PyObject*);
    };
    
    class PotFitEmbAng:public PotFitEmbFunc
    {
    private:
    protected:
    public:
        PotFitEmbAng(const char*);
        ~PotFitEmbAng();
        type0 F(type0);
        type0 dF(type0);
        type0 ddF(type0);
        void DF(type0,type0*);
        void DdF(type0,type0*);
        void find_max_alpha(type0&,type0*,type0*);
        int set_init(PyObject*,type0*&,size_t&);
        int set(PyObject*);
    };
}
#endif
