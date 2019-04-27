#ifndef __MAPP__potfit_funcs__
#define __MAPP__potfit_funcs__
#include "api.h"
#include "global.h"
#include "xmath.h"
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
        std::string A_name;
        std::string dA_name_max;
        std::string A_name_dof;
        
        type0* vars;
        type0* hvars;
        type0* dvars_max;
        type0* dvars_lcl;
        type0* dvars;
        bool* dofs;
        type0 rc;
        size_t nvars;
        virtual type0 F(type0)=0;
        virtual type0 dF(type0)=0;
        virtual type0 ddF(type0)=0;
        virtual void DF(type0,type0,type0*)=0;
        virtual void DdF(type0,type0,type0*)=0;
        virtual void find_max_alpha(type0&)=0;
        virtual void random_neigh(class Random*);
        virtual bool validate()=0;
        virtual int set_init(PyObject*,type0*&,size_t&)=0;
        
        void DF(type0 coef,type0 r){DF(coef,r,dvars_lcl);};
        void DdF(type0 coef,type0 r){DdF(coef,r,dvars_lcl);};
        
        virtual int set_A(PyObject*);
        virtual PyObject* get_A();
        virtual PyObject* get_dA();
        virtual int set_dA_max(PyObject*);
        virtual PyObject* get_dA_max();
        virtual int set_A_dof(PyObject*);
        virtual PyObject* get_A_dof();
    };
    
    class PotFitEmbFunc
    {
    private:
    protected:
    public:
        PotFitEmbFunc(const char*);
        virtual ~PotFitEmbFunc();
        const char* name;
        std::string A_name;
        std::string dA_name_max;
        std::string A_name_dof;
        
        type0* vars;
        type0* hvars;
        type0* dvars_max;
        type0* dvars_lcl;
        type0* dvars;
        bool* dofs;
        size_t nvars;
        virtual type0 F(type0)=0;
        virtual type0 dF(type0)=0;
        virtual type0 ddF(type0)=0;
        virtual void DF(type0,type0,type0*)=0;
        virtual void DdF(type0,type0,type0*)=0;
        virtual void find_max_alpha(type0&)=0;
        virtual void random_neigh(class Random*);
        virtual bool validate()=0;
        virtual int set_init(PyObject*,type0*&,size_t&)=0;
        void DF(type0 rho){DF(1.0,rho,dvars_lcl);};
        void DdF(type0 rho){DdF(1.0,rho,dvars_lcl);};
        void DF(type0 coef,type0 rho){DF(coef,rho,dvars_lcl);};
        void DdF(type0 coef,type0 rho){DdF(coef,rho,dvars_lcl);};
        
        virtual int set_A(PyObject*);
        virtual PyObject* get_A();
        virtual PyObject* get_dA();
        virtual int set_dA_max(PyObject*);
        virtual PyObject* get_dA_max();
        virtual int set_A_dof(PyObject*);
        virtual PyObject* get_A_dof();
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
        void find_max_alpha(type0&);
        void random_neigh(class Random*);
        bool validate();
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
        void find_max_alpha(type0&);
        bool validate();
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
        void find_max_alpha(type0&);
        void random_neigh(class Random*);
        bool validate();
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
        void DF(type0,type0,type0*);
        void DdF(type0,type0,type0*);
        void find_max_alpha(type0&);
        bool validate();
        int set_init(PyObject*,type0*&,size_t&);
        int set(PyObject*);
    };
    
    class PotFitEmbAng:public PotFitEmbFunc
    {
    private:
        type0 lim;
    protected:
    public:
        PotFitEmbAng(const char*);
        ~PotFitEmbAng();
        type0 F(type0);
        type0 dF(type0);
        type0 ddF(type0);
        void DF(type0,type0,type0*);
        void DdF(type0,type0,type0*);
        void find_max_alpha(type0&);
        
        bool validate();
        int set_init(PyObject*,type0*&,size_t&);
        int set(PyObject*);
    };
}
#endif
