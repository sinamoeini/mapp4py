#ifndef __MAPP__thermo_dynamics__
#define __MAPP__thermo_dynamics__
#include "global.h"
#include <string>
#include "MAPP.h"

namespace MAPP_NS
{
    class ThermoQuantity
    {
    private:
    protected:
    public:
        
        ThermoQuantity();
        ThermoQuantity(const char*,type0&);
        ThermoQuantity(const ThermoQuantity&);
        ThermoQuantity(ThermoQuantity&&);
        ~ThermoQuantity();
        ThermoQuantity& operator =(const ThermoQuantity&);
        ThermoQuantity& operator =(ThermoQuantity&&);

        const char* name;
        const int name_lngth;
        type0 const * const ptr;
        void print_header(int);
    };
}

namespace MAPP_NS
{
    class ThermoDynamics
    {
    private:
    protected:
        const int precision;
        const int l_int;
        const int l_double;
        
        ThermoQuantity* qs;
        int nqs;
        
        template<class ...Ts>
        void assign(ThermoQuantity* __qs,const char* name,type0& val,Ts&... vs)
        {
            assign(__qs,name,val);
            assign(__qs+1,vs...);
        }
        
        void assign(ThermoQuantity* __qs,const char* name,type0& val)
        {
            __qs->~ThermoQuantity();
            new (__qs) ThermoQuantity(name,val);
        }
        
        
        
    public:
        ThermoDynamics(const int);
        template<class ...Ts>
        ThermoDynamics(const int __precision,const char* name,Ts&... vs):
        ThermoDynamics(__precision)
        {
            nqs=static_cast<int>((sizeof...(vs)+1)/2);
            qs=new ThermoQuantity[nqs];
            assign(qs,name,vs...);
            
        }
        
        ~ThermoDynamics();

        void print(int);
        void init();
        void fin();
        
        
        template<int N>
        static void print_header(const char (&name)[N],int lngth)
        {
            if(N-1<=lngth)
                MAPP::print_stdout("%*s%*s",(lngth-lngth/2)+(N-1)/2,name,lngth/2-(N-1)/2,"");
            else
                for(int i=0;i<lngth;i++) MAPP::print_stdout("%c",name[i]);
        }
        
    };
}


#endif
