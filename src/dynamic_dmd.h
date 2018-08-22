#ifndef __MAPP__dynamic_dmd__
#define __MAPP__dynamic_dmd__

#include "dynamic.h"
#ifdef OLD_UPDATE
#else
#include "atoms_dmd.h"
#endif

namespace MAPP_NS
{
    class vec;
    template<typename> class Vec;
    class DynamicDMD: public Dynamic
    {
    private:
        const int c_dim;
        const type0 alpha_scale;
        void store_x0();
        bool decide();
        class ForceFieldDMD* ff;
        class AtomsDMD* atoms;
#ifdef OLD_UPDATE
#else
        void neighboring();
        void reset_max_cut();
        template<bool,bool>
        bool __decide();
#endif
    protected:
        
        Vec<type0>* x0;
        Vec<type0>* alpha0;
    
    public:
        DynamicDMD(class AtomsDMD*,class ForceFieldDMD*,bool,
        std::initializer_list<vec*>,std::initializer_list<vec*>,std::initializer_list<vec*> = {});
        ~DynamicDMD();
#ifdef OLD_UPDATE
        void update(vec**,int);
        void update(vec*);
        void update(vec*,type0 (*)[__dim__]);
#else
        template<class...VS>
        void update_wo_x_wo_alpha(VS*&... __vs)
        {
            updt->update_wo_x(__vs...);
        }
        
        template<class...VS>
        void update_w_x_w_alpha(VS*&... __vs)
        {
            if(chng_box)
            {
                updt->update_w_x(atoms->alpha,__vs...);
                
                if(decide()) return;
                
                atoms->x2s_lcl();
                xchng->full_xchng();
                reset_max_cut();
                updt->reset();
                updt->list();
                neighboring();
                store_x0();
            }
            else
            {
                if(decide())
                {
                    updt->update_w_x(atoms->alpha,__vs...);
                    return;
                }
                
                atoms->x2s_lcl();
                xchng->full_xchng();
                reset_max_cut();
                updt->reset();
                updt->list();
                neighboring();
                store_x0();
            }
                
        }
        /*
         this is a dummy class to overcome the fact that we
         cannot do partialy specialize class memebers
         */
        template<bool X,bool ALPHA>
        class Helper
        {
        public:
            template<class...VS>
            static void update(DynamicDMD& dynamic,VS*&... __vs)
            {
                dynamic.update_w_x_w_alpha(__vs...);
            }
            
        };
        
        
        template<bool X=false,bool ALPHA=false,class...VS>
        void update(VS*&... __vs)
        {
            Helper<X,ALPHA>::update(*this,__vs...);
        }
        
        
        template<class...VS>
        void update(type0 (*&__dH)[__dim__],Vec<type0>*& __v,VS*&... __vs)
        {
            updt->update_w_x_w_dH(__dH,__v,__vs...);
        }
#endif
        void init();
        void fin();
    };
}

#ifdef OLD_UPDATE
#else
template<>
class DynamicDMD::Helper<false,false>
{
public:
    template<class...VS>
    static void update(DynamicDMD& dynamic,VS*&... __vs)
    {
        dynamic.update_wo_x_wo_alpha(__vs...);
    }
    
};
template<> bool DynamicDMD::__decide<true,true>();
template<> bool DynamicDMD::__decide<true,false>();
template<> bool DynamicDMD::__decide<false,true>();
#endif

#endif
