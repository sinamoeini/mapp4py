#ifndef __MAPP__dynamic_md__
#define __MAPP__dynamic_md__

#include "dynamic.h"
#ifdef NEW_UPDATE
#include "atoms_md.h"
#endif

namespace MAPP_NS
{
    class vec;
    template<typename> class Vec;
    class DynamicMD: public Dynamic
    {
    private:
        void store_x0();
        bool decide();
        class ForceFieldMD* ff;
        class AtomsMD* atoms;
#ifdef NEW_UPDATE
        void neighboring();
#endif
    protected:
        Vec<type0>* x0;
    public:        
        DynamicMD(class AtomsMD*,class ForceFieldMD*,bool,
                  std::initializer_list<vec*>,std::initializer_list<vec*>,std::initializer_list<vec*> = {});
        ~DynamicMD();
    
        void init_xchng();
        void fin_xchng();
        void init();
        void fin();
#ifdef NEW_UPDATE

        
        
        template<class...VS>
        void update_wo_x(VS*&... __vs)
        {
            updt->update_wo_x(__vs...);
        }
        template<class...VS>
        void update_w_x(VS*&... __vs)
        {
            if(chng_box)
            {
                updt->update_w_x(__vs...);
                
                if(decide()) return;
                
                atoms->x2s_lcl();
                xchng->full_xchng();
                
                updt->reset();
                updt->list();
                neighboring();
                store_x0();
            }
            else
            {
                if(decide())
                {
                    updt->update_w_x(__vs...);
                    return;
                }
                
                atoms->x2s_lcl();
                xchng->full_xchng();
                
                updt->list();
                neighboring();
                
                store_x0();
            }
        }
        
        template<bool X>
        class Helper
        {
        public:
            template<class...VS>
            static void update(DynamicMD& dynamic,VS*&... __vs)
            {
                dynamic.update_w_x(__vs...);
            }
            
        };
        
        
        template<bool X=false,class...VS>
        void update(VS*&... __vs)
        {
            Helper<X>::update(*this,__vs...);
        }
        
        
#else
        void update(vec**,int);
        void update(vec*);
#endif
    };
#ifdef NEW_UPDATE
    template<>
    class DynamicMD::Helper<false>
    {
    public:
        template<class...VS>
        static void update(DynamicMD& dynamic,VS*&... __vs)
        {
            dynamic.update_wo_x(__vs...);
        }
        
    };
#endif
    
}
#endif
