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
        
#ifdef NEW_UPDATE
        bool decide()
        {
            type0 skin_sq=0.25*skin*skin;
            int succ,succ_lcl=1;
            type0* x_vec=atoms->x->begin();
            type0* x0_vec=x0->begin();
            int last_atm=atoms->natms_lcl;
            if(chng_box) last_atm+=atoms->natms_ph;
            for(int iatm=0;succ_lcl && iatm<last_atm;iatm++,x0_vec+=__dim__,x_vec+=__dim__)
                if(Algebra::RSQ<__dim__>(x0_vec,x_vec)>skin_sq)
                    succ_lcl=0;
            
            MPI_Allreduce(&succ_lcl,&succ,1,MPI_INT,MPI_MIN,world);
            if(succ) return true;
            return false;
        }
#else
        bool decide();
#endif
        class ForceFieldMD* ff;
        class AtomsMD* atoms;
#ifdef NEW_UPDATE
        class __Update* __updt;
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
        void update_w_o_x(VS* ...__vs)
        {
            __updt->update_w_o_x(__vs...);
        }
        template<class...VS>
        void update_w_x(VS* ...__vs)
        {
            if(chng_box)
            {
                __updt->update_w_x(atoms->x,__vs...);
                
                if(decide()) return;
                
                atoms->x2s_lcl();
                xchng->full_xchng();
                
                __updt->reset();
                __updt->list();
                neighboring();
                store_x0();
            }
            else
            {
                if(decide())
                {
                    __updt->update_w_x(atoms->x,__vs...);
                    return;
                }
                
                atoms->x2s_lcl();
                xchng->full_xchng();
                
                __updt->list();
                neighboring();
                
                store_x0();
            }
        }
        
        template<bool X>
        class Helper
        {
        public:
            template<class...VS>
            static void update(DynamicMD& dynamic,VS* ...__vs)
            {
                dynamic.update_w_x(__vs...);
            }
            
        };
        
        
        template<bool X=false,class...VS>
        void update(VS* ...__vs)
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
        static void update(DynamicMD& dynamic,VS* ...__vs)
        {
            dynamic.update_w_o_x(__vs...);
        }
        
    };
#endif
    
}
#endif
