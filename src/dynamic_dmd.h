#ifndef __MAPP__dynamic_dmd__
#define __MAPP__dynamic_dmd__

#include "dynamic.h"
#include "atoms_dmd.h"
#include "ff_styles.h"
#include "neighbor_dmd.h"

namespace MAPP_NS
{
    class vec;
    template<typename> class Vec;
    class DynamicDMD: public Dynamic
    {
    friend class DAEOld;
    private:
        const int c_dim;
        const type0 alpha_scale;
        void store_x0();
        bool decide();
        class ForceFieldDMD* ff;
        class AtomsDMD* atoms;
    protected:
        
        Vec<type0>* x0;
        Vec<type0>* alpha0;
    
    public:
        DynamicDMD(class AtomsDMD*,class ForceFieldDMD*,bool,
        std::initializer_list<vec*>,std::initializer_list<vec*>,std::initializer_list<vec*> = {});
        ~DynamicDMD();
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
                atoms->max_cut=ff->max_cut+atoms->comm.skin+alpha_scale*sqrt_2*atoms->max_alpha;
                updt->reset();
                updt->list();
                ff->neighbor->create_list(true);
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
                atoms->max_cut=ff->max_cut+atoms->comm.skin+alpha_scale*sqrt_2*atoms->max_alpha;
                updt->reset();
                updt->list();
                ff->neighbor->create_list(true);
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
        void init();
        void fin();
    };
}

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



namespace MAPP_NS
{
    class vec;
    template<typename> class Vec;
    
    template<bool BC,bool X,bool ALPHA>
    class NewDynamicDMD
    {
    friend class DAEOld;
    private:
        NewDynamic dynamic_sub;
        class ForceFieldDMD* ff;
        class AtomsDMD* atoms;
        const int c_dim;
        const type0 alpha_scale;
        Vec<type0>* x0;
        Vec<type0>* alpha0;
        
        bool decide();
        void store_x0_alpha0();
        void alloc_x0_alpha0();
    protected:
    public:
  
        NewDynamicDMD(AtomsDMD* __atoms,ForceFieldDMD* __ff,
        std::initializer_list<vec*> __updt_vecs,
        std::initializer_list<vec*> __xchng_comp_vecs,
        std::initializer_list<vec*> __arch_vecs):
        dynamic_sub(__atoms,__ff,
        {__atoms->x,__atoms->alpha,__atoms->c,__atoms->elem},__updt_vecs,
        {__atoms->id},__xchng_comp_vecs,
        {},__arch_vecs),
        ff(__ff),
        atoms(__atoms),
        c_dim(__atoms->c->dim),
        alpha_scale(__atoms->xi[__atoms->N-1]),
        x0(NULL),
        alpha0(NULL)
        {
        }
        ~NewDynamicDMD()
        {
        }

        void init()
        {
            
            dynamic_sub.store_arch_vecs();
            dynamic_sub.create_dynamic_vecs();
            alloc_x0_alpha0();
            ff->init();
            ff->neighbor->init();
            atoms->max_cut=ff->max_cut+atoms->comm.skin+alpha_scale*sqrt_2*atoms->max_alpha;
            
            dynamic_sub.create_updt_xchng();
            
            ff->updt=dynamic_sub.updt;
            atoms->x2s_lcl();
            dynamic_sub.xchng->full_xchng();
            dynamic_sub.updt->reset();
            dynamic_sub.updt->list();
            ff->neighbor->create_list(true);
            store_x0_alpha0();
        }
        void fin()
        {
            ff->neighbor->fin();
            ff->fin();
            
            dynamic_sub.destroy_updt_xchng();
            
            delete alpha0;
            alpha0=NULL;
            delete x0;
            x0=NULL;
            
            dynamic_sub.restore_arch_vecs();
            dynamic_sub.destroy_dynamic_vecs();
            dynamic_sub.shrink_to_fit_all();
        }
        
        template<class...VS>
        void update(VS*&... __vs)
        {
            
        }
        void reset(){};
        
        
        void add_xchng(vec* v){dynamic_sub.add_xchng(v);};
        void add_updt(vec* v){dynamic_sub.add_updt(v);};
        
    };

}
using namespace MAPP_NS;
/*--------------------------------------------
 
 --------------------------------------------*/
template<>
template<class...VS>
void NewDynamicDMD<true,true,true>::update(VS*&... __vs)
{
    atoms->update_max_alpha();
    atoms->update_H();
    
    dynamic_sub.updt->update_w_x(atoms->alpha,__vs...);
    
    if(decide()) return;
    reset();
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<>
template<class...VS>
void NewDynamicDMD<false,true,true>::update(VS*&... __vs)
{
    atoms->update_max_alpha();
    
    if(decide())
    {
        dynamic_sub.updt->update_w_x(atoms->alpha,__vs...);
        return;
    }
    reset();
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<>
template<class...VS>
void NewDynamicDMD<true,true,false>::update(VS*&... __vs)
{
    atoms->update_H();
    
    dynamic_sub.updt->update_w_x(__vs...);
    
    if(decide()) return;
    reset();
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<>
template<class...VS>
void NewDynamicDMD<false,true,false>::update(VS*&... __vs)
{
    if(decide())
    {
        dynamic_sub.updt->update_w_x(__vs...);
        return;
    }
    reset();
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<>
template<class...VS>
void NewDynamicDMD<false,false,true>::update(VS*&... __vs)
{
    if(decide())
    {
        dynamic_sub.updt->update_wo_x(atoms->alpha,__vs...);
        return;
    }
    reset();
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<>
template<class...VS>
void NewDynamicDMD<false,false,false>::update(VS*&... __vs)
{
    dynamic_sub.updt->update_wo_x(__vs...);
}




#endif
