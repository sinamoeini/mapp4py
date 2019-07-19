#ifndef __MAPP__dynamic_dmd__
#define __MAPP__dynamic_dmd__

#include "dynamic.h"
#ifdef OLD_UPDATE
#else
#include "atoms_dmd.h"
#include "ff_styles.h"
#include "neighbor_dmd.h"
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

#endif





#ifdef OLD_UPDATE
#else





namespace MAPP_NS
{
    class vec;
    template<typename> class Vec;
    
    template<bool BC,bool X,bool ALPHA>
    class NewDynamicDMD: public NewDynamic
    {
    private:
        const int c_dim;
        const type0 alpha_scale;
        class ForceFieldDMD* ff;
        class AtomsDMD* atoms;


        bool decide();
        void store_x0_alpha0();
        void alloc_x0_alpha0();
    protected:
        Vec<type0>* x0;
        Vec<type0>* alpha0;
    
    public:
  
        NewDynamicDMD(AtomsDMD* __atoms,ForceFieldDMD* __ff,
        std::initializer_list<vec*> __updt_vecs,
        std::initializer_list<vec*> __xchng_comp_vecs,
        std::initializer_list<vec*> __arch_vecs):
        NewDynamic(__atoms,__ff,
        {__atoms->x,__atoms->alpha,__atoms->c,__atoms->elem},__updt_vecs,
        {__atoms->id},__xchng_comp_vecs,
        {},__arch_vecs),
        c_dim(__atoms->c->dim),
        alpha_scale(__atoms->xi[__atoms->N-1]),
        ff(__ff),
        atoms(__atoms),
        x0(NULL),
        alpha0(NULL)
        {
        }
        ~NewDynamicDMD()
        {
        }

        void init()
        {
            
            store_arch_vecs();
            create_dynamic_vecs();
            alloc_x0_alpha0();
            ff->init();
            ff->neighbor->init();
            atoms->max_cut=ff->max_cut+atoms->comm.skin+alpha_scale*sqrt_2*atoms->max_alpha;
            
            xchng=new Exchange(atoms,nxchng_vecs_full);
            updt=new Update(atoms,nupdt_vecs_full,nxchng_vecs_full);
            
            ff->updt=updt;
            atoms->x2s_lcl();
            xchng->full_xchng();
            updt->reset();
            updt->list();
            ff->neighbor->create_list(true);
            store_x0_alpha0();
        }
        void fin()
        {
            ff->neighbor->fin();
            ff->fin();
            
            delete updt;
            updt=NULL;
            delete xchng;
            xchng=NULL;
            delete alpha0;
            alpha0=NULL;
            delete x0;
            x0=NULL;
            
            restore_arch_vecs();
            delete [] atoms->dynamic_vecs;
            atoms->dynamic_vecs=NULL;
            atoms->ndynamic_vecs=0;
            
            destroy_dynamic_vecs();
            restore_arch_vecs();
            
            for(int ivec=0;ivec<atoms->nvecs;ivec++)
                if(!atoms->vecs[ivec]->is_empty())
                {
                    atoms->vecs[ivec]->vec_sz=atoms->natms_lcl;
                    atoms->vecs[ivec]->shrink_to_fit();
                }
            atoms->natms_ph=0;
        }
        
        template<class...VS>
        void update(VS*&... __vs)
        {
            
        }
        
    };

}
using namespace MAPP_NS;


template<>
template<class...VS>
void NewDynamicDMD<true,true,true>::update(VS*&... __vs)
{
    atoms->update_max_alpha();
    atoms->update_H();
    
    updt->update_w_x(atoms->alpha,__vs...);
    
    if(decide()) return;
    
    atoms->x2s_lcl();
    xchng->full_xchng();
    atoms->max_cut=ff->max_cut+atoms->comm.skin+alpha_scale*sqrt_2*atoms->max_alpha;
    updt->reset();
    updt->list();
    ff->neighbor->create_list(true);
    store_x0_alpha0();
}

template<>
template<class...VS>
void NewDynamicDMD<false,true,true>::update(VS*&... __vs)
{
    atoms->update_max_alpha();
    
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
    store_x0_alpha0();
}

template<>
template<class...VS>
void NewDynamicDMD<true,true,false>::update(VS*&... __vs)
{
    atoms->update_H();
    
    updt->update_w_x(__vs...);
    
    if(decide()) return;
    
    
    
    atoms->x2s_lcl();
    xchng->full_xchng();
    updt->reset();
    updt->list();
    ff->neighbor->create_list(true);
    store_x0_alpha0();
}

template<>
template<class...VS>
void NewDynamicDMD<false,true,false>::update(VS*&... __vs)
{
    if(decide())
    {
        updt->update_w_x(__vs...);
        return;
    }
    
    
    atoms->x2s_lcl();
    xchng->full_xchng();
    updt->list();
    ff->neighbor->create_list(false);
    store_x0_alpha0();
}

template<>
template<class...VS>
void NewDynamicDMD<false,false,true>::update(VS*&... __vs)
{
    
    if(decide())
    {
        updt->update_wo_x(atoms->alpha,__vs...);
        return;
    }
    atoms->update_max_alpha();
    
    atoms->x2s_lcl();
    atoms->max_cut=ff->max_cut+atoms->comm.skin+alpha_scale*sqrt_2*atoms->max_alpha;
    xchng->full_xchng_static();
    updt->reset();
    updt->list();
    ff->neighbor->create_list(true);
    store_x0_alpha0();
}

template<>
template<class...VS>
void NewDynamicDMD<false,false,false>::update(VS*&... __vs)
{
    updt->update_wo_x(__vs...);
}
#endif




#endif
