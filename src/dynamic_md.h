#ifndef __MAPP__dynamic_md__
#define __MAPP__dynamic_md__

#include "dynamic.h"
#include "atoms_md.h"
#include "ff_styles.h"
#include "neighbor_md.h"

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
        void neighboring();
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
    };

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
    
}




namespace MAPP_NS
{
    class vec;
    template<typename> class Vec;
    
    template<bool BC,bool X>
    class NewDynamicMD: public NewDynamic
    {
    private:
        class ForceFieldMD* ff;
        class AtomsMD* atoms;


        bool decide();
        void store_x0();
        void alloc_x0();
    protected:
        Vec<type0>* x0;
    
    public:
  
        NewDynamicMD(AtomsMD* __atoms,ForceFieldMD* __ff,
        std::initializer_list<vec*> __updt_vecs,
        std::initializer_list<vec*> __xchng_comp_vecs,
        std::initializer_list<vec*> __arch_vecs):
        NewDynamic(__atoms,__ff,
        {__atoms->x,__atoms->elem},__updt_vecs,
        {__atoms->id},__xchng_comp_vecs,
        {},__arch_vecs),
        ff(__ff),
        atoms(__atoms),
        x0(NULL)
        {
        }
        ~NewDynamicMD()
        {
        }

        void init()
        {
            
            store_arch_vecs();
            create_dynamic_vecs();
            alloc_x0();
            ff->init();
            ff->neighbor->init();
            atoms->max_cut=ff->max_cut+atoms->comm.skin;
            
            xchng=new Exchange(atoms,nxchng_vecs_full);
            updt=new Update(atoms,nupdt_vecs_full,nxchng_vecs_full);
            
            ff->updt=updt;
            atoms->x2s_lcl();
            xchng->full_xchng();
            updt->reset();
            updt->list();
            ff->neighbor->create_list(true);
            store_x0();
        }
        void fin()
        {
            ff->neighbor->fin();
            ff->fin();
            
            delete updt;
            updt=NULL;
            delete xchng;
            xchng=NULL;
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
        void reset(){};
        
    };

}
using namespace MAPP_NS;
/*--------------------------------------------
 
 --------------------------------------------*/
template<>
template<class...VS>
void NewDynamicMD<true,true>::update(VS*&... __vs)
{
    atoms->update_H();
    updt->update_w_x(__vs...);
    if(decide()) return;
    reset();
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<>
template<class...VS>
void NewDynamicMD<false,true>::update(VS*&... __vs)
{
    if(decide())
    {
        updt->update_w_x(__vs...);
        return;
    }
    reset();
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<>
template<class...VS>
void NewDynamicMD<false,false>::update(VS*&... __vs)
{
    updt->update_wo_x(__vs...);
}
#endif
