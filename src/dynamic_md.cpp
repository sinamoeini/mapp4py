#include "comm.h"
#include "dynamic_md.h"
#include "atoms_md.h"
#include "xmath.h"
#include "MAPP.h"
#include "ff_styles.h"
#include "neighbor_md.h"
using namespace MAPP_NS;
/*--------------------------------------------
 
 --------------------------------------------*/
DynamicMD::DynamicMD(AtomsMD* __atoms,ForceFieldMD* __ff,bool __chng_box,
std::initializer_list<vec*> __updt_vecs,
std::initializer_list<vec*> __xchng_comp_vecs,
std::initializer_list<vec*> __arch_vecs):
Dynamic(__atoms,__ff,__chng_box,
{__atoms->x,__atoms->elem},__updt_vecs,
{__atoms->id},__xchng_comp_vecs,
{},__arch_vecs),
ff(__ff),
atoms(__atoms),
x0(NULL)
{
}
/*--------------------------------------------
 
 --------------------------------------------*/
DynamicMD::~DynamicMD()
{
}
/*--------------------------------------------
 
 --------------------------------------------*/
void DynamicMD::init()
{
    store_arch_vecs();
    create_dynamic_vecs();
    x0=new Vec<type0>(atoms,__dim__,"x0");
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
/*--------------------------------------------
 update one vectors
 --------------------------------------------*/
void DynamicMD::fin()
{
    ff->neighbor->fin();
    ff->fin();

    delete updt;
    delete xchng;
    delete x0;
    
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
/*--------------------------------------------
 
 --------------------------------------------*/
void DynamicMD::store_x0()
{
    int last_atm=atoms->natms_lcl;
    if(chng_box) last_atm+=atoms->natms_ph;
    memcpy(x0->begin(),atoms->x->begin(),last_atm*__dim__*sizeof(type0));
}

/*--------------------------------------------
 
 --------------------------------------------*/
bool DynamicMD::decide()
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
/*--------------------------------------------
 
 --------------------------------------------*/
void DynamicMD::init_xchng()
{
    atoms->x2s_lcl();
    xchng->full_xchng();
    updt->list();
}
/*--------------------------------------------
 
 --------------------------------------------*/
void DynamicMD::fin_xchng()
{
    updt->list();
    ff->neighbor->create_list(chng_box);
    store_x0();
}
/*--------------------------------------------
 
 --------------------------------------------*/
void DynamicMD::neighboring()
{
    ff->neighbor->create_list(chng_box);
}
/*------------------------------------------------------------------------------------------------------------------------------------
 _____   _____    _   _   _____        _____   _____    _   _   _____        _____       ___   _       _____   _____
|_   _| |  _  \  | | | | | ____|      |_   _| |  _  \  | | | | | ____|      |  ___|     /   | | |     /  ___/ | ____|
  | |   | |_| |  | | | | | |__          | |   | |_| |  | | | | | |__        | |__      / /| | | |     | |___  | |__
  | |   |  _  /  | | | | |  __|         | |   |  _  /  | | | | |  __|       |  __|    / / | | | |     \___  \ |  __|
  | |   | | \ \  | |_| | | |___         | |   | | \ \  | |_| | | |___       | |      / /  | | | |___   ___| | | |___
  |_|   |_|  \_\ \_____/ |_____|        |_|   |_|  \_\ \_____/ |_____|      |_|     /_/   |_| |_____| /_____/ |_____|
 ------------------------------------------------------------------------------------------------------------------------------------*/
template<>
void NewDynamicMD<true,true>::alloc_x0()
{
    x0=new Vec<type0>(atoms,__dim__,"x0");
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<>
void NewDynamicMD<true,true>::store_x0()
{
    memcpy(x0->begin(),atoms->x->begin(),(atoms->natms_lcl+atoms->natms_ph)*__dim__*sizeof(type0));
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<>
bool NewDynamicMD<true,true>::decide()
{
    type0 skin_sq=0.25*Algebra::pow<2>(dynamic_sub.skin);
    int succ,succ_lcl=1;
    type0* x_vec=atoms->x->begin();
    type0* x0_vec=x0->begin();
    int last_atm=atoms->natms_lcl+atoms->natms_ph;
    for(int iatm=0;succ_lcl && iatm<last_atm;iatm++,x0_vec+=__dim__,x_vec+=__dim__)
        if(Algebra::RSQ<__dim__>(x0_vec,x_vec)>skin_sq)
            succ_lcl=0;
    
    MPI_Allreduce(&succ_lcl,&succ,1,MPI_INT,MPI_MIN,dynamic_sub.world);
    if(succ) return true;
    return false;
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<>
void NewDynamicMD<true,true>::reset()
{
    atoms->x2s_lcl();
    dynamic_sub.xchng->full_xchng();
    dynamic_sub.updt->reset();
    dynamic_sub.updt->list();
    ff->neighbor->create_list(true);
    store_x0();
}
/*------------------------------------------------------------------------------------------------------------------------------------
 _____       ___   _       _____   _____        _____   _____    _   _   _____        _____       ___   _       _____   _____
|  ___|     /   | | |     /  ___/ | ____|      |_   _| |  _  \  | | | | | ____|      |  ___|     /   | | |     /  ___/ | ____|
| |__      / /| | | |     | |___  | |__          | |   | |_| |  | | | | | |__        | |__      / /| | | |     | |___  | |__
|  __|    / / | | | |     \___  \ |  __|         | |   |  _  /  | | | | |  __|       |  __|    / / | | | |     \___  \ |  __|
| |      / /  | | | |___   ___| | | |___         | |   | | \ \  | |_| | | |___       | |      / /  | | | |___   ___| | | |___
|_|     /_/   |_| |_____| /_____/ |_____|        |_|   |_|  \_\ \_____/ |_____|      |_|     /_/   |_| |_____| /_____/ |_____|
 ------------------------------------------------------------------------------------------------------------------------------------*/
template<>
void NewDynamicMD<false,true>::alloc_x0()
{
    x0=new Vec<type0>(atoms,__dim__,"x0");
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<>
void NewDynamicMD<false,true>::store_x0()
{
    memcpy(x0->begin(),atoms->x->begin(),atoms->natms_lcl*__dim__*sizeof(type0));
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<>
bool NewDynamicMD<false,true>::decide()
{
    type0 skin_sq=0.25*Algebra::pow<2>(dynamic_sub.skin);
    int succ,succ_lcl=1;
    type0* x_vec=atoms->x->begin();
    type0* x0_vec=x0->begin();
    int last_atm=atoms->natms_lcl;
    for(int iatm=0;succ_lcl && iatm<last_atm;iatm++,x0_vec+=__dim__,x_vec+=__dim__)
        if(Algebra::RSQ<__dim__>(x0_vec,x_vec)>skin_sq)
            succ_lcl=0;
    
    MPI_Allreduce(&succ_lcl,&succ,1,MPI_INT,MPI_MIN,dynamic_sub.world);
    if(succ) return true;
    return false;
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<>
void NewDynamicMD<false,true>::reset()
{
    atoms->x2s_lcl();
    dynamic_sub.xchng->full_xchng();
    dynamic_sub.updt->list();
    ff->neighbor->create_list(false);
    store_x0();
}
/*------------------------------------------------------------------------------------------------------------------------------------
 _____       ___   _       _____   _____        _____       ___   _       _____   _____        _____       ___   _       _____   _____
|  ___|     /   | | |     /  ___/ | ____|      |  ___|     /   | | |     /  ___/ | ____|      |  ___|     /   | | |     /  ___/ | ____|
| |__      / /| | | |     | |___  | |__        | |__      / /| | | |     | |___  | |__        | |__      / /| | | |     | |___  | |__
|  __|    / / | | | |     \___  \ |  __|       |  __|    / / | | | |     \___  \ |  __|       |  __|    / / | | | |     \___  \ |  __|
| |      / /  | | | |___   ___| | | |___       | |      / /  | | | |___   ___| | | |___       | |      / /  | | | |___   ___| | | |___
|_|     /_/   |_| |_____| /_____/ |_____|      |_|     /_/   |_| |_____| /_____/ |_____|      |_|     /_/   |_| |_____| /_____/ |_____|
 ------------------------------------------------------------------------------------------------------------------------------------*/
template<>
void NewDynamicMD<false,false>::alloc_x0()
{
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<>
void NewDynamicMD<false,false>::store_x0()
{
}
