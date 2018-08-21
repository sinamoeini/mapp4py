#include "comm.h"
#include "dynamic_md.h"
#include "atoms_md.h"
#include "xmath.h"
#include "timer.h"
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
atoms(__atoms),
ff(__ff),
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
#ifdef NEW_UPDATE
#else
    ff->dynamic=this;
#endif
    ff->init();
    ff->neighbor->init();
    atoms->max_cut=ff->max_cut+atoms->comm.skin;
    
    
#ifdef NEW_UPDATE
    xchng=new Exchange(atoms,nxchng_vecs_full);
    updt=new Update(atoms,nupdt_vecs_full,nxchng_vecs_full);
    ff->updt=updt;
#else
    xchng=new OldExchange(atoms,nxchng_vecs_full);
    updt=new OldUpdate(atoms,nupdt_vecs_full,nxchng_vecs_full);
#endif
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
#ifdef NEW_UPDATE
#else
inline
#endif
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
#ifdef NEW_UPDATE
#else
/*--------------------------------------------
 update one vectors
 --------------------------------------------*/
void DynamicMD::update(vec* updt_vec)
{
    update(&updt_vec,1);
}
/*--------------------------------------------
 update a number of vectors
 --------------------------------------------*/
void DynamicMD::update(vec** updt_vecs,int nupdt_vecs)
{
    bool x_xst=false;
    for(int ivec=0;x_xst==false && ivec<nupdt_vecs;ivec++)
        if(updt_vecs[ivec]==atoms->x)
            x_xst=true;
    if(x_xst==false)
    {
        if(nupdt_vecs==1)
            updt->update(updt_vecs[0],false);
        else
            updt->update(updt_vecs,nupdt_vecs,false);
        return;
    }
    
    
    if(chng_box)
    {
        if(nupdt_vecs==1)
            updt->update(atoms->x,true);
        else
            updt->update(updt_vecs,nupdt_vecs,true);

        if(decide())
            return;
        
        atoms->x2s_lcl();
        xchng->full_xchng();
        
        updt->reset();
        updt->list();
        ff->neighbor->create_list(chng_box);
        store_x0();
    }
    else
    {
        if(decide())
        {
            if(nupdt_vecs==1)
                updt->update(atoms->x,true);
            else
                updt->update(updt_vecs,nupdt_vecs,true);
            return;
        }

        atoms->x2s_lcl();
        xchng->full_xchng();
        
        updt->list();
        ff->neighbor->create_list(chng_box);
        
        store_x0();
    }
}
#endif
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
#ifdef NEW_UPDATE
/*--------------------------------------------
 
 --------------------------------------------*/
void DynamicMD::neighboring()
{
    ff->neighbor->create_list(chng_box);
}
#endif
