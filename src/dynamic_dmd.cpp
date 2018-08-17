#include "comm.h"
#include "dynamic_dmd.h"
#include "atoms_dmd.h"
#include "xmath.h"
#include "timer.h"
#include "MAPP.h"
#include "ff_styles.h"
#include "neighbor_dmd.h"
using namespace MAPP_NS;
/*--------------------------------------------
 
 --------------------------------------------*/
DynamicDMD::DynamicDMD(AtomsDMD* __atoms,ForceFieldDMD* __ff,bool __chng_box,
std::initializer_list<vec*> __updt_vecs,
std::initializer_list<vec*> __xchng_comp_vecs,
std::initializer_list<vec*> __arch_vecs):
Dynamic(__atoms,__ff,__chng_box,
{__atoms->x,__atoms->alpha,__atoms->c,__atoms->elem},__updt_vecs,
{__atoms->id},__xchng_comp_vecs,
{},__arch_vecs),
atoms(__atoms),
ff(__ff),
c_dim(__atoms->c->dim),
alpha_scale(__atoms->xi[__atoms->N-1]),
x0(NULL),
alpha0(NULL)
{
}
/*--------------------------------------------
 
 --------------------------------------------*/
DynamicDMD::~DynamicDMD()
{
}
/*--------------------------------------------
 init a simulation
 --------------------------------------------*/
void DynamicDMD::init()
{
    
    store_arch_vecs();
    create_dynamic_vecs();
    x0=new Vec<type0>(atoms,__dim__,"x0");
    alpha0=new Vec<type0>(atoms,c_dim,"alpha0");

    ff->dynamic=this;
    ff->init();
    ff->neighbor->init();
    atoms->max_cut=ff->max_cut+atoms->comm.skin+alpha_scale*sqrt_2*atoms->max_alpha;
    
    xchng=new Exchange(atoms,nxchng_vecs_full);
#ifdef NEW_UPDATE
    updt=new Update(atoms,nupdt_vecs_full,nxchng_vecs_full);
#else
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
void DynamicDMD::fin()
{
    ff->neighbor->fin();
    ff->fin();

    delete updt;
    delete xchng;
    delete alpha0;
    delete x0;
    
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
/*--------------------------------------------
 
 --------------------------------------------*/
void DynamicDMD::store_x0()
{
    int last_atm=atoms->natms_lcl;
    if(chng_box) last_atm+=atoms->natms_ph;
    memcpy(x0->begin(),atoms->x->begin(),last_atm*__dim__*sizeof(type0));
    memcpy(alpha0->begin(),atoms->alpha->begin(),last_atm*c_dim*sizeof(type0));
}
/*--------------------------------------------
 
 --------------------------------------------*/
#ifdef NEW_UPDATE
#else
inline
#endif
bool DynamicDMD::decide()
{
    int succ,succ_lcl=1;
    type0* x_vec=atoms->x->begin();
    type0* x0_vec=x0->begin();
    type0* alpha_vec=atoms->alpha->begin();
    type0* alpha0_vec=alpha0->begin();
    int last_atm=atoms->natms_lcl;
    if(chng_box) last_atm+=atoms->natms_ph;
    
    for(int iatm=0;succ_lcl && iatm<last_atm;iatm++,x0_vec+=__dim__,x_vec+=__dim__,alpha_vec+=c_dim,alpha0_vec+=c_dim)
    {
        type0 dr=sqrt(Algebra::RSQ<__dim__>(x0_vec,x_vec));
        type0 dalpha=alpha_vec[0]-alpha0_vec[0];
        for(int i=0;i<c_dim;i++)
            dalpha=MAX(dalpha,alpha_vec[i]-alpha0_vec[i]);
        
        if(dr+dalpha*alpha_scale>0.5*skin) succ_lcl=0;
    }

    MPI_Allreduce(&succ_lcl,&succ,1,MPI_INT,MPI_MIN,world);
    if(succ) return true;
    return false;
}
#ifdef NEW_UPDATE
#else
/*--------------------------------------------
 update one vectors
 --------------------------------------------*/
void DynamicDMD::update(vec* updt_vec)
{
    update(&updt_vec,1);
}
/*--------------------------------------------
 update one vectors
 --------------------------------------------*/
void DynamicDMD::update(Vec<type0>* dx,type0 (*dH)[__dim__])
{
    updt->update(dx,dH);
}
/*--------------------------------------------
 update a number of vectors
 --------------------------------------------*/
void DynamicDMD::update(vec** updt_vecs,int nupdt_vecs)
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
            if(nupdt_vecs==1)
                updt->update(atoms->x,true);
            else
                updt->update(updt_vecs,nupdt_vecs,true);
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
#endif
#ifdef NEW_UPDATE
/*--------------------------------------------
 
 --------------------------------------------*/
void DynamicDMD::neighboring()
{
    ff->neighbor->create_list(true);
}
#endif
