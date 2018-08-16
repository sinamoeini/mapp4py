#include "atoms.h"
#include "exchange.h"
#include "comm.h"
#include "xmath.h"
/*-------temp_remove-------
#include "ff.h"
#include "neighbor.h"
 */
using namespace MAPP_NS;
/*------------------------------------------------------------------
 _____  __    __  _____   _   _       ___   __   _   _____   _____  
| ____| \ \  / / /  ___| | | | |     /   | |  \ | | /  ___| | ____| 
| |__    \ \/ /  | |     | |_| |    / /| | |   \| | | |     | |__   
|  __|    }  {   | |     |  _  |   / / | | | |\   | | |  _  |  __|  
| |___   / /\ \  | |___  | | | |  / /  | | | | \  | | |_| | | |___  
|_____| /_/  \_\ \_____| |_| |_| /_/   |_| |_|  \_| \_____/ |_____|
 ------------------------------------------------------------------*/
/*--------------------------------------------
 
 --------------------------------------------*/
Exchange::Exchange(Atoms* atoms,int& nxchng_vecs_):
natms_lcl(atoms->natms_lcl),
x(atoms->x),
world(atoms->comm.world),
rank(atoms->comm.rank),
s_lo(atoms->comm.s_lo),
s_hi(atoms->comm.s_hi),
xchng_id(atoms->comm.xchng_id),

vecs(atoms->dynamic_vecs),
nvecs(atoms->ndynamic_vecs),
nxchng_vecs(nxchng_vecs_)
{
    Algebra::V_eq<__dim__*2>(&(atoms->comm.neigh[0][0]), &(neigh[0][0]));
    snd_buff[0]=snd_buff[1]=NULL;
    snd_buff_sz[0]=snd_buff_sz[1]=0;
    snd_buff_cpcty[0]=snd_buff_cpcty[1]=0;
    
    rcv_buff=NULL;
    rcv_buff_sz=0;
    rcv_buff_cpcty=0;
    
    tot_xchng_sz=0;
    for(int ivec=0;ivec<nxchng_vecs;ivec++)
        tot_xchng_sz+=vecs[ivec]->byte_sz;
    
}
/*--------------------------------------------
 destructor
 --------------------------------------------*/
Exchange::~Exchange()
{
    delete [] snd_buff[0];
    delete [] snd_buff[1];
    delete [] rcv_buff;
}
/*--------------------------------------------
 
 --------------------------------------------*/
inline void Exchange::load(int& iatm,int idir)
{
    if(snd_buff_cpcty[idir]<snd_buff_sz[idir]+tot_xchng_sz)
    {
        byte* tmp_buff=new byte[snd_buff_sz[idir]+tot_xchng_sz+buff_grw];
        memcpy(tmp_buff,snd_buff[idir],snd_buff_sz[idir]);
        delete [] snd_buff[idir];
        snd_buff[idir]=tmp_buff;
        snd_buff_cpcty[idir]=snd_buff_sz[idir]+tot_xchng_sz+buff_grw;
    }
    byte* tmp_buff=&snd_buff[idir][snd_buff_sz[idir]];
    
    for(int ivec=0;ivec<nxchng_vecs;ivec++)
        vecs[ivec]->pop_out(tmp_buff,iatm);
    
    snd_buff_sz[idir]+=tot_xchng_sz;
}
/*--------------------------------------------
 
 --------------------------------------------*/
inline void Exchange::load(byte*& buff,int& idir)
{
    if(snd_buff_cpcty[idir]<snd_buff_sz[idir]+tot_xchng_sz)
    {
        byte* tmp_buff=new byte[snd_buff_sz[idir]+tot_xchng_sz+buff_grw];
        memcpy(tmp_buff,snd_buff[idir],snd_buff_sz[idir]);
        delete [] snd_buff[idir];
        snd_buff[idir]=tmp_buff;
        snd_buff_cpcty[idir]=snd_buff_sz[idir]+tot_xchng_sz+buff_grw;
    }
    memcpy(&snd_buff[idir][snd_buff_sz[idir]],buff,tot_xchng_sz);
    snd_buff_sz[idir]+=tot_xchng_sz;
    buff+=tot_xchng_sz;
}
/*--------------------------------------------
 
 --------------------------------------------*/
inline int Exchange::xchng_buff(int idim,int idir)
{
    rcv_buff_sz=0;
    int max_snd_sz;
    MPI_Allreduce(&snd_buff_sz[idir],&max_snd_sz,1,MPI_INT,MPI_MAX,world);
    if(max_snd_sz==0)
        return 0;
    
    MPI_Sendrecv(&snd_buff_sz[idir],1,MPI_INT,neigh[idim][idir],0,
                 &rcv_buff_sz,1,MPI_INT,neigh[idim][1-idir],0,
                 world,MPI_STATUS_IGNORE);
    
    if(rcv_buff_cpcty<rcv_buff_sz)
    {
        delete [] rcv_buff;
        rcv_buff=new byte[rcv_buff_sz+buff_grw];
        rcv_buff_cpcty=rcv_buff_sz+buff_grw;
    }
    
    MPI_Sendrecv(snd_buff[idir],snd_buff_sz[idir],MPI_BYTE,neigh[idim][idir],0,
                 rcv_buff,rcv_buff_sz,MPI_BYTE,neigh[idim][1-idir],0,
                 world,MPI_STATUS_IGNORE);

    
    snd_buff_sz[idir]=0;
    
    /* 
     * here make sure, that we have enough space
     * for the potential atoms to be inserted in
     * xchng vectors
     */
    for(int ivec=0;ivec<nxchng_vecs;ivec++)
        vecs[ivec]->reserve(rcv_buff_sz/tot_xchng_sz);
    
    return 1;
}
/*--------------------------------------------
 
 --------------------------------------------*/
void Exchange::full_xchng()
{
    for(int ivec=0;ivec<nxchng_vecs;ivec++)
        vecs[ivec]->resize(natms_lcl);
    int disp;
    type0 s,ds_lo,ds_hi;
    int iatm;
    int xchng_lcl=0;
    
    for(int idim=0;idim<__dim__;idim++)
    {
        if(rank==neigh[idim][0] && rank==neigh[idim][1])
            continue;
        
        disp=idim*sizeof(type0);
        
        snd_buff_sz[0]=snd_buff_sz[1]=0;
        iatm=0;
        type0* s_ptr=x->begin()+idim;
        unsigned int n=x->vec_sz;
        while(iatm<n)
        {
            //s=(*x)(iatm,idim);
            //s=x->begin()[iatm*__dim__+idim];
            s=*s_ptr;
            if(s_lo[idim]<=s && s<s_hi[idim])
            {
                iatm++;
                s_ptr+=__dim__;
                continue;
            }
            
            if(s>=s_hi[idim])
            {
                ds_hi=s-s_hi[idim];
                ds_lo=1.0+s_lo[idim]-s;
            }
            else
            {
                ds_hi=1.0+s-s_hi[idim];
                ds_lo=s_lo[idim]-s;
            }
            
            if(ds_hi<ds_lo)
                load(iatm,snd_to_frnt);
            else
                load(iatm,snd_to_bhnd);
            n--;
        }
        
        if(snd_buff_sz[0]+snd_buff_sz[1])
            xchng_lcl=1;
        
        byte* buff_tmp;
        for(int idir=0;idir<2;idir++)
        {
            while(xchng_buff(idim,idir))
            {
                buff_tmp=rcv_buff;
                while(rcv_buff_sz)
                {
                    s=*(type0*)(buff_tmp+disp);
                    if(s_lo[idim]<=s && s<s_hi[idim])
                        for(int ivec=0;ivec<nxchng_vecs;ivec++)
                            vecs[ivec]->pop_in(buff_tmp);
                    else
                        load(buff_tmp,idir);
                    
                    rcv_buff_sz-=tot_xchng_sz;
                }
            }
            
        }
    }
    int xchng;
    MPI_Allreduce(&xchng_lcl,&xchng,1,MPI_INT,MPI_MAX,world);
    if(xchng) xchng_id++;
    natms_lcl=x->vec_sz;
}
/*--------------------------------------------
 
 --------------------------------------------*/
void Exchange::full_xchng_all()
{
    if(nvecs==nxchng_vecs) return full_xchng();
    int __nxchng_vecs=nxchng_vecs;
    int __tot_xchng_sz=tot_xchng_sz;
    
    nxchng_vecs=nvecs;
    tot_xchng_sz=0;
    for(int ivec=0;ivec<nxchng_vecs;ivec++)
        tot_xchng_sz+=vecs[ivec]->byte_sz;
    full_xchng();
    nxchng_vecs=__nxchng_vecs;
    tot_xchng_sz=__tot_xchng_sz;
}
/*------------------------------------------------
 _   _   _____   _____       ___   _____   _____  
| | | | |  _  \ |  _  \     /   | |_   _| | ____| 
| | | | | |_| | | | | |    / /| |   | |   | |__   
| | | | |  ___/ | | | |   / / | |   | |   |  __|  
| |_| | | |     | |_| |  / /  | |   | |   | |___  
\_____/ |_|     |_____/ /_/   |_|   |_|   |_____|
 ------------------------------------------------*/
/*--------------------------------------------
 
 --------------------------------------------*/
Update::Update(Atoms* atoms,
int& nupdt_vecs_,int& nxchng_vecs_):
natms_lcl(atoms->natms_lcl),
natms_ph(atoms->natms_ph),
H(atoms->H),
depth_inv(atoms->depth_inv),

rank(atoms->comm.rank),
s_lo(atoms->comm.s_lo),
s_hi(atoms->comm.s_hi),

max_cut(atoms->max_cut),
x(atoms->x),

vecs(atoms->dynamic_vecs),
nvecs(atoms->ndynamic_vecs),
nupdt_vecs(nupdt_vecs_),
nxchng_vecs(nxchng_vecs_)
{

    Algebra::V_eq<__dim__*2>(&(atoms->comm.neigh[0][0]),&(neigh[0][0]));
    snd_buff=NULL;
    snd_buff_cpcty=0;
    
    rcv_buff=NULL;
    rcv_buff_cpcty=0;

    
    tot_updt_vecs_sz=0;
    for(int ivec=0;ivec<nupdt_vecs;ivec++)
        tot_updt_vecs_sz+=vecs[ivec]->byte_sz;
    
    int icurs=0;
    for(int idim=0;idim<__dim__;idim++)
    {
        if(rank==neigh[idim][0] && rank==neigh[idim][1])
            comm_manager[idim]=new LoadUnLoadUpdateSelfComm(this);
        else
            comm_manager[idim]=new LoadUnLoadUpdateComm(this,atoms->comm.world);
        
        // snd_to_bhnd && rcv_fm_frnt
        pbc_correction[idim][0]=(atoms->comm.coords[idim]==atoms->comm.dims[idim]-1);
        icurs++;
        
        // snd_to_frnt && rcv_fm_bhnd
        pbc_correction[idim][1]=(atoms->comm.coords[idim]==0);
        icurs++;
    }

    tot_ncomms=0;
    snd_atms_lst=NULL;
    snd_atms_lst_cpcty=snd_atms_lst_sz=rcv_atms_lst_sz=NULL;
}
/*--------------------------------------------
 
 --------------------------------------------*/
Update::~Update()
{
    delete [] rcv_buff;
    delete [] snd_buff;
    
    for(int i=0;i<tot_ncomms;i++)
        delete [] snd_atms_lst[i];
    
    delete [] snd_atms_lst;
    delete [] snd_atms_lst_cpcty;
    delete [] snd_atms_lst_sz;
    delete [] rcv_atms_lst_sz;
    
    for(int idim=0;idim<__dim__;idim++)
        delete comm_manager[idim];
}
/*--------------------------------------------
 
 --------------------------------------------*/
void Update::reset()
{
    for(int idim=0;idim<__dim__;idim++)
        max_cut_s[idim]=max_cut*depth_inv[idim];
    
    int icurs=0;
    int tot_ncomms_=0;
    for(int idim=0;idim<__dim__;idim++)
    {
        // snd_to_bhnd && rcv_fm_frnt
        s_bnd[idim][0]=s_lo[idim]+max_cut_s[idim];
        tot_ncomms_+=static_cast<int>(max_cut_s[idim]/(s_hi[idim]-s_lo[idim]))+1;
        ncomms[idim][0]=tot_ncomms_;
        icurs++;
        // snd_to_frnt && rcv_fm_bhnd
        s_bnd[idim][1]=s_hi[idim]-max_cut_s[idim];
        tot_ncomms_+=static_cast<int>(max_cut_s[idim]/(s_hi[idim]-s_lo[idim]))+1;
        ncomms[idim][1]=tot_ncomms_;
        icurs++;
    }
    
    if(tot_ncomms_==tot_ncomms)
        return;
    
    for(int i=0;i<tot_ncomms;i++)
        delete [] snd_atms_lst[i];
    delete [] snd_atms_lst;
    delete [] snd_atms_lst_cpcty;
    delete [] snd_atms_lst_sz;
    delete [] rcv_atms_lst_sz;
    
    tot_ncomms=tot_ncomms_;
    snd_atms_lst=new int*[tot_ncomms];
    snd_atms_lst_cpcty=new int[tot_ncomms];
    snd_atms_lst_sz=new int[tot_ncomms];
    rcv_atms_lst_sz=new int[tot_ncomms];
    for(int i=0;i<tot_ncomms;i++)
    {
        snd_atms_lst_cpcty[i]=0;
        snd_atms_lst[i]=NULL;
    }
}
/*--------------------------------------------
 
 --------------------------------------------*/
inline void Update::add_to_snd_lst(int& icomm,int& iatm)
{
    if(snd_atms_lst_sz[icomm]+1>snd_atms_lst_cpcty[icomm])
    {
        int* tmp_lst=new int[snd_atms_lst_sz[icomm]+1+snd_atms_lst_grw];
        memcpy(tmp_lst,snd_atms_lst[icomm],snd_atms_lst_sz[icomm]*sizeof(int));
        delete [] snd_atms_lst[icomm];
        snd_atms_lst[icomm]=tmp_lst;
        snd_atms_lst_cpcty[icomm]=snd_atms_lst_sz[icomm]+1+snd_atms_lst_grw;
    }
    snd_atms_lst[icomm][snd_atms_lst_sz[icomm]]=iatm;
    snd_atms_lst_sz[icomm]++;
}
/*--------------------------------------------
 
 --------------------------------------------*/
inline void Update::reserve_rcv_buff(int xtra)
{
    if(rcv_buff_cpcty<xtra+rcv_buff_sz)
    {
        delete [] rcv_buff;
        rcv_buff=new byte[xtra+rcv_buff_sz+rcv_buff_grw];
        rcv_buff_cpcty=xtra+rcv_buff_sz+rcv_buff_grw;
    }
}
/*--------------------------------------------
 
 --------------------------------------------*/
inline void Update::reserve_snd_buff(int xtra)
{
    if(snd_buff_cpcty<xtra+snd_buff_sz)
    {
        delete [] snd_buff;
        snd_buff=new byte[xtra+snd_buff_sz+snd_buff_grw];
        snd_buff_cpcty=xtra+snd_buff_sz+snd_buff_grw;
    }
}
/*--------------------------------------------
 
 --------------------------------------------*/
void Update::update(vec** updt_vecs,int nupdt_vecs,bool x_xst)
{
    int tot_byte_sz=0;
    for(int ivec=0;ivec<nupdt_vecs;ivec++)
    {
        updt_vecs[ivec]->vec_sz=natms_lcl;
        tot_byte_sz+=updt_vecs[ivec]->byte_sz;
    }
    snd_buff_sz=rcv_buff_sz=0;
    reserve_snd_buff(tot_byte_sz*max_snd_atms_lst_sz);
    reserve_rcv_buff(tot_byte_sz*max_rcv_atms_lst_sz);
    

    int icurs=0;
    int icomm=0;
    for(int idim=0;idim<__dim__;idim++)
    {
        for(int idir=0;idir<2;idir++)
        {
            while(icomm<ncomms[idim][idir])
            {
                comm_manager[idim]->update_mult(icomm
                ,neigh[idim][idir]
                ,neigh[idim][1-idir]
                ,updt_vecs,nupdt_vecs,tot_byte_sz);
                
                if(x_xst && pbc_correction[idim][idir])
                {
                    type0* __x_vec=x->end()-__dim__;
                    if(idir)
                    {
                        for(int iatm=0;iatm<rcv_atms_lst_sz[icomm];iatm++,__x_vec-=__dim__)
                            for(int jdim=0;jdim<idim+1;jdim++)
                                __x_vec[jdim]-=H[idim][jdim];
                        
                    }
                    else
                    {
                        for(int iatm=0;iatm<rcv_atms_lst_sz[icomm];iatm++,__x_vec-=__dim__)
                            for(int jdim=0;jdim<idim+1;jdim++)
                                __x_vec[jdim]+=H[idim][jdim];
                    }
                }
                icomm++;
            }
            icurs++;
        }
    }
}
/*--------------------------------------------
 
 --------------------------------------------*/
void Update::update(vec* updt_vec,bool x_xst)
{
    snd_buff_sz=0;
    reserve_snd_buff(updt_vec->byte_sz*max_snd_atms_lst_sz);
    updt_vec->vec_sz=natms_lcl;
    
    int icurs=0;
    int icomm=0;
    for(int idim=0;idim<__dim__;idim++)
    {
        for(int idir=0;idir<2;idir++)
        {
            while(icomm<ncomms[idim][idir])
            {
                comm_manager[idim]->update_sing(icomm
                ,neigh[idim][idir]
                ,neigh[idim][1-idir]
                ,updt_vec);
                
                if(x_xst && pbc_correction[idim][idir])
                {
                    type0* __x_vec=x->end()-__dim__;
                    if(idir)
                    {
                        for(int iatm=0;iatm<rcv_atms_lst_sz[icomm];iatm++,__x_vec-=__dim__)
                            for(int jdim=0;jdim<idim+1;jdim++)
                                __x_vec[jdim]-=H[idim][jdim];
                        
                    }
                    else
                    {
                        for(int iatm=0;iatm<rcv_atms_lst_sz[icomm];iatm++,__x_vec-=__dim__)
                            for(int jdim=0;jdim<idim+1;jdim++)
                                __x_vec[jdim]+=H[idim][jdim];
                    }
                    
                }
                icomm++;
            }
            icurs++;
        }
    }
}
/*--------------------------------------------
 
 --------------------------------------------*/
void Update::update(vec* updt_vec,type0 (*dH)[__dim__])
{
    
    snd_buff_sz=0;
    reserve_snd_buff(updt_vec->byte_sz*max_snd_atms_lst_sz);
    updt_vec->vec_sz=natms_lcl;
    
    int icurs=0;
    int icomm=0;
    for(int idim=0;idim<__dim__;idim++)
    {
        for(int idir=0;idir<2;idir++)
        {
            while(icomm<ncomms[idim][idir])
            {
                comm_manager[idim]->update_sing(icomm
                ,neigh[idim][idir]
                ,neigh[idim][1-idir]
                ,updt_vec);
                
                if(pbc_correction[idim][idir])
                {
                    
                    type0* __vec=static_cast<type0*>(updt_vec->end())-__dim__;

                    if(idir)
                    {
                        for(int iatm=0;iatm<rcv_atms_lst_sz[icomm];iatm++,__vec-=__dim__)
                            for(int jdim=0;jdim<idim+1;jdim++)
                                __vec[jdim]-=H[idim][jdim];
                    }
                    else
                    {
                        for(int iatm=0;iatm<rcv_atms_lst_sz[icomm];iatm++,__vec-=__dim__)
                            for(int jdim=0;jdim<idim+1;jdim++)
                                __vec[jdim]+=H[idim][jdim];
                    }
                }
                icomm++;
            }
            icurs++;
        }
    }
}
/*--------------------------------------------
 
 --------------------------------------------*/
void Update::list()
{
    natms_ph=0;
    type0* x_vec;
    int x_dim=x->dim;
    int icurs=0;
    int icomm=0;
    int lo_atm,hi_atm;
    int last_atm;
    bool dir;
    type0 inc;
    max_snd_atms_lst_sz=max_rcv_atms_lst_sz=0;
    
    for(int idim=0;idim<__dim__;idim++)
    {
        
        last_atm=x->vec_sz;
        
        inc=1.0;
        dir=true;
        for(int idir=0;idir<2;idir++)
        {
            lo_atm=0;
            hi_atm=last_atm;
            while(icomm<ncomms[idim][idir])
            {
                snd_buff_sz=0;
                rcv_atms_lst_sz[icomm]=snd_atms_lst_sz[icomm]=0;
                x_vec=x->begin();
                for(int iatm=lo_atm;iatm<hi_atm;iatm++)
                    if((x_vec[iatm*x_dim+idim]<s_bnd[idim][idir])==dir)
                        add_to_snd_lst(icomm,iatm);
                
                comm_manager[idim]->load_unload(icomm
                ,neigh[idim][idir]
                ,neigh[idim][1-idir]);
                
                lo_atm=x->vec_sz-rcv_atms_lst_sz[icomm];
                hi_atm=x->vec_sz;
                if(pbc_correction[idim][idir])
                {
                    x_vec=x->begin();
                    for(int iatm=lo_atm;iatm<hi_atm;iatm++)
                        x_vec[iatm*x_dim+idim]+=inc;
                }
                max_snd_atms_lst_sz=MAX(max_snd_atms_lst_sz,snd_atms_lst_sz[icomm]);
                max_rcv_atms_lst_sz=MAX(max_rcv_atms_lst_sz,rcv_atms_lst_sz[icomm]);
                icomm++;
            }
            icurs++;
            inc-=2.0;
            dir=!dir;
        }
    }
    natms_ph=x->vec_sz-natms_lcl;
    for(int ivec=nxchng_vecs;ivec<nvecs;ivec++)
    {
        vecs[ivec]->vec_sz=0;
        vecs[ivec]->resize(x->vec_sz);
    }
}
/*--------------------------------------------
 
 --------------------------------------------*/
void Update::rm_rdndncy()
{
    snd_buff_sz=rcv_buff_sz=0;
    reserve_rcv_buff(max_snd_atms_lst_sz);
    reserve_snd_buff(natms_ph);
    
    byte* mark=snd_buff;
    /*-------temp_remove-------
    forcefield->neighbor->mark_redndnt_ph(mark);
    */
    int rcv_atms_lst_sz_;
    int snd_atms_lst_sz_=0;
    int snd_atms_lst_cpcty_=max_snd_atms_lst_sz;
    int* snd_atms_lst_=NULL;
    if(snd_atms_lst_cpcty_) snd_atms_lst_=new int[snd_atms_lst_cpcty_];
    
    int nlocomm;
    byte* mark_=mark+natms_ph;
    int icurs=2*__dim__-1;
    int jcomm=tot_ncomms-1;
    max_snd_atms_lst_sz=max_rcv_atms_lst_sz=0;
    for(int idim=__dim__-1;idim>-1;idim--)
    {
        for(int idir=1;idir>-1;idir--)
        {
            if(icurs!=0)
                nlocomm=ncomms[idim][idir];
            else
                nlocomm=0;
            
            while(jcomm>nlocomm-1)
            {
                mark_-=rcv_atms_lst_sz[jcomm];
                comm_manager[idim]->xchng_buff(
                neigh[idim][1-idir],rcv_atms_lst_sz[jcomm],mark_,
                neigh[idim][idir],snd_atms_lst_sz[jcomm],rcv_buff);
                
                snd_atms_lst_sz_=0;
                for(int i=0; i<snd_atms_lst_sz[jcomm];i++)
                {
                    if(rcv_buff[i]=='1')
                    {
                        if(snd_atms_lst[jcomm][i]>=natms_lcl)
                            mark[snd_atms_lst[jcomm][i]-natms_lcl]='1';
                        snd_atms_lst_[snd_atms_lst_sz_++]=snd_atms_lst[jcomm][i];
                    }
                }
                memcpy(snd_atms_lst[jcomm],snd_atms_lst_,snd_atms_lst_sz_*sizeof(int));
                snd_atms_lst_sz[jcomm]=snd_atms_lst_sz_;
                
                rcv_atms_lst_sz_=0;
                for(int i=0; i<rcv_atms_lst_sz[jcomm];i++)
                    if(mark_[i]=='1')
                        rcv_atms_lst_sz_++;
                rcv_atms_lst_sz[jcomm]=rcv_atms_lst_sz_;
                
                max_snd_atms_lst_sz=MAX(max_snd_atms_lst_sz,snd_atms_lst_sz[jcomm]);
                max_rcv_atms_lst_sz=MAX(max_rcv_atms_lst_sz,rcv_atms_lst_sz[jcomm]);
                
                jcomm--;
            }
            icurs--;
        }
    }
    
    delete [] snd_atms_lst_;
    
    int old_2_new_cpcty=natms_lcl+natms_ph;
    int* old_2_new=NULL;
    if(old_2_new_cpcty) old_2_new=new int[old_2_new_cpcty];
    
    int list_sz=0;
    int list_cpcty=natms_ph;
    int* list=NULL;
    if(list_cpcty) list=new int[list_cpcty];
    
    for(int iatm=0;iatm<natms_lcl;iatm++)
        old_2_new[iatm]=iatm;
    
    icurs=natms_lcl;
    for(int iatm=natms_lcl;iatm<natms_lcl+natms_ph;iatm++)
        if(mark[iatm-natms_lcl]=='1')
        {
            old_2_new[iatm]=icurs++;
            list[list_sz++]=iatm;
        }
    
    int new_natms_ph=list_sz;

    for(int icomm=0;icomm<tot_ncomms;icomm++)
        for(int i=0; i<snd_atms_lst_sz[icomm];i++)
            snd_atms_lst[icomm][i]=old_2_new[snd_atms_lst[icomm][i]];
    /*-------temp_remove-------
    forcefield->neighbor->rename_atoms(old_2_new);
     */
    delete [] old_2_new;
    
    int* list_=list;
    
    int vec_sz=natms_lcl;
    while(*list_==natms_lcl+icurs)
    {
        list_++;
        vec_sz++;
        list_sz--;
    }
    
    for(int ivec=0;ivec<nupdt_vecs;ivec++)
    {
        vecs[ivec]->vec_sz=vec_sz;
        vecs[ivec]->cpy_pst(list_,list_sz);
    }
    
    delete [] list;
    
    natms_ph=new_natms_ph;
}
/*------------------------------------------------------------------------------------------------------------------------------------
 
 ------------------------------------------------------------------------------------------------------------------------------------*/
Update::LoadUnLoadUpdateComm::
LoadUnLoadUpdateComm(Update* __updt,MPI_Comm& __world):
LoadUnLoadUpdate(),
world(__world),
snd_atms_lst(__updt->snd_atms_lst),
snd_atms_lst_sz(__updt->snd_atms_lst_sz),
rcv_atms_lst_sz(__updt->rcv_atms_lst_sz),
snd_buff(__updt->snd_buff),
snd_buff_sz(__updt->snd_buff_sz),
snd_buff_cpcty(__updt->snd_buff_cpcty),
rcv_buff(__updt->rcv_buff),
rcv_buff_sz(__updt->rcv_buff_sz),
rcv_buff_cpcty(__updt->rcv_buff_cpcty),
vecs(__updt->vecs),
nupdt_vecs(__updt->nupdt_vecs),
tot_updt_vecs_sz(__updt->tot_updt_vecs_sz)
{
}
/*--------------------------------------------
 
 --------------------------------------------*/
inline void Update::LoadUnLoadUpdateComm::load_unload
(int& __icomm,int& __snd_p,int& __rcv_p)
{
    snd_buff_sz=snd_atms_lst_sz[__icomm]*tot_updt_vecs_sz;
    if(snd_buff_cpcty<snd_buff_sz)
    {
        delete [] snd_buff;
        snd_buff=new byte[snd_buff_sz+snd_buff_grw];
        snd_buff_cpcty=snd_buff_sz+snd_buff_grw;
    }
    
    byte* tmp_snd_buff=snd_buff;
    for(int ivec=0;ivec<nupdt_vecs;ivec++)
        vecs[ivec]->cpy(tmp_snd_buff,snd_atms_lst[__icomm],snd_atms_lst_sz[__icomm]);
    
    MPI_Sendrecv(&snd_atms_lst_sz[__icomm],1,MPI_INT,__snd_p,0,
                 &rcv_atms_lst_sz[__icomm],1,MPI_INT,__rcv_p,0,
                 world,MPI_STATUS_IGNORE);
    
    rcv_buff_sz=rcv_atms_lst_sz[__icomm]*tot_updt_vecs_sz;
    if(rcv_buff_cpcty<rcv_buff_sz)
    {
        delete [] rcv_buff;
        rcv_buff=new byte[rcv_buff_sz+rcv_buff_grw];
        rcv_buff_cpcty=rcv_atms_lst_sz[__icomm]*tot_updt_vecs_sz+rcv_buff_grw;
    }

    MPI_Sendrecv(snd_buff,snd_buff_sz,MPI_BYTE,__snd_p,0,
                 rcv_buff,rcv_buff_sz,MPI_BYTE,__rcv_p,0,
                 world,MPI_STATUS_IGNORE);

    
    byte* tmp_rcv_buff=rcv_buff;
    for(int ivec=0;ivec<nupdt_vecs;ivec++)
    {
        vecs[ivec]->reserve(rcv_atms_lst_sz[__icomm]);
        vecs[ivec]->pst(tmp_rcv_buff,rcv_atms_lst_sz[__icomm]);
    }
}
/*--------------------------------------------
 
 --------------------------------------------*/
inline void Update::LoadUnLoadUpdateComm::update_mult
(int& __icomm,int& __snd_p,int& __rcv_p,vec**& __vecs
,int& __nvecs,int& __vecs_byte_sz)
{
    byte* tmp_snd_buff=snd_buff;
    for(int ivec=0;ivec<__nvecs;ivec++)
        __vecs[ivec]->cpy(tmp_snd_buff,snd_atms_lst[__icomm],snd_atms_lst_sz[__icomm]);

    MPI_Sendrecv(snd_buff,snd_atms_lst_sz[__icomm]*__vecs_byte_sz,MPI_BYTE,__snd_p,0,
                 rcv_buff,rcv_atms_lst_sz[__icomm]*__vecs_byte_sz,MPI_BYTE,__rcv_p,0,
                 world,MPI_STATUS_IGNORE);

    byte* tmp_rcv_buff=rcv_buff;
    for(int ivec=0;ivec<__nvecs;ivec++)
        __vecs[ivec]->pst(tmp_rcv_buff,rcv_atms_lst_sz[__icomm]);
}
/*--------------------------------------------
 
 --------------------------------------------*/
inline void Update::LoadUnLoadUpdateComm::update_sing
(int& __icomm,int& __snd_p,int& __rcv_p,vec*& __v)
{
    byte* tmp_snd_buff=snd_buff;
    __v->cpy(tmp_snd_buff,snd_atms_lst[__icomm],snd_atms_lst_sz[__icomm]);

    MPI_Sendrecv(snd_buff,snd_atms_lst_sz[__icomm]*__v->byte_sz,MPI_BYTE,__snd_p,0,
                 __v->end(),rcv_atms_lst_sz[__icomm]*__v->byte_sz,MPI_BYTE,__rcv_p,0,
                 world,MPI_STATUS_IGNORE);
    __v->vec_sz+=rcv_atms_lst_sz[__icomm];
}
/*--------------------------------------------
 
 --------------------------------------------*/
inline void Update::LoadUnLoadUpdateComm::xchng_buff
(int& __snd_p,int& __snd_buff_sz,byte*& __snd_buff
,int& __rcv_p,int& __rcv_buff_sz,byte*& __rcv_buff)
{
    MPI_Sendrecv(__snd_buff,__snd_buff_sz,MPI_BYTE,__snd_p,0,
                 __rcv_buff,__rcv_buff_sz,MPI_BYTE,__rcv_p,0,
                 world,MPI_STATUS_IGNORE);
}
/*------------------------------------------------------------------------------------------------------------------------------------
 
 ------------------------------------------------------------------------------------------------------------------------------------*/
Update::LoadUnLoadUpdateSelfComm::
LoadUnLoadUpdateSelfComm(Update* __updt):
LoadUnLoadUpdate(),
snd_atms_lst(__updt->snd_atms_lst),
snd_atms_lst_sz(__updt->snd_atms_lst_sz),
rcv_atms_lst_sz(__updt->rcv_atms_lst_sz),
vecs(__updt->vecs),
nupdt_vecs(__updt->nupdt_vecs)
{
}
/*--------------------------------------------
 
 --------------------------------------------*/
inline void Update::LoadUnLoadUpdateSelfComm::load_unload
(int& __icomm,int&,int&)
{
    rcv_atms_lst_sz[__icomm]=snd_atms_lst_sz[__icomm];
    
    for(int ivec=0;ivec<nupdt_vecs;ivec++)
    {
        vecs[ivec]->reserve(rcv_atms_lst_sz[__icomm]);
        vecs[ivec]->cpy_pst(snd_atms_lst[__icomm],rcv_atms_lst_sz[__icomm]);
    }
}
/*--------------------------------------------
 
 --------------------------------------------*/
inline void Update::LoadUnLoadUpdateSelfComm::update_mult
(int& __icomm,int&,int&,vec**& vecs
,int& nvecs,int&)
{
    for(int ivec=0;ivec<nvecs;ivec++)
        vecs[ivec]->cpy_pst(snd_atms_lst[__icomm],snd_atms_lst_sz[__icomm]);
}
/*--------------------------------------------
 
 --------------------------------------------*/
inline void Update::LoadUnLoadUpdateSelfComm::update_sing
(int& __icomm,int&,int&,vec*& v)
{
    v->cpy_pst(snd_atms_lst[__icomm],snd_atms_lst_sz[__icomm]);
}
/*--------------------------------------------
 
 --------------------------------------------*/
inline void Update::LoadUnLoadUpdateSelfComm::xchng_buff
(int&,int& __snd_buff_sz,byte*& __snd_buff
,int&,int& __rcv_buff_sz,byte*& __rcv_buff)
{
    memcpy(__rcv_buff,__snd_buff,__rcv_buff_sz);
}
/*------------------------------------------------------------------------------------------------------------------------------------
 
 ------------------------------------------------------------------------------------------------------------------------------------*/
__Update::__Update(Atoms* atoms,
int& nupdt_vecs_,int& nxchng_vecs_):
world(atoms->world),
natms_lcl(atoms->natms_lcl),
natms_ph(atoms->natms_ph),
H(atoms->H),
depth_inv(atoms->depth_inv),

rank(atoms->comm.rank),
s_lo(atoms->comm.s_lo),
s_hi(atoms->comm.s_hi),

max_cut(atoms->max_cut),
x(atoms->x),

vecs(atoms->dynamic_vecs),
nvecs(atoms->ndynamic_vecs),
nupdt_vecs(nupdt_vecs_),
nxchng_vecs(nxchng_vecs_)
{

    Algebra::V_eq<__dim__*2>(&(atoms->comm.neigh[0][0]),&(neigh[0][0]));
    snd_buff=NULL;
    snd_buff_cpcty=0;
    
    rcv_buff=NULL;
    rcv_buff_cpcty=0;

    
    tot_updt_vecs_sz=0;
    for(int ivec=0;ivec<nupdt_vecs;ivec++)
        tot_updt_vecs_sz+=vecs[ivec]->byte_sz;
    
    int icurs=0;
    for(int idim=0;idim<__dim__;idim++)
    {
        self_comm[idim]=(rank==neigh[idim][0] && rank==neigh[idim][1]);
        
        // snd_to_bhnd && rcv_fm_frnt
        pbc_correction[idim][0]=(atoms->comm.coords[idim]==atoms->comm.dims[idim]-1);
        icurs++;
        
        // snd_to_frnt && rcv_fm_bhnd
        pbc_correction[idim][1]=(atoms->comm.coords[idim]==0);
        icurs++;
    }

    tot_ncomms=0;
    snd_atms_lst=NULL;
    snd_atms_lst_cpcty=snd_atms_lst_sz=rcv_atms_lst_sz=NULL;
}
/*--------------------------------------------
 
 --------------------------------------------*/
__Update::~__Update()
{
    delete [] rcv_buff;
    delete [] snd_buff;
    
    for(int i=0;i<tot_ncomms;i++)
        delete [] snd_atms_lst[i];
    
    delete [] snd_atms_lst;
    delete [] snd_atms_lst_cpcty;
    delete [] snd_atms_lst_sz;
    delete [] rcv_atms_lst_sz;
    
    
}
/*--------------------------------------------
 
 --------------------------------------------*/
void __Update::reset()
{
    for(int idim=0;idim<__dim__;idim++)
        max_cut_s[idim]=max_cut*depth_inv[idim];
    
    int icurs=0;
    int tot_ncomms_=0;
    for(int idim=0;idim<__dim__;idim++)
    {
        // snd_to_bhnd && rcv_fm_frnt
        s_bnd[idim][0]=s_lo[idim]+max_cut_s[idim];
        tot_ncomms_+=static_cast<int>(max_cut_s[idim]/(s_hi[idim]-s_lo[idim]))+1;
        ncomms[idim][0]=tot_ncomms_;
        icurs++;
        // snd_to_frnt && rcv_fm_bhnd
        s_bnd[idim][1]=s_hi[idim]-max_cut_s[idim];
        tot_ncomms_+=static_cast<int>(max_cut_s[idim]/(s_hi[idim]-s_lo[idim]))+1;
        ncomms[idim][1]=tot_ncomms_;
        icurs++;
    }
    
    if(tot_ncomms_==tot_ncomms)
        return;
    
    for(int i=0;i<tot_ncomms;i++)
        delete [] snd_atms_lst[i];
    delete [] snd_atms_lst;
    delete [] snd_atms_lst_cpcty;
    delete [] snd_atms_lst_sz;
    delete [] rcv_atms_lst_sz;
    
    tot_ncomms=tot_ncomms_;
    snd_atms_lst=new int*[tot_ncomms];
    snd_atms_lst_cpcty=new int[tot_ncomms];
    snd_atms_lst_sz=new int[tot_ncomms];
    rcv_atms_lst_sz=new int[tot_ncomms];
    for(int i=0;i<tot_ncomms;i++)
    {
        snd_atms_lst_cpcty[i]=0;
        snd_atms_lst[i]=NULL;
    }
}
/*--------------------------------------------
 
 --------------------------------------------*/
void __Update::list()
{
    natms_ph=0;
    max_snd_atms_lst_sz=max_rcv_atms_lst_sz=0;
    int icomm=0;
    __list<0>(icomm);
    natms_ph=x->vec_sz-natms_lcl;
    for(int ivec=nxchng_vecs;ivec<nvecs;ivec++)
    {
        vecs[ivec]->vec_sz=0;
        vecs[ivec]->resize(x->vec_sz);
    }
}
/*--------------------------------------------
 
 --------------------------------------------*/
void __Update::rm_rdndncy()
{
    snd_buff_sz=rcv_buff_sz=0;
    reserve_rcv_buff(max_snd_atms_lst_sz);
    reserve_snd_buff(natms_ph);
    
    byte* mark=snd_buff;
    /*-------temp_remove-------
    forcefield->neighbor->mark_redndnt_ph(mark);
    */
    int rcv_atms_lst_sz_;
    int snd_atms_lst_sz_=0;
    int snd_atms_lst_cpcty_=max_snd_atms_lst_sz;
    int* snd_atms_lst_=NULL;
    if(snd_atms_lst_cpcty_) snd_atms_lst_=new int[snd_atms_lst_cpcty_];
    
    int nlocomm;
    byte* mark_=mark+natms_ph;
    int icurs=2*__dim__-1;
    int jcomm=tot_ncomms-1;
    max_snd_atms_lst_sz=max_rcv_atms_lst_sz=0;
    for(int idim=__dim__-1;idim>-1;idim--)
    {
        for(int idir=1;idir>-1;idir--)
        {
            if(icurs!=0)
                nlocomm=ncomms[idim][idir];
            else
                nlocomm=0;
            
            while(jcomm>nlocomm-1)
            {
                mark_-=rcv_atms_lst_sz[jcomm];
                if(self_comm[idim])
                    self_xchng_buff(neigh[idim][1-idir],rcv_atms_lst_sz[jcomm],mark_,neigh[idim][idir],snd_atms_lst_sz[jcomm],rcv_buff);
                else
                    xchng_buff(neigh[idim][1-idir],rcv_atms_lst_sz[jcomm],mark_,neigh[idim][idir],snd_atms_lst_sz[jcomm],rcv_buff);
                
                snd_atms_lst_sz_=0;
                for(int i=0; i<snd_atms_lst_sz[jcomm];i++)
                {
                    if(rcv_buff[i]=='1')
                    {
                        if(snd_atms_lst[jcomm][i]>=natms_lcl)
                            mark[snd_atms_lst[jcomm][i]-natms_lcl]='1';
                        snd_atms_lst_[snd_atms_lst_sz_++]=snd_atms_lst[jcomm][i];
                    }
                }
                memcpy(snd_atms_lst[jcomm],snd_atms_lst_,snd_atms_lst_sz_*sizeof(int));
                snd_atms_lst_sz[jcomm]=snd_atms_lst_sz_;
                
                rcv_atms_lst_sz_=0;
                for(int i=0; i<rcv_atms_lst_sz[jcomm];i++)
                    if(mark_[i]=='1')
                        rcv_atms_lst_sz_++;
                rcv_atms_lst_sz[jcomm]=rcv_atms_lst_sz_;
                
                max_snd_atms_lst_sz=MAX(max_snd_atms_lst_sz,snd_atms_lst_sz[jcomm]);
                max_rcv_atms_lst_sz=MAX(max_rcv_atms_lst_sz,rcv_atms_lst_sz[jcomm]);
                
                jcomm--;
            }
            icurs--;
        }
    }
    
    delete [] snd_atms_lst_;
    
    int old_2_new_cpcty=natms_lcl+natms_ph;
    int* old_2_new=NULL;
    if(old_2_new_cpcty) old_2_new=new int[old_2_new_cpcty];
    
    int list_sz=0;
    int list_cpcty=natms_ph;
    int* list=NULL;
    if(list_cpcty) list=new int[list_cpcty];
    
    for(int iatm=0;iatm<natms_lcl;iatm++)
        old_2_new[iatm]=iatm;
    
    icurs=natms_lcl;
    for(int iatm=natms_lcl;iatm<natms_lcl+natms_ph;iatm++)
        if(mark[iatm-natms_lcl]=='1')
        {
            old_2_new[iatm]=icurs++;
            list[list_sz++]=iatm;
        }
    
    int new_natms_ph=list_sz;

    for(int icomm=0;icomm<tot_ncomms;icomm++)
        for(int i=0; i<snd_atms_lst_sz[icomm];i++)
            snd_atms_lst[icomm][i]=old_2_new[snd_atms_lst[icomm][i]];
    /*-------temp_remove-------
    forcefield->neighbor->rename_atoms(old_2_new);
     */
    delete [] old_2_new;
    
    int* list_=list;
    
    int vec_sz=natms_lcl;
    while(*list_==natms_lcl+icurs)
    {
        list_++;
        vec_sz++;
        list_sz--;
    }
    
    for(int ivec=0;ivec<nupdt_vecs;ivec++)
    {
        vecs[ivec]->vec_sz=vec_sz;
        vecs[ivec]->cpy_pst(list_,list_sz);
    }
    
    delete [] list;
    
    natms_ph=new_natms_ph;
}
/*------------------------------------------------------------------------------------------------------------------------------------
 
 ------------------------------------------------------------------------------------------------------------------------------------*/
inline void __Update::xchng_buff
(int& __snd_p,int& __snd_buff_sz,byte*& __snd_buff
,int& __rcv_p,int& __rcv_buff_sz,byte*& __rcv_buff)
{
    MPI_Sendrecv(__snd_buff,__snd_buff_sz,MPI_BYTE,__snd_p,0,
                 __rcv_buff,__rcv_buff_sz,MPI_BYTE,__rcv_p,0,
                 world,MPI_STATUS_IGNORE);
}
/*------------------------------------------------------------------------------------------------------------------------------------
 
 ------------------------------------------------------------------------------------------------------------------------------------*/

inline void __Update::self_xchng_buff
(int&,int& __snd_buff_sz,byte*& __snd_buff
,int&,int& __rcv_buff_sz,byte*& __rcv_buff)
{
    memcpy(__rcv_buff,__snd_buff,__rcv_buff_sz);
}
