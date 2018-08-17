#ifndef __MAPP__exchange__
#define __MAPP__exchange__
#include "global.h"
#include <mpi.h>
/*------------------------------------------------------------------
 _____  __    __  _____   _   _       ___   __   _   _____   _____  
| ____| \ \  / / /  ___| | | | |     /   | |  \ | | /  ___| | ____| 
| |__    \ \/ /  | |     | |_| |    / /| | |   \| | | |     | |__   
|  __|    }  {   | |     |  _  |   / / | | | |\   | | |  _  |  __|  
| |___   / /\ \  | |___  | | | |  / /  | | | | \  | | |_| | | |___  
|_____| /_/  \_\ \_____| |_| |_| /_/   |_| |_|  \_| \_____/ |_____|
 ------------------------------------------------------------------*/
namespace MAPP_NS
{
    template<typename> class Vec;
    class vec;
    class Exchange
    {
    private:
        
        /*things that reference cannot be removed*/
        int& natms_lcl;
        Vec<type0>*& x;
        unsigned long& xchng_id;
        
        /*things that reference cannot be removed*/
        const int rank;
        int neigh[__dim__][2];
        type0 (&s_lo)[__dim__];
        type0 (&s_hi)[__dim__];
        
        
        static constexpr int buff_grw=1024;
        byte* snd_buff[2];
        int snd_buff_sz[2];
        int snd_buff_cpcty[2];
        
        byte* rcv_buff;
        int rcv_buff_sz;
        int rcv_buff_cpcty;
        
        void load(int&,int);
        void load(byte*&,int&);
        int xchng_buff(int,int);
        
        vec**& vecs;
        int& nvecs;
        int& nxchng_vecs;
        int tot_xchng_sz;

        MPI_Comm& world;
    protected:
    public:        
        Exchange(class Atoms*,int&);
        ~Exchange();
        void full_xchng();
        void full_xchng_all();
    };
}
using namespace MAPP_NS;
/*------------------------------------------------
 _   _   _____   _____       ___   _____   _____  
| | | | |  _  \ |  _  \     /   | |_   _| | ____| 
| | | | | |_| | | | | |    / /| |   | |   | |__   
| | | | |  ___/ | | | |   / / | |   | |   |  __|  
| |_| | | |     | |_| |  / /  | |   | |   | |___  
\_____/ |_|     |_____/ /_/   |_|   |_|   |_____|
 ------------------------------------------------*/
namespace MAPP_NS
{
    
    class OldUpdate
    {
    private:
        /*things that reference cannot be removed*/
        int& natms_lcl;
        int& natms_ph;
        
        /*things that reference can be removed*/
        type0 (&H)[__dim__][__dim__];
        type0 (&depth_inv)[__dim__];
        type0 (&s_lo)[__dim__];
        type0 (&s_hi)[__dim__];
        const type0& max_cut;
        
        type0 max_cut_s[__dim__];
        const int rank;
        int neigh[__dim__][2];
        
        Vec<type0>*& x;
        
        vec**& vecs;
        int& nvecs;
        int& nxchng_vecs;
        int& nupdt_vecs;
        int tot_updt_vecs_sz;
        
        int tot_ncomms;
        int ncomms[__dim__][2];
        bool pbc_correction[__dim__][2];
        type0 s_bnd[__dim__][2];
        
        int** snd_atms_lst;
        int* snd_atms_lst_sz;
        int* snd_atms_lst_cpcty;
        int max_snd_atms_lst_sz;
        static constexpr int snd_atms_lst_grw=8;

        int* rcv_atms_lst_sz;
        int max_rcv_atms_lst_sz;

        byte* snd_buff;        
        int snd_buff_sz;
        int snd_buff_cpcty;
        static constexpr int snd_buff_grw=1024;

        byte* rcv_buff;
        int rcv_buff_sz;
        int rcv_buff_cpcty;
        static constexpr int rcv_buff_grw=1024;
        
        void add_to_snd_lst(int&,int&);
        void reserve_snd_buff(int);
        void reserve_rcv_buff(int);
        
        class LoadUnLoadUpdate;
        class LoadUnLoadUpdateComm;
        class LoadUnLoadUpdateSelfComm;
        LoadUnLoadUpdate* comm_manager[__dim__];
        
    protected:
    public:
        OldUpdate(Atoms*,int&,int&);
        ~OldUpdate();
        void reset();
        void update(vec*,bool);
        void update(vec*,type0 (*)[__dim__]);
        void update(vec**,int,bool);
        void list();
        void rm_rdndncy();
    };
}
/*----------------------------------------------------------------------------------------
 
 ----------------------------------------------------------------------------------------*/
namespace MAPP_NS
{
    class OldUpdate::LoadUnLoadUpdate
    {
    private:
    protected:
    public:
        LoadUnLoadUpdate(){};
        virtual ~LoadUnLoadUpdate(){};
        virtual void load_unload(int&,int&,int&)=0;
        virtual void update_mult(int&,int&,int&,vec**&,int&,int&)=0;
        virtual void update_sing(int&,int&,int&,vec*&)=0;
        virtual void xchng_buff(int&,int&,byte*&,int&,int&,byte*&)=0;
    };
}
/*----------------------------------------------------------------------------------------
 
 ----------------------------------------------------------------------------------------*/
namespace MAPP_NS
{
    class OldUpdate::LoadUnLoadUpdateComm:
    public OldUpdate::LoadUnLoadUpdate
    {
    private:
        
        int**& snd_atms_lst;
        int*& snd_atms_lst_sz;
        
        int*& rcv_atms_lst_sz;
        
        vec**& vecs;
        int& nupdt_vecs;
        int& tot_updt_vecs_sz;
        MPI_Comm& world;
        
        byte*& snd_buff;
        int& snd_buff_sz;
        int& snd_buff_cpcty;
        
        byte*& rcv_buff;
        int& rcv_buff_sz;
        int& rcv_buff_cpcty;
        
        
        
    protected:
    public:
        LoadUnLoadUpdateComm(OldUpdate*,MPI_Comm&);
        void load_unload(int&,int&,int&);
        void update_mult(int&,int&,int&,vec**&,int&,int&);
        void update_sing(int&,int&,int&,vec*&);
        void xchng_buff(int&,int&,byte*&,int&,int&,byte*&);
    };
}
/*----------------------------------------------------------------------------------------
 
 ----------------------------------------------------------------------------------------*/
namespace MAPP_NS
{
    class OldUpdate::LoadUnLoadUpdateSelfComm:
    public OldUpdate::LoadUnLoadUpdate
    {
    private:
        int**& snd_atms_lst;
        int*& snd_atms_lst_sz;
        
        int*& rcv_atms_lst_sz;
        
        vec**& vecs;
        int& nupdt_vecs;
    protected:
    public:
        LoadUnLoadUpdateSelfComm(OldUpdate*);
        void load_unload(int&,int&,int&);
        void update_mult(int&,int&,int&,vec**&,int&,int&);
        void update_sing(int&,int&,int&,vec*&);
        void xchng_buff(int&,int&,byte*&,int&,int&,byte*&);
    };
}
/*------------------------------------------------
 _   _   _____   _____       ___   _____   _____
| | | | |  _  \ |  _  \     /   | |_   _| | ____|
| | | | | |_| | | | | |    / /| |   | |   | |__
| | | | |  ___/ | | | |   / / | |   | |   |  __|
| |_| | | |     | |_| |  / /  | |   | |   | |___
\_____/ |_|     |_____/ /_/   |_|   |_|   |_____|
 ------------------------------------------------*/
#include "atoms.h"
#include "xmath.h"
namespace MAPP_NS
{
    
    class Update
    {
    private:
        /*things that reference cannot be removed*/
        MPI_Comm& world;
        int& natms_lcl;
        int& natms_ph;
        
        /*things that reference can be removed*/
        type0 (&H)[__dim__][__dim__];
        type0 (&depth_inv)[__dim__];
        type0 (&s_lo)[__dim__];
        type0 (&s_hi)[__dim__];
        const type0& max_cut;
        
        type0 max_cut_s[__dim__];
        const int rank;
        int neigh[__dim__][2];
        
        Vec<type0>*& x;
        
        vec**& vecs;
        int& nvecs;
        int& nxchng_vecs;
        int& nupdt_vecs;
        int tot_updt_vecs_sz;
        
        int tot_ncomms;
        int ncomms[__dim__][2];
        bool pbc_correction[__dim__][2];
        bool self_comm[__dim__];
        type0 s_bnd[__dim__][2];
        
        int** snd_atms_lst[__dim__][2];
        int* snd_atms_lst_sz[__dim__][2];
        int* snd_atms_lst_cpcty[__dim__][2];
        int max_snd_atms_lst_sz;
        static constexpr int snd_atms_lst_grw=8;

        int* rcv_atms_lst_sz[__dim__][2];
        int max_rcv_atms_lst_sz;

        byte* snd_buff;
        int snd_buff_sz;
        int snd_buff_cpcty;
        static constexpr int snd_buff_grw=1024;

        byte* rcv_buff;
        int rcv_buff_sz;
        int rcv_buff_cpcty;
        static constexpr int rcv_buff_grw=1024;
        
        
        template<int idim>
        void start(int (*__neigh)[2],int* __dims,int* __coords)
        {
            self_comm[idim]=(rank==(*__neigh)[0] && rank==(*__neigh)[1]);
            
            // snd_to_bhnd && rcv_fm_frnt
            pbc_correction[idim][0]=(*__coords==*__dims-1);
            neigh[idim][0]=(*__neigh)[0];
            ncomms[idim][0]=0;
            snd_atms_lst[idim][0]=NULL;
            snd_atms_lst_sz[idim][0]=NULL;
            snd_atms_lst_cpcty[idim][0]=NULL;
            rcv_atms_lst_sz[idim][0]=NULL;
            
            // snd_to_frnt && rcv_fm_bhnd
            pbc_correction[idim][1]=(*__coords==0);
            ncomms[idim][1]=0;
            neigh[idim][1]=(*__neigh)[1];
            snd_atms_lst[idim][1]=NULL;
            snd_atms_lst_sz[idim][1]=NULL;
            snd_atms_lst_cpcty[idim][1]=NULL;
            rcv_atms_lst_sz[idim][1]=NULL;
            
            start<idim+1>(__neigh+1,__dims+1,__coords+1);
        }

        template<int idim,int idir>
        void dealloc()
        {
            for(int icomm=0;icomm<ncomms[idim][idir];icomm++)
                delete [] snd_atms_lst[idim][idir][icomm];
            delete [] snd_atms_lst[idim][idir];
            delete [] snd_atms_lst_sz[idim][idir];
            delete [] snd_atms_lst_cpcty[idim][idir];
            delete [] rcv_atms_lst_sz[idim][idir];
        }
        template<int idim,int idir>
        void dealloc_all()
        {
            dealloc<idim,idir>();
            dealloc_all<idim+idir,1-idir>();
        }
        
        
        template<int idim,int idir>
        void realloc(int __ncomms)
        {
            snd_atms_lst[idim][idir]=new int*[__ncomms];
            snd_atms_lst_sz[idim][idir]=new int[__ncomms];
            snd_atms_lst_cpcty[idim][idir]=new int[__ncomms];
            rcv_atms_lst_sz[idim][idir]=new int[__ncomms];
            for(int icomm=0;icomm<__ncomms;icomm++)
            {
                snd_atms_lst[idim][idir][icomm]=NULL;
                snd_atms_lst_cpcty[idim][idir][icomm]=0;
            }
            ncomms[idim][idir]=__ncomms;
        }
        template<int idim>
        void __reset(int& __tot_ncomms)
        {
            max_cut_s[idim]=max_cut*depth_inv[idim];
            
            int __ncomms;
            
            s_bnd[idim][0]=s_lo[idim]+max_cut_s[idim];
            __ncomms=static_cast<int>(max_cut_s[idim]/(s_hi[idim]-s_lo[idim]))+1;
            __tot_ncomms+=__ncomms;
            if(__ncomms!=ncomms[idim][0])
            {
                dealloc<idim,0>();
                realloc<idim,0>(__ncomms);
            }
            
            s_bnd[idim][1]=s_hi[idim]-max_cut_s[idim];
            __ncomms=static_cast<int>(max_cut_s[idim]/(s_hi[idim]-s_lo[idim]))+1;
            __tot_ncomms+=__ncomms;
            if(__ncomms!=ncomms[idim][1])
            {
                dealloc<idim,1>();
                realloc<idim,1>(__ncomms);
            }
            __reset<idim+1>(__tot_ncomms);
        }
        

        void reserve_rcv_buff(int xtra)
        {
            if(rcv_buff_cpcty<xtra+rcv_buff_sz)
            {
                delete [] rcv_buff;
                rcv_buff=new byte[xtra+rcv_buff_sz+rcv_buff_grw];
                rcv_buff_cpcty=xtra+rcv_buff_sz+rcv_buff_grw;
            }
        }
        void reserve_snd_buff(int xtra)
        {
            if(snd_buff_cpcty<xtra+snd_buff_sz)
            {
                delete [] snd_buff;
                snd_buff=new byte[xtra+snd_buff_sz+snd_buff_grw];
                snd_buff_cpcty=xtra+snd_buff_sz+snd_buff_grw;
            }
        }
        
        template<int idim,int idir>
        void add_to_snd_lst(int& __icomm,int& iatm)
        {
            if(snd_atms_lst_sz[idim][idir][__icomm]+1>snd_atms_lst_cpcty[idim][idir][__icomm])
            {
                int* tmp_lst=new int[snd_atms_lst_sz[idim][idir][__icomm]+1+snd_atms_lst_grw];
                memcpy(tmp_lst,snd_atms_lst[idim][idir][__icomm],snd_atms_lst_sz[idim][idir][__icomm]*sizeof(int));
                delete [] snd_atms_lst[idim][idir][__icomm];
                snd_atms_lst[idim][idir][__icomm]=tmp_lst;
                snd_atms_lst_cpcty[idim][idir][__icomm]=snd_atms_lst_sz[idim][idir][__icomm]+1+snd_atms_lst_grw;
            }
            snd_atms_lst[idim][idir][__icomm][snd_atms_lst_sz[idim][idir][__icomm]]=iatm;
            snd_atms_lst_sz[idim][idir][__icomm]++;
        }

        

        template<int idim,int idir>
        void load_unload
        (int& __icomm)
        {
            snd_buff_sz=snd_atms_lst_sz[idim][idir][__icomm]*tot_updt_vecs_sz;
            if(snd_buff_cpcty<snd_buff_sz)
            {
                delete [] snd_buff;
                snd_buff=new byte[snd_buff_sz+snd_buff_grw];
                snd_buff_cpcty=snd_buff_sz+snd_buff_grw;
            }
            
            byte* tmp_snd_buff=snd_buff;
            for(int ivec=0;ivec<nupdt_vecs;ivec++)
                vecs[ivec]->cpy(tmp_snd_buff,snd_atms_lst[idim][idir][__icomm],snd_atms_lst_sz[idim][idir][__icomm]);
            
            MPI_Sendrecv(&snd_atms_lst_sz[idim][idir][__icomm],1,MPI_INT,neigh[idim][idir],0,
                         &rcv_atms_lst_sz[idim][idir][__icomm],1,MPI_INT,neigh[idim][1-idir],0,
                         world,MPI_STATUS_IGNORE);
            
            rcv_buff_sz=rcv_atms_lst_sz[idim][idir][__icomm]*tot_updt_vecs_sz;
            if(rcv_buff_cpcty<rcv_buff_sz)
            {
                delete [] rcv_buff;
                rcv_buff=new byte[rcv_buff_sz+rcv_buff_grw];
                rcv_buff_cpcty=rcv_atms_lst_sz[idim][idir][__icomm]*tot_updt_vecs_sz+rcv_buff_grw;
            }
            
            MPI_Sendrecv(snd_buff,snd_buff_sz,MPI_BYTE,neigh[idim][idir],0,
                         rcv_buff,rcv_buff_sz,MPI_BYTE,neigh[idim][1-idir],0,
                         world,MPI_STATUS_IGNORE);
            
            
            byte* tmp_rcv_buff=rcv_buff;
            for(int ivec=0;ivec<nupdt_vecs;ivec++)
            {
                vecs[ivec]->reserve(rcv_atms_lst_sz[idim][idir][__icomm]);
                vecs[ivec]->pst(tmp_rcv_buff,rcv_atms_lst_sz[idim][idir][__icomm]);
            }
        }

        
        
        template<int idim,int idir>
        void self_load_unload
        (int& __icomm)
        {
            rcv_atms_lst_sz[idim][idir][__icomm]=snd_atms_lst_sz[idim][idir][__icomm];
            
            for(int ivec=0;ivec<nupdt_vecs;ivec++)
            {
                vecs[ivec]->reserve(rcv_atms_lst_sz[idim][idir][__icomm]);
                vecs[ivec]->cpy_pst(snd_atms_lst[idim][idir][__icomm],rcv_atms_lst_sz[idim][idir][__icomm]);
            }
        }
        
        template<int idim,int idir>
        void update_var_cpy(int& __icomm,byte*& __snd_buff,vec* __v)
        {
            __v->cpy(__snd_buff,snd_atms_lst[idim][idir][__icomm],snd_atms_lst_sz[idim][idir][__icomm]);
        }
        template<int idim,int idir,class...VS>
        void update_var_cpy(int& __icomm,byte*& __snd_buff,vec* __v,VS*... __vs)
        {
            update_var_cpy<idim,idir>(__icomm,__snd_buff,__v);
            update_var_cpy<idim,idir>(__icomm,__snd_buff,__vs...);
        }
        
        template<int idim,int idir>
        void update_var_pst(int& __icomm,byte*& __rcv_buff,vec* __v)
        {
            __v->pst(__rcv_buff,rcv_atms_lst_sz[idim][idir][__icomm]);
        }
        
        template<int idim,int idir,class...VS>
        void update_var_pst(int& __icomm,byte*& __rcv_buff,vec* __v,VS* ... __vs)
        {
            update_var_pst<idim,idir>(__icomm,__rcv_buff,__v);
            update_var_pst<idim,idir>(__icomm,__rcv_buff,__vs...);
        }
        
        
        template<int idim,int idir>
        void update_var(int& __icomm,int& __vecs_byte_sz,vec* __v)
        {
            byte* tmp_snd_buff=snd_buff;
            __v->cpy(tmp_snd_buff,snd_atms_lst[idim][idir][__icomm],snd_atms_lst_sz[idim][idir][__icomm]);
            
            MPI_Sendrecv(snd_buff,snd_atms_lst_sz[idim][idir][__icomm]*__v->byte_sz,MPI_BYTE,neigh[idim][idir],0,
                         __v->end(),rcv_atms_lst_sz[idim][idir][__icomm]*__v->byte_sz,MPI_BYTE,neigh[idim][1-idir],0,
                         world,MPI_STATUS_IGNORE);
            __v->vec_sz+=rcv_atms_lst_sz[idim][idir][__icomm];

        }
        template<int idim,int idir,class...VS>
        void update_var(int& __icomm,int& __vecs_byte_sz,vec* __v,VS*... __vs)
        {
            byte* __snd_buff=snd_buff;
            update_var_cpy<idim,idir>(__icomm,__snd_buff,__v,__vs...);
            MPI_Sendrecv(snd_buff,snd_atms_lst_sz[idim][idir][__icomm]*__vecs_byte_sz,MPI_BYTE,neigh[idim][idir],0,
                         rcv_buff,rcv_atms_lst_sz[idim][idir][__icomm]*__vecs_byte_sz,MPI_BYTE,neigh[idim][1-idir],0,
                         world,MPI_STATUS_IGNORE);
            
            byte* __rcv_buff=rcv_buff;
            update_var_pst<idim,idir>(__icomm,__rcv_buff,__v,__vs...);
        }
        
        
        template<int idim,int idir>
        void self_update_var(int& __icomm,int&,vec* __v)
        {
            __v->cpy_pst(snd_atms_lst[idim][idir][__icomm],snd_atms_lst_sz[idim][idir][__icomm]);
        }
        template<int idim,int idir,class...VS>
        void self_update_var(int& __icomm,int& __vecs_byte_sz,vec* __v,VS*... __vs)
        {
            self_update_var<idim,idir>(__icomm,__vecs_byte_sz,__v);
            self_update_var<idim,idir>(__icomm,__vecs_byte_sz,__vs...);
        }
        
        
        template<int idim,int idir,class F,class ...VS>
        void __update(F&& f,int& tot_byte_sz,VS*... vs)
        {
            for(int icomm=0;icomm<ncomms[idim][idir];icomm++)
            {
                if(self_comm[idim])
                    self_update_var<idim,idir>(icomm,tot_byte_sz,vs...);
                else
                    update_var<idim,idir>(icomm,tot_byte_sz,vs...);

                
                if(pbc_correction[idim][idir])
                {
                    type0* x_vec=x->end()-__dim__;
                    for(int iatm=0;iatm<rcv_atms_lst_sz[idim][idir][icomm];iatm++,x_vec-=__dim__)
                        f(H[idim],x_vec);
                }
            }
        }
        
        template<int idim,int idir,class F>
        void __update(F&& f,Vec<type0>* __v,type0* __dH)
        {
            for(int icomm=0;icomm<ncomms[idim][idir];icomm++)
            {
                if(self_comm[idim])
                    self_update_var<idim,idir>(icomm,__v->byte_sz,__v);
                else
                    update_var<idim,idir>(icomm,__v->byte_sz,__v);
                
                
                if(pbc_correction[idim][idir])
                {
                    type0* __vec=__v->end()-__dim__;
                    // I think this is what it is suppose to be but older version says the one after
                    // gotta check later
                    /*
                    for(int iatm=0;iatm<rcv_atms_lst_sz[idim][idir][icomm];iatm++,__vec-=__dim__)
                        f(H[idim],__dH,__vec);
                     */
                    for(int iatm=0;iatm<rcv_atms_lst_sz[idim][idir][icomm];iatm++,__vec-=__dim__)
                        f(H[idim],__vec);
                }
            }
        }
        
        
        template<int idim,int idir,class ...VS>
        void __update(int& tot_byte_sz,VS*... vs)
        {
            for(int icomm=0;icomm<ncomms[idim][idir];icomm++)
            {
                if(self_comm[idim])
                    self_update_var<idim,idir>(icomm,tot_byte_sz,vs...);
                else
                    update_var<idim,idir>(icomm,tot_byte_sz,vs...);
            }
        }
        
        
        
        template<int idim,int idir,class F0,class F1>
        void ___list(int last_atm,F0&& f0,F1&& f1)
        {
            int lo_atm=0;
            int hi_atm=last_atm;
            type0* x_vec;
            for(int icomm=0;icomm<ncomms[idim][idir];icomm++)
            {
                
                snd_buff_sz=0;
                rcv_atms_lst_sz[idim][idir][icomm]=snd_atms_lst_sz[idim][idir][icomm]=0;
                x_vec=x->begin()+__dim__*lo_atm+idim;
                for(int iatm=lo_atm;iatm<hi_atm;iatm++,x_vec+=__dim__)
                    if(f0(*x_vec,s_bnd[idim][idir]))
                        add_to_snd_lst<idim,idir>(icomm,iatm);
                
                if(self_comm[idim])
                    self_load_unload<idim,idir>(icomm);
                else
                    load_unload<idim,idir>(icomm);
                
                
                lo_atm=x->vec_sz-rcv_atms_lst_sz[idim][idir][icomm];
                hi_atm=x->vec_sz;
                if(pbc_correction[idim][idir])
                {
                    x_vec=x->begin()+__dim__*lo_atm+idim;
                    for(int iatm=lo_atm;iatm<hi_atm;iatm++,x_vec+=__dim__)
                        f1(*x_vec);
                }
                max_snd_atms_lst_sz=MAX(max_snd_atms_lst_sz,snd_atms_lst_sz[idim][idir][icomm]);
                max_rcv_atms_lst_sz=MAX(max_rcv_atms_lst_sz,rcv_atms_lst_sz[idim][idir][icomm]);
            }
            
        }
        template<int idim>
        void __list()
        {
            int last_atm=x->vec_sz;
            ___list<idim,0>(last_atm,[](const type0& l,const type0& r)->bool{return (l<r);},[](type0& __x){++__x;});
            ___list<idim,1>(last_atm,[](const type0& l,const type0& r)->bool{return (l>=r);},[](type0& __x){--__x;});
            __list<idim+1>();
        }
        
        
        template<int idim>
        class Helper
        {
        public:
            template<class... VS>
            static void update_w_x(Update& update,int& tot_byte_sz,VS*... vs)
            {
                update.__update<idim,0>(Algebra::V_add<idim+1,type0>,tot_byte_sz,vs...);
                update.__update<idim,1>(Algebra::V_sub<idim+1,type0>,tot_byte_sz,vs...);
                Helper<idim+1>::update_w_x(update,tot_byte_sz,vs...);
            }
            template<class... VS>
            static void update_wo_x(Update& update,int& tot_byte_sz,VS*... vs)
            {
                update.__update<idim,0>(tot_byte_sz,vs...);
                update.__update<idim,1>(tot_byte_sz,vs...);
                Helper<idim+1>::update_wo_x(update,tot_byte_sz,vs...);
            }
            
            
            static void update_w_x_w_dH(Update& update,Vec<type0>* __v,type0 (*__dH)[__dim__])
            {
                // I think this is what it is suppose to be but older version says the one after
                // gotta check later
                /*
                update.__update<idim,0>(Algebra::V_add2<idim+1,type0>,__v,*__dH);
                update.__update<idim,1>(Algebra::V_sub2<idim+1,type0>,__v,*__dH);
                */
                update.__update<idim,0>(Algebra::V_add<idim+1,type0>,__v,*__dH);
                update.__update<idim,1>(Algebra::V_sub<idim+1,type0>,__v,*__dH);
                
                Helper<idim+1>::update_w_x_w_dH(update,__v,__dH+1);
            }
        };
        
        
        
        
        int reset_vs(vec* __v)
        {
            __v->vec_sz=natms_lcl;
            return __v->byte_sz;
        }
        template<class ...VS>
        int reset_vs(vec* __v,VS*... __vs)
        {
            return reset_vs(__v)+reset_vs(__vs...);
        }

        template<int idim,int idir>
        void xchng_buff
        (int& __icomm,byte*& __snd_buff,byte*& __rcv_buff)
        {
            MPI_Sendrecv(__snd_buff,rcv_atms_lst_sz[idim][idir][__icomm],MPI_BYTE,neigh[idim][1-idir],0,
                         __rcv_buff,snd_atms_lst_sz[idim][idir][__icomm],MPI_BYTE,neigh[idim][idir],0,
                         world,MPI_STATUS_IGNORE);
        }
        template<int idim,int idir>
        void self_xchng_buff
        (int& __icomm,byte*& __snd_buff,byte*& __rcv_buff)
        {
           memcpy(__rcv_buff,__snd_buff,snd_atms_lst_sz[idim][idir][__icomm]);
        }
        
        
        template<int idim,int idir>
        void __rm_rdndncy(byte* mark,byte*& __mark,int* __snd_atms_lst)
        {
            int __snd_atms_lst_sz,__rcv_atms_lst_sz;
            for(int icomm=ncomms[idim][idir];icomm>-1;icomm--)
            {
                __mark-=rcv_atms_lst_sz[idim][idir][icomm];
                if(self_comm[idim])
                    self_xchng_buff<idim,idir>(icomm,__mark,rcv_buff);
                else
                    xchng_buff<idim,idir>(icomm,__mark,rcv_buff);
                __snd_atms_lst_sz=0;
                for(int i=0; i<snd_atms_lst_sz[idim][idir][icomm];i++)
                {
                    if(rcv_buff[i]=='1')
                    {
                        if(snd_atms_lst[idim][idir][icomm][i]>=natms_lcl)
                            mark[snd_atms_lst[idim][idir][icomm][i]-natms_lcl]='1';
                        __snd_atms_lst[__snd_atms_lst_sz++]=snd_atms_lst[idim][idir][icomm][i];
                    }
                }
                memcpy(snd_atms_lst[idim][idir][icomm],__snd_atms_lst,__snd_atms_lst_sz*sizeof(int));
                snd_atms_lst_sz[idim][idir][icomm]=__snd_atms_lst_sz;
                
                __rcv_atms_lst_sz=0;
                for(int i=0; i<rcv_atms_lst_sz[idim][idir][icomm];i++)
                    if(__mark[i]=='1')
                        __rcv_atms_lst_sz++;
                rcv_atms_lst_sz[idim][idir][icomm]=__rcv_atms_lst_sz;
                
                max_snd_atms_lst_sz=MAX(max_snd_atms_lst_sz,snd_atms_lst_sz[idim][idir][icomm]);
                max_rcv_atms_lst_sz=MAX(max_rcv_atms_lst_sz,rcv_atms_lst_sz[idim][idir][icomm]);
            }
            
            __rm_rdndncy<idim+idir-1,1-idir>(mark,__mark,__snd_atms_lst);
        }

        template<int idim,int idir>
        void __rm_rdndncy_old_2_new(int* old_2_new)
        {
            for(int icomm=0;icomm<ncomms[idim][idir];icomm++)
                for(int i=0; i<snd_atms_lst_sz[idim][idir][icomm];i++)
                    snd_atms_lst[idim][idir][icomm][i]=old_2_new[snd_atms_lst[idim][idir][icomm][i]];
            __rm_rdndncy_old_2_new<idim+idir,1-idir>(old_2_new);
        }
        
    protected:
    public:
        Update(Atoms* atoms,
        int& __nupdt_vecs,int& __nxchng_vecs):
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
        nupdt_vecs(__nupdt_vecs),
        nxchng_vecs(__nxchng_vecs)
        {
            start<0>(atoms->comm.neigh,atoms->comm.dims,atoms->comm.coords);
            snd_buff=NULL;
            snd_buff_cpcty=0;
            
            rcv_buff=NULL;
            rcv_buff_cpcty=0;
            
            tot_ncomms=0;
            tot_updt_vecs_sz=0;
            for(int ivec=0;ivec<nupdt_vecs;ivec++)
                tot_updt_vecs_sz+=vecs[ivec]->byte_sz;
            
        }
        
        ~Update()
        {
            delete [] rcv_buff;
            delete [] snd_buff;
            dealloc_all<0,0>();
        }
        void reset()
        {
            int __tot_ncomms=0;
            __reset<0>(__tot_ncomms);
            tot_ncomms=__tot_ncomms;
        }
        void list()
        {
            natms_ph=0;
            max_snd_atms_lst_sz=max_rcv_atms_lst_sz=0;
            __list<0>();
            natms_ph=x->vec_sz-natms_lcl;
            for(int ivec=nxchng_vecs;ivec<nvecs;ivec++)
            {
                vecs[ivec]->vec_sz=0;
                vecs[ivec]->resize(x->vec_sz);
            }
        }
        
        template<class ...VS>
        void update_w_x(VS*... __vs)
        {
            int tot_byte_sz=reset_vs(__vs...);
            snd_buff_sz=rcv_buff_sz=0;
            reserve_snd_buff(tot_byte_sz*max_snd_atms_lst_sz);
            reserve_rcv_buff(tot_byte_sz*max_rcv_atms_lst_sz);
            Helper<0>::update_w_x(*this,tot_byte_sz,__vs...);
            
        }
        
        template<class ...VS>
        void update_wo_x(VS*... __vs)
        {
            int tot_byte_sz=reset_vs(__vs...);
            snd_buff_sz=rcv_buff_sz=0;
            reserve_snd_buff(tot_byte_sz*max_snd_atms_lst_sz);
            reserve_rcv_buff(tot_byte_sz*max_rcv_atms_lst_sz);
            Helper<0>::update_wo_x(*this,tot_byte_sz,__vs...);
            
        }
        
        void update_w_x_w_dH(Vec<type0>*__v,type0 (*__dH)[__dim__])
        {
            Helper<0>::update_w_x_w_dH(*this,__v,__dH);
        }
        
        void rm_rdndncy();
    };
    
    template<>
    inline void Update::start<__dim__>(int (*)[2],int*,int*){};
    template<>
    inline void Update::__list<__dim__>(){};
    template<>
    inline void Update::__reset<__dim__>(int&){};
    template<>
    inline void Update::dealloc_all<__dim__,0>(){};
    template<>
    inline void Update::__rm_rdndncy<-1,1>(byte*,byte*&,int*){}
    template<>
    inline void Update::__rm_rdndncy_old_2_new<__dim__,0>(int*){}
    template<>
    class Update::Helper<__dim__>
    {
    public:
        template<class... VS>
        static void update_w_x(Update&,int&,VS*...){}
        template<class... VS>
        static void update_wo_x(Update&,int&,VS*...){}
        static void update_w_x_w_dH(Update&,Vec<type0>*,type0 (*)[__dim__]){}
    };
}
#endif 
