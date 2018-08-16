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
    
    class Update
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
        Update(Atoms*,int&,int&);
        ~Update();
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
    class Update::LoadUnLoadUpdate
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
    class Update::LoadUnLoadUpdateComm:
    public Update::LoadUnLoadUpdate
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
        LoadUnLoadUpdateComm(Update*,MPI_Comm&);
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
    class Update::LoadUnLoadUpdateSelfComm:
    public Update::LoadUnLoadUpdate
    {
    private:
        int**& snd_atms_lst;
        int*& snd_atms_lst_sz;
        
        int*& rcv_atms_lst_sz;
        
        vec**& vecs;
        int& nupdt_vecs;
    protected:
    public:
        LoadUnLoadUpdateSelfComm(Update*);
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
    
    class __Update
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
        
        void add_to_snd_lst(int& icomm,int& iatm)
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
        
        
        void xchng_buff(int&,int&,byte*&,int&,int&,byte*&);
        void self_xchng_buff(int&,int&,byte*&,int&,int&,byte*&);
        
        
        void load_unload
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
        
        void update_var(int& __icomm,int& __snd_p,int& __rcv_p,int& __vecs_byte_sz,vec* __v)
        {
            byte* tmp_snd_buff=snd_buff;
            __v->cpy(tmp_snd_buff,snd_atms_lst[__icomm],snd_atms_lst_sz[__icomm]);
            
            MPI_Sendrecv(snd_buff,snd_atms_lst_sz[__icomm]*__vecs_byte_sz,MPI_BYTE,__snd_p,0,
                         __v->end(),rcv_atms_lst_sz[__icomm]*__vecs_byte_sz,MPI_BYTE,__rcv_p,0,
                         world,MPI_STATUS_IGNORE);
            __v->vec_sz+=rcv_atms_lst_sz[__icomm];
        }
        
        void self_load_unload
        (int& __icomm,int&,int&)
        {
            rcv_atms_lst_sz[__icomm]=snd_atms_lst_sz[__icomm];
            
            for(int ivec=0;ivec<nupdt_vecs;ivec++)
            {
                vecs[ivec]->reserve(rcv_atms_lst_sz[__icomm]);
                vecs[ivec]->cpy_pst(snd_atms_lst[__icomm],rcv_atms_lst_sz[__icomm]);
            }
        }
        
        void update_var_cpy(int& __icomm,byte*& __snd_buff,vec* __v)
        {
            __v->cpy(__snd_buff,snd_atms_lst[__icomm],snd_atms_lst_sz[__icomm]);
        }
        template<class...VS>
        void update_var_cpy(int& __icomm,byte*& __snd_buff,vec* __v,VS*... __vs)
        {
            update_var_cpy(__icomm,__snd_buff,__v);
            update_var_cpy(__icomm,__snd_buff,__vs...);
        }
        
        void update_var_pst(int& __icomm,byte*& __rcv_buff,vec* __v)
        {
            __v->pst(__rcv_buff,rcv_atms_lst_sz[__icomm]);
        }
        
        template<class...VS>
        void update_var_pst(int& __icomm,byte*& __rcv_buff,vec* __v,VS* ... __vs)
        {
            update_var_pst(__icomm,__rcv_buff,__v);
            update_var_pst(__icomm,__rcv_buff,__vs...);
        }
        
        template<class...VS>
        void update_var(int& __icomm,int& __snd_p,int& __rcv_p,int& __vecs_byte_sz,vec* __v,VS*... __vs)
        {
            byte* __snd_buff=snd_buff;
            update_var_cpy(__icomm,__snd_buff,__v,__vs...);
            MPI_Sendrecv(snd_buff,snd_atms_lst_sz[__icomm]*__vecs_byte_sz,MPI_BYTE,__snd_p,0,
                         rcv_buff,rcv_atms_lst_sz[__icomm]*__vecs_byte_sz,MPI_BYTE,__rcv_p,0,
                         world,MPI_STATUS_IGNORE);
            
            byte* __rcv_buff=rcv_buff;
            update_var_pst(__icomm,__rcv_buff,__v,__vs...);
        }
        
        

        void self_update_var(int& __icomm,int&,int&,int&,vec* __v)
        {
            __v->cpy_pst(snd_atms_lst[__icomm],snd_atms_lst_sz[__icomm]);
        }
        template<class...VS>
        void self_update_var(int& __icomm,int& __snd_p,int& __rcv_p,int& __vecs_byte_sz,vec* __v,VS*... __vs)
        {
            self_update_var(__icomm,__snd_p,__rcv_p,__vecs_byte_sz,__v);
            self_update_var(__icomm,__snd_p,__rcv_p,__vecs_byte_sz,__vs...);
        }
        
        
        template<int idim,int idir,class F,class ...VS>
        void __update_w_x(int& icomm,int& tot_byte_sz,F&& f,VS*... vs)
        {
            while(icomm<ncomms[idim][idir])
            {
                if(self_comm[idim])
                    self_update_var(icomm,neigh[idim][idir],neigh[idim][1-idir],tot_byte_sz,vs...);
                else
                    update_var(icomm,neigh[idim][idir],neigh[idim][1-idir],tot_byte_sz,vs...);

                
                if(pbc_correction[idim][idir])
                {
                    type0* x_vec=x->end()-__dim__;
                    for(int iatm=0;iatm<rcv_atms_lst_sz[icomm];iatm++,x_vec-=__dim__)
                        f(H[idim],x_vec);
                }
                icomm++;
            }
        }
        
        template<int idim,int idir,class ...VS>
        void __update_w_o_x(int& icomm,int& tot_byte_sz,VS*... vs)
        {
            while(icomm<ncomms[idim][idir])
            {
                if(self_comm[idim])
                    self_update_var(icomm,neigh[idim][idir],neigh[idim][1-idir],tot_byte_sz,vs...);
                else
                    update_var(icomm,neigh[idim][idir],neigh[idim][1-idir],tot_byte_sz,vs...);
                icomm++;
            }
        }
        
        
        
        template<int idim,int idir,class F0,class F1>
        void ___list(int last_atm,int& icomm,F0&& f0,F1&& f1)
        {
            int lo_atm=0;
            int hi_atm=last_atm;
            type0* x_vec;
            while(icomm<ncomms[idim][idir])
            {
                
                snd_buff_sz=0;
                rcv_atms_lst_sz[icomm]=snd_atms_lst_sz[icomm]=0;
                x_vec=x->begin()+__dim__*lo_atm+idim;
                for(int iatm=lo_atm;iatm<hi_atm;iatm++,x_vec+=__dim__)
                    if(f0(*x_vec,s_bnd[idim][idir]))
                        add_to_snd_lst(icomm,iatm);
                
                if(self_comm[idim])
                    self_load_unload(icomm,neigh[idim][idir],neigh[idim][1-idir]);
                else
                    load_unload(icomm,neigh[idim][idir],neigh[idim][1-idir]);
                
                
                lo_atm=x->vec_sz-rcv_atms_lst_sz[icomm];
                hi_atm=x->vec_sz;
                if(pbc_correction[idim][idir])
                {
                    x_vec=x->begin()+__dim__*lo_atm+idim;
                    for(int iatm=lo_atm;iatm<hi_atm;iatm++,x_vec+=__dim__)
                        f1(*x_vec);
                }
                max_snd_atms_lst_sz=MAX(max_snd_atms_lst_sz,snd_atms_lst_sz[icomm]);
                max_rcv_atms_lst_sz=MAX(max_rcv_atms_lst_sz,rcv_atms_lst_sz[icomm]);
                icomm++;
            }
            
        }
        template<int idim>
        void __list(int& icomm)
        {
            int last_atm=x->vec_sz;
            ___list<idim,0>(last_atm,icomm,[](const type0& l,const type0& r)->bool{return (l<r);},[](type0& __x){++__x;});
            ___list<idim,1>(last_atm,icomm,[](const type0& l,const type0& r)->bool{return (l>=r);},[](type0& __x){--__x;});
            __list<idim+1>(icomm);
        }
        
        
        template<int idim>
        class UpdateHelper
        {
        public:
            template<class... VS>
            static void update_w_x(__Update& update,int& icomm,int& tot_byte_sz,VS*... vs)
            {
                update.__update_w_x<idim,0>(icomm,tot_byte_sz,Algebra::V_add<idim+1,type0>,vs...);
                update.__update_w_x<idim,1>(icomm,tot_byte_sz,Algebra::V_sub<idim+1,type0>,vs...);
                UpdateHelper<idim+1>::update_w_x(update,icomm,tot_byte_sz,vs...);
            }
            template<class... VS>
            static void update_w_o_x(__Update& update,int& icomm,int& tot_byte_sz,VS*... vs)
            {
                update.__update_w_o_x<idim,0>(icomm,tot_byte_sz,vs...);
                update.__update_w_o_x<idim,1>(icomm,tot_byte_sz,vs...);
                UpdateHelper<idim+1>::update_w_o_x(update,icomm,tot_byte_sz,vs...);
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

        
    protected:
    public:
        __Update(Atoms*,int&,int&);
        ~__Update();
        void reset();
        void list();
        void rm_rdndncy();
        
        template<class ...VS>
        void update_w_x(VS*... __vs)
        {
            int tot_byte_sz=reset_vs(__vs...);
            snd_buff_sz=rcv_buff_sz=0;
            reserve_snd_buff(tot_byte_sz*max_snd_atms_lst_sz);
            reserve_rcv_buff(tot_byte_sz*max_rcv_atms_lst_sz);
            int icomm=0;
            UpdateHelper<0>::update_w_x(*this,icomm,tot_byte_sz,__vs...);
            
        }
        
        template<class ...VS>
        void update_w_o_x(VS*... __vs)
        {
            int tot_byte_sz=reset_vs(__vs...);
            snd_buff_sz=rcv_buff_sz=0;
            reserve_snd_buff(tot_byte_sz*max_snd_atms_lst_sz);
            reserve_rcv_buff(tot_byte_sz*max_rcv_atms_lst_sz);
            int icomm=0;
            UpdateHelper<0>::update_w_o_x(*this,icomm,tot_byte_sz,__vs...);
            
        }
    };
    
    template<>
    inline void __Update::__list<__dim__>(int&){};
    template<>
    class __Update::UpdateHelper<__dim__>
    {
    public:
        template<class... VS>
        static void update_w_x(__Update&,int&,int&,VS*...)
        {
        }
        template<class... VS>
        static void update_w_o_x(__Update&,int&,int&,VS*...)
        {}
    };
}
#endif 
