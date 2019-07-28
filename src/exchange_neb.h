#ifndef __MAPP__exchange_neb__
#define __MAPP__exchange_neb__
#include "atoms.h"
#include "global.h"
#include "mpi_compat.h"
#include "xmath.h"
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
    class UpdateNEB
    {
    private:
        MPI_Comm& world;
        
        int neigh[__dim__][2];
        int ncomms[__dim__][2];

        
        int comm_dims[__dim__];
        int comm_coords[__dim__];
        
        static constexpr int buff_grw=1024;
        
        
        static constexpr int id_sz=static_cast<int>(sizeof(id_type));
        static constexpr int i_sz=static_cast<int>(sizeof(int));
        static constexpr int ch_sz=static_cast<int>(sizeof(char));
        
        
        
        byte* rcv_buff;
        int rcv_buff_sz;
        int rcv_buff_cpcty;
        int rcv_max_sz;
        
        
        
        byte* snd_buff[2*__dim__+1];
        int snd_buff_strd[2*__dim__+1];
        int snd_buff_sz[2*__dim__+1];
        int snd_buff_cpcty[2*__dim__+1];
        int snd_max_sz[2*__dim__+1];
        int tot_snd_max_sz;
        int nsttc_atms;
        
        
        int** lsts;
        int* rcv_szs;
        int* snd_szs;
        int lsts_sz;
        int lsts_cpcty;
        static constexpr int lsts_grw=8;
        
        bool* same;
        int* vec0_lst;
        int vec0_lst_sz;
        int* vec1_lst;
        int vec1_lst_sz;
        
        
        byte* updt_snd_buff[2*__dim__+1];
        int updt_snd_buff_sz[2*__dim__+1];
        byte* updt_rcv_buff;
        int updt_strd;
        
        template<int n>
        void __rsrv_buffs(byte** buff,int* sz,int* cpcty,int* strd,const int* xtra)
        {
            if(*cpcty<*sz+*xtra**strd)
            {
                int __cpcty=*sz+*xtra**strd+buff_grw;
                byte* tmp_buff=new byte[__cpcty];
                memcpy(tmp_buff,*buff,*sz);
                delete [] *buff;
                *buff=tmp_buff;
                *cpcty=__cpcty;
            }
            
            __rsrv_buffs<n-1>(buff+1,sz+1,cpcty+1,strd+1,xtra+1);
        }
        
        template<int offset,int n>
        void rsrv_buffs(byte** buff,int* sz,int* cpcty,int* strd,const int* xtra)
        {
            __rsrv_buffs<n-offset>(buff+offset,sz+offset,cpcty+offset,strd+offset,xtra+offset);
        }
        
        
        void add2lsts(int* __lst,int& __snd_sz,int& __rcv_sz)
        {
            if(lsts_cpcty<lsts_sz+1)
            {
                int __lsts_cpcty=lsts_sz+1+lsts_grw;
                
                int** __lsts=new int*[__lsts_cpcty];
                memcpy(__lsts,lsts,lsts_sz*sizeof(int*));
                delete [] lsts;
                lsts=__lsts;
                
                int* __rcv_szs=new int[__lsts_cpcty];
                memcpy(__rcv_szs,rcv_szs,lsts_sz*sizeof(int));
                delete [] rcv_szs;
                rcv_szs=__rcv_szs;
                
                int* __snd_szs=new int[__lsts_cpcty];
                memcpy(__snd_szs,snd_szs,lsts_sz*sizeof(int));
                delete [] snd_szs;
                snd_szs=__snd_szs;
                
                lsts_cpcty=__lsts_cpcty;
            }
            
            lsts[lsts_sz]=__lst;
            rcv_szs[lsts_sz]=__rcv_sz;
            snd_szs[lsts_sz]=__snd_sz;
            ++lsts_sz;
        }
        
        template<int idim,int idir>
        void update_snd_buff(int*& lst,byte* __rcv_buff,const int& __rcv_sz,byte** __snd_buff,int* __snd_buff_sz,int* __snd_buff_strd)
        {
            int j;
            const int src_strd=(__dim__-idim)*i_sz+id_sz+ch_sz;
            for(int i=0;i<__rcv_sz;i++,__rcv_buff+=src_strd)
            {
                j=lst[i];
                memcpy(__snd_buff[j]+__snd_buff_sz[j]*__snd_buff_strd[j],
                       __rcv_buff+(src_strd-__snd_buff_strd[j]),__snd_buff_strd[j]);
                __snd_buff_sz[j]+=__snd_buff_strd[j];
            }
        }
        
        template<int idim>
        void ret_dsp(const type0* s,const int* coord,const int* dim,int* icoord)
        {
            *icoord=static_cast<int>(*s*static_cast<type0>(*dim));
            if(*s<static_cast<type0>(*icoord)/static_cast<type0>(*dim))
                ++(*icoord);
            
            *icoord-=*coord;
            if(*icoord>*dim/2)  *icoord-=*dim;
            if(*icoord<-*dim/2) *icoord+=*dim;
            ret_dsp<idim-1>(s+1,coord+1,dim+1,icoord+1);
        }
        
        template<int idim>
        int ret_lst(const int* n)
        {
            if(*n>0) return snd_to_frnt;
            if(*n<0) return snd_to_bhnd;
            return 2+ret_lst<idim-1>(n+1);
            
        }
        
        
        
        
        
        template<int idim,int idir>
        int* get_lst(byte* __rcv_buff,int __rcv_sz,int* xtra)
        {
            Algebra::zero<2*__dim__+1-idim*2-idir>(xtra+ idim*2+idir);
            int* lst=new int[__rcv_sz];
            
            const int strd=(__dim__-idim)*i_sz+id_sz+ch_sz;
            int* pos;
            int j;
            for(int i=0;i<__rcv_sz;i++,__rcv_buff+=strd)
            {
                pos=reinterpret_cast<int*>(__rcv_buff);
                // be careful here let me double check
                *pos+=1-2*idir;
                j=2*idim+ret_lst<idim>(pos);
                lst[i]=j;
                ++xtra[j];
            }
            return lst;
        }
        
        
        template<int idim,int idir>
        int* deal_w_rcv_buff(byte* __rcv_buff,int __rcv_sz)
        {
            if(__rcv_sz==0) return NULL;
            int xtra[2*__dim__+1];
            int* lst=get_lst<idim,idir>(__rcv_buff,__rcv_sz,xtra);
            rsrv_buffs<idim*2-idir,2*__dim__+1>(snd_buff,snd_buff_sz,snd_buff_cpcty,snd_buff_strd,xtra);
            update_snd_buff<idim,idir>(lst,__rcv_buff,__rcv_sz,snd_buff,snd_buff_sz,snd_buff_strd);
            return lst;
        }
        
        
        void init_lst(type0* b0,type0* x0,type0* b1,type0* x1,id_type* __id,int __natms_lcl)
        {
            
            Algebra::zero<2*__dim__+1>(snd_szs);
            lsts_sz=0;
            
            // consider the worst case scenario
            constexpr int strd=__dim__*i_sz+id_sz+ch_sz;
            if(rcv_buff_cpcty<__natms_lcl*2*strd)
            {
                delete [] rcv_buff;
                rcv_buff_cpcty=__natms_lcl*2*strd+buff_grw;
                rcv_buff=new byte[rcv_buff_cpcty];
            }
            rcv_buff_sz=0;

            
            // this is the worst case scenrio as well
            int* lst= __natms_lcl==0 ? NULL:new int[__natms_lcl*2];
            
            delete [] same;
            same= __natms_lcl==0 ? NULL:new bool[__natms_lcl];
            
            
            int xtra[2*__dim__+1];
            Algebra::zero<2*__dim__+1>(xtra);
            type0 s[__dim__];

            
            byte ibuff[2*__dim__*i_sz+2*id_sz+2*ch_sz];
            int* ibuff_dst0=reinterpret_cast<int*>(ibuff);
            int* ibuff_dst1=reinterpret_cast<int*>(ibuff+strd);
            id_type* ibuff_id0=reinterpret_cast<id_type*>(ibuff+__dim__*i_sz);
            id_type* ibuff_id1=reinterpret_cast<id_type*>(ibuff+__dim__*i_sz+strd);
            char* ibuff_cat0=reinterpret_cast<char*>(ibuff+__dim__*i_sz+id_sz);
            char* ibuff_cat1=reinterpret_cast<char*>(ibuff+__dim__*i_sz+id_sz+strd);
            
            int* __lst=lst;
            bool* __same=same;
            type0* __x0=x0;
            type0* __x1=x1;
            byte* __rcv_buff=rcv_buff;
            for(int i=0;i<__natms_lcl;i++)
            {
                Algebra::X2S<__dim__>(b0,__x0,s);
                ret_dsp<__dim__>(s,comm_coords,comm_dims,ibuff_dst0);
                
                Algebra::X2S<__dim__>(b1,__x1,s);
                ret_dsp<__dim__>(s,comm_coords,comm_dims,ibuff_dst1);
                
                
                if(Algebra::is_same<__dim__>(ibuff_dst0,ibuff_dst1))
                {
                    *__same=true;
                    *ibuff_id0=__id[i];
                    *ibuff_cat0=0;
                    *__lst=ret_lst<__dim__>(ibuff_dst0);
                    ++xtra[*__lst];
                    ++__lst;
                    
                    
                    memcpy(__rcv_buff,ibuff,strd);
                    __rcv_buff+=strd;
                    rcv_buff_sz+=strd;
                    
                }
                else
                {
                    *__same=false;
                    *ibuff_id0=__id[i];
                    *ibuff_cat0=-1;
                    *__lst=ret_lst<__dim__>(ibuff_dst0);
                    ++xtra[*__lst];
                    ++__lst;
                    
                    *ibuff_id1=__id[i];
                    *ibuff_cat1=1;
                    *__lst=ret_lst<__dim__>(ibuff_dst1);
                    ++xtra[*__lst];
                    ++__lst;
                    
                    memcpy(__rcv_buff,ibuff,strd*2);
                    __rcv_buff+=2*strd;
                    rcv_buff_sz+=2*strd;
                }
                
                

                ++__same;
                __x0+=__dim__;
                __x1+=__dim__;
            }

            
            

            
            rsrv_buffs<0,2*__dim__+1>(snd_buff,snd_buff_sz,snd_buff_cpcty,snd_buff_strd,xtra);
            int __rcv_sz=rcv_buff_sz/strd;

            update_snd_buff<__dim__,0>(lst,rcv_buff,__rcv_sz,snd_buff,snd_buff_sz,snd_buff_strd);
            int __snd_sz=0;
            add2lsts(__lst,__snd_sz,__rcv_sz);
            nsttc_atms=xtra[2*__dim__];
        }
        
        
        template<int idim,int idir>
        void __xchng_lst()
        {
            int max_snd_sz=1;
            const int strd=(__dim__-idim)*i_sz+id_sz+ch_sz;
            ncomms[idim][idir]=0;
            
            int __snd_sz,__rcv_sz;
            
            while(max_snd_sz!=0)
            {
                __snd_sz=snd_buff_sz[idim*2+idir]/strd;
                
                // figure out maximum number of components to be sent in this stage
                MPI_Allreduce(&__snd_sz,&max_snd_sz,1,MPI_INT,MPI_MAX,world);
                // if all of the sent numbers are zero we have nothing to do
                if(max_snd_sz==0) continue;

                MPI_Sendrecv(&__snd_sz,1,MPI_INT,neigh[idim][idir],0,
                             &__rcv_sz,1,MPI_INT,neigh[idim][1-idir],0,
                             world,MPI_STATUS_IGNORE);
                
                
                if(rcv_buff_cpcty<__rcv_sz*strd)
                {
                    delete [] rcv_buff;
                    rcv_buff_cpcty=__rcv_sz*strd+buff_grw;
                    rcv_buff=new byte[rcv_buff_cpcty];
                }
                
                
                MPI_Sendrecv(snd_buff[idim*2+idir],__snd_sz*strd,MPI_BYTE,neigh[idim][idir],0,
                             rcv_buff,__rcv_sz*strd,MPI_BYTE,neigh[idim][1-idir],0,
                             world,MPI_STATUS_IGNORE);
                
                
                //set the maximums
                snd_max_sz[idim*2+idir]=MAX(snd_max_sz[idim*2+idir],__snd_sz);
                rcv_max_sz=MAX(rcv_max_sz,__rcv_sz);
                
                //reset the active snd_buffer
                snd_buff_sz[idim*2+idir]=0;
               
                // deal with the rcv_buff and add lst to lsts
                add2lsts(deal_w_rcv_buff<idim,idir>(rcv_buff,__rcv_sz),__snd_sz,__rcv_sz);
                ++ncomms[idim][idir];
            }
        }
        
        template<int idim>
        void xchng_lst()
        {
            __xchng_lst<__dim__-idim,0>();
            __xchng_lst<__dim__-idim,1>();
            return xchng_lst<idim-1>();
        }
        
        
        void fin_lst()
        {
            const int strd=id_sz+ch_sz;
            byte* buff=snd_buff[2*__dim__];
            int n0=nsttc_atms;
            int n=snd_buff_sz[2*__dim__]/strd;
            int n1=n-n0;
            
            int* key_1= n1==0 ? NULL:new int[n1];
            for(int i=0;i<n1;i++) key_1[i]=n0+i;
            XMath::quicksort(key_1,key_1+n1,
            [&buff](int* ikey,int* jkey)
            {return (*reinterpret_cast<id_type*>(buff+*ikey*strd) < *reinterpret_cast<id_type*>(buff+*jkey*strd));},
            [](int* ikey,int* jkey){std::swap(*ikey,*jkey);});
            
            
            vec0_lst_sz=vec1_lst_sz=0;
            byte* ch_buff=buff+i_sz;
            char mrkr;
            for(int i=0;i<n;i++)
            {
                mrkr=*reinterpret_cast<char*>(ch_buff);
                if(mrkr==-1)
                {
                    ++vec0_lst_sz;
                }
                else if(mrkr==1)
                {
                    ++vec1_lst_sz;
                }
                else
                {
                    ++vec0_lst_sz;
                    ++vec1_lst_sz;
                }
                
                ch_buff+=strd;
            }
            
            delete [] vec0_lst;
            vec0_lst=vec0_lst_sz==0 ? NULL:new int[vec0_lst_sz];
            vec0_lst_sz=0;
            
            delete [] vec1_lst;
            vec1_lst=vec1_lst_sz==0 ? NULL:new int[vec1_lst_sz];
            vec1_lst_sz=0;
            
            int i0=0;
            int i1=0;
            int i=0;

            ch_buff=buff+i_sz;
            while(i0<n0 && i1<n1)
            {
                if(*reinterpret_cast<id_type*>(buff+i0*strd)<*reinterpret_cast<id_type*>(buff+key_1[i1]*strd))
                {
                    i=i0;
                    i0++;
                }
                else
                {
                    i=key_1[i1];
                    i1++;
                }
                mrkr=*reinterpret_cast<char*>(ch_buff+i*strd);
                if(mrkr==-1)
                {
                    vec0_lst[vec0_lst_sz++]=i;
                }
                else if(mrkr==1)
                {
                    vec1_lst[vec1_lst_sz++]=i;
                }
                else
                {
                    vec0_lst[vec0_lst_sz++]=i;
                    vec1_lst[vec1_lst_sz++]=i;
                }
                
            }
            
            while(i0<n0)
            {
                i=i0;
                i0++;
                
                mrkr=*reinterpret_cast<char*>(ch_buff+i*strd);
                if(mrkr==-1)
                {
                    vec0_lst[vec0_lst_sz++]=i;
                }
                else if(mrkr==1)
                {
                    vec1_lst[vec1_lst_sz++]=i;
                }
                else
                {
                    vec0_lst[vec0_lst_sz++]=i;
                    vec1_lst[vec1_lst_sz++]=i;
                }
            }
            
            while(i1<n1)
            {
                i=key_1[i1];
                i1++;
                
                mrkr=*reinterpret_cast<char*>(ch_buff+i*strd);
                if(mrkr==-1)
                {
                    vec0_lst[vec0_lst_sz++]=i;
                }
                else if(mrkr==1)
                {
                    vec1_lst[vec1_lst_sz++]=i;
                }
                else
                {
                    vec0_lst[vec0_lst_sz++]=i;
                    vec1_lst[vec1_lst_sz++]=i;
                }
            }
            
            delete [] key_1;
            
        }
        
        void mk_lsts(type0* b0,type0* x0,type0* b1,type0* x1,id_type* __id,int __natms_lcl)
        {
            init_lst(b0,x0,b1,x1,__id,__natms_lcl);
            xchng_lst<__dim__>();
            fin_lst();
            tot_snd_max_sz=Algebra::sum<2*__dim__+1>(snd_max_sz);
            
            delete [] updt_snd_buff[0];
            updt_snd_buff[0]=tot_snd_max_sz==0 ? NULL:new byte[updt_strd*tot_snd_max_sz];
            Algebra::Do<2*__dim__>::func([this](int i)
            {updt_snd_buff[1+i]=updt_snd_buff[i]+snd_max_sz[i]*updt_strd;});
            
            delete [] updt_rcv_buff;
            updt_rcv_buff=rcv_max_sz==0 ? NULL:new byte[updt_strd*rcv_max_sz];
            
        }
        
        void updt(byte* vec,const int strd)
        {
            if(strd>updt_strd)
            {
                updt_strd=strd;
                delete [] updt_snd_buff[0];
                updt_snd_buff[0]=tot_snd_max_sz==0 ? NULL:new byte[updt_strd*tot_snd_max_sz];
                Algebra::Do<2*__dim__>::func([this](int i)
                {updt_snd_buff[1+i]=updt_snd_buff[i]+snd_max_sz[i]*updt_strd;});
                
                delete [] updt_rcv_buff;
                updt_rcv_buff=rcv_max_sz==0 ? NULL:new byte[updt_strd*rcv_max_sz];
            }
            
            Algebra::zero<2*__dim__+1>(updt_snd_buff_sz);


            
            int iupdt_lst;
            byte* buff=vec;
            for(int i=0,j=0;j<rcv_szs[0];j++,i++)
            {
                if(same[i])
                {
                    iupdt_lst=lsts[0][j];
                    memcpy(buff,updt_snd_buff[iupdt_lst]+updt_snd_buff_sz[iupdt_lst]*strd,strd);
                    ++updt_snd_buff_sz[iupdt_lst];
                }
                else
                {
                    iupdt_lst=lsts[0][j];
                    memcpy(buff,updt_snd_buff[iupdt_lst]+updt_snd_buff_sz[iupdt_lst]*strd,strd);
                    ++updt_snd_buff_sz[iupdt_lst];
                    j++;
                    iupdt_lst=lsts[0][j];
                    memcpy(buff,updt_snd_buff[iupdt_lst]+updt_snd_buff_sz[iupdt_lst]*strd,strd);
                    ++updt_snd_buff_sz[iupdt_lst];
                }
                buff+=strd;
            }
            int ilst=1;
            __updt<__dim__>(strd,ilst);
            
        }
        
        
        
        
        template<int idim,int idir>
        void __updt__(const int strd,int& ilst)
        {
            int iupdt_lst;
            for(int icomm=0;icomm<ncomms[idim][idir];icomm++)
            {
                
                MPI_Sendrecv(updt_snd_buff[idim*2+idir],snd_szs[ilst]*strd,MPI_BYTE,neigh[idim][idir],0,
                             updt_rcv_buff,rcv_szs[ilst]*strd,MPI_BYTE,neigh[idim][1-idir],0,
                             world,MPI_STATUS_IGNORE);
                
                updt_snd_buff_sz[idim*2+idir]=0;
                for(int i=0;i<rcv_szs[ilst];i++)
                {
                    iupdt_lst=lsts[ilst][i];
                    memcpy(updt_rcv_buff+i*strd,updt_snd_buff[iupdt_lst]+updt_snd_buff_sz[iupdt_lst]*strd,strd);
                    ++updt_snd_buff_sz[iupdt_lst];
                }
                ++ilst;
            }
        }
        
        template<int idim>
        void __updt(const int strd,int& ilst)
        {
            __updt__<__dim__-idim,0>(strd,ilst);
            __updt__<__dim__-idim,1>(strd,ilst);
            __updt<idim-1>(strd,ilst);
        }
        
        void fin_updt(const int strd,type0* vec0,type0* vec1)
        {
            for(int i=0;i<vec0_lst_sz;i++)
                memcpy(updt_snd_buff[2*__dim__]+vec0_lst[i]*strd,vec0+i*strd,strd);
            
            for(int i=0;i<vec1_lst_sz;i++)
                memcpy(updt_snd_buff[2*__dim__]+vec1_lst[i]*strd,vec1+i*strd,strd);
                
        }
        
        
    protected:
    public:
        
        UpdateNEB(Atoms* atoms):
        world(atoms->comm.world)
        {
            
        }
        
        
        
    };
    
    template<>
    void UpdateNEB::xchng_lst<0>(){}
    template<>
    void UpdateNEB::__rsrv_buffs<0>(byte**,int*,int*,int*,const int*){}
    template<>
    int UpdateNEB::ret_lst<0>(const int*){return 0;}
    template<>
    void UpdateNEB::ret_dsp<0>(const type0*,const int*,const int*,int*){}
    template<>
    void UpdateNEB::__updt<0>(const int,int&){}
    
}
#endif
