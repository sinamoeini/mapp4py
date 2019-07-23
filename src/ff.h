/*--------------------------------------------
 Created by Sina on 07/16/13.
 Copyright (c) 2013 MIT. All rights reserved.
 --------------------------------------------*/
#ifndef __MAPP__ff__
#define __MAPP__ff__
#include "global.h"
#include "comm.h"
#include "atoms.h"
#include "exchange.h"
namespace MAPP_NS
{
    class ForceField
    {
    private:
    protected:
        Update* updt;

        template<class...VS>
        void update(VS*&... __vs)
        {
            updt->update_wo_x(__vs...);
        }
        const size_t nelems;
        virtual void __force_calc()=0;
        virtual void __energy_calc()=0;
        
        type0 __vec[__nvoigt__+1];
        type0 __vec_lcl[__nvoigt__+1];
        
        MPI_Comm& world;        
    public:
        ForceField(Atoms*);
        virtual ~ForceField();

        
        virtual void pre_init()=0;
        virtual void post_fin()=0;
        virtual void init()=0;
        virtual void fin()=0;
        type0* rsq_crd;
        type0** cut;
        type0** cut_sq;
        type0** cut_sk_sq;
        

        virtual type0 value()=0;
        virtual type0* derivative()=0;
        
    };
}
#endif

