/*--------------------------------------------
 Created by Sina on 07/16/13.
 Copyright (c) 2013 MIT. All rights reserved.
 --------------------------------------------*/
#ifndef __MAPP__ff__
#define __MAPP__ff__
#include "global.h"
#include "comm.h"
#include "atoms.h"

namespace MAPP_NS
{
    template<typename> class Vec;
    class ForceField
    {
    private:
    protected:
        const size_t nelems;
        type0 nrgy_strss_lcl[1+__dim__*(__dim__+1)/2];
        virtual void force_calc()=0;
        virtual void energy_calc()=0;
        
        MPI_Comm& world;        
    public:
        ForceField(Atoms*);
        virtual ~ForceField();

        
        virtual void setup()=0;
        virtual void init()=0;
        virtual void fin()=0;
        type0* rsq_crd;
        type0** cut;
        type0** cut_sq;
        type0** cut_sk_sq;

        virtual void reset();
        

        virtual type0 value_timer()=0;
        virtual void derivative_timer()=0;
        virtual void derivative_timer(type0(*&)[__dim__])=0;
        
        
        Vec<type0>* f;
        type0 nrgy_strss[1+__dim__*(__dim__+1)/2];
    };
}
#endif

