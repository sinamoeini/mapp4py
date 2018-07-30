#include "global.h"
#include <stddef.h>
#include <string>
#ifndef __MAPP__mapp_func__
#define __MAPP__mapp_func__
namespace MAPP_NS
{
    class EAMFunc
    {
    private:
    protected:
    public:
        elem_type nelems;
        std::string* names;
        type0** rc;

        void alloc_ptrs()
        {
            names=new std::string[nelems];
            rc=new type0*[nelems];
            *rc=new type0[nelems*nelems];

            for(int i=1;i<nelems;i++)
                rc[i]=rc[i-1]+nelems;
        }
        void dealloc_ptrs()
        {
            delete [] *rc;
            delete [] rc;
            delete [] names;
        }
        elem_type* remap(const std::string* __names,elem_type __nelems)
        {
            if(__nelems==0) return NULL;
            elem_type* map=new elem_type[__nelems];
            for(elem_type i=0;i<__nelems;i++)
            {
                bool found=false;
                for(elem_type j=0;j<nelems&&!found;j++)
                    if(strcmp(__names[i].c_str(),names[j].c_str())==0)
                    {
                        found=true;
                        map[i]=j;
                    }
                if(!found)
                    throw std::string("parameters for element ")+__names[i]+std::string("was not found in the forcefield");
            }
            
            return map;
        }
        
        EAMFunc(elem_type __nelems):nelems(__nelems){alloc_ptrs();};
        virtual ~EAMFunc(){dealloc_ptrs();};
        virtual type0 F(const elem_type&,const type0&)=0;
        virtual type0 rho(const elem_type&,const elem_type&,const type0&)=0;
        virtual type0 phi(const elem_type&,const elem_type&,const type0&)=0;
        virtual type0 fpair(const elem_type&,const elem_type&,const type0&,const type0&,const type0&)=0;
    };
}
typedef MAPP_NS::EAMFunc* create_eam_ff_t();

/*--------------------------------------------

 --------------------------------------------*/
namespace MAPP_NS
{
    class ADPFunc
    {
    private:
    protected:
    public:
        elem_type nelems;
        std::string* names;
        type0** rc;

        void alloc_ptrs()
        {
            names=new std::string[nelems];
            rc=new type0*[nelems];
            *rc=new type0[nelems*nelems];

            for(int i=1;i<nelems;i++)
                rc[i]=rc[i-1]+nelems;
        }
        void dealloc_ptrs()
        {
            delete [] *rc;
            delete [] rc;
            delete [] names;
        }
        elem_type* remap(const std::string* __names,elem_type __nelems)
        {
            if(__nelems==0) return NULL;
            elem_type* map=new elem_type[__nelems];
            for(elem_type i=0;i<__nelems;i++)
            {
                bool found=false;
                for(elem_type j=0;j<nelems&&!found;j++)
                    if(strcmp(__names[i].c_str(),names[j].c_str())==0)
                    {
                        found=true;
                        map[i]=j;
                    }
                if(!found)
                    throw std::string("parameters for element ")+__names[i]+std::string("was not found in the forcefield");
            }
            
            return map;
        }
        
        ADPFunc(elem_type __nelems):nelems(__nelems){alloc_ptrs();};
        virtual ~ADPFunc(){dealloc_ptrs();};
        virtual type0 F(const elem_type&,const type0&)=0;
        virtual type0 rho(const elem_type&,const elem_type&,const type0&)=0;
        virtual type0 phi(const elem_type&,const elem_type&,const type0&)=0;
        virtual type0 fpair(const elem_type&,const elem_type&,const type0&,const type0&,const type0&)=0;
    };
}

#endif
