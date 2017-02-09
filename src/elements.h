/*--------------------------------------------
 Created by Sina on 06/27/14.
 Copyright (c) 2014 MIT. All rights reserved.
 --------------------------------------------*/
#ifndef __MAPP__elements__
#define __MAPP__elements__
#include "global.h"
#include "api.h"
#include <stdio.h>
struct _object;
typedef _object PyObject;
namespace MAPP_NS
{
    class Elements
    {
    private:
    protected:
    public:
        var<size_t> __nelems__;
        size_t nelems;
        elem_type __nelems;
        char** names;
        type0* masses;
        
        Elements();
        Elements(const Elements&);
        Elements(Elements&&);
        Elements& operator=(const Elements&);
        Elements& operator=(Elements&&);
        bool operator==(const Elements&);
        ~Elements();
        elem_type add_type(const type0,const char*);
        elem_type find(const char*);
        void assign_color_rad(const char*,type0(&)[4]);
        PyObject* get_dict();
    };
}

#endif 
