#include "min_cg3_dmd.h"
#include <stdlib.h>
/*
this is the code that was used to generate all this crap

#include <stdio.h>
#include <string>
const char* newvecargs []=
{"__dim__","c_dim","c_dim"};
const char* true_str="true";
const char* false_str="false";
const char* __dim___str="__dim__";
const char* c_dim_str="c_dim";

std::string VECTENS1_new_args(bool BC,bool X,bool A,bool C)
{
    std::string ans="VECTENS1(atoms,";
    ans+=BC?"true":"false";
    if(X) ans+=",__dim__";
    if(A) ans+=",c_dim";
    if(C) ans+=",c_dim";
    ans+=")";
    return ans;
}

std::string VECTENS1_f_args(bool BC,bool X,bool A,bool C)
{
    std::string ans="VECTENS1(atoms,";
    ans+=BC?"true":"false";
    ans+=",ff->F_H";
    if(X) ans+=",ff->f";
    if(A) ans+=",ff->f_alpha";
    if(C) ans+=",ff->mu";
    ans+=")";
    return ans;
}

std::string VECTENS0_new_args(bool BC,bool X,bool A,bool C)
{
    std::string ans="VECTENS0(atoms,";
    ans+=BC?"true":"false";
    if(X || BC) ans+=",__dim__";
    if(A) ans+=",c_dim";
    if(C) ans+=",c_dim";
    ans+=")";
    return ans;
}


std::string VECTENS0_x_args(bool BC,bool X,bool A,bool C)
{
    std::string ans="VECTENS0(atoms,";
    ans+=BC?"true":"false";
    ans+=",atoms->H";
    if(X || BC) ans+=",atoms->x";
    if(A) ans+=",atoms->alpha";
    if(C) ans+=",atoms->c";
    ans+=")";
    return ans;
}
std::string VECTENS0_x_d_args(bool BC,bool X,bool A,bool C)
{
    std::string ans="VECTENS0(atoms,";
    ans+=BC?"true":"false";
    int n=0;
    if(BC==true && X==true)
    {
        ans+=",__dim__";
        n++;
    }
    else if(BC==false && X==true)
    {
        ans+=",h.vecs["+std::to_string(n++)+"]";
    }
    else if(BC==true && X==false)
    {
        ans+=",__dim__";
    }
 
    if(A) ans+=",h.vecs["+std::to_string(n++)+"]";
    if(C) ans+=",h.vecs["+std::to_string(n++)+"]";
 
    ans+=")";
    return ans;
}
void print_fsep()
{
}
void print_fsep_large()
{

}
std::string header(const char* c0,bool BC,bool X,bool A,bool C,const char* c1)
{
    std::string ans="template<>\n";
    ans+=c0;
    ans+=" MinDMDHandler<";
    ans+=BC?true_str:false_str;
    ans+=",";
    ans+=X?true_str:false_str;
    ans+=",";
    ans+=A?true_str:false_str;
    ans+=",";
    ans+=C?true_str:false_str;
    ans+=">::";
    ans+=c1;
    return ans;
}

std::string header_dec(const char* c0,bool BC,bool X,bool A,bool C,const char* c1)
{
    std::string ans="template<>";
    ans+=c0;
    ans+=" MinDMDHandler<";
    ans+=BC?true_str:false_str;
    ans+=",";
    ans+=X?true_str:false_str;
    ans+=",";
    ans+=A?true_str:false_str;
    ans+=",";
    ans+=C?true_str:false_str;
    ans+=">::";
    ans+=c1;
    ans+=";";
    return ans;
}
void print_init(bool BC,bool X,bool A,bool C)
{
    printf("%s\n",header("void",BC,X,A,C,"init()").c_str());

    printf("{\n");
 
    printf("    f.~VECTENS1();\n");
    printf("    new (&f) %s;\n",VECTENS1_f_args(BC,X,A,C).c_str());
    printf("    f0.~VECTENS1();\n");
    printf("    new (&f0) %s;\n",VECTENS1_new_args(BC,X,A,C).c_str());
    printf("    h.~VECTENS1();\n");
    printf("    new (&h) %s;\n",VECTENS1_new_args(BC,X,A,C).c_str());
    printf("\n");
 
    printf("    x.~VECTENS0();\n");
    printf("    new (&x) %s;\n",VECTENS0_x_args(BC,X,A,C).c_str());
    printf("    x0.~VECTENS0();\n");
    printf("    new (&x0) %s;\n",VECTENS0_new_args(BC,X,A,C).c_str());
    printf("    x_d.~VECTENS0();\n");
    printf("    new (&x_d) %s;\n",VECTENS0_x_d_args(BC,X,A,C).c_str());
 
    printf("    dynamic=new NewDynamicDMD<%s,%s,%s>(atoms,ff,{},\n",
           BC?true_str:false_str,
           BC||X?true_str:false_str,
           A?true_str:false_str);

    int n0=(X?1:0)+(A?1:0)+(C?1:0);
    int n1=(BC||X?1:0)+(A?1:0)+(C?1:0);
    int n2=(BC?1:0);
    printf("    {\n        ");
    if(n0+n1+n2==0)
        printf("atoms->x_dof,atoms->alpha_dof,atoms->c_dof");
    else
        printf("atoms->x_dof,atoms->alpha_dof,atoms->c_dof,");
    if(n0)
    {
        printf("\n        ");
        for(int i=0;i<n0;i++)
            printf("h.vecs[%d],",i);
        printf("\n        ");
        for(int i=0;i<n0;i++)
            printf("f0.vecs[%d],",i);
    }
    if(n1)
    {
        printf("\n        ");
        for(int i=0;i<n1-1;i++)
            printf("x0.vecs[%d],",i);
        if(n2)
            printf("x0.vecs[%d],",n1-1);
        else
            printf("x0.vecs[%d]",n1-1);
    }
    if(n2)
        printf("\n        x_d.vecs[0]");
    printf("\n    },{});");
 
 
 
    printf("\n}\n");
}
void print_add_extra_vec(bool BC,bool X,bool A,bool C)
{
    printf("%s\n",header("void",BC,X,A,C,"add_extra_vec(VECTENS1& __v)").c_str());
    printf("{\n");
    int n0=(X?1:0)+(A?1:0)+(C?1:0);
    printf("    __v.~VECTENS1();\n");
    printf("    new (&__v) %s;\n",VECTENS1_new_args(BC,X,A,C).c_str());
    for(int i=0;i<n0;i++)
        printf("    dynamic->add_xchng(__v.vecs[%d]);\n",i);
    printf("}\n");
}



void print_get_max_alpha(bool BC,bool X,bool A,bool C)
{
    printf("%s\n",header("type0",BC,X,A,C,"get_max_alpha()").c_str());
    printf("{\n");
    printf("    type0 max_a=0.0,max_a_lcl=std::numeric_limits<type0>::infinity();\n");
    int n=0;
    if(X || BC) printf("    max_alpha_lcl_x(max_a_lcl,x_d.vecs[%d]->begin());\n",n++);
    if(A) printf("    max_alpha_lcl_alpha(max_a_lcl,x_d.vecs[%d]->begin());\n",n++);
    if(C) printf("    max_alpha_lcl_c(max_a_lcl,x_d.vecs[%d]->begin());\n",n++);
    printf("    MPI_Allreduce(&max_a_lcl,&max_a,1,Vec<type0>::MPI_T,MPI_MIN,atoms->world);\n");
    printf("    return max_a;\n");
    printf("}\n");
}
void print_force_calc(bool BC,bool X,bool A,bool C)
{
    printf("%s\n",header("void",BC,X,A,C,"force_calc()").c_str());
    printf("{\n");
    printf("    ff->derivative();\n");
    if(BC) printf("    post_force_calc_box();\n");
    if(C) printf("    post_force_calc_c();\n");
    printf("}\n");
}
void print_prep(bool BC,bool X,bool A,bool C)
{
    printf("%s\n",header("void",BC,X,A,C,"prep()").c_str());
    printf("{\n");
    if(BC && X) printf("    prep_x_d();\n");
    if(BC && !X) printf("    prep_x_d_affine();\n");
    printf("}\n");
}

void print_F(bool BC,bool X,bool A,bool C)
{
    printf("%s\n",header("type0",BC,X,A,C,"F(type0 __alpha)").c_str());
    printf("{\n");
    printf("    x=x0+__alpha*x_d;\n");
    if(C) printf("    aprop_c();\n");
    if(C) printf("    dynamic->update(atoms->c);\n");
    else printf("    dynamic->update();\n");
    printf("    return ff->value();\n");
    printf("}\n");
}
void print_F_reset(bool BC,bool X,bool A,bool C)
{
    printf("%s\n",header("void",BC,X,A,C,"F_reset()").c_str());
    printf("{\n");
    printf("    x=x0;\n");
    if(C) printf("    dynamic->update(atoms->c);\n");
    else printf("    dynamic->update();\n");
    printf("}\n");
}

void print_all_funcs(bool BC,bool X,bool A,bool C)
{
    print_fsep_large();
    print_init(BC,X,A,C);
    print_fsep();
    print_add_extra_vec(BC,X,A,C);
 
//    print_fsep();
//    print_prep(BC,X,A,C);
//    print_fsep();
//    print_force_calc(BC,X,A,C);
//    print_fsep();
//    print_F(BC,X,A,C);
//    print_fsep();
//    print_F_reset(BC,X,A,C);
 
}

void print_all_decs(bool BC,bool X,bool A,bool C)
{
    printf("        %s\n",header_dec("void",BC,X,A,C,"init()").c_str());
    printf("        %s\n",header_dec("void",BC,X,A,C,"add_extra_vec(VECTENS1&)").c_str());
//    printf("        %s\n",header_dec("type0",BC,X,A,C,"get_max_alpha()").c_str());
//    printf("        %s\n",header_dec("void",BC,X,A,C,"prep()").c_str());
//    printf("        %s\n",header_dec("void",BC,X,A,C,"force_calc()").c_str());
//    printf("        %s\n",header_dec("type0",BC,X,A,C,"F(type0)").c_str());
//    printf("        %s\n",header_dec("void",BC,X,A,C,"F_reset()").c_str());
    printf("\n\n");
}

void print_MinDMDHandler_t(bool BC,bool X,bool A,bool C)
{
 
    printf("    template<>\n");
    printf("    class MinDMDHandler_t<%s,%s,%s,%s>\n",
           BC?true_str:false_str,
           X?true_str:false_str,
           A?true_str:false_str,
           C?true_str:false_str
           );
    printf("    {\n");
    printf("    public:\n");
 
    int n0=(BC||X?1:0)+(A?1:0)+(C?1:0);
    printf("        typedef VecTens<type0,%d> VECTENS0;\n",n0);
    int n1=(X?1:0)+(A?1:0)+(C?1:0);
    printf("        typedef VecTens<type0,%d> VECTENS1;\n",n1);
 
    printf("    };\n");
}


void print_options(const char* ls,bool BC,bool X,bool A,bool C)
{
    printf("    if(dynamic_cast<%s*>(ls)!=NULL && B_DOF==%s &&  X_DOF==%s &&  ALPHA_DOF==%s &&  C_DOF==%s)\n",ls,
           BC?true_str:false_str,
           X?true_str:false_str,
           A?true_str:false_str,
           C?true_str:false_str);
    printf("        return run<%s,%s,%s,%s>(dynamic_cast<%s*>(ls),nsteps);\n",
          BC?true_str:false_str,
          X?true_str:false_str,
          A?true_str:false_str,
          C?true_str:false_str,
          ls);
}

int main()
{
    for(int i0=0;i0<2;i0++)
        for(int i1=0;i1<2;i1++)
            for(int i2=0;i2<2;i2++)
                for(int i3=0;i3<2;i3++)
                {
                    if(i0+i1+i2+i3==0) continue;
                    //print_MinDMDHandler_t(i0,i1,i2,i3);
                    //print_all_decs(i0,i1,i2,i3);
                    print_all_funcs(i0,i1,i2,i3);
                    //print_options("LineSearchBrent",i0,i1,i2,i3);
                    //print_options("LineSearchGoldenSection",i0,i1,i2,i3);
                    //print_options("LineSearchBackTrack",i0,i1,i2,i3);
                }
    return 0;
}

*/



using namespace MAPP_NS;
/*----------------------------------------------------------------------------------------
 
 ----------------------------------------------------------------------------------------*/
template<>
void MinDMDHandler<false,false,false,true>::init()
{
    f.~VECTENS1();
    new (&f) VECTENS1(atoms,false,ff->F_H,ff->mu);
    f0.~VECTENS1();
    new (&f0) VECTENS1(atoms,false,c_dim);
    h.~VECTENS1();
    new (&h) VECTENS1(atoms,false,c_dim);

    x.~VECTENS0();
    new (&x) VECTENS0(atoms,false,atoms->H,atoms->c);
    x0.~VECTENS0();
    new (&x0) VECTENS0(atoms,false,c_dim);
    x_d.~VECTENS0();
    new (&x_d) VECTENS0(atoms,false,h.vecs[0]);
    dynamic=new NewDynamicDMD<false,false,false>(atoms,ff,{},
    {
        atoms->x_dof,atoms->alpha_dof,atoms->c_dof,
        h.vecs[0],
        f0.vecs[0],
        x0.vecs[0]
    },{});
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<>
void MinDMDHandler<false,false,false,true>::add_extra_vec(VECTENS1& __v)
{
    __v.~VECTENS1();
    new (&__v) VECTENS1(atoms,false,c_dim);
    dynamic->add_xchng(__v.vecs[0]);
}
/*----------------------------------------------------------------------------------------
 
 ----------------------------------------------------------------------------------------*/
template<>
void MinDMDHandler<false,false,true,false>::init()
{
    f.~VECTENS1();
    new (&f) VECTENS1(atoms,false,ff->F_H,ff->f_alpha);
    f0.~VECTENS1();
    new (&f0) VECTENS1(atoms,false,c_dim);
    h.~VECTENS1();
    new (&h) VECTENS1(atoms,false,c_dim);

    x.~VECTENS0();
    new (&x) VECTENS0(atoms,false,atoms->H,atoms->alpha);
    x0.~VECTENS0();
    new (&x0) VECTENS0(atoms,false,c_dim);
    x_d.~VECTENS0();
    new (&x_d) VECTENS0(atoms,false,h.vecs[0]);
    dynamic=new NewDynamicDMD<false,false,true>(atoms,ff,{},
    {
        atoms->x_dof,atoms->alpha_dof,atoms->c_dof,
        h.vecs[0],
        f0.vecs[0],
        x0.vecs[0]
    },{});
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<>
void MinDMDHandler<false,false,true,false>::add_extra_vec(VECTENS1& __v)
{
    __v.~VECTENS1();
    new (&__v) VECTENS1(atoms,false,c_dim);
    dynamic->add_xchng(__v.vecs[0]);
}
/*----------------------------------------------------------------------------------------
 
 ----------------------------------------------------------------------------------------*/
template<>
void MinDMDHandler<false,false,true,true>::init()
{
    f.~VECTENS1();
    new (&f) VECTENS1(atoms,false,ff->F_H,ff->f_alpha,ff->mu);
    f0.~VECTENS1();
    new (&f0) VECTENS1(atoms,false,c_dim,c_dim);
    h.~VECTENS1();
    new (&h) VECTENS1(atoms,false,c_dim,c_dim);

    x.~VECTENS0();
    new (&x) VECTENS0(atoms,false,atoms->H,atoms->alpha,atoms->c);
    x0.~VECTENS0();
    new (&x0) VECTENS0(atoms,false,c_dim,c_dim);
    x_d.~VECTENS0();
    new (&x_d) VECTENS0(atoms,false,h.vecs[0],h.vecs[1]);
    dynamic=new NewDynamicDMD<false,false,true>(atoms,ff,{},
    {
        atoms->x_dof,atoms->alpha_dof,atoms->c_dof,
        h.vecs[0],h.vecs[1],
        f0.vecs[0],f0.vecs[1],
        x0.vecs[0],x0.vecs[1]
    },{});
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<>
void MinDMDHandler<false,false,true,true>::add_extra_vec(VECTENS1& __v)
{
    __v.~VECTENS1();
    new (&__v) VECTENS1(atoms,false,c_dim,c_dim);
    dynamic->add_xchng(__v.vecs[0]);
    dynamic->add_xchng(__v.vecs[1]);
}
/*----------------------------------------------------------------------------------------
 
 ----------------------------------------------------------------------------------------*/
template<>
void MinDMDHandler<false,true,false,false>::init()
{
    f.~VECTENS1();
    new (&f) VECTENS1(atoms,false,ff->F_H,ff->f);
    f0.~VECTENS1();
    new (&f0) VECTENS1(atoms,false,__dim__);
    h.~VECTENS1();
    new (&h) VECTENS1(atoms,false,__dim__);

    x.~VECTENS0();
    new (&x) VECTENS0(atoms,false,atoms->H,atoms->x);
    x0.~VECTENS0();
    new (&x0) VECTENS0(atoms,false,__dim__);
    x_d.~VECTENS0();
    new (&x_d) VECTENS0(atoms,false,h.vecs[0]);
    dynamic=new NewDynamicDMD<false,true,false>(atoms,ff,{},
    {
        atoms->x_dof,atoms->alpha_dof,atoms->c_dof,
        h.vecs[0],
        f0.vecs[0],
        x0.vecs[0]
    },{});
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<>
void MinDMDHandler<false,true,false,false>::add_extra_vec(VECTENS1& __v)
{
    __v.~VECTENS1();
    new (&__v) VECTENS1(atoms,false,__dim__);
    dynamic->add_xchng(__v.vecs[0]);
}
/*----------------------------------------------------------------------------------------
 
 ----------------------------------------------------------------------------------------*/
template<>
void MinDMDHandler<false,true,false,true>::init()
{
    f.~VECTENS1();
    new (&f) VECTENS1(atoms,false,ff->F_H,ff->f,ff->mu);
    f0.~VECTENS1();
    new (&f0) VECTENS1(atoms,false,__dim__,c_dim);
    h.~VECTENS1();
    new (&h) VECTENS1(atoms,false,__dim__,c_dim);

    x.~VECTENS0();
    new (&x) VECTENS0(atoms,false,atoms->H,atoms->x,atoms->c);
    x0.~VECTENS0();
    new (&x0) VECTENS0(atoms,false,__dim__,c_dim);
    x_d.~VECTENS0();
    new (&x_d) VECTENS0(atoms,false,h.vecs[0],h.vecs[1]);
    dynamic=new NewDynamicDMD<false,true,false>(atoms,ff,{},
    {
        atoms->x_dof,atoms->alpha_dof,atoms->c_dof,
        h.vecs[0],h.vecs[1],
        f0.vecs[0],f0.vecs[1],
        x0.vecs[0],x0.vecs[1]
    },{});
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<>
void MinDMDHandler<false,true,false,true>::add_extra_vec(VECTENS1& __v)
{
    __v.~VECTENS1();
    new (&__v) VECTENS1(atoms,false,__dim__,c_dim);
    dynamic->add_xchng(__v.vecs[0]);
    dynamic->add_xchng(__v.vecs[1]);
}
/*----------------------------------------------------------------------------------------
 
 ----------------------------------------------------------------------------------------*/
template<>
void MinDMDHandler<false,true,true,false>::init()
{
    f.~VECTENS1();
    new (&f) VECTENS1(atoms,false,ff->F_H,ff->f,ff->f_alpha);
    f0.~VECTENS1();
    new (&f0) VECTENS1(atoms,false,__dim__,c_dim);
    h.~VECTENS1();
    new (&h) VECTENS1(atoms,false,__dim__,c_dim);

    x.~VECTENS0();
    new (&x) VECTENS0(atoms,false,atoms->H,atoms->x,atoms->alpha);
    x0.~VECTENS0();
    new (&x0) VECTENS0(atoms,false,__dim__,c_dim);
    x_d.~VECTENS0();
    new (&x_d) VECTENS0(atoms,false,h.vecs[0],h.vecs[1]);
    dynamic=new NewDynamicDMD<false,true,true>(atoms,ff,{},
    {
        atoms->x_dof,atoms->alpha_dof,atoms->c_dof,
        h.vecs[0],h.vecs[1],
        f0.vecs[0],f0.vecs[1],
        x0.vecs[0],x0.vecs[1]
    },{});
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<>
void MinDMDHandler<false,true,true,false>::add_extra_vec(VECTENS1& __v)
{
    __v.~VECTENS1();
    new (&__v) VECTENS1(atoms,false,__dim__,c_dim);
    dynamic->add_xchng(__v.vecs[0]);
    dynamic->add_xchng(__v.vecs[1]);
}
/*----------------------------------------------------------------------------------------
 
 ----------------------------------------------------------------------------------------*/
template<>
void MinDMDHandler<false,true,true,true>::init()
{
    f.~VECTENS1();
    new (&f) VECTENS1(atoms,false,ff->F_H,ff->f,ff->f_alpha,ff->mu);
    f0.~VECTENS1();
    new (&f0) VECTENS1(atoms,false,__dim__,c_dim,c_dim);
    h.~VECTENS1();
    new (&h) VECTENS1(atoms,false,__dim__,c_dim,c_dim);

    x.~VECTENS0();
    new (&x) VECTENS0(atoms,false,atoms->H,atoms->x,atoms->alpha,atoms->c);
    x0.~VECTENS0();
    new (&x0) VECTENS0(atoms,false,__dim__,c_dim,c_dim);
    x_d.~VECTENS0();
    new (&x_d) VECTENS0(atoms,false,h.vecs[0],h.vecs[1],h.vecs[2]);
    dynamic=new NewDynamicDMD<false,true,true>(atoms,ff,{},
    {
        atoms->x_dof,atoms->alpha_dof,atoms->c_dof,
        h.vecs[0],h.vecs[1],h.vecs[2],
        f0.vecs[0],f0.vecs[1],f0.vecs[2],
        x0.vecs[0],x0.vecs[1],x0.vecs[2]
    },{});
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<>
void MinDMDHandler<false,true,true,true>::add_extra_vec(VECTENS1& __v)
{
    __v.~VECTENS1();
    new (&__v) VECTENS1(atoms,false,__dim__,c_dim,c_dim);
    dynamic->add_xchng(__v.vecs[0]);
    dynamic->add_xchng(__v.vecs[1]);
    dynamic->add_xchng(__v.vecs[2]);
}
/*----------------------------------------------------------------------------------------
 
 ----------------------------------------------------------------------------------------*/
template<>
void MinDMDHandler<true,false,false,false>::init()
{
    f.~VECTENS1();
    new (&f) VECTENS1(atoms,true,ff->F_H);
    f0.~VECTENS1();
    new (&f0) VECTENS1(atoms,true);
    h.~VECTENS1();
    new (&h) VECTENS1(atoms,true);

    x.~VECTENS0();
    new (&x) VECTENS0(atoms,true,atoms->H,atoms->x);
    x0.~VECTENS0();
    new (&x0) VECTENS0(atoms,true,__dim__);
    x_d.~VECTENS0();
    new (&x_d) VECTENS0(atoms,true,__dim__);
    dynamic=new NewDynamicDMD<true,true,false>(atoms,ff,{},
    {
        atoms->x_dof,atoms->alpha_dof,atoms->c_dof,
        x0.vecs[0],
        x_d.vecs[0]
    },{});
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<>
void MinDMDHandler<true,false,false,false>::add_extra_vec(VECTENS1& __v)
{
    __v.~VECTENS1();
    new (&__v) VECTENS1(atoms,true);
}
/*----------------------------------------------------------------------------------------
 
 ----------------------------------------------------------------------------------------*/
template<>
void MinDMDHandler<true,false,false,true>::init()
{
    f.~VECTENS1();
    new (&f) VECTENS1(atoms,true,ff->F_H,ff->mu);
    f0.~VECTENS1();
    new (&f0) VECTENS1(atoms,true,c_dim);
    h.~VECTENS1();
    new (&h) VECTENS1(atoms,true,c_dim);

    x.~VECTENS0();
    new (&x) VECTENS0(atoms,true,atoms->H,atoms->x,atoms->c);
    x0.~VECTENS0();
    new (&x0) VECTENS0(atoms,true,__dim__,c_dim);
    x_d.~VECTENS0();
    new (&x_d) VECTENS0(atoms,true,__dim__,h.vecs[0]);
    dynamic=new NewDynamicDMD<true,true,false>(atoms,ff,{},
    {
        atoms->x_dof,atoms->alpha_dof,atoms->c_dof,
        h.vecs[0],
        f0.vecs[0],
        x0.vecs[0],x0.vecs[1],
        x_d.vecs[0]
    },{});
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<>
void MinDMDHandler<true,false,false,true>::add_extra_vec(VECTENS1& __v)
{
    __v.~VECTENS1();
    new (&__v) VECTENS1(atoms,true,c_dim);
    dynamic->add_xchng(__v.vecs[0]);
}
/*----------------------------------------------------------------------------------------
 
 ----------------------------------------------------------------------------------------*/
template<>
void MinDMDHandler<true,false,true,false>::init()
{
    f.~VECTENS1();
    new (&f) VECTENS1(atoms,true,ff->F_H,ff->f_alpha);
    f0.~VECTENS1();
    new (&f0) VECTENS1(atoms,true,c_dim);
    h.~VECTENS1();
    new (&h) VECTENS1(atoms,true,c_dim);

    x.~VECTENS0();
    new (&x) VECTENS0(atoms,true,atoms->H,atoms->x,atoms->alpha);
    x0.~VECTENS0();
    new (&x0) VECTENS0(atoms,true,__dim__,c_dim);
    x_d.~VECTENS0();
    new (&x_d) VECTENS0(atoms,true,__dim__,h.vecs[0]);
    dynamic=new NewDynamicDMD<true,true,true>(atoms,ff,{},
    {
        atoms->x_dof,atoms->alpha_dof,atoms->c_dof,
        h.vecs[0],
        f0.vecs[0],
        x0.vecs[0],x0.vecs[1],
        x_d.vecs[0]
    },{});
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<>
void MinDMDHandler<true,false,true,false>::add_extra_vec(VECTENS1& __v)
{
    __v.~VECTENS1();
    new (&__v) VECTENS1(atoms,true,c_dim);
    dynamic->add_xchng(__v.vecs[0]);
}
/*----------------------------------------------------------------------------------------
 
 ----------------------------------------------------------------------------------------*/
template<>
void MinDMDHandler<true,false,true,true>::init()
{
    f.~VECTENS1();
    new (&f) VECTENS1(atoms,true,ff->F_H,ff->f_alpha,ff->mu);
    f0.~VECTENS1();
    new (&f0) VECTENS1(atoms,true,c_dim,c_dim);
    h.~VECTENS1();
    new (&h) VECTENS1(atoms,true,c_dim,c_dim);

    x.~VECTENS0();
    new (&x) VECTENS0(atoms,true,atoms->H,atoms->x,atoms->alpha,atoms->c);
    x0.~VECTENS0();
    new (&x0) VECTENS0(atoms,true,__dim__,c_dim,c_dim);
    x_d.~VECTENS0();
    new (&x_d) VECTENS0(atoms,true,__dim__,h.vecs[0],h.vecs[1]);
    dynamic=new NewDynamicDMD<true,true,true>(atoms,ff,{},
    {
        atoms->x_dof,atoms->alpha_dof,atoms->c_dof,
        h.vecs[0],h.vecs[1],
        f0.vecs[0],f0.vecs[1],
        x0.vecs[0],x0.vecs[1],x0.vecs[2],
        x_d.vecs[0]
    },{});
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<>
void MinDMDHandler<true,false,true,true>::add_extra_vec(VECTENS1& __v)
{
    __v.~VECTENS1();
    new (&__v) VECTENS1(atoms,true,c_dim,c_dim);
    dynamic->add_xchng(__v.vecs[0]);
    dynamic->add_xchng(__v.vecs[1]);
}
/*----------------------------------------------------------------------------------------
 
 ----------------------------------------------------------------------------------------*/
template<>
void MinDMDHandler<true,true,false,false>::init()
{
    f.~VECTENS1();
    new (&f) VECTENS1(atoms,true,ff->F_H,ff->f);
    f0.~VECTENS1();
    new (&f0) VECTENS1(atoms,true,__dim__);
    h.~VECTENS1();
    new (&h) VECTENS1(atoms,true,__dim__);

    x.~VECTENS0();
    new (&x) VECTENS0(atoms,true,atoms->H,atoms->x);
    x0.~VECTENS0();
    new (&x0) VECTENS0(atoms,true,__dim__);
    x_d.~VECTENS0();
    new (&x_d) VECTENS0(atoms,true,__dim__);
    dynamic=new NewDynamicDMD<true,true,false>(atoms,ff,{},
    {
        atoms->x_dof,atoms->alpha_dof,atoms->c_dof,
        h.vecs[0],
        f0.vecs[0],
        x0.vecs[0],
        x_d.vecs[0]
    },{});
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<>
void MinDMDHandler<true,true,false,false>::add_extra_vec(VECTENS1& __v)
{
    __v.~VECTENS1();
    new (&__v) VECTENS1(atoms,true,__dim__);
    dynamic->add_xchng(__v.vecs[0]);
}
/*----------------------------------------------------------------------------------------
 
 ----------------------------------------------------------------------------------------*/
template<>
void MinDMDHandler<true,true,false,true>::init()
{
    f.~VECTENS1();
    new (&f) VECTENS1(atoms,true,ff->F_H,ff->f,ff->mu);
    f0.~VECTENS1();
    new (&f0) VECTENS1(atoms,true,__dim__,c_dim);
    h.~VECTENS1();
    new (&h) VECTENS1(atoms,true,__dim__,c_dim);

    x.~VECTENS0();
    new (&x) VECTENS0(atoms,true,atoms->H,atoms->x,atoms->c);
    x0.~VECTENS0();
    new (&x0) VECTENS0(atoms,true,__dim__,c_dim);
    x_d.~VECTENS0();
    new (&x_d) VECTENS0(atoms,true,__dim__,h.vecs[1]);
    dynamic=new NewDynamicDMD<true,true,false>(atoms,ff,{},
    {
        atoms->x_dof,atoms->alpha_dof,atoms->c_dof,
        h.vecs[0],h.vecs[1],
        f0.vecs[0],f0.vecs[1],
        x0.vecs[0],x0.vecs[1],
        x_d.vecs[0]
    },{});
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<>
void MinDMDHandler<true,true,false,true>::add_extra_vec(VECTENS1& __v)
{
    __v.~VECTENS1();
    new (&__v) VECTENS1(atoms,true,__dim__,c_dim);
    dynamic->add_xchng(__v.vecs[0]);
    dynamic->add_xchng(__v.vecs[1]);
}
/*----------------------------------------------------------------------------------------
 
 ----------------------------------------------------------------------------------------*/
template<>
void MinDMDHandler<true,true,true,false>::init()
{
    f.~VECTENS1();
    new (&f) VECTENS1(atoms,true,ff->F_H,ff->f,ff->f_alpha);
    f0.~VECTENS1();
    new (&f0) VECTENS1(atoms,true,__dim__,c_dim);
    h.~VECTENS1();
    new (&h) VECTENS1(atoms,true,__dim__,c_dim);

    x.~VECTENS0();
    new (&x) VECTENS0(atoms,true,atoms->H,atoms->x,atoms->alpha);
    x0.~VECTENS0();
    new (&x0) VECTENS0(atoms,true,__dim__,c_dim);
    x_d.~VECTENS0();
    new (&x_d) VECTENS0(atoms,true,__dim__,h.vecs[1]);
    dynamic=new NewDynamicDMD<true,true,true>(atoms,ff,{},
    {
        atoms->x_dof,atoms->alpha_dof,atoms->c_dof,
        h.vecs[0],h.vecs[1],
        f0.vecs[0],f0.vecs[1],
        x0.vecs[0],x0.vecs[1],
        x_d.vecs[0]
    },{});
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<>
void MinDMDHandler<true,true,true,false>::add_extra_vec(VECTENS1& __v)
{
    __v.~VECTENS1();
    new (&__v) VECTENS1(atoms,true,__dim__,c_dim);
    dynamic->add_xchng(__v.vecs[0]);
    dynamic->add_xchng(__v.vecs[1]);
}
/*----------------------------------------------------------------------------------------
 
 ----------------------------------------------------------------------------------------*/
template<>
void MinDMDHandler<true,true,true,true>::init()
{
    f.~VECTENS1();
    new (&f) VECTENS1(atoms,true,ff->F_H,ff->f,ff->f_alpha,ff->mu);
    f0.~VECTENS1();
    new (&f0) VECTENS1(atoms,true,__dim__,c_dim,c_dim);
    h.~VECTENS1();
    new (&h) VECTENS1(atoms,true,__dim__,c_dim,c_dim);

    x.~VECTENS0();
    new (&x) VECTENS0(atoms,true,atoms->H,atoms->x,atoms->alpha,atoms->c);
    x0.~VECTENS0();
    new (&x0) VECTENS0(atoms,true,__dim__,c_dim,c_dim);
    x_d.~VECTENS0();
    new (&x_d) VECTENS0(atoms,true,__dim__,h.vecs[1],h.vecs[2]);
    dynamic=new NewDynamicDMD<true,true,true>(atoms,ff,{},
    {
        atoms->x_dof,atoms->alpha_dof,atoms->c_dof,
        h.vecs[0],h.vecs[1],h.vecs[2],
        f0.vecs[0],f0.vecs[1],f0.vecs[2],
        x0.vecs[0],x0.vecs[1],x0.vecs[2],
        x_d.vecs[0]
    },{});
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<>
void MinDMDHandler<true,true,true,true>::add_extra_vec(VECTENS1& __v)
{
    __v.~VECTENS1();
    new (&__v) VECTENS1(atoms,true,__dim__,c_dim,c_dim);
    dynamic->add_xchng(__v.vecs[0]);
    dynamic->add_xchng(__v.vecs[1]);
    dynamic->add_xchng(__v.vecs[2]);
}










/*--------------------------------------------
 constructor
 --------------------------------------------*/
MinCG3DMD::MinCG3DMD(type0 __e_tol,
bool(&__H_dof)[__dim__][__dim__],bool __affine,type0 __max_dx,type0 __max_dalpha,LineSearch* __ls):
Min(__e_tol,__H_dof,__affine,__max_dx,__ls),
atoms(NULL),
ff(NULL),
max_dalpha(__max_dalpha),
xprt(NULL),
X_DOF(true),
ALPHA_DOF(true),
C_DOF(false)
{
}
/*--------------------------------------------
 destructor
 --------------------------------------------*/
MinCG3DMD::~MinCG3DMD()
{
    atoms=NULL;
    ff=NULL;
}
/*--------------------------------------------
 pre run check it throw excepctions
 --------------------------------------------*/
void MinCG3DMD::pre_run_chk(AtomsDMD* __atoms,ForceFieldDMD* __ff)
{
    //check if configuration is loaded
    if(!__atoms)
        throw std::string("cannot start minimization without initial conditions");
    
    //check if force field is loaded
    if(!__ff)
        throw std::string("cannot start minimization without governing equations (force field)");
    
    if(std::isnan(__atoms->kB))
        throw std::string("boltzmann constant should be set prior to minimizatiom");
    
    if(std::isnan(__atoms->hP))
        throw std::string("planck constant should be set prior to minimizatiom");

    
    if(std::isnan(__atoms->temp))
        throw std::string("temperature should be set prior to minimizatiom");
    
    if(C_DOF)
    {
        const int c_dim=__atoms->c_dim;
        type0* c=__atoms->c->begin();
        int natms_lcl=__atoms->natms_lcl;
        bool chk=true;
        type0 cv;
        for(int i=0;i<natms_lcl && chk;i++,c+=c_dim)
        {
            cv=1.0;
            for(int j=0;j<c_dim;j++)
            {
                if(c[j]<0.0) continue;
                if(c[j]==0.0) chk=false;
                cv-=c[j];
            }
            if(cv<=0.0) chk=false;
        }
        
        int err_lcl=chk ? 0:1;
        int err;
        MPI_Allreduce(&err_lcl,&err,1,Vec<int>::MPI_T,MPI_MAX,__atoms->world);
        if(err)
            throw std::string("for energy minimization no c can be either 0.0 or 1.0");
    }
}
/*--------------------------------------------
 init before a run
 --------------------------------------------*/
void MinCG3DMD::init()
{
    if(xprt)
    {
        try
        {
            xprt->atoms=atoms;
            xprt->init();
        }
        catch(std::string& err_msg)
        {
            fin();
            throw err_msg;
        }
    }
}
/*--------------------------------------------
 finishing minimization
 --------------------------------------------*/
void MinCG3DMD::fin()
{
    if(xprt)
    {
        xprt->fin();
        xprt->atoms=NULL;
    }
}
/*--------------------------------------------
 min
 --------------------------------------------*/
void MinCG3DMD::run(int nsteps)
{
    if(dynamic_cast<LineSearchBrent*>(ls)!=NULL && chng_box==false &&  X_DOF==false &&  ALPHA_DOF==false &&  C_DOF==true)
        return run<false,false,false,true>(dynamic_cast<LineSearchBrent*>(ls),nsteps);
    if(dynamic_cast<LineSearchGoldenSection*>(ls)!=NULL && chng_box==false &&  X_DOF==false &&  ALPHA_DOF==false &&  C_DOF==true)
        return run<false,false,false,true>(dynamic_cast<LineSearchGoldenSection*>(ls),nsteps);
    if(dynamic_cast<LineSearchBackTrack*>(ls)!=NULL && chng_box==false &&  X_DOF==false &&  ALPHA_DOF==false &&  C_DOF==true)
        return run<false,false,false,true>(dynamic_cast<LineSearchBackTrack*>(ls),nsteps);
    if(dynamic_cast<LineSearchBrent*>(ls)!=NULL && chng_box==false &&  X_DOF==false &&  ALPHA_DOF==true &&  C_DOF==false)
        return run<false,false,true,false>(dynamic_cast<LineSearchBrent*>(ls),nsteps);
    if(dynamic_cast<LineSearchGoldenSection*>(ls)!=NULL && chng_box==false &&  X_DOF==false &&  ALPHA_DOF==true &&  C_DOF==false)
        return run<false,false,true,false>(dynamic_cast<LineSearchGoldenSection*>(ls),nsteps);
    if(dynamic_cast<LineSearchBackTrack*>(ls)!=NULL && chng_box==false &&  X_DOF==false &&  ALPHA_DOF==true &&  C_DOF==false)
        return run<false,false,true,false>(dynamic_cast<LineSearchBackTrack*>(ls),nsteps);
    if(dynamic_cast<LineSearchBrent*>(ls)!=NULL && chng_box==false &&  X_DOF==false &&  ALPHA_DOF==true &&  C_DOF==true)
        return run<false,false,true,true>(dynamic_cast<LineSearchBrent*>(ls),nsteps);
    if(dynamic_cast<LineSearchGoldenSection*>(ls)!=NULL && chng_box==false &&  X_DOF==false &&  ALPHA_DOF==true &&  C_DOF==true)
        return run<false,false,true,true>(dynamic_cast<LineSearchGoldenSection*>(ls),nsteps);
    if(dynamic_cast<LineSearchBackTrack*>(ls)!=NULL && chng_box==false &&  X_DOF==false &&  ALPHA_DOF==true &&  C_DOF==true)
        return run<false,false,true,true>(dynamic_cast<LineSearchBackTrack*>(ls),nsteps);
    if(dynamic_cast<LineSearchBrent*>(ls)!=NULL && chng_box==false &&  X_DOF==true &&  ALPHA_DOF==false &&  C_DOF==false)
        return run<false,true,false,false>(dynamic_cast<LineSearchBrent*>(ls),nsteps);
    if(dynamic_cast<LineSearchGoldenSection*>(ls)!=NULL && chng_box==false &&  X_DOF==true &&  ALPHA_DOF==false &&  C_DOF==false)
        return run<false,true,false,false>(dynamic_cast<LineSearchGoldenSection*>(ls),nsteps);
    if(dynamic_cast<LineSearchBackTrack*>(ls)!=NULL && chng_box==false &&  X_DOF==true &&  ALPHA_DOF==false &&  C_DOF==false)
        return run<false,true,false,false>(dynamic_cast<LineSearchBackTrack*>(ls),nsteps);
    if(dynamic_cast<LineSearchBrent*>(ls)!=NULL && chng_box==false &&  X_DOF==true &&  ALPHA_DOF==false &&  C_DOF==true)
        return run<false,true,false,true>(dynamic_cast<LineSearchBrent*>(ls),nsteps);
    if(dynamic_cast<LineSearchGoldenSection*>(ls)!=NULL && chng_box==false &&  X_DOF==true &&  ALPHA_DOF==false &&  C_DOF==true)
        return run<false,true,false,true>(dynamic_cast<LineSearchGoldenSection*>(ls),nsteps);
    if(dynamic_cast<LineSearchBackTrack*>(ls)!=NULL && chng_box==false &&  X_DOF==true &&  ALPHA_DOF==false &&  C_DOF==true)
        return run<false,true,false,true>(dynamic_cast<LineSearchBackTrack*>(ls),nsteps);
    if(dynamic_cast<LineSearchBrent*>(ls)!=NULL && chng_box==false &&  X_DOF==true &&  ALPHA_DOF==true &&  C_DOF==false)
        return run<false,true,true,false>(dynamic_cast<LineSearchBrent*>(ls),nsteps);
    if(dynamic_cast<LineSearchGoldenSection*>(ls)!=NULL && chng_box==false &&  X_DOF==true &&  ALPHA_DOF==true &&  C_DOF==false)
        return run<false,true,true,false>(dynamic_cast<LineSearchGoldenSection*>(ls),nsteps);
    if(dynamic_cast<LineSearchBackTrack*>(ls)!=NULL && chng_box==false &&  X_DOF==true &&  ALPHA_DOF==true &&  C_DOF==false)
        return run<false,true,true,false>(dynamic_cast<LineSearchBackTrack*>(ls),nsteps);
    if(dynamic_cast<LineSearchBrent*>(ls)!=NULL && chng_box==false &&  X_DOF==true &&  ALPHA_DOF==true &&  C_DOF==true)
        return run<false,true,true,true>(dynamic_cast<LineSearchBrent*>(ls),nsteps);
    if(dynamic_cast<LineSearchGoldenSection*>(ls)!=NULL && chng_box==false &&  X_DOF==true &&  ALPHA_DOF==true &&  C_DOF==true)
        return run<false,true,true,true>(dynamic_cast<LineSearchGoldenSection*>(ls),nsteps);
    if(dynamic_cast<LineSearchBackTrack*>(ls)!=NULL && chng_box==false &&  X_DOF==true &&  ALPHA_DOF==true &&  C_DOF==true)
        return run<false,true,true,true>(dynamic_cast<LineSearchBackTrack*>(ls),nsteps);
    if(dynamic_cast<LineSearchBrent*>(ls)!=NULL && chng_box==true &&  X_DOF==false &&  ALPHA_DOF==false &&  C_DOF==false)
        return run<true,false,false,false>(dynamic_cast<LineSearchBrent*>(ls),nsteps);
    if(dynamic_cast<LineSearchGoldenSection*>(ls)!=NULL && chng_box==true &&  X_DOF==false &&  ALPHA_DOF==false &&  C_DOF==false)
        return run<true,false,false,false>(dynamic_cast<LineSearchGoldenSection*>(ls),nsteps);
    if(dynamic_cast<LineSearchBackTrack*>(ls)!=NULL && chng_box==true &&  X_DOF==false &&  ALPHA_DOF==false &&  C_DOF==false)
        return run<true,false,false,false>(dynamic_cast<LineSearchBackTrack*>(ls),nsteps);
    if(dynamic_cast<LineSearchBrent*>(ls)!=NULL && chng_box==true &&  X_DOF==false &&  ALPHA_DOF==false &&  C_DOF==true)
        return run<true,false,false,true>(dynamic_cast<LineSearchBrent*>(ls),nsteps);
    if(dynamic_cast<LineSearchGoldenSection*>(ls)!=NULL && chng_box==true &&  X_DOF==false &&  ALPHA_DOF==false &&  C_DOF==true)
        return run<true,false,false,true>(dynamic_cast<LineSearchGoldenSection*>(ls),nsteps);
    if(dynamic_cast<LineSearchBackTrack*>(ls)!=NULL && chng_box==true &&  X_DOF==false &&  ALPHA_DOF==false &&  C_DOF==true)
        return run<true,false,false,true>(dynamic_cast<LineSearchBackTrack*>(ls),nsteps);
    if(dynamic_cast<LineSearchBrent*>(ls)!=NULL && chng_box==true &&  X_DOF==false &&  ALPHA_DOF==true &&  C_DOF==false)
        return run<true,false,true,false>(dynamic_cast<LineSearchBrent*>(ls),nsteps);
    if(dynamic_cast<LineSearchGoldenSection*>(ls)!=NULL && chng_box==true &&  X_DOF==false &&  ALPHA_DOF==true &&  C_DOF==false)
        return run<true,false,true,false>(dynamic_cast<LineSearchGoldenSection*>(ls),nsteps);
    if(dynamic_cast<LineSearchBackTrack*>(ls)!=NULL && chng_box==true &&  X_DOF==false &&  ALPHA_DOF==true &&  C_DOF==false)
        return run<true,false,true,false>(dynamic_cast<LineSearchBackTrack*>(ls),nsteps);
    if(dynamic_cast<LineSearchBrent*>(ls)!=NULL && chng_box==true &&  X_DOF==false &&  ALPHA_DOF==true &&  C_DOF==true)
        return run<true,false,true,true>(dynamic_cast<LineSearchBrent*>(ls),nsteps);
    if(dynamic_cast<LineSearchGoldenSection*>(ls)!=NULL && chng_box==true &&  X_DOF==false &&  ALPHA_DOF==true &&  C_DOF==true)
        return run<true,false,true,true>(dynamic_cast<LineSearchGoldenSection*>(ls),nsteps);
    if(dynamic_cast<LineSearchBackTrack*>(ls)!=NULL && chng_box==true &&  X_DOF==false &&  ALPHA_DOF==true &&  C_DOF==true)
        return run<true,false,true,true>(dynamic_cast<LineSearchBackTrack*>(ls),nsteps);
    if(dynamic_cast<LineSearchBrent*>(ls)!=NULL && chng_box==true &&  X_DOF==true &&  ALPHA_DOF==false &&  C_DOF==false)
        return run<true,true,false,false>(dynamic_cast<LineSearchBrent*>(ls),nsteps);
    if(dynamic_cast<LineSearchGoldenSection*>(ls)!=NULL && chng_box==true &&  X_DOF==true &&  ALPHA_DOF==false &&  C_DOF==false)
        return run<true,true,false,false>(dynamic_cast<LineSearchGoldenSection*>(ls),nsteps);
    if(dynamic_cast<LineSearchBackTrack*>(ls)!=NULL && chng_box==true &&  X_DOF==true &&  ALPHA_DOF==false &&  C_DOF==false)
        return run<true,true,false,false>(dynamic_cast<LineSearchBackTrack*>(ls),nsteps);
    if(dynamic_cast<LineSearchBrent*>(ls)!=NULL && chng_box==true &&  X_DOF==true &&  ALPHA_DOF==false &&  C_DOF==true)
        return run<true,true,false,true>(dynamic_cast<LineSearchBrent*>(ls),nsteps);
    if(dynamic_cast<LineSearchGoldenSection*>(ls)!=NULL && chng_box==true &&  X_DOF==true &&  ALPHA_DOF==false &&  C_DOF==true)
        return run<true,true,false,true>(dynamic_cast<LineSearchGoldenSection*>(ls),nsteps);
    if(dynamic_cast<LineSearchBackTrack*>(ls)!=NULL && chng_box==true &&  X_DOF==true &&  ALPHA_DOF==false &&  C_DOF==true)
        return run<true,true,false,true>(dynamic_cast<LineSearchBackTrack*>(ls),nsteps);
    if(dynamic_cast<LineSearchBrent*>(ls)!=NULL && chng_box==true &&  X_DOF==true &&  ALPHA_DOF==true &&  C_DOF==false)
        return run<true,true,true,false>(dynamic_cast<LineSearchBrent*>(ls),nsteps);
    if(dynamic_cast<LineSearchGoldenSection*>(ls)!=NULL && chng_box==true &&  X_DOF==true &&  ALPHA_DOF==true &&  C_DOF==false)
        return run<true,true,true,false>(dynamic_cast<LineSearchGoldenSection*>(ls),nsteps);
    if(dynamic_cast<LineSearchBackTrack*>(ls)!=NULL && chng_box==true &&  X_DOF==true &&  ALPHA_DOF==true &&  C_DOF==false)
        return run<true,true,true,false>(dynamic_cast<LineSearchBackTrack*>(ls),nsteps);
    if(dynamic_cast<LineSearchBrent*>(ls)!=NULL && chng_box==true &&  X_DOF==true &&  ALPHA_DOF==true &&  C_DOF==true)
        return run<true,true,true,true>(dynamic_cast<LineSearchBrent*>(ls),nsteps);
    if(dynamic_cast<LineSearchGoldenSection*>(ls)!=NULL && chng_box==true &&  X_DOF==true &&  ALPHA_DOF==true &&  C_DOF==true)
        return run<true,true,true,true>(dynamic_cast<LineSearchGoldenSection*>(ls),nsteps);
    if(dynamic_cast<LineSearchBackTrack*>(ls)!=NULL && chng_box==true &&  X_DOF==true &&  ALPHA_DOF==true &&  C_DOF==true)
        return run<true,true,true,true>(dynamic_cast<LineSearchBackTrack*>(ls),nsteps);
}
/*------------------------------------------------------------------------------------------------------------------------------------
 
 ------------------------------------------------------------------------------------------------------------------------------------*/
PyObject* MinCG3DMD::__new__(PyTypeObject* type,PyObject* args,PyObject* kwds)
{
    Object* __self=reinterpret_cast<Object*>(type->tp_alloc(type,0));
    PyObject* self=reinterpret_cast<PyObject*>(__self);
    return self;
}
/*--------------------------------------------
 
 --------------------------------------------*/
int MinCG3DMD::__init__(PyObject* self,PyObject* args,PyObject* kwds)
{
    FuncAPI<type0,symm<bool[__dim__][__dim__]>,bool,type0,type0,OP<LineSearch>> f("__init__",{"e_tol","H_dof","affine","max_dx","max_dalpha","ls"});
    f.noptionals=6;
    f.logics<0>()[0]=VLogics("ge",0.0);
    f.logics<3>()[0]=VLogics("gt",0.0);
    f.logics<4>()[0]=VLogics("gt",0.0);
    
    //set the defualts
    f.val<0>()=sqrt(std::numeric_limits<type0>::epsilon());
    for(int i=0;i<__dim__;i++) for(int j=0;j<__dim__;j++)f.val<1>()[i][j]=false;
    f.val<2>()=false;
    f.val<3>()=1.0;
    f.val<4>()=0.1;
    PyObject* empty_tuple=PyTuple_New(0);
    PyObject* empty_dict=PyDict_New();
    PyObject* __ls=LineSearchBackTrack::__new__(&LineSearchBackTrack::TypeObject,empty_tuple,empty_dict);
    LineSearchBackTrack::__init__(__ls,empty_tuple,empty_dict);
    Py_DECREF(empty_dict);
    Py_DECREF(empty_tuple);
    f.val<5>().ob=__ls;
    
    
    if(f(args,kwds)==-1) return -1;
    
    
    
    Object* __self=reinterpret_cast<Object*>(self);
    Py_INCREF(f.val<5>().ob);
    __self->ls=reinterpret_cast<LineSearch::Object*>(f.val<5>().ob);
    __self->min=new MinCG3DMD(f.val<0>(),f.val<1>(),f.val<2>(),f.val<3>(),f.val<4>(),&(__self->ls->ls));
    __self->xprt=NULL;
    
    return 0;
}
/*--------------------------------------------
 
 --------------------------------------------*/
PyObject* MinCG3DMD::__alloc__(PyTypeObject* type,Py_ssize_t)
{
    Object* __self=new Object;
    Py_TYPE(__self)=type;
    Py_REFCNT(__self)=1;
    __self->min=NULL;
    __self->ls=NULL;
    __self->xprt=NULL;
    return reinterpret_cast<PyObject*>(__self);
}
/*--------------------------------------------
 
 --------------------------------------------*/
void MinCG3DMD::__dealloc__(PyObject* self)
{
    Object* __self=reinterpret_cast<Object*>(self);
    delete __self->min;
    __self->min=NULL;
    if(__self->ls) Py_DECREF(__self->ls);
    __self->ls=NULL;
    if(__self->xprt) Py_DECREF(__self->xprt);
    __self->xprt=NULL;
    delete __self;
}
/*--------------------------------------------*/
PyTypeObject MinCG3DMD::TypeObject={PyObject_HEAD_INIT(NULL)};
/*--------------------------------------------*/
int MinCG3DMD::setup_tp()
{
    TypeObject.tp_name="mapp.dmd.min_cg";
    TypeObject.tp_doc=R"---(
    __init__(e_tol=1.0e-8,H_dof=[[False],[False,False],[False,False,False]],affine=False,max_dx=1.0,max_dalpha=0.1,ls=mapp.dmd.ls_bt())
    
    CG minimization algorithm
        
    Parameters
    ----------
    e_tol : double
       Energy tolerance criterion for stopping minimization
    H_dof : symm<bool[dim][dim]>
       Unitcell degrees of freedom during minimization, here dim is the dimension of simulation
    affine : bool
       If set to True atomic displacements would be affine
    max_dx : double
       Maximum displacement of any atom in one step of minimization
    max_dalpha : double
       Maximum change in alpha component of any atom in one step of minimization
    ls : mapp.ls
       Line search method

    Notes
    -----
    Cojugate Gradient (CG) algorithm for minimization.    

    )---";
    
    TypeObject.tp_flags=Py_TPFLAGS_DEFAULT;
    TypeObject.tp_basicsize=sizeof(Object);
    
    TypeObject.tp_new=__new__;
    TypeObject.tp_init=__init__;
    TypeObject.tp_alloc=__alloc__;
    TypeObject.tp_dealloc=__dealloc__;
    setup_tp_methods();
    TypeObject.tp_methods=methods;
    setup_tp_getset();
    TypeObject.tp_getset=getset;
    
    int ichk=PyType_Ready(&TypeObject);
    if(ichk<0) return ichk;
    Py_INCREF(&TypeObject);
    GET_WRAPPER_DOC(TypeObject,__init__)=NULL;
    return ichk;
}
/*--------------------------------------------*/
PyGetSetDef MinCG3DMD::getset[]=EmptyPyGetSetDef(11);
/*--------------------------------------------*/
void MinCG3DMD::setup_tp_getset()
{
    getset_e_tol(getset[0]);
    getset_H_dof(getset[1]);
    getset_max_dx(getset[2]);
    getset_max_dalpha(getset[3]);
    getset_ls(getset[4]);
    getset_ntally(getset[5]);
    getset_export(getset[6]);
    getset_x_dof(getset[7]);
    getset_alpha_dof(getset[8]);
    getset_c_dof(getset[9]);
}
/*--------------------------------------------*/
PyMethodDef MinCG3DMD::methods[]=EmptyPyMethodDef(2);
/*--------------------------------------------*/
void MinCG3DMD::setup_tp_methods()
{
    ml_run(methods[0]);
}
/*--------------------------------------------
 
 --------------------------------------------*/
void MinCG3DMD::getset_max_dalpha(PyGetSetDef& getset)
{
    getset.name=(char*)"max_dalpha";
    getset.doc=(char*)R"---(
    (double) mximum alpha change
    
    Maximum change in alpha component of any atom in one step of minimization
    )---";
    getset.get=[](PyObject* self,void*)->PyObject*
    {
        return var<type0>::build(reinterpret_cast<Object*>(self)->min->max_dx);
    };
    getset.set=[](PyObject* self,PyObject* op,void*)->int
    {
        VarAPI<type0> max_dalpha("max_dalpha");
        max_dalpha.logics[0]=VLogics("gt",0.0);
        int ichk=max_dalpha.set(op);
        if(ichk==-1) return -1;
        reinterpret_cast<Object*>(self)->min->max_dalpha=max_dalpha.val;
        return 0;
    };
}
/*--------------------------------------------
 
 --------------------------------------------*/
void MinCG3DMD::getset_export(PyGetSetDef& getset)
{
    getset.name=(char*)"export";
    getset.doc=(char*)R"---(
    (mapp.dmd.export) export object
    
    Export object to record the snapshots of the system while minimizing
    )---";
    getset.get=[](PyObject* self,void*)->PyObject*
    {
        ExportDMD::Object* xprt=reinterpret_cast<Object*>(self)->xprt;
        if(!xprt) Py_RETURN_NONE;
        Py_INCREF(xprt);
        return reinterpret_cast<PyObject*>(xprt);
    };
    getset.set=[](PyObject* self,PyObject* op,void*)->int
    {
        VarAPI<OP<ExportDMD>> xprt("export");
        int ichk=xprt.set(op);
        if(ichk==-1) return -1;
        if(reinterpret_cast<Object*>(self)->xprt) Py_DECREF(reinterpret_cast<Object*>(self)->xprt);
        Py_INCREF(xprt.val.ob);
        reinterpret_cast<Object*>(self)->xprt=reinterpret_cast<ExportDMD::Object*>(xprt.val.ob);
        return 0;
    };
}
/*--------------------------------------------
 
 --------------------------------------------*/
void MinCG3DMD::getset_x_dof(PyGetSetDef& getset)
{
    getset.name=(char*)"x_dof";
    getset.doc=(char*)R"---(
    (bool) if set true position of atoms will considered as degrees of freedom
    
    If set true position of atoms will considered as degrees of freedom
    )---";
    getset.get=[](PyObject* self,void*)->PyObject*
    {
        return var<bool>::build(reinterpret_cast<Object*>(self)->min->X_DOF);
    };
    getset.set=[](PyObject* self,PyObject* op,void*)->int
    {
        VarAPI<bool> x_dof("x_dof");
        if(x_dof.set(op)==-1) return -1;
        reinterpret_cast<Object*>(self)->min->X_DOF=x_dof.val;
        return 0;
    };
}
/*--------------------------------------------
 
 --------------------------------------------*/
void MinCG3DMD::getset_alpha_dof(PyGetSetDef& getset)
{
    getset.name=(char*)"alpha_dof";
    getset.doc=(char*)R"---(
    (bool) if set true alpha of atoms will considered as degrees of freedom
    
    If set true alpha of atoms will considered as degrees of freedom
    )---";
    getset.get=[](PyObject* self,void*)->PyObject*
    {
        return var<bool>::build(reinterpret_cast<Object*>(self)->min->ALPHA_DOF);
    };
    getset.set=[](PyObject* self,PyObject* op,void*)->int
    {
        VarAPI<bool> alpha_dof("alpha_dof");
        if(alpha_dof.set(op)==-1) return -1;
        reinterpret_cast<Object*>(self)->min->ALPHA_DOF=alpha_dof.val;
        return 0;
    };
}
/*--------------------------------------------
 
 --------------------------------------------*/
void MinCG3DMD::getset_c_dof(PyGetSetDef& getset)
{
    getset.name=(char*)"c_dof";
    getset.doc=(char*)R"---(
    (bool) if set true c of atoms will considered as degrees of freedom
        
        If set true c of atoms will considered as degrees of freedom
        )---";
        getset.get=[](PyObject* self,void*)->PyObject*
    {
        return var<bool>::build(reinterpret_cast<Object*>(self)->min->C_DOF);
    };
    getset.set=[](PyObject* self,PyObject* op,void*)->int
    {
        VarAPI<bool> c_dof("c_dof");
        if(c_dof.set(op)==-1) return -1;
        reinterpret_cast<Object*>(self)->min->C_DOF=c_dof.val;
        return 0;
    };
}
/*--------------------------------------------
 
 --------------------------------------------*/
void MinCG3DMD::ml_run(PyMethodDef& tp_methods)
{
    tp_methods.ml_flags=METH_VARARGS | METH_KEYWORDS;
    tp_methods.ml_name="run";
    tp_methods.ml_meth=(PyCFunction)(PyCFunctionWithKeywords)(
    [](PyObject* self,PyObject* args,PyObject* kwds)->PyObject*
    {
        Object* __self=reinterpret_cast<Object*>(self);
        FuncAPI<OP<AtomsDMD>,int> f("run",{"atoms","max_nsteps"});
        f.logics<1>()[0]=VLogics("ge",0);
        
        
        
        if(f(args,kwds)) return NULL;
        
        
        AtomsDMD* __atoms=reinterpret_cast<AtomsDMD::Object*>(f.val<0>().ob)->atoms;
        ForceFieldDMD* __ff=reinterpret_cast<AtomsDMD::Object*>(f.val<0>().ob)->ff;
        ExportDMD* __xprt=__self->xprt==NULL ? NULL:__self->xprt->xprt;
        

        
        try
        {
            __self->min->pre_run_chk(__atoms,__ff);
        }
        catch(std::string& err_msg)
        {
            PyErr_SetString(PyExc_TypeError,err_msg.c_str());
            return NULL;
        }
        
        __self->min->atoms=__atoms;
        __self->min->ff=__ff;
        __self->min->xprt=__xprt;
        
        try
        {
            __self->min->run(f.val<1>());
        }
        catch(std::string& err_msg)
        {
            PyErr_SetString(PyExc_TypeError,err_msg.c_str());
            return NULL;
        }
        
        __self->min->xprt=NULL;
        __self->min->ff=NULL;
        __self->min->atoms=NULL;
        
        Py_RETURN_NONE;
    });
    
    tp_methods.ml_doc=(char*)R"---(
    run(atoms,max_nsteps)
   
    Execute minimization
    
    This method starts the energy minimization for a given atoms object and maximum number of steps.
    
    Parameters
    ----------
    atoms : mapp.md.atoms
        System of interest
    max_nsteps : int
        Maximum number of steps to achieve energy minimization
        
    Returns
    -------
    None

    )---";
}



