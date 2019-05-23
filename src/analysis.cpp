#include "analysis.h"
//#include "api.h"
#include "atoms_md.h"
#include "dynamic_md.h"
#include "neighbor_md.h"
#include "ff_md.h"
#include "memory.h"
using namespace MAPP_NS;
//using namespace Analysis;
/*--------------------------------------------
 
 --------------------------------------------*/
void BCTPolarity::ml_bct_polarity(PyMethodDef& tp_methods)
{
    tp_methods.ml_flags=METH_VARARGS | METH_KEYWORDS;
    tp_methods.ml_name="bct_polarity";
    tp_methods.ml_meth=(PyCFunction)(PyCFunctionWithKeywords)(
    [](PyObject* self,PyObject* args,PyObject* kwds)->PyObject*
    {
        AtomsMD::Object* __self=reinterpret_cast<AtomsMD::Object*>(self);
        FuncAPI<int,symm<type0**>,std::string*> f("bct_polarity",{"n","r_c","elems"});
        f.noptionals=1;
        f.logics<0>()[0]=VLogics("gt",0);
        f.logics<1>()[0]=VLogics("ge",0.0);
        AtomsMD* __atoms=__self->atoms;
        const std::string* names=__atoms->elements.names;
        const size_t nelems=__atoms->elements.nelems;
        if(f(args,kwds)) return NULL;
        if(f.remap<2,1>("elements present in system",names,nelems))
            return NULL;
        type0 max_cut=0.0;
        for(size_t i=0;i<nelems;i++)
            for(size_t j=0;j<i+1;j++)
                max_cut=MAX(max_cut,f.val<1>()[i][j]);
        if(max_cut==0.0)
        {
            PyErr_Format(PyExc_TypeError,"at least one cutoff radius should be greater than 0.0");
            return NULL;
        }
        
        
        type0* crl=polarity(__self->atoms,f.mov<1>(), f.val<0>());
        
        size_t __n=static_cast<size_t>(f.val<0>());
        size_t* __np=&__n;
        PyObject* op=var<type0*>::build(crl,&__np);
        Memory::dealloc(crl);
        return op;
    });

    tp_methods.ml_doc=R"---(
        
    )---";

}
/*--------------------------------------------
 
 --------------------------------------------*/
type0* BCTPolarity::polarity(AtomsMD* __atoms,type0**&& __cut,int n)
{
    
    ForceFieldMD* __ff=new ForceFieldZero(__atoms,std::move(__cut));
    __ff->neighbor->pair_wise=false;
    
    DynamicMD* dynamic=new DynamicMD(__atoms,__ff,false,{},{__atoms->x_d,__atoms->x_dof},{});
    dynamic->init();
    
    Vec<type0>* bct_pol=new Vec<type0>(__atoms,3,"bct_pol");
    int** neighbor_list=__ff->neighbor->neighbor_list;
    int* neighbor_list_size=__ff->neighbor->neighbor_list_size;
    const int natms_lcl=__atoms->natms_lcl;
    type0* xvec=__atoms->x->begin();
    type0* pvec=bct_pol->begin();
    for(int i=0;i<natms_lcl*__dim__;i++)
        pvec[i]=0.0;
    
    int* dummy=NULL;
    int* __neigh_list=NULL;
    type0* rsqs=NULL;
    Memory::alloc(dummy,__atoms->natms_lcl+__atoms->natms_ph);
    Memory::alloc(__neigh_list,__atoms->natms_lcl+__atoms->natms_ph);
    Memory::alloc(rsqs,__atoms->natms_lcl+__atoms->natms_ph);
    for(int iatm=0;iatm<natms_lcl;iatm++)
    {
        memcpy(__neigh_list,neighbor_list[iatm],neighbor_list_size[iatm]*sizeof(int));
        for(int jatm=0;jatm<neighbor_list_size[iatm];jatm++)
        {
            dummy[jatm]=jatm;
            rsqs[jatm]=Algebra::RSQ<__dim__>(xvec+iatm*__dim__,xvec+neighbor_list[iatm][jatm]*__dim__);
        }
        
         XMath::quicksort(dummy,dummy+neighbor_list_size[iatm]
        ,[&rsqs](int* i,int* j){return (rsqs[*i]<rsqs[*j]);}
        ,[&neighbor_list,&iatm](int* i,int* j){
            std::swap(*i,*j);
            
        });
        
        for(int jatm=0;jatm<neighbor_list_size[iatm];jatm++)
            neighbor_list[iatm][jatm]=__neigh_list[dummy[jatm]];
        
        calc_polarity(xvec, iatm,neighbor_list[iatm],neighbor_list_size[iatm],pvec+iatm*__dim__);
    }
    Memory::dealloc(__neigh_list);
    Memory::dealloc(rsqs);
    Memory::dealloc(dummy);
    
    dynamic->update_wo_x(bct_pol);
    
    type0* crl_lcl=NULL;
    Memory::alloc(crl_lcl,n);
    int* ns_lcl=NULL;
    Memory::alloc(ns_lcl,n);
    type0* crl=NULL;
    Memory::alloc(crl,n);
    int* ns=NULL;
    Memory::alloc(ns,n);
    for(int i=0;i<n;i++)
    {
        crl[i]=crl_lcl[i]=0.0;
        ns[i]=ns_lcl[i]=0;
    }
    
    
    type0 dr=__ff->max_cut/static_cast<type0>(n);
    type0 rsq,tmp;
    int p;
    type0 max_cut_sq=Algebra::pow<2>(__ff->max_cut);
    xvec=__atoms->x->begin();
    pvec=bct_pol->begin();
    for(int iatm=0;iatm<natms_lcl;iatm++)
    {
        const int list_size=neighbor_list_size[iatm];
        tmp=Algebra::V_mul_V<__dim__>(pvec+iatm*__dim__,pvec+iatm*__dim__);
        if(tmp==0.0) continue;
        for(int j=0,jatm;j<list_size;j++)
        {
            jatm=neighbor_list[iatm][j];
            rsq=Algebra::RSQ<__dim__>(xvec+iatm*__dim__,xvec+jatm*__dim__);
            if(rsq>=max_cut_sq) continue;
            tmp=Algebra::V_mul_V<__dim__>(pvec+jatm*__dim__,pvec+jatm*__dim__);
            if(tmp==0.0) continue;
            p=static_cast<int>(sqrt(rsq)/dr);
            //crl_lcl[p]+=fabs(Algebra::V_mul_V<__dim__>(pvec+iatm*__dim__,pvec+jatm*__dim__));
            crl_lcl[p]+=asin(fabs(Algebra::V_mul_V<__dim__>(pvec+iatm*__dim__,pvec+jatm*__dim__)));
            ns_lcl[p]++;
        }
    }
    
    MPI_Allreduce(crl_lcl,crl,n,Vec<type0>::MPI_T, MPI_SUM,__atoms->world);
    MPI_Allreduce(ns_lcl,ns,n,Vec<int>::MPI_T, MPI_SUM,__atoms->world);
    for(int i=0;i<n;i++)
        if(ns[i]!=0)
            crl[i]/=static_cast<type0>(ns[i]);
    
    Memory::dealloc(ns);
    Memory::dealloc(ns_lcl);
    Memory::dealloc(crl_lcl);
    
    delete  bct_pol;
    
    dynamic->fin();
    delete dynamic;
    delete __ff;
    
    return crl;
}
/*--------------------------------------------
 
 --------------------------------------------*/
int BCTPolarity::calc_polarity(type0* x,int iatm,int* lst,int sz,type0* polarity)
{
    
    if(sz<14) return 0;
    
    auto calc_rsq=
    [&x,&lst](int& i,int& j)->type0
    {
        return Algebra::RSQ<__dim__>(x+__dim__*lst[i],x+__dim__*lst[j]);
    };

    
    int lcl_neigh [14][6];
    type0 rsqs[14];
    type0 max_rsq[14];
    int key[14];
    
    
    for(int i=0;i<14;i++)
    {
        for(int j=0;j<14;j++)
        {
            key[j]=j;
            rsqs[j]=calc_rsq(i,j);
        }
        
        XMath::quicksort(key,key+14,
        [&rsqs](int* ikey,int* jkey){return (rsqs[*ikey] < rsqs[*jkey]); },
        [](int* ikey,int* jkey){std::swap(*ikey,*jkey);});
        
        memcpy(lcl_neigh[i],key+1,6*sizeof(int));
        max_rsq[i]=rsqs[key[6]];
    }
    
    for(int i=0;i<14;i++)
        key[i]=i;
    
    XMath::quicksort(key,key+14,
    [&max_rsq](int* ikey,int* jkey){return (max_rsq[*ikey] < max_rsq[*jkey]); },
    [](int* ikey,int* jkey){std::swap(*ikey,*jkey);});
    
    
    enum{H,S};
    int type[14];
    for(int i=0;i<8;i++)
        type[key[i]]=H;
    
    for(int i=8;i<14;i++)
        type[key[i]]=S;
    

    
    for(int i=0;i<14;i++)
    {
        if(type[i]==H)
        {
            int H_count=0;
            int S_count=0;
            for(int j=0;j<6;j++)
            {
                if(type[lcl_neigh[i][j]]==H)
                    H_count++;
                else
                    S_count++;
            }
            if(H_count!=3 || S_count!=3)
                return 0;
        }
        else
        {
            int H_count=0;
            for(int j=0;j<4;j++)
                if(type[lcl_neigh[i][j]]==H)
                    H_count++;
            
            if(H_count!=4)
                return 0;
        }
    }
    /*-------------------------------
     
     
         yz----------xyz
        /.            /|
       / .           / |
      /  .          /  |
     z------------zx   |
     |   .         |   |
     |   .         |   |
     |   y.........|..xy
     |  .          |  /
     | .           | /
     |.            |/
     o-------------x
     
     0:o    [-1,-1,-1]
     1:x    [+1,-1,-1]
     2:y    [-1,+1,-1]
     3:z    [-1,-1,+1]
     4:xyz  [+1,+1,+1]
     5:yz   [-1,+1,+1]
     6:zx   [+1,-1,+1]
     7:xy   [+1,+1,-1]
     -------------------------------*/
    
    
    
    int new_arrang[14];
    
    new_arrang[0]=key[0];
    int ic=1;
    for(int i=0,j;i<6;i++)
    {
        j=lcl_neigh[new_arrang[0]][i];
        if(type[j]==H)
            new_arrang[ic++]=j;
    }
    
    type0 HH[3][3];
    for(int i=0;i<3;i++)
        for(int j=0;j<3;j++)
            HH[i][j]=x[lst[new_arrang[i+1]]*__dim__+j]-x[lst[new_arrang[0]]*__dim__+j];
    
    if(M3DET(HH)<0.0)
        std::swap(new_arrang[2],new_arrang[3]);
    
    ic=6;
    for(int i=0,j;i<6;i++)
    {
        j=lcl_neigh[new_arrang[1]][i];
        if(type[j]==H && j!=new_arrang[0])
            new_arrang[ic++]=j;
    }
    
    for(int i=0,j;i<6;i++)
    {
        j=lcl_neigh[new_arrang[2]][i];
        if(type[j]!=H)
            continue;
        
        if(new_arrang[6]==j)
            std::swap(new_arrang[6],new_arrang[7]);
        if(new_arrang[0]!=j && new_arrang[6]!=j && new_arrang[7]!=j)
            new_arrang[5]=j;
    }
    
    new_arrang[4]=-1;
    for(int i=0;i<14 && new_arrang[4]==-1;i++)
    {
        if(type[i]!=H) continue;
        if(i==new_arrang[0]) continue;
        if(i==new_arrang[1]) continue;
        if(i==new_arrang[2]) continue;
        if(i==new_arrang[3]) continue;
        if(i==new_arrang[5]) continue;
        if(i==new_arrang[6]) continue;
        if(i==new_arrang[7]) continue;
        new_arrang[4]=i;
    }
    
    
    
    type0 F[3][3];
    type0 a[4];
    for(int i=0;i<3;i++)
    {
        for(int j=0;j<4;j++)
            a[j]=x[lst[new_arrang[j+4]]*__dim__+i]-x[lst[new_arrang[j]]*__dim__+i];
        
        F[i][0]=0.125*(a[0]-a[1]+a[2]+a[3]);
        F[i][1]=0.125*(a[0]+a[1]-a[2]+a[3]);
        F[i][2]=0.125*(a[0]+a[1]+a[2]-a[3]);
        
    }
    
    
    for(int i=0;i<3;i++)
        for(int j=0;j<3;j++)
        {
            HH[i][j]=0.0;
            for(int k=0;k<3;k++)
                HH[i][j]+=F[i][k]*F[j][k];
        }
    type0 p[3];
    
    /*
    if(atoms->my_p==0 && iatm<10)
    {
        printf("{");
        for(int i=0;i<3;i++)
        {
            printf("{%lf",HH[i][0]);
            for(int j=1;j<3;j++)
                printf(",%lf",HH[i][j]);
            printf("}");
            if(i!=2)
                printf(",");
        }
        printf("}\n");
    }
    */
    test(HH,p);
    polarity[0]=p[0];
    polarity[1]=p[1];
    polarity[2]=p[2];
    
    return 1;
    
}
/*--------------------------------------------
 
 --------------------------------------------*/
type0 BCTPolarity::test(type0 (&A)[3][3],type0 (&polar)[3])
{
    type0 p1=A[1][0]*A[1][0]+A[2][0]*A[2][0]+A[2][1]*A[2][1];
    
    
    if(p1==0.0)
    {
        polar[0]=polar[1]=polar[2]=0.0;
        if(A[0][0]>A[1][1] && A[0][0]>A[2][2])
        {
            polar[0]=1.0;
            return A[0][0];
        }
        else if(A[1][1]>A[0][0] && A[1][1]>A[2][2])
        {
            polar[1]=1.0;
            return A[1][1];
        }
        else if(A[2][2]>A[0][0] && A[2][2]>A[1][1])
        {
            polar[2]=1.0;
            return A[2][2];
        }
        
        
        
    }
    type0 eig[3];
    
    type0 q=(A[0][0]+A[1][1]+A[2][2])/3.0;
    type0 p2=(A[0][0]-q)*(A[0][0]-q)+(A[1][1]-q)*(A[1][1]-q)+(A[2][2]-q)*(A[2][2]-q)+2.0*p1;
    type0 p=sqrt(p2/6.0);
    type0 r=0.5*(
                 (A[0][0]-q)*(A[1][1]-q)*(A[2][2]-q)
                 +2.0*A[2][0]*A[1][0]*A[2][1]
                 -(A[0][0]-q)*A[2][1]*A[2][1]
                 -(A[1][1]-q)*A[2][0]*A[2][0]
                 -(A[2][2]-q)*A[1][0]*A[1][0])/(p*p*p);
    
    type0 phi;
    if(r<=-1.0)
        phi=M_PI/3.0;
    else if(r>=1.0)
        phi=0.0;
    else
        phi=acos(r)/3.0;
    
    eig[0]=q+2.0*p*cos(phi);
    eig[2]=q+2.0*p*cos(phi+2.0*M_PI/3.0);
    eig[1]=3.0*q-eig[0]-eig[2];
    
    
    auto inner=
    [] (type0(&x)[3],type0(&y)[3])->type0
    {
        
        return (x[0]*y[0]+x[1]*y[1]+x[2]*y[2]);
        
    };
    
    auto cross=
    [] (type0(&a)[3],type0(&b)[3],type0(&c)[3])->void
    {
        c[0]=a[1]*b[2]-a[2]*b[1];
        c[1]=a[2]*b[0]-a[0]*b[2];
        c[2]=a[0]*b[1]-a[1]*b[0];
    };
    
    auto normalize=
    [] (type0(&x)[3])->void
    {
        type0 r=sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]);
        x[0]/=r;
        x[1]/=r;
        x[2]/=r;
    };
    
    
    
    

    A[0][0]-=eig[0];
    A[1][1]-=eig[0];
    A[2][2]-=eig[0];
    
    type0 prods[3];
    normalize(A[0]);
    normalize(A[1]);
    normalize(A[2]);
    prods[0]=fabs(inner(A[1],A[2]));
    prods[1]=fabs(inner(A[2],A[0]));
    prods[2]=fabs(inner(A[0],A[1]));
    
    if(prods[0]<prods[1] && prods[0]<prods[2])
    {
        cross(A[1],A[2],polar);
    }
    else if(prods[1]<prods[0] && prods[1]<prods[2])
    {
        cross(A[2],A[0],polar);
    }
    else
    {
        cross(A[0],A[1],polar);
    }
    
    normalize(polar);
    
    return eig[0];
}
/*--------------------------------------------
 
 --------------------------------------------*/
void BCTPolarity::QR(type0(&A)[3][3],type0(&Q)[3][3],type0(&R)[3][3])
{
    auto normalize=
    [] (type0(&x)[3])->void
    {
        type0 r=0.0;
        r=sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]);
        x[0]/=r;
        x[1]/=r;
        x[2]/=r;
    };
    
    auto inner=
    [] (type0(&x)[3],type0(&y)[3])->type0
    {
        
        return (x[0]*y[0]+x[1]*y[1]+x[2]*y[2]);
        
    };

    
    memcpy(Q[0],A[0],3*sizeof(type0));
    normalize(Q[0]);
    
    
    type0 m,n=inner(A[1],Q[0]);
    for(int i=0;i<3;i++)
        Q[1][i]=A[1][i]-n*Q[0][i];
    normalize(Q[1]);
    
    n=inner(A[2],Q[0]);
    m=inner(A[2],Q[1]);
    for(int i=0;i<3;i++)
        Q[2][i]=A[2][i]-n*Q[0][i]-m*Q[1][i];
    normalize(Q[2]);
    
    for(int i=0;i<3;i++)
        for(int j=i;j<3;j++)
        {
            R[i][j]=0.0;
            for(int k=0;k<3;k++)
                R[i][j]+=Q[k][i]*A[k][j];
        }
    
    R[1][0]=R[2][0]=R[2][1]=0.0;
    
    
}
