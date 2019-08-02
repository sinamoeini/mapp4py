#ifndef __MAPP__cell__
#define __MAPP__cell__
namespace MAPP_NS
{
    template<const int M>
    class Cell
    {
    private:
        class Atoms* atoms;
        type0 cell_size[__dim__];
        type0 s_lo[__dim__];
        type0 s_hi[__dim__];
        type0 cut_s[__dim__];
        int cell_denom[__dim__];
        int ncells_per_dim[__dim__];
        int ncells;
        
        const int nneighs;
        int* rel_neigh_lst;
        int (*rel_neigh_lst_coord)[__dim__];
        
        int* head_atm;
    public:
        template<class F0,class F1,class F2>
        void DO(bool,F0,F1,F2);
        
        Cell(class Atoms*);
        ~Cell();
        void setup();
        static int neigh_size();
        
    };
}
using namespace MAPP_NS;
/*--------------------------------------------
 constructor
 --------------------------------------------*/
template<const int M>
Cell<M>::Cell(Atoms* __atoms):
atoms(__atoms),
nneighs(neigh_size()),
head_atm(NULL)
{
    ncells=0;
    rel_neigh_lst=new int[nneighs];
    rel_neigh_lst_coord=new int[nneighs][__dim__];
    
    int countr[__dim__]={DESIG(__dim__,-M)};
    constexpr int max_nneighs=Algebra::pow<2*M+1,__dim__>();
    
    int sum,icurs=0;
    constexpr int rc_sq=M*M;
    for(int i=0;i<max_nneighs;i++)
    {
        sum=0;
        for(int j=0;j<__dim__;j++)
        {
            sum+=countr[j]*countr[j];
            if(countr[j]!=0)
                sum+=1-2*std::abs(countr[j]);
        }
        
        if(sum<rc_sq)
        {
            Algebra::V_eq<__dim__>(countr,rel_neigh_lst_coord[icurs]);
            icurs++;
        }
        
        countr[0]++;
        for(int j=0;j<__dim__-1;j++)
            if(countr[j]==M+1)
            {
                countr[j]=-M;
                countr[j+1]++;
            }
    }
    
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<const int M>
Cell<M>::~Cell()
{
    delete [] rel_neigh_lst_coord;
    delete [] rel_neigh_lst;
    delete [] head_atm;
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<const int M>
int Cell<M>::neigh_size()
{
    int countr[__dim__]={DESIG(__dim__,-M)};
    constexpr int max_nneighs=Algebra::pow<2*M+1,__dim__>();
    constexpr int rc_sq=M*M;
    int ans=0;
    for(int i=0;i<max_nneighs;i++)
    {
        int sum=0;
        for(int j=0;j<__dim__;j++)
        {
            sum+=countr[j]*countr[j];
            if(countr[j]!=0)
                sum+=1-2*std::abs(countr[j]);
        }
        
        if(sum<rc_sq) ++ans;
        
        countr[0]++;
        for(int j=0;j<__dim__-1;j++)
            if(countr[j]==M+1)
            {
                countr[j]=-M;
                countr[j+1]++;
            }
    }
    
    return ans;
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<const int M>
void Cell<M>::setup()
{
    ncells=1;
    Algebra::V_eq_x_mul_V<__dim__>(atoms->max_cut,atoms->depth_inv,cut_s);
    Algebra::V_eq<__dim__>(atoms->s_hi,s_hi);
    Algebra::V_eq<__dim__>(atoms->s_lo,s_lo);
    for(int i=0;i<__dim__;i++)
    {
        cell_size[i]=cut_s[i]/static_cast<type0>(M);
        ncells_per_dim[i]=static_cast<int>
        ((s_hi[i]-s_lo[i])/cell_size[i])+1+2*M;
        cell_denom[i]=ncells;
        ncells*=ncells_per_dim[i];
    }
    
    delete [] head_atm;
    head_atm=new int[ncells];
    
    for(int i=0;i<nneighs;i++)
        rel_neigh_lst[i]=Algebra::V_mul_V<__dim__>(cell_denom,rel_neigh_lst_coord[i]);
}
/*--------------------------------------------
 construct the cell list
 --------------------------------------------*/
template<const int M>template<class F0,class F1,class F2>
void Cell<M>::DO(bool chng_box,F0 pre,F1 mid,F2 post)
{
    if(chng_box) setup();
    
    for(int i=0;i<ncells;i++)
        head_atm[i]=-1;
    
    const int natms_lcl=atoms->natms_lcl;
    const int nall_lcl=atoms->natms_ph+natms_lcl;
    
    int* next_vec=NULL;
    int* cell_vec=NULL;
    if(nall_lcl)
    {
        next_vec=new int[nall_lcl];
        cell_vec=new int[nall_lcl];
    }
    
    const int s_dim=atoms->x->dim;
    type0* s=atoms->x->begin()+(nall_lcl-1)*s_dim;

    for(int i=nall_lcl-1,cell_no;i>natms_lcl-1;i--,s-=s_dim)
    {
        cell_no=0;
        Algebra::Do<__dim__>::func([&s,&cell_no,this](int i){
        if(s[i]<s_lo[i])
        {
            //0<=y y<m
            cell_no+=cell_denom[i]*(MIN(static_cast<int>((s[i]-(s_lo[i]-cut_s[i]))/cell_size[i]),M-1));
        }
        else if(s_hi[i]<=s[i])
        {
            //y<m+n+1
            cell_no+=cell_denom[i]*(MIN(MAX(static_cast<int>((s[i]-s_lo[i])/cell_size[i]),ncells_per_dim[i]-2*M-1),ncells_per_dim[i]-1-M)+M);
        }
        else
        {
            //m<=y y<=m+n+1
            cell_no+=cell_denom[i]*(MIN(static_cast<int>((s[i]-s_lo[i])/cell_size[i]),ncells_per_dim[i]-2*M-1)+M);
        }});

        cell_vec[i]=cell_no;
        next_vec[i]=head_atm[cell_no];
        head_atm[cell_no]=i;
    }
    
    for(int i=natms_lcl-1,cell_no;i>-1;i--,s-=s_dim)
    {
        cell_no=0;
        Algebra::Do<__dim__>::func([&s,&cell_no,this](int i){
        cell_no+=cell_denom[i]*(MIN(static_cast<int>((s[i]-s_lo[i])/cell_size[i]),ncells_per_dim[i]-2*M-1)+M);});
        
        cell_vec[i]=cell_no;
        next_vec[i]=head_atm[cell_no];
        head_atm[cell_no]=i;
    }
    
    atoms->s2x_all();
    
    
    for(int i=0,j,icell,jcell,ineigh;i<natms_lcl;i++)
    {
        pre(i);
        icell=cell_vec[i];
        j=-1;
        ineigh=-1;
        while(j==-1 && ++ineigh<nneighs)
        {
            jcell=icell+rel_neigh_lst[ineigh];
            j=head_atm[jcell];
            while(j!=-1)
            {
                if(i!=j) mid(i,j);
                j=next_vec[j];
            }
        }
        post(i);
    }
    
    delete [] cell_vec;
    delete [] next_vec;
    
}

#endif
