#include "var.h"
#include <mpi.h>
/*--------------------------------------------*/
Var** Var::var_stack=NULL;
size_t Var::stack_size=0;
size_t Var::stack_capacity=0;
size_t Var::stack_grow=10;
const std::string type_attr<bool>::name(){return std::string("bool");};
const bool type_attr<bool>::zero=false;
const std::string type_attr<char>::name(){return std::string("char");};
const char type_attr<char>::zero=0;
const std::string type_attr<unsigned char>::name(){return std::string("unsigned char");};
const unsigned char type_attr<unsigned char>::zero=0;
const std::string type_attr<short>::name(){return std::string("short");};
const short type_attr<short>::zero=0;
const std::string type_attr<unsigned short>::name(){return std::string("unsigned short");};
const unsigned short type_attr<unsigned short>::zero=0;
const std::string type_attr<int>::name(){return std::string("int");};
const int type_attr<int>::zero=0;
const std::string type_attr<unsigned int>::name(){return std::string("unsigned int");};
const unsigned int type_attr<unsigned int>::zero=0;
const std::string type_attr<long>::name(){return std::string("long");};
const long type_attr<long>::zero=0;
const std::string type_attr<unsigned long>::name(){return std::string("unsigned long");};
const unsigned long type_attr<unsigned long>::zero=0;
const std::string type_attr<long long>::name(){return std::string("long long");};
const long long type_attr<long long>::zero=0;
const std::string type_attr<unsigned long long>::name(){return std::string("unsigned long long");};
const unsigned long long type_attr<unsigned long long>::zero=0;
const std::string type_attr<float>::name(){return std::string("float");};
const float type_attr<float>::zero=0;
const std::string type_attr<double>::name(){return std::string("double");};
const double type_attr<double>::zero=0;
const std::string type_attr<long double>::name(){return std::string("long double");};
const long double type_attr<long double>::zero=0;
const std::string type_attr<std::string>::name(){return std::string("std::string");};
const std::string type_attr<std::string>::zero=std::string();

static_assert(sizeof(bool)==sizeof(char),"sizeof(bool)!=sizeof(char)");
static_assert(sizeof(bool)==sizeof(npy_bool),"sizeof(bool)!=sizeof(npy_bool)");
static_assert(static_cast<npy_bool>(true)==NPY_TRUE,"true!=NPY_TRUE");
static_assert(static_cast<npy_bool>(false)==NPY_FALSE,"true!=NPY_TRUE");
/*--------------------------------------------*/
/*--------------------------------------------
 
 --------------------------------------------*/
void Var::push(Var* variable)
{
    if(stack_size+1>stack_capacity)
    {
        Var** var_stack_=new Var*[stack_capacity+1+stack_grow];
        memcpy(var_stack_,var_stack,stack_size*sizeof(Var*));
        delete [] var_stack;
        var_stack=var_stack_;
        stack_capacity+=1+stack_grow;
    }
    
    var_stack[stack_size]=variable;
    stack_size++;
}
/*--------------------------------------------
 
 --------------------------------------------*/
void Var::pop(Var* variable)
{
    size_t ivar=stack_size-1;
    for(;var_stack[ivar]!= variable;ivar--){}
    var_stack[ivar]=var_stack[stack_size-1];
    
    stack_size--;
    if(stack_capacity-stack_size>stack_grow)
    {
        Var** var_stack_=new Var*[stack_size];
        memcpy(var_stack_,var_stack,stack_size*sizeof(Var*));
        delete [] var_stack;
        var_stack=var_stack_;
        stack_capacity=stack_size;
    }
}
/*--------------------------------------------
 
 --------------------------------------------*/
void Var::replace(Var* key,Var* value)
{
    size_t ivar=stack_size-1;
    for(;var_stack[ivar]!= key;ivar--){}
    var_stack[ivar]=value;
}
/*----------------------------------------------------------------------------------------
 
 ----------------------------------------------------------------------------------------*/
Var::Var(int rank_,size_t hash_code_,const char* __name):
rank(rank_),
hash_code(hash_code_),
name(__name),
obj_set(false),
pushed(true),
root(true)
{
    push(this);
}
/*--------------------------------------------
 
 --------------------------------------------*/
Var::Var(int rank_,size_t hash_code_,std::string&& name_):
rank(rank_),
hash_code(hash_code_),
name(name_),
obj_set(false),
pushed(true),
root(false)
{
    push(this);
}
/*--------------------------------------------
 
 --------------------------------------------*/
Var::Var(int __rank,size_t hash_code_,const std::string& __name):
rank(__rank),
hash_code(hash_code_),
name(__name),
obj_set(false),
pushed(true),
root(false)
{
    push(this);
}
/*--------------------------------------------
 
 --------------------------------------------*/
Var::Var(int rank_,size_t hash_code_):
rank(rank_),
hash_code(hash_code_),
obj_set(false),
pushed(false),
root(false)
{
}
/*--------------------------------------------
 
 --------------------------------------------*/
Var::~Var()
{
    if(pushed) pop(this);
}
/*--------------------------------------------
 
 --------------------------------------------*/
Var::Var(const Var& r):
rank(r.rank),
hash_code(r.hash_code),
obj_set(r.obj_set),
size(r.size),
pushed(r.pushed),
root(r.root)
{
    if(pushed) push(this);
}
/*--------------------------------------------
 
 --------------------------------------------*/
Var::Var(Var&& r):
rank(r.rank),
hash_code(r.hash_code),
name(std::move(r.name)),
obj_set(r.obj_set),
size(r.size),
pushed(r.pushed),
root(r.root)
{
    if(pushed) replace(&r,this);
    
    r.obj_set=false;
    r.size=0;
    r.pushed=false;
    r.root=false;
}
/*--------------------------------------------
 
 --------------------------------------------*/
size_t* Var::is_array(const int r)
{
    if(!rank) return NULL;
    size_t** sz_lst=size_list();
    size_t sz=1;
    for(int i=0;i<r;i++)
    {
        for(int j=1;j<sz;j++)
        {
            if(sz_lst[i][j]!=sz_lst[i][0])
            {
                delete [] *sz_lst;
                delete [] sz_lst;
                return NULL;
            }
        }
        sz*=sz_lst[i][0];
    }
    
    size_t* shape=new size_t[r];
    for(int i=0;i<r;i++)
        shape[i]=sz_lst[i][0];
    
    delete [] *sz_lst;
    delete [] sz_lst;
    
    
    return shape;
}
/*--------------------------------------------
 
 --------------------------------------------*/
size_t Var::is_square()
{
    size_t* shape=is_array(2);
    if(shape && shape[0]==shape[1])
    {
        size_t ans=shape[0];
        delete [] shape;
        return ans;
    }
    
    delete [] shape;
    return 0;
}
/*--------------------------------------------
 
 --------------------------------------------*/
size_t Var::is_triangular()
{
    if(rank<2) return 0;
    
    size_t** sz_lst=size_list();
    for(size_t i=0;i<sz_lst[0][0];i++)
        if(sz_lst[1][i]!=i+1)
        {
            delete [] *sz_lst;
            delete [] sz_lst;
            return 0;
        }
    
    size_t ans=sz_lst[0][0];
    delete [] *sz_lst;
    delete [] sz_lst;
    return ans;
}
/*--------------------------------------------
 
 --------------------------------------------*/
size_t** Var::size_list()
{
    if(!rank) return NULL;
    size_t* nelem=(rank==0) ? NULL:new size_t[rank];
    for(int i=0;i<rank;i++) nelem[i]=0;
    size_list(nelem);
    size_t tot=0;
    for(int i=0;i<rank;i++) tot+=nelem[i];
    
    size_t** sz_lst=new size_t*[rank];
    *sz_lst=new size_t[tot];
    for(int i=1;i<rank;i++)
        sz_lst[i]=sz_lst[i-1]+nelem[i-1];
    size_t** sh=(rank==0) ? NULL:new size_t*[rank];
    for(int i=1;i<rank;i++) sh[i]=sz_lst[i];
    size_list(sh);
    delete [] sh;
    delete [] nelem;
    return sz_lst;
}
/*--------------------------------------------
 
 --------------------------------------------*/
void Var::size_list(size_t* sz)
{
    *sz+=size;
    if(rank==1) return;
    for(size_t i=0;i<size;i++)
        (*this)[i].size_list(sz+1);
}
/*--------------------------------------------
 
 --------------------------------------------*/
void Var::size_list(size_t** sz)
{
    **sz=size;
    *sz+=1;
    if(rank==1) return;
    for(size_t i=0;i<size;i++)
        (*this)[i].size_list(sz+1);
}
/*--------------------------------------------
 
 --------------------------------------------*/
bool Var::is_set()
{
    return obj_set;
}
/*--------------------------------------------
 
 --------------------------------------------*/
bool Var::belongs_to(Var& r)
{
    if(this->hash_code!= r.hash_code)
        return false;
    
    if(this->rank==r.rank)
        return *this==r;
    
    if(rank==r.rank+1)
        for(size_t i=0;i<size;i++)
            if(*this==r[i])
                return true;
    
    return false;
}
/*--------------------------------------------
 
 --------------------------------------------*/
bool Var::operator == (Var& r)
{
    if(this->hash_code!= r.hash_code || this->rank!= r.rank || this->size!= r.size)
        return false;
    for(size_t i=0;i<size;i++)
        if(!((*this)[i]==r[i]))
            return false;
    return true;
}
/*--------------------------------------------
 
 --------------------------------------------*/
bool Var::operator != (Var& r)
{
    return !operator==(r);
}
/*--------------------------------------------
 
 --------------------------------------------*/
bool Var::operator > (Var&)
{
    return false;
}
/*--------------------------------------------
 
 --------------------------------------------*/
bool Var::operator >= (Var&)
{
    return false;
}
/*--------------------------------------------
 
 --------------------------------------------*/
bool Var::operator <= (Var&)
{
    return false;
}
/*--------------------------------------------
 
 --------------------------------------------*/
bool Var::operator < (Var&)
{
    return false;
}
