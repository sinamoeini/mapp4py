#ifndef __MAPP__print__
#define __MAPP__print__
#include <string>
#include <limits.h>
class Print
{
private:
    template<int>
    static constexpr int ndigits();
public:
    
    template <class... Ts>
    inline static char* vprintf(char*,Ts...);
    template <class... Ts>
    inline static char* vprintf(const char*,Ts...);
    
    inline static char* vprintf(const char*);
    
    inline static char* vprintf(char*);
   
    template<class>
    inline static constexpr int min_length();
    template<class>
    inline static constexpr int max_length();
    
    template<class>
    inline static constexpr int delta_length();
    
    /*
    template<class T>
    inline static char* to_str(T& v)
    {
        if(std::is_same<T,int>::value)
            return vprintf("%d",v);
        if(std::is_same<T,size_t>::value)
            return vprintf("%zu",v);
        else if(std::is_same<T,double>::value)
            return vprintf("%g",v);
        return NULL;
    }
    
    template <typename... T0>
    inline static void append(char*& buff,const char* format,T0... vals)
    {
        size_t len=snprintf(NULL,0,format,vals...)+1;
        if(buff!=NULL)
        {
            size_t old_len=strlen(buff);
            char* buff_=new char[old_len+len];
            memcpy(buff_,buff,old_len);
            sprintf(buff_+old_len,format,vals...);
            delete [] buff;
            buff=buff_;
        }
        else
        {
            buff=new char[len];
            sprintf(buff,format,vals...);
        }
    }
    
    static inline void append(char*& buff,const char* format)
    {
        size_t len=strlen(format)+1;
        if(buff)
        {
            size_t old_len=strlen(buff);
            char* buff_=new char[old_len+len];
            memcpy(buff_,buff,old_len);
            memcpy(buff_+old_len,format,len);
            delete [] buff;
            buff=buff_;
        }
        else
        {
            buff=new char[len];
            memcpy(buff,format,len);
        }
    }
    
    static inline void prepend(const char* format,char*& buff)
    {
        size_t len=strlen(format)+1;
        if(buff)
        {
            size_t old_len=strlen(buff);
            char* buff_=new char[old_len+len];
            memcpy(buff_,format,len-1);
            memcpy(buff_+len-1,buff,old_len+1);
            delete [] buff;
            buff=buff_;
        }
        else
        {
            buff=new char[len];
            memcpy(buff,format,len);
        }
    }
    
    template<class T0,class T1,class ... Ts>
    inline static char* concatenate(T0 l,T1 r, Ts ... args)
    {
        return concatenate(concatenate(l,r),args...);
    }
    
    inline static char* concatenate(char* str)
    {
        return str;
    }
    
    inline static char* concatenate(const char* str)
    {
        size_t str_sz=strlen(str);
        char* __str=new char [str_sz+1];
        memcpy(__str,str,str_sz+1);
        return __str;
    }
    
    inline static char* concatenate(char* l,char* r)
    {
        if(l && r)
        {
            size_t l_sz=strlen(l);
            size_t r_sz=strlen(r);
            char* ans=new char[r_sz+l_sz+1];
            memcpy(ans,l,l_sz);
            memcpy(ans+l_sz,r,r_sz);
            ans[r_sz+l_sz]='\0';
            delete [] l;
            delete [] r;
            return ans;
        }
        if(r) return r;
        return l;
    }
    
    inline static char* concatenate(const char* l,char* r)
    {
        if(r)
        {
            size_t l_sz=strlen(l);
            size_t r_sz=strlen(r);
            char* ans=new char[r_sz+l_sz+1];
            memcpy(ans,l,l_sz);
            memcpy(ans+l_sz,r,r_sz);
            ans[r_sz+l_sz]='\0';
            delete [] r;
            return ans;
        }

        return concatenate(l);
    }
    
    inline static char* concatenate(char* l,const char* r)
    {
        if(l)
        {
            size_t l_sz=strlen(l);
            size_t r_sz=strlen(r);
            char* ans=new char[r_sz+l_sz+1];
            memcpy(ans,l,l_sz);
            memcpy(ans+l_sz,r,r_sz);
            ans[r_sz+l_sz]='\0';
            delete [] l;
            return ans;
        }
        
        return concatenate(r);
    }
    inline static char* concatenate(const char* l,const char* r)
    {
        size_t l_sz=strlen(l);
        size_t r_sz=strlen(r);
        char* ans=new char[r_sz+l_sz+1];
        memcpy(ans,l,l_sz);
        memcpy(ans+l_sz,r,r_sz);
        ans[r_sz+l_sz]='\0';
        return ans;
    }*/
};
/*--------------------------------------------
 
 --------------------------------------------*/
template <class... Ts>
inline char* Print::vprintf(char* format,Ts... vals)
{
    size_t len=snprintf(NULL,0,format,vals...)+1;
    char* buff=new char[len];
    sprintf(buff,format,vals...);
    delete [] format;
    format=NULL;
    return buff;
}
/*--------------------------------------------
 
 --------------------------------------------*/
template <class... Ts>
inline char* Print::vprintf(const char* format,Ts... vals)
{
    size_t len=snprintf(NULL,0,format,vals...)+1;
    char* buff=new char[len];
    sprintf(buff,format,vals...);
    
    return buff;
}
/*--------------------------------------------
 
 --------------------------------------------*/
inline char* Print::vprintf(const char* format)
{
    size_t len=strlen(format)+1;
    char* buff=new char[len];
    memcpy(buff,format,len*sizeof(char));
    return buff;
}
/*--------------------------------------------
 
 --------------------------------------------*/
inline char* Print::vprintf(char* format)
{
    return format;
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<>
inline constexpr int Print::ndigits<0>()
{
    return 0;
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<int N>
constexpr int Print::ndigits()
{
    return 1+ndigits<N/10>();
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<>
inline constexpr int Print::min_length<float>()
{
    return 6;
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<>
inline constexpr int Print::min_length<double>()
{
    return 6;
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<>
inline constexpr int Print::min_length<long double>()
{
    return 6;
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<class T>
inline constexpr int Print::min_length()
{
    return 1;
}
#define MAX(A,B) (A>B?A:B)
/*--------------------------------------------
 
 --------------------------------------------*/
template<>
inline constexpr int Print::max_length<float>()
{
    return ndigits<std::numeric_limits<float>::max_exponent10>()+5;
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<>
inline constexpr int Print::max_length<double>()
{
    return ndigits<std::numeric_limits<double>::max_exponent10>()+5;
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<>
inline constexpr int Print::max_length<long double>()
{
    return ndigits<std::numeric_limits<long double>::max_exponent10>()+5;
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<class T>
inline constexpr int Print::max_length()
{
    return static_cast<int>(MAX(std::to_string(std::numeric_limits<T>::max()).length(),std::to_string(std::numeric_limits<T>::min()).length()));
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<class T>
inline constexpr int Print::delta_length()
{
    return max_length<T>()-min_length<T>();
}
#endif


