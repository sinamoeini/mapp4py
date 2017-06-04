#ifndef __MAPP__print__
#define __MAPP__print__
#include <string>
#include <cstring>
#include <iostream>
#include <limits>
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
    
    template<class T>
    inline static std::string to_string(T __val)
    {
        if(std::is_integral<T>::value)
        {
            if(std::is_unsigned<T>::value)
                return std::to_string(static_cast<unsigned long long>(__val));
            else
                return std::to_string(static_cast<long long>(__val));
        }
        
        return std::to_string(static_cast<long double>(__val));
    };
    
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
    return static_cast<int>(MAX(Print::to_string(std::numeric_limits<T>::max()).length(),Print::to_string(std::numeric_limits<T>::min()).length()));
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<class T>
inline constexpr int Print::delta_length()
{
    return max_length<T>()-min_length<T>();
}

#endif


