
/* min & max macro */
#define MIN(A,B) (A<B?A:B)
#define MAX(A,B) (A>B?A:B)

#define RESTRICT __restrict__


#ifndef __MAPP__type_def__
#define __MAPP__type_def__

#define TOSTRING(i) std::to_string(static_cast<long long>(i))

typedef double type0;
typedef unsigned char elem_type;
typedef unsigned char byte;

namespace MAPP_NS
{
    template<int N, elem_type ... Rest>
    class dof_mask
    {
    public:
    static constexpr const elem_type (&value)[sizeof...(Rest)+N] = dof_mask<N-1,1<<(static_cast<int>(sizeof(elem_type)*8)-N), Rest...>::value;
    };
        
    template<elem_type... Rest>
    class dof_mask<0, Rest...>
    {
    public:
        static constexpr elem_type value[sizeof...(Rest)] = {Rest... };
    };
        
    template<elem_type... Rest>
    constexpr elem_type dof_mask<0, Rest...>::value[];
    
    
    constexpr int __dim__=3;
    constexpr type0 ___dim___=static_cast<type0>(__dim__);
    constexpr int __nvoigt__=__dim__*(__dim__+1)/2;
    constexpr elem_type elem_mask=1<<(8*sizeof(elem_type)-__dim__);
    constexpr elem_type dofs_mask=((1<<__dim__)-1)*elem_mask;
    constexpr type0 sqrt_2=1.4142135623730951454746218587388284504413604736328125;
}
#endif /* type_def_h */
