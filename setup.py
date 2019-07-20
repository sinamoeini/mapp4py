try:
    import setuptools
except ImportError:
    setuptools = None
if setuptools:
    from setuptools import Extension
    from setuptools import setup
    from setuptools import find_packages
else:
    from distutils.core import Extension
    from distutils.core import setup


def has_cxx_11():
    import tempfile,os,sys
    fd, fname = tempfile.mkstemp('.cpp','conftest', text=True)
    f = os.fdopen(fd, 'w')
    try:
        f.write('#ifndef __cplusplus\n');
        f.write('#error \"This is not a C++ compiler\"\n');
        f.write('#else\n');
        f.write('namespace cxx11\n');
        f.write('{\n');
        f.write('  namespace test_static_assert\n');
        f.write('  {\n');
        f.write('    template <typename T>\n');
        f.write('    struct check\n');
        f.write('    {\n');
        f.write('      static_assert(sizeof(int) <= sizeof(T), \"not big enough\");\n');
        f.write('    };\n');
        f.write('  }\n');
        f.write('  namespace test_final_override\n');
        f.write('  {\n');
        f.write('    struct Base\n');
        f.write('    {\n');
        f.write('      virtual ~Base() {}\n');
        f.write('      virtual void f() {}\n');
        f.write('    };\n');
        f.write('    struct Derived : public Base\n');
        f.write('    {\n');
        f.write('      virtual ~Derived() override {}\n');
        f.write('      virtual void f() override {}\n');
        f.write('    };\n');
        f.write('\n');
        f.write('  }\n');
        f.write('  namespace test_double_right_angle_brackets\n');
        f.write('  {\n');
        f.write('    template < typename T >\n');
        f.write('    struct check {};\n');
        f.write('    typedef check<void> single_type;\n');
        f.write('    typedef check<check<void>> double_type;\n');
        f.write('    typedef check<check<check<void>>> triple_type;\n');
        f.write('    typedef check<check<check<check<void>>>> quadruple_type;\n');
        f.write('\n');
        f.write('  }\n');
        f.write('  namespace test_decltype\n');
        f.write('  {\n');
        f.write('    int\n');
        f.write('    f()\n');
        f.write('    {\n');
        f.write('      int a = 1;\n');
        f.write('      decltype(a) b = 2;\n');
        f.write('      return a + b;\n');
        f.write('    }\n');
        f.write('  }\n');
        f.write('  namespace test_type_deduction\n');
        f.write('  {\n');
        f.write('    template < typename T1, typename T2 >\n');
        f.write('    struct is_same\n');
        f.write('    {\n');
        f.write('      static const bool value = false;\n');
        f.write('    };\n');
        f.write('\n');
        f.write('    template < typename T >\n');
        f.write('    struct is_same<T, T>\n');
        f.write('    {\n');
        f.write('      static const bool value = true;\n');
        f.write('    };\n');
        f.write('    template < typename T1, typename T2 >\n');
        f.write('    auto\n');
        f.write('    add(T1 a1, T2 a2) -> decltype(a1 + a2)\n');
        f.write('    {\n');
        f.write('      return a1 + a2;\n');
        f.write('    }\n');
        f.write('    int\n');
        f.write('    test(const int c, volatile int v)\n');
        f.write('    {\n');
        f.write('      static_assert(is_same<int, decltype(0)>::value == true, \"\");\n');
        f.write('      static_assert(is_same<int, decltype(c)>::value == false, \"\");\n');
        f.write('      static_assert(is_same<int, decltype(v)>::value == false, \"\");\n');
        f.write('      auto ac = c;\n');
        f.write('      auto av = v;\n');
        f.write('      auto sumi = ac + av + \'x\';\n');
        f.write('      auto sumf = ac + av + 1.0;\n');
        f.write('      static_assert(is_same<int, decltype(ac)>::value == true, \"\");\n');
        f.write('      static_assert(is_same<int, decltype(av)>::value == true, \"\");\n');
        f.write('      static_assert(is_same<int, decltype(sumi)>::value == true, \"\");\n');
        f.write('      static_assert(is_same<int, decltype(sumf)>::value == false, \"\");\n');
        f.write('      static_assert(is_same<int, decltype(add(c, v))>::value == true, \"\");\n');
        f.write('      return (sumf > 0.0) ? sumi : add(c, v);\n');
        f.write('    }\n');
        f.write('\n');
        f.write('  }\n');
        f.write('  namespace test_noexcept\n');
        f.write('  {\n');
        f.write('\n');
        f.write('    int f() { return 0; }\n');
        f.write('    int g() noexcept { return 0; }\n');
        f.write('\n');
        f.write('    static_assert(noexcept(f()) == false, \"\");\n');
        f.write('    static_assert(noexcept(g()) == true, \"\");\n');
        f.write('\n');
        f.write('  }\n');
        f.write('  namespace test_constexpr\n');
        f.write('  {\n');
        f.write('    template < typename CharT >\n');
        f.write('    unsigned long constexpr\n');
        f.write('    strlen_c_r(const CharT *const s, const unsigned long acc) noexcept\n');
        f.write('    {\n');
        f.write('      return *s ? strlen_c_r(s + 1, acc + 1) : acc;\n');
        f.write('    }\n');
        f.write('    template < typename CharT >\n');
        f.write('    unsigned long constexpr\n');
        f.write('    strlen_c(const CharT *const s) noexcept\n');
        f.write('    {\n');
        f.write('      return strlen_c_r(s, 0UL);\n');
        f.write('    }\n');
        f.write('    static_assert(strlen_c(\"\") == 0UL, \"\");\n');
        f.write('    static_assert(strlen_c(\"1\") == 1UL, \"\");\n');
        f.write('    static_assert(strlen_c(\"example\") == 7UL, \"\");\n');
        f.write('    static_assert(strlen_c(\"another\\0example\") == 7UL, \"\");\n');
        f.write('  }\n');
        f.write('  namespace test_rvalue_references\n');
        f.write('  {\n');
        f.write('    template < int N >\n');
        f.write('    struct answer\n');
        f.write('    {\n');
        f.write('      static constexpr int value = N;\n');
        f.write('    };\n');
        f.write('    answer<1> f(int&)       { return answer<1>(); }\n');
        f.write('    answer<2> f(const int&) { return answer<2>(); }\n');
        f.write('    answer<3> f(int&&)      { return answer<3>(); }\n');
        f.write('    void\n');
        f.write('    test()\n');
        f.write('    {\n');
        f.write('      int i = 0;\n');
        f.write('      const int c = 0;\n');
        f.write('      static_assert(decltype(f(i))::value == 1, \"\");\n');
        f.write('      static_assert(decltype(f(c))::value == 2, \"\");\n');
        f.write('      static_assert(decltype(f(0))::value == 3, \"\");\n');
        f.write('    }\n');
        f.write('\n');
        f.write('  }\n');
        f.write('  namespace test_uniform_initialization\n');
        f.write('  {\n');
        f.write('    struct test\n');
        f.write('    {\n');
        f.write('      static const int zero {};\n');
        f.write('      static const int one {1};\n');
        f.write('    };\n');
        f.write('    static_assert(test::zero == 0, \"\");\n');
        f.write('    static_assert(test::one == 1, \"\");\n');
        f.write('  }\n');
        f.write('  namespace test_lambdas\n');
        f.write('  {\n');
        f.write('    void\n');
        f.write('    test1()\n');
        f.write('    {\n');
        f.write('      auto lambda1 = [](){};\n');
        f.write('      auto lambda2 = lambda1;\n');
        f.write('      lambda1();\n');
        f.write('      lambda2();\n');
        f.write('    }\n');
        f.write('    int\n');
        f.write('    test2()\n');
        f.write('    {\n');
        f.write('      auto a = [](int i, int j){ return i + j; }(1, 2);\n');
        f.write('      auto b = []() -> int { return \'0\'; }();\n');
        f.write('      auto c = [=](){ return a + b; }();\n');
        f.write('      auto d = [&](){ return c; }();\n');
        f.write('      auto e = [a, &b](int x) mutable {\n');
        f.write('        const auto identity = [](int y){ return y; };\n');
        f.write('        for (auto i = 0; i < a; ++i)\n');
        f.write('          a += b--;\n');
        f.write('        return x + identity(a + b);\n');
        f.write('      }(0);\n');
        f.write('      return a + b + c + d + e;\n');
        f.write('    }\n');
        f.write('    int\n');
        f.write('    test3()\n');
        f.write('    {\n');
        f.write('      const auto nullary = [](){ return 0; };\n');
        f.write('      const auto unary = [](int x){ return x; };\n');
        f.write('      using nullary_t = decltype(nullary);\n');
        f.write('      using unary_t = decltype(unary);\n');
        f.write('      const auto higher1st = [](nullary_t f){ return f(); };\n');
        f.write('      const auto higher2nd = [unary](nullary_t f1){\n');
        f.write('        return [unary, f1](unary_t f2){ return f2(unary(f1())); };\n');
        f.write('      };\n');
        f.write('      return higher1st(nullary) + higher2nd(nullary)(unary);\n');
        f.write('    }\n');
        f.write('\n');
        f.write('  }\n');
        f.write('  namespace test_variadic_templates\n');
        f.write('  {\n');
        f.write('    template <int...>\n');
        f.write('    struct sum;\n');
        f.write('    template <int N0, int... N1toN>\n');
        f.write('    struct sum<N0, N1toN...>\n');
        f.write('    {\n');
        f.write('      static constexpr auto value = N0 + sum<N1toN...>::value;\n');
        f.write('    };\n');
        f.write('    template <>\n');
        f.write('    struct sum<>\n');
        f.write('    {\n');
        f.write('      static constexpr auto value = 0;\n');
        f.write('    };\n');
        f.write('    static_assert(sum<>::value == 0, \"\");\n');
        f.write('    static_assert(sum<1>::value == 1, \"\");\n');
        f.write('    static_assert(sum<23>::value == 23, \"\");\n');
        f.write('    static_assert(sum<1, 2>::value == 3, \"\");\n');
        f.write('    static_assert(sum<5, 5, 11>::value == 21, \"\");\n');
        f.write('    static_assert(sum<2, 3, 5, 7, 11, 13>::value == 41, \"\");\n');
        f.write('  }\n');
        f.write('  namespace test_template_alias_sfinae\n');
        f.write('  {\n');
        f.write('    struct foo {};\n');
        f.write('    template<typename T>\n');
        f.write('    using member = typename T::member_type;\n');
        f.write('    template<typename T>\n');
        f.write('    void func(...) {}\n');
        f.write('    template<typename T>\n');
        f.write('    void func(member<T>*) {}\n');
        f.write('    void test();\n');
        f.write('    void test() { func<foo>(0); }\n');
        f.write('  }\n');
        f.write('  template<class T>\n');
        f.write('  class class_test_spc\n');
        f.write('  {\n');
        f.write('  public:\n');
        f.write('    int test_func(){return 0;}\n');
        f.write('  };\n');
        f.write('}\n');
        f.write('using namespace cxx11;\n');
        f.write('template<>\n');
        f.write('int class_test_spc<int>::test_func()\n');
        f.write('{return 1;}\n');
        f.write('#endif\n');
    finally:
        f.close();

    from distutils.sysconfig import customize_compiler
    from distutils.ccompiler import new_compiler
    compiler=new_compiler();
    customize_compiler(compiler);

    success=True;
    devnull=open(os.devnull,'w')
    oldstderr=os.dup(sys.stderr.fileno())
    os.dup2(devnull.fileno(),sys.stderr.fileno())
    try:
            obj_files=compiler.compile([fname],extra_postargs=['-std=c++11']);
    except:
        pass;
        success=False;

    os.dup2(oldstderr, sys.stderr.fileno())
    devnull.close()
    return success;


def find_mpi_cxx(cxx_11):
    from os import environ
    mpi_cxx_lst=['mpic++','mpicxx','mpiCC','hcp','mpxlC_r','mpxlC','mpCC','cmpic++'];
    mpi_cxx=None;
    if cxx_11==False:
        if environ.get('MPICXX') is None:
            i=0;
            mpi_cxx=None;
            from distutils.spawn import find_executable
            found_any=False;
            while i<mpi_cxx_lst.__len__() and mpi_cxx==None:
                if find_executable(mpi_cxx_lst[i]) is not None:
                    found_any=True;
                    mpi_cxx=mpi_cxx_lst[i];
                    environ['CC']=mpi_cxx;
                    if has_cxx_11()==False:
                        mpi_cxx=None;
                i+=1;
            if mpi_cxx==None and found_any==True:
                raise Exception('default compiler does not supprt c++11, neither does any of c++ mpi compilers in path');
            if mpi_cxx==None and found_any==False:
                raise Exception('default compiler does not supprt c++11, and could not find any c++ mpi compilers in path');
        else:
            environ['CC']=environ.get('MPICXX');
            if has_cxx_11()==False:
                raise Exception('default compiler does not supprt c++11, neither does any of provided c++ mpi compiler (MPICXX)');
            mpi_cxx=environ.get('MPICXX');
    else:
        if environ.get('MPICXX') is None:
            i=0;
            mpi_cxx=None;
            from distutils.spawn import find_executable
            while i<mpi_cxx_lst.__len__() and mpi_cxx==None:
                if find_executable(mpi_cxx_lst[i]) is not None:
                    mpi_cxx=mpi_cxx_lst[i];
                i+=1;
            if mpi_cxx==None:
                raise Exception('could not find any c++ mpi compilers in path')
        else:
            mpi_cxx=environ.get('MPICXX');
    return mpi_cxx;

def mpi_cxx_params(mpi_cxx,ext):
    line='';
    cmd=mpi_cxx+' -show';
    from subprocess import check_output
    import sys
    try:
        if sys.version_info < (3,0):
            line=check_output(cmd.split());
        else:
            line=check_output(cmd.split()).decode('utf-8');
    except:
        pass;
        err_msg='could not excute \"'+cmd+'\"';
        raise Exception(err_msg)

    from distutils.util import split_quoted

    words=split_quoted(line);
    append_next_word = None
    for word in words[1:]:
        if append_next_word is not None:
            append_next_word.append(word)
            append_next_word = None
            continue

        suffix = os.path.splitext(word)[1]
        switch = word[0:2] ; value = word[2:]

        if suffix in (".c", ".cc", ".cpp", ".cxx", ".c++", ".m", ".mm"):
            # hmm, should we do something about C vs. C++ sources?
            # or leave it up to the CCompiler implementation to
            # worry about?
            ext.sources.append(word)
        elif switch == "-I":
            ext.include_dirs.append(value)
        elif switch == "-D":
            equals = value.find("=")
            if equals == -1:        # bare "-DFOO" -- no value
                ext.define_macros.append((value, None))
            else:                   # "-DFOO=blah"
                ext.define_macros.append((value[0:equals],
                                          value[equals+2:]))
        elif switch == "-U":
            ext.undef_macros.append(value)
        elif switch == "-C":        # only here 'cause makesetup has it!
            ext.extra_compile_args.append(word)
        elif switch == "-l":
            ext.libraries.append(value)
        elif switch == "-L":
            ext.library_dirs.append(value)
        elif switch == "-R":
            ext.runtime_library_dirs.append(value)
        elif word == "-rpath":
            append_next_word = ext.runtime_library_dirs
        elif word == "-Xlinker":
            ext.extra_link_args.append(word)
            append_next_word = ext.extra_link_args
        elif word.startswith("-Wl,"):
            ext.extra_link_args.append(word)
        elif word == "-Xcompiler":
            append_next_word = ext.extra_compile_args
        elif switch == "-u":
            ext.extra_link_args.append(word)
            if not value:
                append_next_word = ext.extra_link_args
        elif suffix in (".a", ".so", ".sl", ".o", ".dylib"):
            # NB. a really faithful emulation of makesetup would
            # append a .o file to extra_objects only if it
            # had a slash in it; otherwise, it would s/.o/.c/
            # and append it to sources.  Hmmmm.
            ext.extra_objects.append(word)
        else:
            ext.extra_compile_args.append(word)
def get_libmpi_so(mpi_cxx):
    f = open('helloworld.cpp','w');
    
    f.write('#include <stdio.h>\n');
    f.write('#include <cstdlib>\n');
    f.write('#include <iostream>\n');
    f.write('#include <mpi.h>\n');
    f.write('int main(int nargs,char* args[])\n');
    f.write('{\n');
    f.write('    MPI_Init(&nargs,&args);\n');
    f.write('    int rank,size;\n');
    f.write('    MPI_Comm_rank(MPI_COMM_WORLD,&rank);\n');
    f.write('    MPI_Comm_size(MPI_COMM_WORLD,&size);\n');
    f.write('    printf(\"Hello world I am proc %d out of %d\\n\",rank,size);\n');
    f.write('    MPI_Finalize();\n');
    f.write('    return EXIT_SUCCESS;\n');
    f.write('}\n');
    f.close();
    

    
    returncode=0;
    cmd=mpi_cxx+' helloworld.cpp -o helloworld';
    try:
        returncode=subprocess.call(cmd.split());
    except:
        pass;
        raise Exception('could not excute \"'+cmd+'\"');
    if returncode!=0:
        raise Exception('could not excute \"'+cmd+'\"');

    if platform=='darwin':
        cmd='otool -L helloworld';
        ptrn='libmpi.*.dylib'
        regex_ptrn='libmpi.[0-9]*.dylib$'
    else:
        cmd='ldd helloworld';
        ptrn='libmpi.so.*'
        regex_ptrn='libmpi.so.[0-9]*$'


    try:
        if sys.version_info < (3,0):
            result=subprocess.check_output(cmd.split());
        else:
            result=subprocess.check_output(cmd.split()).decode('utf-8');
    except:
        pass;
        raise Exception('could not excute \"'+cmd+'\"')

    def find_ptrn(result,ptrn):
        for line in result.split():
            if re.search(ptrn,line):
                return line;
        return None;

    so_file=find_ptrn(result,regex_ptrn);
    if(so_file==None):
        raise Exception('could not find file \"'+ptrn+'\"');

    cmd='rm -f helloworld.cpp helloworld';
    try:
        returncode=subprocess.call(cmd.split());
    except:
        pass;
        raise Exception('could not excute \"'+cmd+'\"');
    if returncode!=0:
        raise Exception('could not excute \"'+cmd+'\"');
    return so_file.split('/').pop()
def numpy_params(ext):
    try:
        from numpy import __path__ as npy_path
    except:
        pass;
        raise Exception('mapp4py requires numpy');
    
    ext.libraries.append('npymath')
    ext.library_dirs.append(npy_path[0]+'/core/lib')
    ext.include_dirs.append(npy_path[0]+'/core/include');


cxx_11=has_cxx_11()
mpi_cxx='';
try:
    mpi_cxx=find_mpi_cxx(cxx_11);
except:
    raise;




import os
cpp_files=[]
cpp_files+=['src/'+ each for each in os.listdir('src') if each.endswith('.cpp')]
h_files=[]
h_files+=['src/'+ each for each in os.listdir('src') if each.endswith('.h')]

mapp_ext=Extension('mapp4py',sources=cpp_files);
mapp_ext.extra_compile_args.append('-std=c++11');
numpy_params(mapp_ext);
if cxx_11==True:
    mpi_cxx_params(mpi_cxx,mapp_ext);


setup(name ='mapp4py',
      version = '0.0.0',
      headers = h_files,
      description = 'MIT Atomistic Parallel Package',
      author = 'Sina Moeini',
      author_email = 'sinam@mit.edu',
      url = 'https://github.com/sinamoeini/mapp4py',
      ext_modules = [mapp_ext])







