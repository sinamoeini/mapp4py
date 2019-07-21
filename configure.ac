AC_INIT([mapp4py],[0.0.0])
AC_PREREQ([2.68])



AC_CONFIG_MACRO_DIR([m4])

AC_CONFIG_SRCDIR([src/main.cpp])
AC_CONFIG_HEADERS([config.h])

AC_CONFIG_AUX_DIR([build-aux])


AC_CONFIG_FILES([Makefile])

AC_PROG_CXX

AC_LANG([C++])

AX_MPICXX([],[AC_MSG_ERROR(could not find mpi library)])
CXX=$MPICXX
LIBS=$MPILIBS" "$LIBS
AX_LIBMPI_SONAME_FUNC

AX_CXX_COMPILE_STDCXX_11(noext,mandatory)


AX_PYTHON([AC_MSG_ERROR(could not find a bonafide version of python with numpy)])


bindir=$PYTHON_SITE
CXXFLAGS="-fPIC "$CXXFLAGS
CPPFLAGS=-I$PYTHON_INC" "-I$NPY_INCLUDE" "$CPPFLAGS
LDFLAGS=-L$NPY_LIB_DIR" "$LDFLAGS
LIBS="-lnpymath "$LIBS


AM_INIT_AUTOMAKE([1.11 -Wall -Wno-extra-portability foreign subdir-objects])
AM_MAINTAINER_MODE([enable])

#LT_PREREQ([2.2])
#LT_INIT([dlopen shared])




AC_OUTPUT
