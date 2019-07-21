# ===========================================================================
#          https://www.gnu.org/software/autoconf-archive/ax_mpi.html
# ===========================================================================
#
# SYNOPSIS
#
#   AX_MPI([ACTION-IF-FOUND[, ACTION-IF-NOT-FOUND]])
#
# DESCRIPTION
#
#   This macro tries to find out how to compile programs that use MPI
#   (Message Passing Interface), a standard API for parallel process
#   communication (see http://www-unix.mcs.anl.gov/mpi/)
#
#   On success, it sets the MPICC, MPICXX, MPIF77, or MPIFC output variable
#   to the name of the MPI compiler, depending upon the current language.
#   (This may just be $CC/$CXX/$F77/$FC, but is more often something like
#   mpicc/mpiCC/mpif77/mpif90.) It also sets MPILIBS to any libraries that
#   are needed for linking MPI (e.g. -lmpi or -lfmpi, if a special
#   MPICC/MPICXX/MPIF77/MPIFC was not found).
#
#   Note that this macro should be used only if you just have a few source
#   files that need to be compiled using MPI. In particular, you should
#   neither overwrite CC/CXX/F77/FC with the values of
#   MPICC/MPICXX/MPIF77/MPIFC, nor assume that you can use the same flags
#   etc. as the standard compilers. If you want to compile a whole program
#   using the MPI compiler commands, use one of the macros
#   AX_PROG_{CC,CXX,FC}_MPI.
#
#   ACTION-IF-FOUND is a list of shell commands to run if an MPI library is
#   found, and ACTION-IF-NOT-FOUND is a list of commands to run if it is not
#   found. If ACTION-IF-FOUND is not specified, the default action will
#   define HAVE_MPI.
#
# LICENSE
#
#   Copyright (c) 2008 Steven G. Johnson <stevenj@alum.mit.edu>
#   Copyright (c) 2008 Julian C. Cummings <cummings@cacr.caltech.edu>
#
#   This program is free software: you can redistribute it and/or modify it
#   under the terms of the GNU General Public License as published by the
#   Free Software Foundation, either version 3 of the License, or (at your
#   option) any later version.
#
#   This program is distributed in the hope that it will be useful, but
#   WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
#   Public License for more details.
#
#   You should have received a copy of the GNU General Public License along
#   with this program. If not, see <https://www.gnu.org/licenses/>.
#
#   As a special exception, the respective Autoconf Macro's copyright owner
#   gives unlimited permission to copy, distribute and modify the configure
#   scripts that are the output of Autoconf when processing the Macro. You
#   need not follow the terms of the GNU General Public License when using
#   or distributing such scripts, even though portions of the text of the
#   Macro appear in them. The GNU General Public License (GPL) does govern
#   all other use of the material that constitutes the Autoconf Macro.
#
#   This special exception to the GPL applies to versions of the Autoconf
#   Macro released by the Autoconf Archive. When you make and distribute a
#   modified version of the Autoconf Macro, you may extend this special
#   exception to the GPL to apply to your modified version as well.

#serial 978

AC_DEFUN([AX_LIBMPI_SONAME_FUNC], [
AC_PREREQ(2.50) dnl for AC_LANG_CASE



AC_MSG_CHECKING(whether dlopen is needed for mpi)
AC_COMPILE_IFELSE([AC_LANG_SOURCE([_AX_LIBMPI_SONAME_FUNC_test])],
[
    AC_MSG_RESULT(no)
],
[
    AC_MSG_RESULT(yes)
    ax_libmpi_dlopen=yes








AC_CANONICAL_HOST
AC_PROG_GREP
AC_PROG_AWK
AC_CHECK_PROGS(AX_CAT,cat,AC_MSG_ERROR([could not find cat]))

AS_CASE([$host_os],
[darwin*],
AC_CHECK_PROGS(AX_OTOOL, otool otool64,AC_MSG_ERROR([could not find otool or otool64])),

[linux*],
AC_CHECK_PROGS(AX_OBJDUMP, objdump,AC_MSG_ERROR([could not find objdump]))

)




AC_MSG_CHECKING(for path to mpi library)
AX_LIBMPI_SONAME=

$AX_CAT > conftest.cpp <<EOF
#include <cstdlib>
#include <mpi.h>
int main(int ,char* [[]])
{
    MPI_Init(NULL,NULL);
    MPI_Finalize();
    return EXIT_SUCCESS;
}
EOF

$CXX $LIBS -o conftest conftest.cpp

AS_CASE([$host_os],
[darwin*],
[AX_LIBMPI_SONAME=`$AX_OTOOL -L conftest | $AWK '{if (NR > 1) {print [$]1}}' | $GREP libmpi.*.dylib | $AWK  '{if(NR==1){print [$]1; exit}}'`],



[linux*],
[AX_LIBMPI_SONAME=`$AX_OBJDUMP -p conftest | $GREP NEEDED | $GREP libmpi.so.* |  $AWK  '{if(NR==1){print [$]NF; exit}}'`]


)
rm -rf conftest conftest.cpp


    if test x = x"$AX_LIBMPI_SONAME"
    then
        AC_MSG_RESULT([no])
        AC_MSG_ERROR([could not find soname of mpi library])
    else
        AC_MSG_RESULT($AX_LIBMPI_SONAME)
    fi

    AC_DEFINE_UNQUOTED(LIBMPI_SONAME,["$AX_LIBMPI_SONAME"],[dec])

])




])



m4_define([_AX_LIBMPI_SONAME_FUNC_test], [[

#include <mpi.h>
#ifdef OMPI_MAJOR_VERSION
#if OMPI_MAJOR_VERSION < 3 || OMPI_MAJOR_VERSION >= 10
#error "needs dlopen"
#endif
#endif

]])

