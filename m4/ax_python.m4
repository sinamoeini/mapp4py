# ===========================================================================
#        https://www.gnu.org/software/autoconf-archive/ax_python.html
# ===========================================================================
#
# SYNOPSIS
#
#   AX_PYTHON
#
# DESCRIPTION
#
#   This macro does a complete Python development environment check.
#
#   It checks for all known versions. When it finds an executable, it looks
#   to find the header files and library.
#
#   It sets PYTHON_BIN to the name of the python executable,
#   PYTHON_INCLUDE_DIR to the directory holding the header files, and
#   PYTHON_LIB to the name of the Python library.
#
#   This macro calls AC_SUBST on PYTHON_BIN (via AC_CHECK_PROG),
#   PYTHON_INCLUDE_DIR and PYTHON_LIB.
#
# LICENSE
#
#   Copyright (c) 2008 Michael Tindal
#
#   This program is free software; you can redistribute it and/or modify it
#   under the terms of the GNU General Public License as published by the
#   Free Software Foundation; either version 2 of the License, or (at your
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

#serial 18

AC_DEFUN([AX_PYTHON],
[AC_MSG_CHECKING(for python build information)
AC_MSG_RESULT([])
AC_ARG_VAR(PYTHON,[path to python executable])
AC_ARG_VAR(PYTHON_SITE,[path to site-packages directory of python])

ax_python_bin_list="python3 python"

if test -z $PYTHON;
then
    __ax_python=
else
    ax_python_bin_list=`basename $PYTHON`" "$ax_python_bin_list
    __ax_python=$PYTHON
fi

for python in $ax_python_bin_list
do
    AC_CHECK_PROGS(__ax_python, [$python])
    ax_python=$__ax_python
    __ax_python=

    if test x$ax_python != x
    then
        ax_python_inc=`$ax_python -c "from distutils.sysconfig import get_config_var; print(get_config_var('CONFINCLUDEPY'))"`
        ax_python_site_packages=
        if test -z $PYTHON_SITE
        then
            ax_python_site_packages=`$ax_python -c "[import site; print(site.getsitepackages()[0]])"`
        else
            AC_MSG_CHECKING(whether provided site-packages is consistent with $python)
            ax_python_site_packages_test=`$ax_python -c "[import site; print('$PYTHON_SITE' in site.getsitepackages()])"`

            if test x$ax_python_site_packages_test = xTrue
            then
                AC_MSG_RESULT(yes)
                ax_python_site_packages=$PYTHON_SITE
            else
                AC_MSG_RESULT(no)
            fi
        fi
        if test x$ax_python_site_packages != x
        then
            AC_MSG_CHECKING(for numpy)
            $ax_python -c "import numpy" 2>/dev/null
            if test $? -eq 0;
            then
                AC_MSG_RESULT(yes)
                ax_npy_core=`$ax_python -c "[from numpy import __path__ as npy_path; print(npy_path[0]+'/core');]"`;
                ax_saved_ldflags=$LDFLAGS
                LDFLAGS=-L"$ax_npy_core""/lib"
                ax_npy_lib=no
                AC_CHECK_LIB(npymath,npy_log,[ax_npy_lib=$ax_npy_core"/lib";])
                LDFLAGS=$ax_saved_ldflags
                if test x$ax_npy_lib != xno
                then
                    ax_saved_cppflags=$CPPFLAGS
                    CPPFLAGS=-I"$ax_npy_core""/include"
                    ax_npy_inc=
                    AC_CHECK_HEADER([numpy/numpyconfig.h],[ax_npy_inc=$ax_npy_core"/include";])
                    CPPFLAGS=$ax_saved_cppflags
                    if test x$ax_npy_inc != x
                    then
                        break;
                    fi
                fi
            fi

        fi
    fi
done


if test x$ax_npy_inc = x
then
    $1
    :
else
    PYTHON=$ax_python;
    PYTHON_SITE=$ax_python_site_packages;
    PYTHON_INC=$ax_python_inc;
    NPY_LIB_DIR=$ax_npy_lib;
    NPY_INCLUDE=$ax_npy_inc;
fi
])dnl
