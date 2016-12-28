#                                               -*- Autoconf -*-
#--------------------------------------------------------------------
# Process this file with autoconf to produce a configure script.
#--------------------------------------------------------------------

AC_PREREQ([2.63])
AC_INIT([Structured], [1.0], [GRAD510], [Structured.tar.bz2])
AC_CONFIG_SRCDIR([lib/src/BoxFramework/IntVect.cpp])

#--------------------------------------------------------------------
# Created files
#--------------------------------------------------------------------

AC_CONFIG_HEADERS([lib/src/BoxFramework/Config.H])
AC_CONFIG_FILES([Common/mk/Make.defs])

#--------------------------------------------------------------------
# System information
#--------------------------------------------------------------------

AC_CANONICAL_HOST
AC_CHECK_SIZEOF([void *])
if test "$ac_cv_sizeof_void_p" -eq "8"; then
  ch_int_systemptrsize="64"
elif test "$ac_cv_sizeof_void_p" -eq "4"; then
  ch_int_systemptrsize="32"
fi

#--------------------------------------------------------------------
# Top directory
#--------------------------------------------------------------------

#AC_SUBST([TOPDIR],[`pwd | sed '1{s/ /\\ /}'`])
#AC_SUBST([TOPDIR],[/Users/Mike/Desktop/cppPrograms/Grad510/Structured_cappsm/])
AC_SUBST([TOPDIR],[`pwd`])

#--------------------------------------------------------------------
# Check for MPI
#--------------------------------------------------------------------

AC_MSG_CHECKING([for whether to build with MPI])
AC_ARG_WITH([mpi],
   [AC_HELP_STRING([--with-mpi=prefix], [use MPI (prefix to compilers)])],
   [], [with_mpi=no])
if test "$with_mpi" != "no"; then
   use_mpi=yes
   CPPFLAGS="$CPPFLAGS-DUSE_MPI "
   mpicxxname=mpicxx
   if test -z $with_mpi || test "$with_mpi" = "yes"; then
      case "$host_os" in
         linux*) FC=mpif90
                 CC=mpicc
                 CXX=mpicxx ;;
      esac
   else
      FC=$with_mpi/mpif90
      CC=$with_mpi/mpicc
      CXX=$with_mpi/mpicxx
   fi
   case "$host_os" in
         linux*)
# Try to determine the C++ compiler (-V usually gives version.  For g++ it
# reports an error but states its name)
                SCXX=`$CXX -V 2>&1 | sed '1{s/.*Compaq.*/cxx/;s/.*Intel.*/icpc/;s/.*Sun.*/sunCC/;s/.*Power.*/xlc++/;s/.*g++.*/g++/;q;}'`
# Try to determine the C compiler (-V usually gives version.  For gcc it reports
# an error but states its name)
                SCC=`$CC -V 2>&1 | sed '1{s/.*Compaq.*/ccc/;s/.*Intel.*/icc/;s/.*Sun.*/suncc/;s/.*Power.*/xlc/;s/.*gcc.*/gcc/;q;}'` ;;
   esac      
else
   use_mpi=no
fi
AC_MSG_RESULT($use_mpi)

#--------------------------------------------------------------------
# Check for OpenMP
#--------------------------------------------------------------------

PARALLELCXXFLAGS=""
AC_MSG_CHECKING([whether to build with OpenMP])
AC_ARG_ENABLE(openmp,
   AC_HELP_STRING([--enable-openmp], [build with OpenMP (default=no)]),,
   enableval=no)
if test "$enableval" = "yes"; then
   use_openmp=yes
   # we really should check for compiler support   
   PARALLELCXXFLAGS="-fopenmp"
else
   use_openmp=no
fi
AC_MSG_RESULT($use_openmp)

#--------------------------------------------------------------------
# Get the C++ compiler
#--------------------------------------------------------------------

if test "$use_mpi" = "yes"; then
   PCXX=$CXX
   CXX=$SCXX
   AC_MSG_NOTICE([Using underlying C++ compiler $SCXX for configure tests])
fi
AC_PROG_CXX([$CXX])
if test "$use_mpi" = "no"; then
   SCXX=$CXX
fi

# Check for debug
AC_MSG_CHECKING([if release version is enabled])
AC_ARG_ENABLE(release,
   AC_HELP_STRING([--enable-release], [build released version (default=no)]),,
   enableval=no)
if test "$enableval" = "yes"; then
   DEBUG=no
else
   DEBUG=yes
fi
AC_MSG_RESULT($enableval)

# Set debug, optimization, and regular flags for known compilers
AC_MSG_NOTICE([Using C++ compiler: $CXX])
AC_MSG_NOTICE([Using underlying C++ compiler: $SCXX])
case "$SCXX" in
   g++)             DEBUGFLAGS="-g"
                    OPTFLAGS="-O3 -funroll-loops"
                    REGFLAGS="-march=native -std=c++14 -pedantic -Wall -Wno-unused-local-typedefs -Wno-unknown-pragmas -Wno-vla" ;;
   clang++)         DEBUGFLAGS="-g"
                    OPTFLAGS="-O3 -funroll-loops"
                    REGFLAGS="-march=native -std=c++14 -pedantic -Wall -Wno-unused-local-typedefs -Wno-unknown-pragmas -Wno-vla-extension" ;;
   icpc)            DEBUGFLAGS="-g"
                    OPTFLAGS="-O3 -funroll-loops"
                    REGFLAGS="-march=native -std=c++14 -pedantic -Wall -Wno-unknown-pragmas" ;;
esac

# Set compilation flags
if test "$DEBUG" = "yes"; then
   CXXFLAGS="$DEBUGFLAGS $REGFLAGS $PARALLELCXXFLAGS"
else
   CXXFLAGS="$OPTFLAGS $REGFLAGS $PARALLELCXXFLAGS"
   CPPFLAGS="$CPPFLAGS-DRELEASE "
fi

#--------------------------------------------------------------------
# Get the C compiler
#--------------------------------------------------------------------

# Use the serial compiler for these tests
if test "$use_mpi" = "yes"; then
   PCC=$CC
   CC=$SCC
   AC_MSG_NOTICE([Using underlying C compiler $SCC for configure tests])
fi
AC_PROG_CC()
if test "$use_mpi" = "no"; then
   SCC=$CC
fi

#--------------------------------------------------------------------
# Extensions
#--------------------------------------------------------------------

AC_SUBST([LIBEXT],[a])

#--------------------------------------------------------------------
# Program characteristics
#--------------------------------------------------------------------

# Check for single precision
AC_MSG_CHECKING([for precision])
AC_ARG_ENABLE(single_precision,
   AC_HELP_STRING([--enable-single_precision],
   [use single precision (default=no)]),, enableval=no)
if test "$enableval" = "yes"; then
   CPPFLAGS="$CPPFLAGS-DUSE_SINGLE_PRECISION "
   use_precision="single"
else
   use_precision="double"
fi
AC_MSG_RESULT($use_precision)

# Check for space dimension
AC_MSG_CHECKING([for space dim])
AC_ARG_WITH([dim2],
  [AC_HELP_STRING([--with-dim2=prefix], [use dimension 2])],
  [], [with_dim2=no])
if test "$with_dim2" != "no"; then
 CPPFLAGS="$CPPFLAGS-DSPACEDIM2 "
 use_dim="2"
else
  use_dim="3"
fi
AC_MSG_RESULT($use_dim)

#--------------------------------------------------------------------
# Required libraries
#--------------------------------------------------------------------

LIBS="-lm $LIBS"
LDD_SEARCH=""

# Z (may be required for hdf5)

AC_CHECK_LIB([z], [compress])

# CGNS

AC_MSG_CHECKING([if CGNS directory was specified])
AC_ARG_WITH([cgns],
   [AC_HELP_STRING([--with-cgns=CGNSdir],
      [directory containing CGNS installation])],
   [cgnsdir=$withval], [cgnsdir=no])
if test -z "$cgnsdir" || test "$cgnsdir" = "no"; then
   AC_MSG_RESULT(no)
   AC_CHECK_HEADER([cgnslib.h], [havecgnsh="yes"], [havecgnsh="no"])
   if test "$havecgnsh" = "no"; then
   # Search this directory and home directory
      AC_MSG_CHECKING([local, /usr/local/cgns, and home directories for cgnslib.h])
      for d in \
         `ls -dr * 2>/dev/null` \
         `ls -dr /usr/local/cgns/* 2>/dev/null` \
         `ls -d $HOME/local/cgns/* 2>/dev/null` \
         `ls -dr $HOME 2>/dev/null` ; do
         if test -f "$d/cgnslib.h"; then
            CGNSHDIR=`(cd $d; pwd)`
         elif test -f "$d/include/cgnslib.h"; then
            CGNSHDIR=`(cd $d/include; pwd)`
         fi
         if test ! -z "$CGNSHDIR"; then
            break
         fi
      done
      AC_MSG_RESULT($CGNSHDIR)
      if test -z "$CGNSHDIR"; then
         AC_MSG_WARN([cannot find CGNS header:  use --with-cgns=CGNSdir])
      else
         CPPFLAGS="$CPPFLAGS-I$CGNSHDIR "
      fi
   fi
   AC_CHECK_LIB([cgns], [cg_open], [havecgnslib="yes"], [havecgnslib="no"])
   if test "$havecgnslib" = "no"; then
   # First try adding hdf5 and sz to libs (often they are added by libcgns)
      AC_MSG_NOTICE([Adding -lhdf5 and -lsz to LIBS and rechecking for libcgns])
      cgns_str_savelibs="$LIBS"
      LIBS="-lhdf5 -lsz $LIBS"
      AS_UNSET([ac_cv_lib_cgns_cg_open])
      AC_CHECK_LIB([cgns], [cg_open], [havecgnslib="yes"], [havecgnslib="no"])
      if test "$havecgnslib" = "no"; then
         # Revert libs
         LIBS="$cgns_str_savelibs"
         # Second, search this directory and home directory
         AC_MSG_CHECKING([local, /usr/local/cgns, and home directories for libcgns.so])
         for d in \
            `ls -dr * 2>/dev/null` \
            `ls -d /usr/local/cgns/* 2>/dev/null` \
            `ls -d $HOME/local/cgns/* 2>/dev/null` \
            `ls -dr $HOME 2>/dev/null`; do
            if test -f "$d/libcgns.so"; then
               CGNSLIBDIR=`(cd $d; pwd)`
            elif test -f "$d/lib/libcgns.so"; then
               CGNSLIBDIR=`(cd $d/lib; pwd)`
            fi
            if test ! -z "$CGNSLIBDIR"; then
               break
            fi
         done
         AC_MSG_RESULT($CGNSLIBDIR)
         if test -z "$CGNSLIBDIR"; then
            AC_MSG_WARN([cannot find CGNS library:  use --with-cgns=CGNSdir])
         else
            LIBS="-lcgns $LIBS"
            LDFLAGS="$LDFLAGS-L$CGNSLIBDIR "
            LDD_SEARCH="$LDD_SEARCH,-rpath,$CGNSLIBDIR"
         fi
      else
         LIBS="-lcgns $LIBS"
      fi
   else
      LIBS="-lcgns $LIBS"
   fi
elif test "$cgnsdir" = "none"; then
   AC_MSG_RESULT(disabled)
   CPPFLAGS="$CPPFLAGS-DNO_CGNS "
else
   AC_MSG_RESULT($cgnsdir)
   if test -f "$cgnsdir/include/cgnslib.h"; then
      CGNSHDIR=`(cd $cgnsdir/include; pwd)`
      CPPFLAGS="$CPPFLAGS-I$CGNSHDIR "
   else
      AC_MSG_WARN([$cgnsdir/include directory does not contain cgnslib.h])
   fi
   if test -f "$cgnsdir/lib/libcgns.so"; then
      LIBS="-lcgns $LIBS"
      CGNSLIBDIR=`(cd $cgnsdir/lib; pwd)`
      LDFLAGS="$LDFLAGS-L$CGNSLIBDIR "
      LDD_SEARCH="$LDD_SEARCH,-rpath,$CGNSLIBDIR"
   else
      AC_MSG_WARN([$cgnsdir/lib directory does not contain libcgns.so])
   fi
fi

# Add search path to linker flags
if test -z "$LDD_SEARCH"; then
   LDFLAGS="$LDFLAGS$LIBS "
else
   LDFLAGS="$LDFLAGS$LIBS-Wl$LDD_SEARCH "
fi

#--------------------------------------------------------------------
# Standard checks for required programs
#--------------------------------------------------------------------

AC_PROG_INSTALL

#--------------------------------------------------------------------
# Additional programs
#--------------------------------------------------------------------

# Library archiver
AC_ARG_VAR([AR], [Library archiver])
AC_CHECK_PROG([AR], [ar], [ar rcs])

# Stream editor
AC_ARG_VAR([SED], [Stream editor])
AC_CHECK_PROG([SED], [sed], [sed], [])

# Perl
AC_ARG_VAR([PERL], [Perl runtime])
AC_CHECK_PROG([PERL], [perl], [perl], [])

# C shell
AC_ARG_VAR([CSHELLCMD], [C shell])
AC_CHECK_PROG([CSHELLCMD], [csh], [csh -f -c], [])

# Doxygen
AC_ARG_VAR([DOXYGEN], [Doxygen compiler])
AC_CHECK_PROG([DOXYGEN], [doxygen], [doxygen], [])


#--------------------------------------------------------------------
# Install directories
#--------------------------------------------------------------------

if test -z $prefix || test $prefix = NONE; then
  prefix=$ac_default_prefix
fi
if test -z $exec_prefix || test $exec_prefix = NONE; then
  exec_prefix=$prefix
fi

AC_MSG_CHECKING([where to install exectutables])
eval "INSTBINDIR=$bindir"
AC_MSG_RESULT([$INSTBINDIR])

#--------------------------------------------------------------------
# Custom tests use C++
#--------------------------------------------------------------------

AC_LANG([C++])

#--------------------------------------------------------------------
# Check for x86 vector extensions and intrinsics
#--------------------------------------------------------------------

ch_str_veccpu_feature=unknown
ch_str_veccxx_feature=unknown
ch_int_veccpu_featurebits=0
ch_int_veccxx_featurebits=0
ch_str_sse_h=unknown
ch_str_sse2_h=unknown
ch_str_sse3_h=unknown
ch_str_ssse3_h=unknown
ch_str_sse4d1_h=unknown
ch_str_sse4d2_h=unknown
ch_str_sse4a_h=unknown
ch_str_avx_h=unknown
ch_str_fma_h=unknown
ch_str_avx2_h=unknown
ch_str_intelvec_h=unknown
# Vector operators are tested for sse, sse2, avx, and avx2
ch_stat_vecopp=unknown

AC_MSG_CHECKING([for x86 vector extensions])

case $host_cpu in
  *86*)

#--Test hardware

    ch_stat_cpuid=no
    # Avoid using an optimization flags -- the assembly can cause problems.
    # We will also gradually add sse flags
    ch_str_savecxxflags="$CXXFLAGS"
    CXXFLAGS=""
    AC_RUN_IFELSE([AC_LANG_PROGRAM(
[[#include <fstream>
#include <string>]],
[[unsigned bx, cx, dx;
unsigned volatile res = 0;
std::string str("");
__asm__ __volatile__ ("cpuid" : "=c"(cx), "=d"(dx) : "a"(0x1));
if (dx & (1 << 25)) { res |= (1 <<  0); str += "sse"; }
if (dx & (1 << 26)) { res |= (1 <<  1); str += " sse2"; }
if (cx & (1 <<  0)) { res |= (1 <<  2); str += " sse3"; }
if (cx & (1 <<  9)) { res |= (1 <<  3); str += " ssse3"; }
if (cx & (1 << 19)) { res |= (1 <<  4); str += " sse4.1"; }
if (cx & (1 << 20)) { res |= (1 <<  5); str += " sse4.2"; }
if (cx & (1 << 28)) { res |= (1 <<  7); str += " avx1"; }
if (cx & (1 << 12)) { res |= (1 <<  8); str += " fma"; }
__asm__ __volatile__ ("cpuid" : "=c"(cx) : "a"(0x80000001));
if (cx & (1 <<  6)) { res |= (1 <<  6); str += " sse4a"; }
__asm__ __volatile__ ("cpuid" : "=b"(bx) : "a"(0x7), "c"(0x0));
if (bx & (1 <<  5)) { res |= (1 <<  9); str += " avx2"; }
std::ofstream ofile("conftest.data");
ofile << str << " featurebits=" << res;
ofile.close();]])],
                  [ch_stat_cpuid=yes])
    if test "x$ch_stat_cpuid" = "xyes" && test -r "conftest.data"; then
      ch_str_veccpu_feature=`cat conftest.data`
      ch_str_have_cpusse=`expr "x$ch_str_veccpu_feature" : 'x.*\(sse\) '`
      ch_str_have_cpusse2=`expr "x$ch_str_veccpu_feature" : 'x.*\(sse2\) '`
      ch_str_have_cpusse3=`expr "x$ch_str_veccpu_feature" : 'x.* \(sse3\) '`
      ch_str_have_cpussse3=`expr "x$ch_str_veccpu_feature" : 'x.*\(ssse3\) '`
      ch_str_have_cpusse4d1=`expr "x$ch_str_veccpu_feature" : 'x.*\(sse4\.1\) '`
      ch_str_have_cpusse4d2=`expr "x$ch_str_veccpu_feature" : 'x.*\(sse4\.2\) '`
      ch_str_have_cpusse4a=`expr "x$ch_str_veccpu_feature" : 'x.*\(sse4a\) '`
      ch_str_have_cpuavx1=`expr "x$ch_str_veccpu_feature" : 'x.*\(avx1\) '`
      ch_str_have_cpufma=`expr "x$ch_str_veccpu_feature" : 'x.*\(fma\) '`
      ch_str_have_cpuavx2=`expr "x$ch_str_veccpu_feature" : 'x.*\(avx2\) '`
      [ch_int_veccpu_featurebits=`expr "x$ch_str_veccpu_feature" : 'x.*featurebits=\([0-9]*\)'`]
#DBG      echo "$ch_str_have_cpusse x $ch_str_have_cpusse2 x $ch_str_have_cpusse3 x $ch_str_have_cpussse3 x $ch_str_have_cpusse4d1 x $ch_str_have_cpusse4d2 x $ch_str_have_cpusse4a x $ch_str_have_cpuavx1 x $ch_str_have_cpufma x $ch_str_have_cpuavx2 x $ch_int_veccpu_featurebits"
      ch_str_veccpu_feature=`expr "x$ch_str_veccpu_feature" : 'x\(.*\) featurebits=.*'`

#--Test compiler
#--SSE

      case $ch_str_cxxmake in
        portland) CXXFLAGS="-Mvect=sse" ;;
        *       ) CXXFLAGS="-msse" ;;
      esac
      if test "x$ch_str_have_cpusse" = "xsse"; then
        AC_COMPILE_IFELSE([AC_LANG_PROGRAM(
[[#include <xmmintrin.h>  // SSE]],
[[__m128 test_v4sf = _mm_set1_ps(1.f) + _mm_setr_ps(1.f, 2.f, 3.f, 4.f);]])],
                          [ch_str_sse_h="xmmintrin.h"])
        if test "x$ch_str_sse_h" = "xxmmintrin.h"; then
          if test "x$ch_stat_vecopp" != "xno"; then
            ch_stat_vecopp="yes"
          fi
        else
          # Vector ops do not work
          ch_stat_vecopp="no"
          AC_COMPILE_IFELSE([AC_LANG_PROGRAM(
[[#include <xmmintrin.h>  // SSE]],
[[__m128 test_v4sf = _mm_add_ps(_mm_set1_ps(1.f),
                   _mm_setr_ps(1.f, 2.f, 3.f, 4.f));]])],
                            [ch_str_sse_h="xmmintrin.h"])
        fi
        if test "x$ch_str_sse_h" = "xxmmintrin.h"; then
          ch_str_veccxx_feature="sse"
          ch_int_veccxx_featurebits=1
          AC_COMPILE_IFELSE([AC_LANG_PROGRAM(
[[#include <immintrin.h>  // General Vec]],
[[__m128 test_v4sf = _mm_add_ps(_mm_set1_ps(1.f),
                              _mm_setr_ps(1.f, 2.f, 3.f, 4.f));]])],
                            [ch_str_intelvec_h="immintrin.h"],
                            [ch_str_intelvec_h="unknown"])
        fi
      fi
#DBG      echo "Values: $ch_str_veccxx_feature $ch_int_veccxx_featurebits $ch_str_sse_h $ch_str_intelvec_h"

#--SSE2

      case $ch_str_cxxmake in
        portland) CXXFLAGS="-Mvect=sse" ;;
        *       ) CXXFLAGS="-msse2" ;;
      esac
      if test "x$ch_str_have_cpusse2" = "xsse2"; then
        AC_COMPILE_IFELSE([AC_LANG_PROGRAM(
[[#include <emmintrin.h>  // SSE2]],
[[__m128d test_v2df = _mm_set1_pd(1.) + _mm_setr_pd(1., 2.);]])],
                          [ch_str_sse2_h="emmintrin.h"])
        if test "x$ch_str_sse2_h" = "xemmintrin.h"; then
          if test "x$ch_stat_vecopp" != "xno"; then
            ch_stat_vecopp="yes"
          fi
        else
          # Vector ops do not work
          ch_stat_vecopp="no"
          AC_COMPILE_IFELSE([AC_LANG_PROGRAM(
[[#include <emmintrin.h>  // SSE2]],
[[__m128d test_v2df = _mm_add_pd(_mm_set1_pd(1.), _mm_setr_pd(1., 2.));]])],
                            [ch_str_sse2_h="emmintrin.h"])
        fi
        if test "x$ch_str_sse2_h" = "xemmintrin.h"; then
          ch_str_veccxx_feature="$ch_str_veccxx_feature sse2"
          ch_int_veccxx_featurebits=`expr "$ch_int_veccxx_featurebits" + '2'`
          AC_COMPILE_IFELSE([AC_LANG_PROGRAM(
[[#include <immintrin.h>  // General Vec]],
[[__m128d test_v2df = _mm_add_pd(_mm_set1_pd(1.), _mm_setr_pd(1., 2.));]])],
                            [ch_str_intelvec_h="immintrin.h"],
                            [ch_str_intelvec_h="unknown"])
        fi
      fi
#DBG      echo "Values: $ch_str_veccxx_feature $ch_int_veccxx_featurebits $ch_str_sse2_h $ch_str_intelvec_h"

#--SSE3

      case $ch_str_cxxmake in
        portland) CXXFLAGS="-Mvect=sse" ;;
        *       ) CXXFLAGS="-msse3" ;;
      esac
      if test "x$ch_str_have_cpusse3" = "xsse3"; then
        AC_COMPILE_IFELSE([AC_LANG_PROGRAM(
[[#include <pmmintrin.h>  // SSE3]],
[[__m128 test_v4sf = _mm_hadd_ps(_mm_set1_ps(1.f),
                               _mm_setr_ps(1.f, 2.f, 3.f, 4.f));]])],
                          [ch_str_sse3_h="pmmintrin.h"])
        if test "x$ch_str_sse3_h" = "xpmmintrin.h"; then
          ch_str_veccxx_feature="$ch_str_veccxx_feature sse3"
          ch_int_veccxx_featurebits=`expr "$ch_int_veccxx_featurebits" + '4'`
          AC_COMPILE_IFELSE([AC_LANG_PROGRAM(
[[#include <immintrin.h>  // General Vec]],
[[__m128 test_v4sf = _mm_hadd_ps(_mm_set1_ps(1.f),
                               _mm_setr_ps(1.f, 2.f, 3.f, 4.f));]])],
                            [ch_str_intelvec_h="immintrin.h"],
                            [ch_str_intelvec_h="unknown"])
        fi
      fi
#DBG      echo "Values: $ch_str_veccxx_feature $ch_int_veccxx_featurebits $ch_str_sse3_h $ch_str_intelvec_h"

#--SSSE3

      case $ch_str_cxxmake in
        portland) CXXFLAGS="-Mvect=sse" ;;
        *       ) CXXFLAGS="-mssse3" ;;
      esac
      if test "x$ch_str_have_cpussse3" = "xssse3"; then
        AC_COMPILE_IFELSE([AC_LANG_PROGRAM(
[[#include <tmmintrin.h>  // SSSE3]],
[[__m128i test_v4si = _mm_sign_epi32(_mm_set1_epi32(2),
                                   _mm_setr_epi32(1, -1, 1, -1));]])],
                          [ch_str_ssse3_h="tmmintrin.h"])
        if test "x$ch_str_ssse3_h" = "xtmmintrin.h"; then
          ch_str_veccxx_feature="$ch_str_veccxx_feature ssse3"
          ch_int_veccxx_featurebits=`expr "$ch_int_veccxx_featurebits" + '8'`
          AC_COMPILE_IFELSE([AC_LANG_PROGRAM(
[[#include <immintrin.h>  // General Vec]],
[[__m128i test_v4si = _mm_sign_epi32(_mm_set1_epi32(2),
                                   _mm_setr_epi32(1, -1, 1, -1));]])],
                            [ch_str_intelvec_h="immintrin.h"],
                            [ch_str_intelvec_h="unknown"])
        fi
      fi
#DBG      echo "Values: $ch_str_veccxx_feature $ch_int_veccxx_featurebits $ch_str_ssse3_h $ch_str_intelvec_h"

#--SSE4.1

      case $ch_str_cxxmake in
        portland) CXXFLAGS="-Mvect=sse" ;;
        *       ) CXXFLAGS="-msse4.1" ;;
      esac
      if test "x$ch_str_have_cpusse4d1" = "xsse4.1"; then
        AC_COMPILE_IFELSE([AC_LANG_PROGRAM(
[[#include <smmintrin.h>  // SSE4.1]],
[[__m128 test_v4sf = _mm_dp_ps(_mm_set1_ps(1.f),
                             _mm_setr_ps(1.f, 2.f, 3.f, 0.f),
                             113);]])],
                          [ch_str_sse4d1_h="smmintrin.h"])
        if test "x$ch_str_sse4d1_h" = "xsmmintrin.h"; then
          ch_str_veccxx_feature="$ch_str_veccxx_feature sse4.1"
          ch_int_veccxx_featurebits=`expr "$ch_int_veccxx_featurebits" + '16'`
          AC_COMPILE_IFELSE([AC_LANG_PROGRAM(
[[#include <immintrin.h>  // General Vec]],
[[__m128 test_v4sf = _mm_dp_ps(_mm_set1_ps(1.f),
                             _mm_setr_ps(1.f, 2.f, 3.f, 0.f),
                             113);]])],
                            [ch_str_intelvec_h="immintrin.h"],
                            [ch_str_intelvec_h="unknown"])
        fi
      fi
#DBG      echo "Values: $ch_str_veccxx_feature $ch_int_veccxx_featurebits $ch_str_sse4d1_h $ch_str_intelvec_h"

#--SSE4.2

      case $ch_str_cxxmake in
        portland) CXXFLAGS="-Mvect=sse" ;;
        *       ) CXXFLAGS="-msse4.2" ;;
      esac
      if test "x$ch_str_have_cpusse4d2" = "xsse4.2"; then
        AC_COMPILE_IFELSE([AC_LANG_PROGRAM(
[[#include <nmmintrin.h>  // SSE4.2]],
[[__m128i test_v2di = _mm_cmpgt_epi64(_mm_set1_epi32(1),
                                    _mm_setr_epi32(1, 2, 3, 4));]])],
                          [ch_str_sse4d2_h="nmmintrin.h"])
        if test "x$ch_str_sse4d2_h" = "xnmmintrin.h"; then
          ch_str_veccxx_feature="$ch_str_veccxx_feature sse4.2"
          ch_int_veccxx_featurebits=`expr "$ch_int_veccxx_featurebits" + '32'`
          AC_COMPILE_IFELSE([AC_LANG_PROGRAM(
[[#include <immintrin.h>  // General Vec]],
[[__m128i test_v2di = _mm_cmpgt_epi64(_mm_set1_epi32(1),
                                    _mm_setr_epi32(1, 2, 3, 4));]])],
                            [ch_str_intelvec_h="immintrin.h"],
                            [ch_str_intelvec_h="unknown"])
        fi
      fi
#DBG      echo "Values: $ch_str_veccxx_feature $ch_int_veccxx_featurebits $ch_str_sse4d2_h $ch_str_intelvec_h"

#--SSE4a

      case $ch_str_cxxmake in
        # Intel likes to pretend that AMD extensions don't exist :P
        intel   ) CXXFLAGS="" ;;
        portland) CXXFLAGS="-Mvect=sse" ;;
        *       ) CXXFLAGS="-msse4a" ;;
      esac
      if test "x$ch_str_have_cpusse4a" = "xsse4a"; then
        AC_COMPILE_IFELSE([AC_LANG_PROGRAM(
[[#include <ammintrin.h>  // SSE4a]],
[[__m128i test_v2di = _mm_extracti_si64(_mm_set1_epi32(6), 1, 2);]])],
                          [ch_str_sse4a_h="ammintrin.h"])
        if test "x$ch_str_sse4a_h" = "xammintrin.h"; then
          ch_str_veccxx_feature="$ch_str_veccxx_feature sse4a"
          ch_int_veccxx_featurebits=`expr "$ch_int_veccxx_featurebits" + '64'`
          # The general "immintrin.h" doesn't seem to work with SSE4a
        fi
      fi
#DBG      echo "Values: $ch_str_veccxx_feature $ch_int_veccxx_featurebits $ch_str_sse4a_h $ch_str_intelvec_h"

#--AVX

      case $ch_str_cxxmake in
        portland) CXXFLAGS="-Mvect=sse" ;;
        *       ) CXXFLAGS="-mavx" ;;
      esac
      if test "x$ch_str_have_cpuavx1" = "xavx1"; then
        AC_COMPILE_IFELSE([AC_LANG_PROGRAM(
[[#include <immintrin.h>  // AVX]],
[[__m256d test_v4df = _mm256_set1_pd(1.) + _mm256_setr_pd(1., 2., 3., 4.);]])],
                          [ch_str_avx_h="immintrin.h"])
        if test "x$ch_str_avx_h" = "ximmintrin.h"; then
          if test "x$ch_stat_vecopp" != "xno"; then
            ch_stat_vecopp="yes"
          fi
        else
          # Vector ops do not work
          ch_stat_vecopp="no"
          AC_COMPILE_IFELSE([AC_LANG_PROGRAM(
[[#include <immintrin.h>  // AVX]],
[[__m256d test_v4df = _mm256_add_pd(_mm256_set1_pd(1.),
                                  _mm256_setr_pd(1., 2., 3., 4.));]])],
                            [ch_str_avx_h="immintrin.h"])
        fi
        if test "x$ch_str_avx_h" = "ximmintrin.h"; then
          ch_str_veccxx_feature="$ch_str_veccxx_feature avx"
          ch_int_veccxx_featurebits=`expr "$ch_int_veccxx_featurebits" + '128'`
          # Seems we can only use immintrin.h reliably
          ch_str_intelvec_h="immintrin.h"
        fi
      fi
#DBG      echo "Values: $ch_str_veccxx_feature $ch_int_veccxx_featurebits $ch_str_avx_h $ch_str_intelvec_h"

#--FMA

      case $ch_str_cxxmake in
        portland) CXXFLAGS="-Mvect=sse" ;;
        *       ) CXXFLAGS="-mfma" ;;
      esac
      if test "x$ch_str_have_cpufma" = "xfma"; then
        AC_COMPILE_IFELSE([AC_LANG_PROGRAM(
[[#include <immintrin.h>  // FMA]],
[[__m256d test_v4df = _mm256_fmadd_pd(_mm256_set1_pd(2.), _mm256_set1_pd(3.),
                                    _mm256_setr_pd(1., 2., 3., 4.));]])],
                          [ch_str_fma_h="immintrin.h"])
        if test "x$ch_str_fma_h" = "ximmintrin.h"; then
          ch_str_veccxx_feature="$ch_str_veccxx_feature fma"
          ch_int_veccxx_featurebits=`expr "$ch_int_veccxx_featurebits" + '256'`
          # Seems we can only use immintrin.h reliably
          ch_str_intelvec_h="immintrin.h"
        fi
      fi
#DBG      echo "Values: $ch_str_veccxx_feature $ch_int_veccxx_featurebits $ch_str_fma_h $ch_str_intelvec_h"

#--AVX2

      case $ch_str_cxxmake in
        portland) CXXFLAGS="-Mvect=sse" ;;
        *       ) CXXFLAGS="-mavx2" ;;
      esac
      if test "x$ch_str_have_cpuavx2" = "xavx2"; then
        AC_COMPILE_IFELSE([AC_LANG_PROGRAM(
[[#include <immintrin.h>  // AVX2]],
[[__m256i test_v8si =
  _mm256_set1_epi32(1) + _mm256_setr_epi32(1, 2, 3, 4, 5, 6, 7, 8);]])],
                          [ch_str_avx2_h="immintrin.h"])
        if test "x$ch_str_avx2_h" = "ximmintrin.h"; then
          if test "x$ch_stat_vecopp" != "xno"; then
            ch_stat_vecopp="yes"
          fi
        else
          # Vector ops do not work
          ch_stat_vecopp="no"
          AC_COMPILE_IFELSE([AC_LANG_PROGRAM(
[[#include <immintrin.h>  // AVX2]],
[[__m256i test_v8si = _mm256_add_epi32(
  _mm256_set1_epi32(1), _mm256_setr_epi32(1, 2, 3, 4, 5, 6, 7, 8));]])],
                            [ch_str_avx2_h="immintrin.h"])
        fi
        if test "x$ch_str_avx2_h" = "ximmintrin.h"; then
          ch_str_veccxx_feature="$ch_str_veccxx_feature avx2"
          ch_int_veccxx_featurebits=`expr "$ch_int_veccxx_featurebits" + '512'`
          # Seems we can only use immintrin.h reliably
          ch_str_intelvec_h="immintrin.h"
        fi
      fi
#DBG      echo "Values: $ch_str_veccxx_feature $ch_int_veccxx_featurebits $ch_str_avx2_h $ch_str_intelvec_h"
    fi

# Reset flags we changed
    CXXFLAGS="$ch_str_savecxxflags"

#--Results

    AC_MSG_RESULT([$ch_str_veccxx_feature])
    if test "$ch_int_veccpu_featurebits" -ne "$ch_int_veccxx_featurebits"; then
      if test "$ch_int_veccxx_featurebits" -eq "0"; then
        AC_MSG_NOTICE([The CPU supports the vector extensions: $ch_str_veccpu_feature, but compiler intrinsics do not work])
      else
        AC_MSG_NOTICE([The CPU supports the vector extensions: $ch_str_veccpu_feature, but compiler intrinsics only work for: $ch_str_veccxx_feature])
      fi
    fi
    if test "$ch_int_veccpu_featurebits" -ne "0"; then
      AC_DEFINE_UNQUOTED([CHDEF_SYSTEM_X86VECEXT_CPU_BITS],
                         [$ch_int_veccpu_featurebits])
    fi
    if test "$ch_int_veccxx_featurebits" -ne "0"; then
      AC_DEFINE_UNQUOTED([CHDEF_SYSTEM_X86VECEXT_COMPILER_BITS],
                         [$ch_int_veccxx_featurebits])
    fi
    if test "x$ch_str_sse_h" != "xunknown"; then
      AC_DEFINE_UNQUOTED([CHDEF_SYSTEM_X86VECEXT_SSE_H], [<$ch_str_sse_h>])
    fi
    if test "x$ch_str_sse2_h" != "xunknown"; then
      AC_DEFINE_UNQUOTED([CHDEF_SYSTEM_X86VECEXT_SSE2_H], [<$ch_str_sse2_h>])
    fi
    if test "x$ch_str_sse3_h" != "xunknown"; then
      AC_DEFINE_UNQUOTED([CHDEF_SYSTEM_X86VECEXT_SSE3_H], [<$ch_str_sse3_h>])
    fi
    if test "x$ch_str_ssse3_h" != "xunknown"; then
      AC_DEFINE_UNQUOTED([CHDEF_SYSTEM_X86VECEXT_SSSE3_H], [<$ch_str_ssse3_h>])
    fi
    if test "x$ch_str_sse4d1_h" != "xunknown"; then
      AC_DEFINE_UNQUOTED([CHDEF_SYSTEM_X86VECEXT_SSE4D1_H], [<$ch_str_sse4d1_h>])
    fi
    if test "x$ch_str_sse4d2_h" != "xunknown"; then
      AC_DEFINE_UNQUOTED([CHDEF_SYSTEM_X86VECEXT_SSE4D2_H], [<$ch_str_sse4d2_h>])
    fi
    if test "x$ch_str_sse4a_h" != "xunknown"; then
      AC_DEFINE_UNQUOTED([CHDEF_SYSTEM_X86VECEXT_SSE4A_H], [<$ch_str_sse4a_h>])
    fi
    if test "x$ch_str_avx_h" != "xunknown"; then
      AC_DEFINE_UNQUOTED([CHDEF_SYSTEM_X86VECEXT_AVX_H], [<$ch_str_avx_h>])
    fi
    if test "x$ch_str_fma_h" != "xunknown"; then
      AC_DEFINE_UNQUOTED([CHDEF_SYSTEM_X86VECEXT_FMA_H], [<$ch_str_fma_h>])
    fi
    if test "x$ch_str_avx2_h" != "xunknown"; then
      AC_DEFINE_UNQUOTED([CHDEF_SYSTEM_X86VECEXT_AVX2_H], [<$ch_str_avx2_h>])
    fi
    if test "x$ch_str_intelvec_h" != "xunknown"; then
      AC_DEFINE_UNQUOTED([CHDEF_SYSTEM_X86VECEXT_INTEL_H], [<$ch_str_intelvec_h>])
    fi
    if test "x$ch_stat_vecopp" = "xyes"; then
      AC_DEFINE([CHDEF_SYSTEM_X86VECEXT_HAVE_OPERATOR])
    fi
    ;;

  *    )
    AC_MSG_RESULT([not an x86 processor])
    ;;
esac

#--------------------------------------------------------------------
# Check for CUDA
#--------------------------------------------------------------------

# This sets the define USE_GPU compiler option
use_gpu=no
AC_ARG_ENABLE([gpu],
              [AC_HELP_STRING([--enable-gpu],
                              [use gpu for computation (default=no)])],
              [use_gpu=yes])

ch_stat_cudart=unknown
ch_stat_cuda=unknown
ch_str_nvccpath=system
ch_str_cudatkdir=unknown
# Driver API
ch_str_cudahdir=unknown
ch_str_cudalibdir=unknown
# Runtime API
ch_str_cudarthdir=unknown
ch_str_cudartlibdir=unknown
# Device properties
ch_int_cudadevicecount=0
ch_str_cudaattributes=unknown
ch_str_cudaattributesprint=""
ch_int_cudaattr_computemajor=0
ch_int_cudaattr_computeminor=0

AC_ARG_VAR([NVCC], [Nvidia CUDA compiler])
AC_ARG_VAR([NVCCFLAGS], [Nvidia CUDA compiler flags])
AC_ARG_VAR([PTXASFLAGS], [Options for Nvidia ptxas tool])

# nvcc should normally be in the path
AC_ARG_WITH([nvcc],
            [AC_HELP_STRING([--with-nvcc=<path>],
                            [path to Nvidia nvcc CUDA compiler])],
            [ch_str_nvccpath=$withval])
if test "x$ch_str_nvccpath" = "xsystem"; then
  AC_CHECK_PROG([ch_str_nvcc], [nvcc], [nvcc], [unknown])
else  
  AC_PATH_PROG([ch_str_nvcc], [nvcc], [unknown], [$ch_str_nvccpath])
fi
if test "x$use_gpu" = "xyes" && test "x$ch_str_nvcc" != "xunknown"; then
  ch_stat_cuda="havenvcc"

  # First check for the runtime API, normally installed to /usr/local/cuda.  We
  # also check for the driver API header which is usually installed with the
  # toolkit
  AC_ARG_WITH([cudatk],
              [AC_HELP_STRING([--with-cudatk=<path>],
                              [path to CUDA toolkit installation])],
              [ch_str_cudatkdir=$withval])
  if test "x$ch_str_cudatkdir" = "xunknown"; then
    # Runtime header
    AC_CHECK_HEADER([cuda_runtime.h],
                    [ch_str_cudarthdir="system"])
    if test "x$ch_str_cudarthdir" = "xunknown"; then
      AC_MSG_CHECKING([standard locations for cuda_runtime.h])
      if test -f "/usr/local/cuda/include/cuda_runtime.h"; then
        ch_str_cudarthdir="/usr/local/cuda/include"
        AC_MSG_RESULT([$ch_str_cudarthdir])
      else
        AC_MSG_RESULT([no])
      fi
    fi
    # Driver header
    AC_CHECK_HEADER([cuda.h],
                    [ch_str_cudahdir="system"])
    if test "x$ch_str_cudahdir" = "xunknown"; then
      AC_MSG_CHECKING([standard locations for cuda.h])
      if test -f "/usr/local/cuda/include/cuda.h"; then
        ch_str_cudahdir="/usr/local/cuda/include"
        AC_MSG_RESULT([$ch_str_cudahdir])
      else
        AC_MSG_RESULT([no])
      fi
    fi
    # On CFDP systems, libcudart.so is set in the runtime library path but
    # not in the compile time library path.  So this will usually fail.
    AC_CHECK_LIB([cudart],
                 [cudaMemcpy],
                 [ch_str_cudartlibdir="system"])
    if test "x$ch_str_cudartlibdir" = "xunknown"; then
      AC_MSG_CHECKING([standard locations for libcudart.so])
      if test "$ch_int_systemptrsize" -eq "64"; then
        if test -f "/usr/local/cuda/lib64/libcudart.so"; then
          ch_str_cudartlibdir="/usr/local/cuda/lib64"
        fi
      else
        if test -f "/usr/local/cuda/lib/libcudart.so"; then
          ch_str_cudartlibdir="/usr/local/cuda/lib"
        fi
      fi
      if test "x$ch_str_cudartlibdir" = "xunknown"; then
        AC_MSG_RESULT([no])
      else
        AC_MSG_RESULT([$ch_str_cudartlibdir])
      fi
    fi
  else
    if test -f "$ch_str_cudatkdir/include/cuda_runtime.h"; then
      ch_str_cudarthdir=$ch_str_cudatkdir/include
    fi
    if test -f "$ch_str_cudatkdir/include/cuda.h"; then
      ch_str_cudahdir=$ch_str_cudatkdir/include
    fi
    if test "$ch_int_systemptrsize" -eq "64"; then
      if test -f "$ch_str_cudatkdir/lib64/libcudart.so"; then
        ch_str_cudartlibdir="$ch_str_cudatkdir/lib64"
      fi
    else
      if test -f "$ch_str_cudatkdir/lib/libcudart.so"; then
        ch_str_cudartlibdir="$ch_str_cudatkdir/lib"
      fi
    fi
  fi

  # Cuda is normally installed in the system at /usr/lib/
  AC_CHECK_LIB([cuda], [cuInit], [ch_stat_cuda="havedriver"])

  # Print status
  if test "x$ch_str_cudarthdir" = "xunknown"; then
    AC_MSG_WARN([CUDA runtime header cuda_runtime.h not found])
  else
    AC_MSG_NOTICE([CUDA runtime header location: $ch_str_cudarthdir])
  fi
  if test "x$ch_str_cudartlibdir" = "xunknown"; then
    AC_MSG_WARN([CUDA runtime library libcudart.so not found])
  else
    AC_MSG_NOTICE([CUDA runtime library location: $ch_str_cudartlibdir])
  fi
  if test "x$ch_str_cudahdir" = "xunknown"; then
    AC_MSG_WARN([CUDA driver header cuda.h not found])
  else
    AC_MSG_NOTICE([CUDA driver header location: $ch_str_cudahdir])
  fi
  if test "x$ch_stat_cuda" = "xhavedriver"; then
    AC_MSG_NOTICE([CUDA driver library location: system])
  else
    AC_MSG_WARN([CUDA driver library libcuda.so not found])
  fi

  if test "x$ch_str_cudartlibdir" != "xunknown"; then

#--Have nvcc and runtime

    if test "x$ch_stat_cuda" = "xhavedriver"; then

#--Have nvcc, runtime, and driver

      ch_stat_cuda="yes"
      # Avoid using any optimization flags in the following tests -- the
      # assembly can cause problems.
      ch_str_savecxxflags="$CXXFLAGS"
      ch_str_savecppflags="$CPPFLAGS"
      CXXFLAGS="-lcuda"
      if test "x$ch_str_cudahdir" != "xsystem"; then
        CPPFLAGS="-I$ch_str_cudahdir"
      fi

#--Just run init

      AC_RUN_IFELSE([AC_LANG_PROGRAM(
[[#include <cuda.h>]],
[[CUresult cuer = cuInit(0);
return cuer;]])],
                    [ch_stat_cuda="yes"],
                    [ch_stat_cuda="failed"])
      if test "x$ch_stat_cuda" = "xfailed"; then
        AC_MSG_WARN([The CUDA toolkit and driver were found but a simple program failed.])
      else

#--Look for CUDA devices

        AC_MSG_CHECKING([number of CUDA devices])
        AC_RUN_IFELSE([AC_LANG_PROGRAM(
[[#include <fstream>
#include <cuda.h>]],
[[CUresult cuer = cuInit(0);
if (cuer) return 1;
int deviceCount;
cuer = cuDeviceGetCount(&deviceCount);
if (cuer) return 1;
std::ofstream ofile("conftest.data");
ofile << deviceCount;
ofile.close();
return 0;]])],
                      [ch_stat_cuda="yes"],
                      [ch_stat_cuda="failed"])
        if test "x$ch_stat_cuda" = "xyes"; then
          ch_int_cudadevicecount=`cat conftest.data`
        fi
        AC_MSG_RESULT([$ch_int_cudadevicecount])
        if test "$ch_int_cudadevicecount" -eq "0"; then
          AC_MSG_NOTICE([No CUDA devices found])
        else
          if test "$ch_int_cudadevicecount" -gt "1"; then
            AC_MSG_NOTICE([More than 1 CUDA device found.  Properties taken from first device])
          fi

#--Name and other attributes.  This test program reproduces in a similar form
#--output from Nvidia's deviceQueryDrv sample

          AC_RUN_IFELSE([AC_LANG_PROGRAM(
[[#include <fstream>
#include <iomanip>
#include <cuda.h>]],
[[CUresult cuer = cuInit(0);
if (cuer) return 1;
int dev = 0;

int driverVersion;
char deviceName[256];
int major;
int minor;
size_t totalGlobalMem;
int multiProcessorCount;
int totalConstantMem;
int sharedMemPerMP;
int sharedMemPerBlock;
int registersPerMP;
int registersPerBlock;
int warpSize;
int maxThreadsPerMP;
int maxThreadsPerBlock;
int maxBlockDimX;
int maxBlockDimY;
int maxBlockDimZ;
int maxGridDimX;
int maxGridDimY;
int maxGridDimZ;

cuer = cuDeviceComputeCapability(&major, &minor, dev);
if (cuer) return 1;
cuer = cuDeviceGetName(deviceName, 256, dev);
if (cuer) return 1;
cuer = cuDriverGetVersion(&driverVersion);
if (cuer) return 1;
cuer = cuDeviceTotalMem(&totalGlobalMem, dev);
if (cuer) return 1;

cuer = cuDeviceGetAttribute(&multiProcessorCount,
                            CU_DEVICE_ATTRIBUTE_MULTIPROCESSOR_COUNT, dev);
if (cuer) return 1;
cuer = cuDeviceGetAttribute(&totalConstantMem,
                            CU_DEVICE_ATTRIBUTE_TOTAL_CONSTANT_MEMORY, dev);
if (cuer) return 1;
cuer = cuDeviceGetAttribute(&sharedMemPerMP,
                            CU_DEVICE_ATTRIBUTE_MAX_SHARED_MEMORY_PER_MULTIPROCESSOR,
                            dev);
if (cuer) return 1;
cuer = cuDeviceGetAttribute(&sharedMemPerBlock,
                            CU_DEVICE_ATTRIBUTE_MAX_SHARED_MEMORY_PER_BLOCK,
                            dev);
if (cuer) return 1;
cuer = cuDeviceGetAttribute(&registersPerMP,
                            CU_DEVICE_ATTRIBUTE_MAX_REGISTERS_PER_MULTIPROCESSOR, dev);
if (cuer) return 1;
cuer = cuDeviceGetAttribute(&registersPerBlock,
                            CU_DEVICE_ATTRIBUTE_MAX_REGISTERS_PER_BLOCK, dev);
if (cuer) return 1;
cuer = cuDeviceGetAttribute(&warpSize,
                            CU_DEVICE_ATTRIBUTE_WARP_SIZE, dev);
if (cuer) return 1;
cuer = cuDeviceGetAttribute(&maxThreadsPerMP,
                            CU_DEVICE_ATTRIBUTE_MAX_THREADS_PER_MULTIPROCESSOR,
                            dev);
if (cuer) return 1;
cuer = cuDeviceGetAttribute(&maxThreadsPerBlock,
                            CU_DEVICE_ATTRIBUTE_MAX_THREADS_PER_BLOCK, dev);
if (cuer) return 1;
cuer = cuDeviceGetAttribute(&maxBlockDimX,
                            CU_DEVICE_ATTRIBUTE_MAX_BLOCK_DIM_X, dev);
if (cuer) return 1;
cuer = cuDeviceGetAttribute(&maxBlockDimY,
                            CU_DEVICE_ATTRIBUTE_MAX_BLOCK_DIM_Y, dev);
if (cuer) return 1;
cuer = cuDeviceGetAttribute(&maxBlockDimZ,
                            CU_DEVICE_ATTRIBUTE_MAX_BLOCK_DIM_Z, dev);
if (cuer) return 1;
cuer = cuDeviceGetAttribute(&maxGridDimX,
                            CU_DEVICE_ATTRIBUTE_MAX_GRID_DIM_X, dev);
if (cuer) return 1;
cuer = cuDeviceGetAttribute(&maxGridDimY,
                            CU_DEVICE_ATTRIBUTE_MAX_GRID_DIM_Y, dev);
if (cuer) return 1;
cuer = cuDeviceGetAttribute(&maxGridDimZ,
                            CU_DEVICE_ATTRIBUTE_MAX_GRID_DIM_Z, dev);
if (cuer) return 1;

std::ofstream ofile("conftest.data");
ofile << "DRIVER_VERSION=" << driverVersion << '=';
ofile << " DEVICE_NAME=" << deviceName << '=';
ofile << " COMPUTE_MAJOR=" << major << '=';
ofile << " COMPUTE_MINOR=" << minor << '=';
ofile << " TOTAL_GLOBAL_MEMORY=" << totalGlobalMem << '=';
ofile << " MULTIPROCESSOR_COUNT=" << multiProcessorCount << '=';
ofile << " TOTAL_CONSTANT_MEMORY=" << totalConstantMem << '=';
ofile << " TOTAL_SHARED_MEMORY_PER_MULTIPROCESSOR=" << sharedMemPerMP << '=';
ofile << " MAX_SHARED_MEMORY_PER_BLOCK=" << sharedMemPerBlock << '=';
ofile << " TOTAL_REGISTERS_PER_MULTIPROCESSOR=" << registersPerMP << '=';
ofile << " MAX_REGISTERS_PER_BLOCK=" << registersPerBlock << '=';
ofile << " WARP_SIZE=" << warpSize << '=';
ofile << " MAX_THREADS_PER_MULTIPROCESSOR=" << maxThreadsPerMP << '=';
ofile << " MAX_THREADS_PER_BLOCK=" << maxThreadsPerBlock << '=';
ofile << " MAX_BLOCK_DIM_X=" << maxBlockDimX << '=';
ofile << " MAX_BLOCK_DIM_Y=" << maxBlockDimY << '=';
ofile << " MAX_BLOCK_DIM_Z=" << maxBlockDimZ << '=';
ofile << " MAX_GRID_DIM_X=" << maxGridDimX << '=';
ofile << " MAX_GRID_DIM_Y=" << maxGridDimY << '=';
ofile << " MAX_GRID_DIM_Z=" << maxGridDimZ << '=';
ofile.close();
const char *const indent = "  ";
const int labelWidth = 50;
std::ofstream pfile("conftest.print");
pfile << indent << std::left << std::setw(labelWidth)
      << "CUDA driver version:"
      << driverVersion/1000 << '.' << (driverVersion%100)/10
      << std::endl;
pfile << indent << std::left << std::setw(labelWidth)
      << "Device name:"
      << deviceName
      << std::endl;
pfile << indent << std::left << std::setw(labelWidth)
      << "CUDA compute capability:"
      << major << '.' << minor
      << std::endl;
pfile << indent << std::left << std::setw(labelWidth)
      << "Total amount of global memory:"
      << totalGlobalMem << " bytes"
      << std::endl;
pfile << indent << std::left << std::setw(labelWidth)
      << "Number of multiprocessors:"
      << multiProcessorCount
      << std::endl;
pfile << indent << std::left << std::setw(labelWidth)
      << "Total amount of constant memory:"
      << totalConstantMem << " bytes"
      << std::endl;
pfile << indent << std::left << std::setw(labelWidth)
      << "Total amount of shared memory per multiprocessor:"
      << sharedMemPerMP << " bytes"
      << std::endl;
pfile << indent << std::left << std::setw(labelWidth)
      << "Maximum amount of shared memory per block:"
      << sharedMemPerBlock << " bytes"
      << std::endl;
pfile << indent << std::left << std::setw(labelWidth)
      << "Total number of registers per multiprocessor:"
      << registersPerMP
      << std::endl;
pfile << indent << std::left << std::setw(labelWidth)
      << "Max number of registers available per block:"
      << registersPerBlock
      << std::endl;
pfile << indent << std::left << std::setw(labelWidth)
      << "Warp size:"
      << warpSize
      << std::endl;
pfile << indent << std::left << std::setw(labelWidth)
      << "Maximum number of threads per multiprocessor:"
      << maxThreadsPerMP
      << std::endl;
pfile << indent << std::left << std::setw(labelWidth)
      << "Maximum number of threads per block:"
      << maxThreadsPerBlock
      << std::endl;
pfile << indent << std::left << std::setw(labelWidth)
      << "Maximum sizes of each dimension of a block:"
      << maxBlockDimX << " x " << maxBlockDimY << " x " << maxBlockDimZ
      << std::endl;
pfile << indent << std::left << std::setw(labelWidth)
      << "Maximum sizes of each dimension of a grid:"
      << maxGridDimX << " x " << maxGridDimY << " x " << maxGridDimZ
      << std::endl;
pfile.close();
return 0;]])],
                        [ch_stat_cuda="yes"],
                        [ch_stat_cuda="failed"])
          if test "x$ch_stat_cuda" = "xyes"; then
            ch_str_cudaattributes=`cat conftest.data`
            ch_str_cudaattributesprint=`cat conftest.print`
            # Compute codes
            [ch_int_cudaattr_computemajor=`expr "x$ch_str_cudaattributes" : 'x.*COMPUTE_MAJOR=\([^=]*\)'`]
            [ch_int_cudaattr_computeminor=`expr "x$ch_str_cudaattributes" : 'x.*COMPUTE_MINOR=\([^=]*\)'`]
            AC_DEFINE_UNQUOTED([CHDEF_SYSTEM_CUDA_DRIVER_VERSION],
                               [[`expr "x$ch_str_cudaattributes" : 'x.*DRIVER_VERSION=\([^=]*\)'`]])
            AC_DEFINE_UNQUOTED([CHDEF_SYSTEM_CUDA_DEVICE_NAME],
                               [["`expr "x$ch_str_cudaattributes" : 'x.*DEVICE_NAME=\([^=]*\)'`"]])
            AC_DEFINE_UNQUOTED([CHDEF_SYSTEM_CUDA_COMPUTE_CAPABILITY],
                               [$ch_int_cudaattr_computemajor.$ch_int_cudaattr_computeminor])
            AC_DEFINE_UNQUOTED([CHDEF_SYSTEM_CUDAATTR_TOTAL_GLOBAL_MEMORY],
                               [[`expr "x$ch_str_cudaattributes" : 'x.*TOTAL_GLOBAL_MEMORY=\([^=]*\)'`]])
            AC_DEFINE_UNQUOTED([CHDEF_SYSTEM_CUDAATTR_TOTAL_CONSTANT_MEMORY],
                               [[`expr "x$ch_str_cudaattributes" : 'x.*TOTAL_CONSTANT_MEMORY=\([^=]*\)'`]])
            AC_DEFINE_UNQUOTED([CHDEF_SYSTEM_CUDAATTR_TOTAL_SHARED_MEMORY_PER_MULTIPROCESSOR],
                               [[`expr "x$ch_str_cudaattributes" : 'x.*TOTAL_SHARED_MEMORY_PER_MULTIPROCESSOR=\([^=]*\)'`]])
            AC_DEFINE_UNQUOTED([CHDEF_SYSTEM_CUDAATTR_MAX_SHARED_MEMORY_PER_BLOCK],
                               [[`expr "x$ch_str_cudaattributes" : 'x.*MAX_SHARED_MEMORY_PER_BLOCK=\([^=]*\)'`]])
            AC_DEFINE_UNQUOTED([CHDEF_SYSTEM_CUDAATTR_TOTAL_REGISTERS_PER_MULTIPROCESSOR],
                               [[`expr "x$ch_str_cudaattributes" : 'x.*TOTAL_REGISTERS_PER_MULTIPROCESSOR=\([^=]*\)'`]])
            AC_DEFINE_UNQUOTED([CHDEF_SYSTEM_CUDAATTR_MAX_REGISTERS_PER_BLOCK],
                               [[`expr "x$ch_str_cudaattributes" : 'x.*MAX_REGISTERS_PER_BLOCK=\([^=]*\)'`]])
            AC_DEFINE_UNQUOTED([CHDEF_SYSTEM_CUDAATTR_WARP_SIZE],
                               [[`expr "x$ch_str_cudaattributes" : 'x.*WARP_SIZE=\([^=]*\)'`]])
            AC_DEFINE_UNQUOTED([CHDEF_SYSTEM_CUDAATTR_MAX_THREADS_PER_MULTIPROCESSOR],
                               [[`expr "x$ch_str_cudaattributes" : 'x.*MAX_THREADS_PER_MULTIPROCESSOR=\([^=]*\)'`]])
            AC_DEFINE_UNQUOTED([CHDEF_SYSTEM_CUDAATTR_MAX_THREADS_PER_BLOCK],
                               [[`expr "x$ch_str_cudaattributes" : 'x.*MAX_THREADS_PER_BLOCK=\([^=]*\)'`]])
            AC_DEFINE_UNQUOTED([CHDEF_SYSTEM_CUDAATTR_MAX_BLOCK_DIM_X],
                               [[`expr "x$ch_str_cudaattributes" : 'x.*MAX_BLOCK_DIM_X=\([^=]*\)'`]])
            AC_DEFINE_UNQUOTED([CHDEF_SYSTEM_CUDAATTR_MAX_BLOCK_DIM_Y],
                               [[`expr "x$ch_str_cudaattributes" : 'x.*MAX_BLOCK_DIM_Y=\([^=]*\)'`]])
            AC_DEFINE_UNQUOTED([CHDEF_SYSTEM_CUDAATTR_MAX_BLOCK_DIM_Z],
                               [[`expr "x$ch_str_cudaattributes" : 'x.*MAX_BLOCK_DIM_Z=\([^=]*\)'`]])
            AC_DEFINE_UNQUOTED([CHDEF_SYSTEM_CUDAATTR_MAX_GRID_DIM_X],
                               [[`expr "x$ch_str_cudaattributes" : 'x.*MAX_GRID_DIM_X=\([^=]*\)'`]])
            AC_DEFINE_UNQUOTED([CHDEF_SYSTEM_CUDAATTR_MAX_GRID_DIM_Y],
                               [[`expr "x$ch_str_cudaattributes" : 'x.*MAX_GRID_DIM_Y=\([^=]*\)'`]])
            AC_DEFINE_UNQUOTED([CHDEF_SYSTEM_CUDAATTR_MAX_GRID_DIM_Z],
                               [[`expr "x$ch_str_cudaattributes" : 'x.*MAX_GRID_DIM_Z=\([^=]*\)'`]])
          fi
# end have device
        fi
# end init passed
      fi
      CXXFLAGS="$ch_str_savecxxflags"
      CPPFLAGS="$ch_str_savecppflags"
      # capability 2.1 maps to code 2.1 and arch 2.0.  As of 3.5, all others
      # map directly.  Note that arch must be virtual (compute_) if code is
      # specified
      if test "$ch_int_cudaattr_computemajor" -eq "2" && test "$ch_int_cudaattr_computeminor" -eq "1"; then
        ch_str_gpuarch="compute_20"
      else
        ch_str_gpuarch="compute_$ch_int_cudaattr_computemajor$ch_int_cudaattr_computeminor"
      fi
      ch_str_gpucode="sm_$ch_int_cudaattr_computemajor$ch_int_cudaattr_computeminor"
# else have no driver
    else

#--Have nvcc and runtime but no driver

      ch_stat_cuda="yes"
      AC_MSG_WARN([The CUDA toolkit and runtime library were found but not the driver.  You can compile, but not execute, CUDA code on this system.])
      AC_MSG_WARN([Device properties are set to reasonable defaults for compute capability 3.0.])
      # Outut string
      ch_str_cudaattributesprint=$'  CUDA driver version:                              8.0\n  Device name:                                      Unknown\n  CUDA compute capability:                          3.0\n  Total amount of global memory:                    2091974656 bytes\n  Number of multiprocessors:                        1\n  Total amount of constant memory:                  65536 bytes\n  Total amount of shared memory per multiprocessor: 49152 bytes\n  Maximum amount of shared memory per block:        49152 bytes\n  Total number of registers per multiprocessor:     65536\n  Max number of registers available per block:      65536\n  Warp size:                                        32\n  Maximum number of threads per multiprocessor:     2048\n  Maximum number of threads per block:              1024\n  Maximum sizes of each dimension of a block:       1024 x 1024 x 64\n  Maximum sizes of each dimension of a grid:        2147483647 x 65535 x 65535'
      # Compute codes
      AC_DEFINE_UNQUOTED([CHDEF_SYSTEM_CUDA_DRIVER_VERSION],
                         [8.0])
      AC_DEFINE_UNQUOTED([CHDEF_SYSTEM_CUDA_DEVICE_NAME],
                         [Unknown])
      AC_DEFINE_UNQUOTED([CHDEF_SYSTEM_CUDA_COMPUTE_CAPABILITY],
                         [3.0])
      AC_DEFINE_UNQUOTED([CHDEF_SYSTEM_CUDAATTR_TOTAL_GLOBAL_MEMORY],
                         [2091974656])
      AC_DEFINE_UNQUOTED([CHDEF_SYSTEM_CUDAATTR_TOTAL_CONSTANT_MEMORY],
                         [65536])
      AC_DEFINE_UNQUOTED([CHDEF_SYSTEM_CUDAATTR_TOTAL_SHARED_MEMORY_PER_MULTIPROCESSOR],
                         [49152])
      AC_DEFINE_UNQUOTED([CHDEF_SYSTEM_CUDAATTR_MAX_SHARED_MEMORY_PER_BLOCK],
                         [49152])
      AC_DEFINE_UNQUOTED([CHDEF_SYSTEM_CUDAATTR_TOTAL_REGISTERS_PER_MULTIPROCESSOR],
                         [65536])
      AC_DEFINE_UNQUOTED([CHDEF_SYSTEM_CUDAATTR_MAX_REGISTERS_PER_BLOCK],
                         [65536])
      AC_DEFINE_UNQUOTED([CHDEF_SYSTEM_CUDAATTR_WARP_SIZE],
                         [32])
      AC_DEFINE_UNQUOTED([CHDEF_SYSTEM_CUDAATTR_MAX_THREADS_PER_MULTIPROCESSOR],
                         [2048])
      AC_DEFINE_UNQUOTED([CHDEF_SYSTEM_CUDAATTR_MAX_THREADS_PER_BLOCK],
                         [1024])
      AC_DEFINE_UNQUOTED([CHDEF_SYSTEM_CUDAATTR_MAX_BLOCK_DIM_X],
                         [1024])
      AC_DEFINE_UNQUOTED([CHDEF_SYSTEM_CUDAATTR_MAX_BLOCK_DIM_Y],
                         [1024])
      AC_DEFINE_UNQUOTED([CHDEF_SYSTEM_CUDAATTR_MAX_BLOCK_DIM_Z],
                         [64])
      AC_DEFINE_UNQUOTED([CHDEF_SYSTEM_CUDAATTR_MAX_GRID_DIM_X],
                         [2147483647])
      AC_DEFINE_UNQUOTED([CHDEF_SYSTEM_CUDAATTR_MAX_GRID_DIM_Y],
                         [65535])
      AC_DEFINE_UNQUOTED([CHDEF_SYSTEM_CUDAATTR_MAX_GRID_DIM_Z],
                         [65535])
####################Assuming default values for device properties])
      ch_str_gpuarch="compute_30"
      ch_str_gpucode="compute_30"
# end driver
    fi
  
    # Make settings for Cuda
    NVCC="$ch_str_nvcc"
    if test "x$ch_str_cudahdir" != "xsystem"; then
      CPPFLAGS="$CPPFLAGS-I$ch_str_cudahdir "
    fi
    if test "x$ch_str_cudartlibdir" != "xsystem"; then
      LDFLAGS="$LDFLAGS-L$ch_str_cudartlibdir -lcudart -Wl,-rpath,$ch_str_cudartlibdir "
    else
      LDFLAGS="$LDFLAGS -lcudart "
    fi
    # use c++11 and disable pedantic on host
    NVCCFLAGS="--std=c++11 "
    CXXFLAGS=`echo $CXXFLAGS | sed '1{s/ -pedantic//}'`
    # pointer size
    if test "x$ch_int_systemptrsize" != "xunknown"; then
      NVCCFLAGS="$NVCCFLAGS-m$ch_int_systemptrsize "
    fi
    NVCCFLAGS="$NVCCFLAGS-gencode=arch=$ch_str_gpuarch,code=$ch_str_gpucode "
    PTXASFLAGS="--ptxas-options=-v"
# end have runtime
  fi
# end have nvcc
fi

if test "x$use_gpu" = "xyes"; then
  if test "x$ch_stat_cuda" = "xyes"; then
    CPPFLAGS="$CPPFLAGS-DUSE_GPU "
  else
    AC_MSG_WARN([GPU requested but CUDA installation was not found!])
  fi
fi

#--------------------------------------------------------------------
# Reset the C++ and C compilers to mpi wrapper if using MPI
#--------------------------------------------------------------------

if test "$use_mpi" = "yes"; then
   CXX=$PCXX
   CC=$PCC
fi

AC_SUBST([SCXX])
AC_SUBST([SCC])
AC_SUBST([INSTBINDIR])

AC_OUTPUT

#--------------------------------------------------------------------
# Print useful information
#--------------------------------------------------------------------

# Cuda attributes
if test "x$ch_stat_cuda" = "xyes"; then
  echo ""
  echo "CUDA Information and Attributes for Device 0"
  echo ""
  echo "$ch_str_cudaattributesprint"
fi

echo ""
echo "Configuration summary"
echo   "  Precision       : $use_precision"
echo   "  Space Dimension : $use_dim"
if test "$DEBUG" = "yes"; then
  echo "  Optimization    : debug"
else
  echo "  Optimization    : release"
fi
if test "x$use_mpi" = "xyes"; then
  echo "  MPI             : enabled"
else
  echo "  MPI             : disabled"
fi
if test "x$use_openmp" = "xyes"; then
  echo "  OpenMP          : enabled"
else
  echo "  OpenMP          : disabled"
fi
if test "x$use_gpu" = "xyes"; then
  echo "  GPU             : enabled"
else
  echo "  GPU             : disabled"
fi
