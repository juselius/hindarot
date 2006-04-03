#
# Process this file with autoconf to produce a configure script.
#
# Written by Jonas Juselius <jonas@iki.fi>
#
# $Id$

AC_INIT([HINDAROT],[1.0],[jonas@iki.fi])
AC_CONFIG_SRCDIR(hindarot.f90)
AC_CONFIG_HEADER(config.h)
AC_CONFIG_AUX_DIR(config)
AC_PREFIX_DEFAULT([/usr/local])
AC_PROG_MAKE_SET
AC_CANONICAL_HOST

# Check for programs.

CHECK_GNU_MAKE
if test "$ifGNUmake" = "#"; then
	AC_MSG_ERROR([No working version of GNU make found!])
fi
AC_PROG_RANLIB
AC_PROG_INSTALL
AC_CHECK_PROGS(LD,ld,false)
AC_CHECK_PROGS(AR,ar,false)
AC_CHECK_PROGS(STRIP,strip,true)
AC_CHECK_PROGS(ETAGS,etags,true)
AC_CHECK_PROGS(CTAGS,ctags,true)
AC_PROG_LN_S

AC_CACHE_SAVE

# ACX_ macros are defined in aclocal.m4
# Setup compilers for different architechtures
case $host in 
	*i*86*-linux*)
		try_cc="gcc icc cc"
		try_fort="ifc pgf90"
		build_type="linux_i386"
		;;
	*alpha*-linux*)
		try_cc="ccc cc gcc"
		try_fort="fort f90"
		build_type="linux_alpha"
		;;
	*dec*-osf*)
		try_cc="ccc cc gcc"
		try_fort="fort f90"
		build_type="dec_osf"
		;;
	*ibm*-aix*)
		try_cc="cc gcc"
		try_fort="xlf90 f90"
		build_type="ibm_aix"
		;;
	*sgi*-irix*)
		try_cc="cc gcc"
		try_fort="f90"
		build_type="sgi_irix"
		;;
	*hppa*-hpux*)
		try_cc="cc gcc"
		try_fort="f90"
		build_type="hp_hpux"
		;;
	*sparc*-solaris*)
		try_cc="cc gcc"
		try_fort="f90"
		build_type="sun_solaris"
		;;
	*) 	
		AC_MSG_WARN([Unknown architecture! This might not work...])
		try_cc="cc gcc ccc cl"
		try_fort="f90 xlf90 ifc pgf90 epcf90 f95 xlf95 lf95 fort"
		build_type="unknown"
		;;
esac

AC_PROG_CC([$try_cc])
AC_PROG_F90([$try_fort])

ACX_BUILD_FLAGS([$build_type])
FFLAGS="$FFLAGS $finclude"
AC_CACHE_SAVE

AC_F90_LIBRARY_LDFLAGS
AC_F90_WRAPPERS

AC_PROG_YACC
AC_PROG_LEX
test "$YACC" = "bison -y" && YACC="$YACC -d"

ACX_GETKW([], AC_MSG_ERROR([libgetkw.a not found!]))

dnl AC_CACHE_SAVE
dnl ACX_BLAS(,AC_MSG_WARN([No working blas library found!]))
dnl ACX_LAPACK(,AC_MSG_WARN([No working lapack library found!]))

#AC_CACHE_SAVE
#AC_MSG_NOTICE([Checking if MPI is available])
#AC_LANG_PUSH([Fortran 90])
#ACX_MPI([F90=$MPIF90])
#AC_LANG_POP
#fdefs="$fdefs -DHAVE_MPI"
#AC_CACHE_SAVE

AC_CACHE_SAVE

ACX_SUBST_BUILD_FLAGS 
AC_SUBST(top_bindir,['$(topdir)'])
AC_SUBST(top_includedir,['$(topdir)/include'])
AC_SUBST(top_libdir,['$(topdir)/lib'])

# Checks for header files.
AC_HEADER_STDC
#AC_CHECK_HEADERS([alloca.h malloc.h stdlib.h string.h unistd.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_C_CONST

# Checks for library functions.
AC_FUNC_MALLOC
AC_CHECK_FUNCS([strdup])

AC_CONFIG_FILES([Makefile Config.mk Rules.mk])
AC_OUTPUT
