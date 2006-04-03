dnl
dnl @synopsis ACX_GETKW([ACTION-IF-FOUND[, ACTION-IF-NOT-FOUND]])
dnl

AC_DEFUN([ACX_GETKW], [
AC_PREREQ(2.57)
AC_REQUIRE([AC_PROG_F90])

acx_getkw_ok=no
acx_getkwf_h_ok=no
acx_getkw_save_FFLAGS="$FFLAGS"
acx_getkw_save_LIBS="$LIBS"
acx_getkw_lib="-lgetkw"
acx_getkw_include_paths="-I. -I./include -I/usr/include -I/usr/local/include"
acx_getkw_lib_paths="-L. -L./lib -L/usr/lib -L/usr/local/lib"

AC_ARG_WITH([libgetkw], AC_HELP_STRING([--with-libgetkw=/usr/local],
	[base path to the libgetkw installation]))

case $with_libgetkw in
        yes | "") GETKW_LIBS="$acx_getkw_lib"
				  GETKW_INCLUDES="" ;;
        no) acx_getkw_ok=disabled 
			GETKW_LIBS=""
			GETKW_INCLUDES="" ;;
        */*) GETKW_LIBS="-L$with_libgetkw $acx_getkw_lib"
		     GETKW_INCLUDES="-I$with_libgetkw"
			 acx_getkw_include_paths="-I$with_libgetkw $acx_getkw_include_paths"
			 acx_getkw_lib_paths="-L$with_libgetkw $acx_getkw_lib_paths" ;;
        *) AC_MSG_WARN([Invalid path: $with_libgetkw]) ;;
esac

test "x$with_libgetkw" = "x" && GETKW_LIBS="$acx_getkw_lib"

if test "x$acx_getkw_ok" != "xdisabled"; then
LIBS="$LIBS $GETKW_LIBS"
FFLAGS="$GETKW_INCLUDES $FFLAGS"

AC_LANG_PUSH([Fortran 90])
AC_MSG_CHECKING([for getkw.h])
AC_COMPILE_IFELSE(AC_LANG_PROGRAM([[include 'string_m.h']],
	[[include 'getkwf.h']]), 
   [acx_getkwf_h_ok=yes], [GETKW_INCLUDES=""])
if test $acx_getkwf_h_ok = no; then
for i in $acx_getkw_include_paths; do
  FFLAGS="$i $acx_getkw_save_FFLAGS"
  AC_COMPILE_IFELSE(AC_LANG_PROGRAM([[include 'string_m.h']],[[include 'getkwf.h']]), 
    [acx_getkwf_h_ok=yes], [GETKW_INCLUDES=""])
   if test $acx_getkwf_h_ok = yes; then
	GETKW_INCLUDES="$i"
	break
   fi
done
fi
AC_MSG_RESULT($acx_getkwf_h_ok)

if test $acx_getkwf_h_ok = yes; then
  AC_MSG_CHECKING([for end_parse in $acx_getkw_lib])
AC_LINK_IFELSE(AC_LANG_PROGRAM([],[[call end_parse]]), 
  [acx_getkw_ok=yes], [GETKW_LIBS=""])
if test $acx_getkw_ok = no; then
for i in $acx_getkw_lib_paths; do
  FFLAGS="$i $GETKW_INCLUDES $acx_getkw_save_FFLAGS"
AC_LINK_IFELSE(AC_LANG_PROGRAM([],[[call end_parse]]),
  [acx_getkw_ok=yes], [GETKW_LIBS=""])
  if test $acx_getkw_ok = yes; then
	GETKW_LIBS="$i $acx_getkw_lib"
    break
  fi
done
fi
  AC_MSG_RESULT($acx_getkw_ok)
fi
AC_LANG_POP()

LIBS="$acx_getkw_save_LIBS"
FFLAFS="$acx_getkw_save_FFLAGS"
fi # disabled

AC_SUBST(GETKW_LIBS)
AC_SUBST(GETKW_INCLUDES)

# Finally, execute ACTION-IF-FOUND/ACTION-IF-NOT-FOUND:
if test x"$acx_getkw_ok" = xyes; then
        ifelse([$1],,AC_DEFINE(HAVE_LIBGETKW,1,
		[Define if you have the getkw library.]),[$1])
        :
else
        if test "x$acx_getkw_ok" = "xdisabled"; then
        	acx_getkw_ok=no
		else
        	acx_getkw_ok=no
        	$2
		fi
fi
])dnl ACX_GETKW
