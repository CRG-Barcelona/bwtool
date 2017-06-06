#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_INIT([bwtool],[1.0])
AC_CANONICAL_SYSTEM
AM_INIT_AUTOMAKE
AC_CONFIG_HEADER([config.h])
AC_CONFIG_SRCDIR([bwtool.c])

# Checks for programs.
AC_PROG_CC
AC_PROG_RANLIB
AC_LANG([C])
AM_PROG_CC_C_O
AM_MAINTAINER_MODE([disable])

# Stuff for warnings
case $host in
     *ppc*)	AC_SUBST([MACHTYPE_WARN],["-Wall -Wno-unused-variable"]) ;;
     *)         AC_SUBST([MACHTYPE_WARN],["-Wall -Wformat -Wimplicit -Wreturn-type -Wuninitialized"]) ;;
esac

# Checks for header files.
AC_CHECK_HEADERS([jkweb/common.h beato/bigs.h arpa/inet.h fcntl.h float.h limits.h math.h netdb.h netinet/in.h pthread.h stdlib.h string.h strings.h sys/file.h sys/socket.h sys/statvfs.h sys/time.h syslog.h termios.h unistd.h utime.h])

# REQUIRED libraries

AC_SEARCH_LIBS([pthread_create], [pthread], [has_pthread=yes], [has_pthread=no])
if test "x${has_pthread}" != xyes; then
   AC_MSG_ERROR([
   ---------------------------------------------
   Unable to find libpthread required for compiling
   libjkweb.a. Install this and run ./configure
   with the appropriate LDFLAGS/CFLAGS/etc.
   See the INSTALL for examples.
   ---------------------------------------------])
fi
AC_SEARCH_LIBS([pow], [m], [has_math=yes], [has_math=no])
if test "x${has_math}" != xyes; then
   AC_MSG_ERROR([
   ---------------------------------------------
   Unable to find libm (math) required for compiling
   bwtool. Install this and run ./configure
   with the appropriate LDFLAGS/CFLAGS/etc.
   See the INSTALL for examples.
   ---------------------------------------------])
fi
AC_SEARCH_LIBS([zlibVersion], [z], [has_zlib=yes], [has_zlib=no])
if test "x${has_zlib}" != xyes; then
   AC_MSG_ERROR([
   ---------------------------------------------
   Unable to find zlib required for compiling
   bwtool. Install this and run ./configure
   with the appropriate LDFLAGS/CFLAGS/etc.
   See the INSTALL for examples.
   ---------------------------------------------])
fi

AC_SEARCH_LIBS([slSort],[jkweb],[has_jkweb=yes], [has_jkweb=no])
if test "x${has_jkweb}" != xyes; then
   AC_MSG_ERROR([
   ---------------------------------------------
   Unable to find libjkweb (libbeato)
   required for compiling
   bwtool. Install this and run ./configure
   with the appropriate LDFLAGS/CFLAGS/etc.
   See the INSTALL for examples.
   ---------------------------------------------])
fi

dnl This is a manual test instead of AC_SEARCH_LIBS because the
dnl include file is in a subdir of the include paths.
OLDLIBS=$LIBS;
LIBS="-lbeato $LIBS";
AC_MSG_CHECKING([checking for beato library containing perBaseWigFree])
    AC_LINK_IFELSE([
      AC_LANG_PROGRAM([[
        char perBaseWigFree ();
        ]],
        [[
        return perBaseWigFree();
        ]])],
      [AC_MSG_RESULT([-lbeato]); has_beato=yes;],      
      [AC_MSG_RESULT([no])
       LIBS=$OLDLIBS; LIBS="-lbeato -lhts $LIBS";
       AC_MSG_CHECKING([checking for beato library (compiled with htslib) containing perBaseWigFree])
          AC_LINK_IFELSE([
            AC_LANG_PROGRAM([[
            char perBaseWigFree ();
            ]],
            [[
            return perBaseWigFree();
            ]])],
           [AC_MSG_RESULT([-lbeato -lhts]); has_beato=yes;],
           [AC_MSG_RESULT([no]); has_beato=no;
         ])
     ])
if test "x${has_beato}" != xyes; then
   AC_MSG_ERROR([
   ---------------------------------------------
   libbeato, required for bwtool does not seem to
   be installed.  Install this and run ./configure
   with the appropriate LDFLAGS/CFLAGS/etc.
   See the INSTALL for examples.
   ---------------------------------------------])
fi

# Large file support
AC_SYS_LARGEFILE

# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_C_INLINE
AC_TYPE_INT16_T
AC_TYPE_INT32_T
AC_TYPE_INT64_T
AC_TYPE_INT8_T
AC_TYPE_OFF_T
AC_TYPE_PID_T
AC_TYPE_SIZE_T
AC_TYPE_SSIZE_T
AC_TYPE_UID_T
AC_TYPE_UINT16_T
AC_TYPE_UINT32_T
AC_TYPE_UINT8_T

# Checks for library functions.
AC_FUNC_FORK
AC_FUNC_FSEEKO
AC_FUNC_MALLOC
AC_FUNC_MKTIME
AC_FUNC_MMAP
AC_FUNC_REALLOC
AC_FUNC_STRTOD
AC_CHECK_FUNCS([bzero dup2 floor getcwd gettimeofday memmove memset mkdir munmap pow putenv regcomp rint select setenv socket sqrt strcasecmp strchr strcspn strdup strerror strrchr strspn strstr strtol tzset uname utime])

AC_CONFIG_FILES([Makefile tests/Makefile])

AC_OUTPUT
