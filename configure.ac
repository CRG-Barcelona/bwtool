#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_INIT([bwtool],[1.0])
AC_CANONICAL_SYSTEM
AM_INIT_AUTOMAKE
AC_CONFIG_HEADER([config.h])
AC_CONFIG_SRCDIR([lib/metaBig.c])

# Checks for programs.
AC_PROG_CC
AC_PROG_RANLIB
AC_LANG([C])
AM_PROG_CC_C_O

# Stuff for the -DCOLOR32 define. This is library-independent so just having the define
# is enough
AC_ARG_ENABLE([color-32],
	[AS_HELP_STRING([--disable-color-32],[disable 32-bit color in PNG etc (default enabled)])], 
	[color32=${enableval}], 
	[color32=yes])

   # to catch that one scenario where some does --disable-color-32=no which makes no sense
   if test "x${color32}" = xno; then
      color32=yes
   fi
   # finally define the variable to put in config.h
   if test "x${color32}" = xyes; then 
      AC_DEFINE([COLOR32],[1],[Use 32-bit color in PNG etc])
   fi

# Stuff for MACHTYPE_ 
case $host in
     *x86_64*)  AC_DEFINE([MACHTYPE_x86_64],[],[Used in some rare cases.]) ;;
     *i386*)    AC_DEFINE([MACHTYPE_i386],[],[Used in some rare cases.]) ;;
     *ppc*)     AC_DEFINE([MACHTYPE_ppc],[],[Used in some rare cases.]) ;;
     *alpha*)   AC_DEFINE([MACHTYPE_alpha],[],[Used in some rare cases.]) ;;
     *sparc*)   AC_DEFINE([MACHTYPE_sparc],[],[Used in some rare cases.]) ;;
     *)         AC_DEFINE([MACHTYPE_unknown],[],[Used in some rare cases.]) ;;
esac

# Stuff for warnings
case $host in
     *ppc*)	AC_SUBST([MACHTYPE_WARN],["-Wall -Wno-unused-variable"]) ;;
     *)         AC_SUBST([MACHTYPE_WARN],["-Wall -Wformat -Wimplicit -Wreturn-type -Wuninitialized"]) ;;
esac

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
   libjkweb.a. Install this and run ./configure
   with the appropriate LDFLAGS/CFLAGS/etc.
   See the INSTALL for examples.
   ---------------------------------------------])
fi
AC_SEARCH_LIBS([zlibVersion], [z], [has_zlib=yes], [has_zlib=no])
if test "x${has_zlib}" != xyes; then
   AC_MSG_ERROR([
   ---------------------------------------------
   Unable to find zlib required for compiling
   libjkweb.a. Install this and run ./configure
   with the appropriate LDFLAGS/CFLAGS/etc.
   See the INSTALL for examples.
   ---------------------------------------------])
fi
AC_SEARCH_LIBS([samopen], [bam], [has_bam=yes], [has_bam=no])
if test "x${has_bam}" = xyes; then
   AC_DEFINE([USE_BAM],[1],[Turn on BAM support.])
fi
AC_SEARCH_LIBS([SSL_accept], [ssl], [has_ssl=yes], [has_ssl=no])
AC_SEARCH_LIBS([SHA1_Init], [crypto], [has_ssl=yes], [has_ssl=no])
if test "x${has_ssl}" = xyes; then
   AC_DEFINE([USE_SSL],[1],[Turn on SSL support in https.c])
fi
AC_SEARCH_LIBS([png_read_png], [png], [has_png=yes], [has_png=no])
if test "x${has_png}" = xyes; then
   AC_DEFINE([USE_PNG],[1],[Turn on PNG support in pngwrite.c])
fi

# Checks for header files.
AC_CHECK_HEADERS([arpa/inet.h fcntl.h float.h limits.h netdb.h netinet/in.h stdlib.h string.h strings.h sys/file.h sys/socket.h sys/statvfs.h sys/time.h syslog.h termios.h unistd.h utime.h])

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

AC_CONFIG_FILES([Makefile lib/Makefile libjkweb/Makefile bwtool/Makefile])
AC_OUTPUT
