dnl Hmmm ... not sure what this does
AC_INIT

dnl Check for compiler
AC_PROG_CC

dnl Checks for libraries.
dnl Replace `main' with a function in -lm:
AC_CHECK_LIB(m, isnan)

dnl Checks for header files.
AC_HEADER_STDC

dnl Checks for typedefs, structures, and compiler characteristics.
AC_C_CONST
AC_STRUCT_TM

dnl Checks for library functions.
AC_FUNC_VPRINTF
AC_CHECK_FUNCS(strcspn)
AC_CHECK_FUNCS(qsort)

dnl Figure out what host we're on -- needed for ${host} below.
AC_CANONICAL_HOST

dnl
dnl C flags
dnl
if test "${GCC}" != "yes" ; then
    case "${host}" in
          *-*-hpux* )
	       CFLAGS="${CFLAGS} -Wp,-H30000"
	       ;;
          *-dec-osf* )
               CFLAGS="${CFLAGS} -ieee"
               ;;
    esac
else
    CFLAGS="-O3 -Wall -D __unix"
    #CFLAGS="-ggdb -Wall -D __unix"
    case "${host}" in
          alpha* )
	       CFLAGS="${CFLAGS} -mieee"
	       ;;
    esac
fi


AC_OUTPUT(Makefile)
