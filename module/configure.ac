AC_INIT(dedxqa,[1.00])
AC_CONFIG_SRCDIR([configure.ac])

AM_INIT_AUTOMAKE
AC_PROG_CXX(CC g++)
LT_INIT([disable-static])

if test $ac_cv_prog_gxx = yes; then
   CXXFLAGS="$CXXFLAGS -Wall -Werror -Wextra -Wshadow -Wno-unused-variable -Wno-return-type -Wno-unused-parameter -Wno-unused-function"
fi

AC_CONFIG_FILES([Makefile])
AC_OUTPUT
