# this is example-file: configure.ac

# initial information about the project
AC_INIT([spaced_seed],[0.1],[mchen@cs.stonybrook.edu])

# check if the source folder is correct
AC_CONFIG_SRCDIR([src/spaced_seed.cpp])

# Checks for programs

# check for C++ compiler
AC_PROG_CXX

# automake initialisation (mandatory) and check for minimal automake API version 1.9
AM_INIT_AUTOMAKE([1.9])
AM_INIT_AUTOMAKE([-Wall -Werror foreign])

# use the C++ compiler for the following checks
AC_LANG([C++])

# Checks for header files.
AC_HEADER_STDC
AC_CHECK_HEADERS([string])
AC_CHECK_HEADERS([iostream])
AC_CHECK_HEADERS([limits])
AC_CHECK_FUNCS([memset])
AC_CHECK_HEADERS([fcntl.h])
AC_CHECK_HEADERS([limits.h]) 

AC_C_INLINE 
AC_FUNC_MALLOC
AC_FUNC_MMAP
AC_HEADER_STDBOOL
AC_PROG_CC

# Checks for typedefs, structures, and compiler characteristics.
AC_TYPE_SIZE_T

# distribute additional compiler and linker flags
# --> set these variables instead of CXXFLAGS or LDFLAGS
AC_SUBST([AM_CXXFLAGS])
AC_SUBST([AM_LDFLAGS])
AC_SUBST([LIBS])

# files to generate via autotools (.am or .in source files)
AC_CONFIG_FILES([Makefile])
AC_CONFIG_FILES([src/Makefile])

# generate the final Makefile etc.
AC_OUTPUT
