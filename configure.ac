AC_INIT()

dnl Switch to a C++ compiler, and check if it works.
AC_LANG(C++)
AC_PROG_CXX

dnl Neet to check std=c++11
dnl Check for HDF5 support
AC_CONFIG_MACRO_DIR([m4])
m4_include([m4/ax_lib_hdf5.m4])
AX_LIB_HDF5([serial])	 
if test "$with_hdf5" = "no"; then
	AC_MSG_ERROR([Unable to find HDF5, we need serial HDF5.])
fi

dnl TEST_AND_SET_CXXFLAG(flag, [program])
dnl
dnl This attempts to compile a program with a certain compiler flag.
dnl If no program is given, then the minimal C++ program is compiled, and 
dnl this tests just the validity of the compiler flag. 
dnl
define([TEST_AND_SET_CXXFLAG],[
	AC_MSG_CHECKING([if compiler flag $1 works])	
	
	dnl Store the current CXXFLAGS
	save_CXXFLAGS="$CXXFLAGS"

	dnl Append the flag of interest
	CXXFLAGS="$CXXFLAGS $1"
	
	dnl Create an M4 macro, "prog", which expands to a C++ program.
	dnl This should either be a default one or the one specified.
	dnl Note that macros are not local, but there is a stack so push
	dnl the definition on to the stack to prevent clobbering a definition
	dnl that might already exist.
	m4_if([$2],[],[pushdef(prog, [int main(){}])], [pushdef(prog, [$2])])
	
	flag_test=0
	
	dnl See if the compiler runs
	AC_COMPILE_IFELSE([AC_LANG_SOURCE([prog])], [flag_test=1],[flag_test=0])
	
	dnl De-clobber the "prog" macro
	popdef([prog])

	if test $flag_test = 1
	then
		AC_MSG_RESULT(yes)
	else
		AC_MSG_RESULT(no)
		dnl The flag doesn't work, so restore the old CXXFLAGS
		CXXFLAGS="$save_CXXFLAGS"
	fi
])

TEST_AND_SET_CXXFLAG(-Wall)
TEST_AND_SET_CXXFLAG(-Wextra)
TEST_AND_SET_CXXFLAG(-W)
TEST_AND_SET_CXXFLAG(-O3)
TEST_AND_SET_CXXFLAG(-g)

dnl Process Makefile.in to create Makefile
AC_CONFIG_FILES([makefile])
AC_OUTPUT

