dnl Copyright (C) 2018-2025: Andrew D. Smith
dnl
dnl Authors: Andrew D. Smith
dnl
dnl This is free software: you can redistribute it and/or modify it under the
dnl terms of the GNU General Public License as published by the Free Software
dnl Foundation, either version 3 of the License, or (at your option) any later
dnl version.
dnl
dnl This software is distributed in the hope that it will be useful, but
dnl WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
dnl or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
dnl for more details.

AC_INIT([preseq],
        [3.2.0],
        [andrewds@usc.edu],
        [preseq],
        [https://github.com/smithlabcode/preseq])

dnl the config.h is not currently #included in the source, and only used to
dnl keep command lines short.
AC_CONFIG_HEADERS([config.h])
AM_INIT_AUTOMAKE([subdir-objects foreign])

AC_CONFIG_MACRO_DIR([m4])
AC_LANG(C++)
AC_PROG_CXX
AX_CXX_COMPILE_STDCXX(17, [noext], [mandatory])
AC_PROG_RANLIB

hts_fail_msg="

Failed to locate the HTSLib library on your system. Please use the LDFLAGS and
CPPFLAGS variables to specify the directories where the HTSLib library and
headers can be found.

If you want to disable HTSLib when building preseq, use the --disable-htslib
argument to the configure script.
"

libdeflate_fail_msg="

Failed to locate libdeflate on your system. This is a common dependency of
HTSLib. If your HTSLib was built without libdeflate, you can give the
--disable-libdeflate argument to configure. If you do not need to directly
process SAM/BAM files with preseq, consider using the argument
--without-htslib to build preseq without support for reading SAM/BAM.

If you need libdeflate, please use the LDFLAGS and CPPFLAGS variables to
specify the directories where the library and headers can be found on your
system.
"

AC_ARG_ENABLE([htslib],
  [AS_HELP_STRING([--enable-htslib], [use HTSLib to allow BAM/SAM input])],
  [enable_htslib=$enableval], [enable_htslib=no]
)

AC_ARG_ENABLE([libdeflate],
  [AS_HELP_STRING([--enable-libdeflate],
                  [use libdeflate which might be needed by HTSLib])],
  [enable_libdeflate=$enableval], [enable_libdeflate=no]
)

AS_IF([test "x$enable_htslib" = xyes], [
  AX_PTHREAD([],
    [AC_MSG_ERROR([Threads required for HTSLib but threads library not found])]
  )
  AC_SEARCH_LIBS(
    [gzopen], [z], [],
    [AC_MSG_FAILURE(["ZLib required for HTSLib but not found"])]
  )
  AS_IF([test "x$enable_libdeflate" = xyes], [
    AC_SEARCH_LIBS([libdeflate_deflate_compress], [deflate],
      [AC_CHECK_LIB([hts], [hts_version], [], [
        AC_MSG_FAILURE([$hts_fail_msg])], [-pthread -lz -ldeflate])],
        [AC_MSG_ERROR([$libdeflate_fail_msg])])
    ],
    [AC_CHECK_LIB([hts], [hts_version], [],
      [AC_MSG_FAILURE([$hts_fail_msg])], [-pthread -lz])])
])
AM_CONDITIONAL([ENABLE_HTSLIB], [test "x$enable_htslib" = "xyes"])

AC_CONFIG_FILES([Makefile])

dnl make the test data files available in the build tree
AC_CONFIG_LINKS([
tests/md5sum.txt:tests/md5sum.txt
tests/c_curve_input.hist:tests/data/c_curve_input.hist
tests/lc_extrap_input.vals:tests/data/lc_extrap_input.vals
tests/gc_extrap_input.bed:tests/data/gc_extrap_input.bed
])

AC_OUTPUT
