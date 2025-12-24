dnl This file is part of preseq
dnl
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

dnl check for HTSLib if requested
hts_fail_msg="

Failed to locate HTSLib on your system. Please use the LDFLAGS and
CPPFLAGS variables to specify the directories where the HTSLib library
and headers can be found.

If you want to disable HTSLib when building preseq, use the --disable-htslib
argument to the configure script.
"
AC_ARG_ENABLE([hts],
  [AS_HELP_STRING([--disable-hts], [Enable HTSLib @<:@yes@:>@])],
  [enable_hts=no], [enable_hts=yes])
AS_IF([test "x$enable_hts" = "xyes"],
  [AC_CHECK_LIB([hts], [hts_version], [],
  [AC_MSG_FAILURE([$hts_fail_msg])])]
)
AM_CONDITIONAL([ENABLE_HTS], [test "x$enable_hts" = "xyes"])

dnl arg for using libdeflate, which might happen by default anyway
AC_ARG_WITH([libdeflate],
  [AS_HELP_STRING([--with-libdeflate], [use libdeflate for reading BAM files])],
  [with_libdeflate=yes], [with_libdeflate=no])
AS_IF([test "x$with_libdeflate" = "xyes"],
  [
    AC_SEARCH_LIBS([libdeflate_deflate_compress], [deflate], [],
                   [AC_MSG_ERROR([--with-libdeflate specified but libdeflate not found])])
  ])

AC_CONFIG_FILES([Makefile])

dnl make the test data files available in the build tree
AC_CONFIG_LINKS([
tests/md5sum.txt:tests/md5sum.txt
tests/c_curve_input.hist:tests/data/c_curve_input.hist
tests/lc_extrap_input.vals:tests/data/lc_extrap_input.vals
tests/gc_extrap_input.bed:tests/data/gc_extrap_input.bed
])

AC_OUTPUT
