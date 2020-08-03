dnl This file is part of preseq
dnl
dnl Copyright (C) 2020: Andrew D. Smith
dnl
dnl Authors: Andrew D. Smith
dnl
dnl This is free software: you can redistribute it and/or modify it
dnl under the terms of the GNU General Public License as published by
dnl the Free Software Foundation, either version 3 of the License, or
dnl (at your option) any later version.
dnl
dnl This software is distributed in the hope that it will be useful,
dnl but WITHOUT ANY WARRANTY; without even the implied warranty of
dnl MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
dnl General Public License for more details.

AC_INIT([preseq], [3.0.1], [andrewds@usc.edu],
        [preseq], [https://github.com/smithlabcode/preseq])
dnl the config.h is not currently #included in the source, and only
dnl used to keep command lines short.
AC_CONFIG_HEADERS([config.h])
AM_INIT_AUTOMAKE([subdir-objects foreign])

dnl Less output to the terminal. Use "--disable-silent-rules" when
dnl running ./configure to reverse this.
AM_SILENT_RULES([yes])

AC_CONFIG_MACRO_DIR([m4])
AC_LANG(C++)
AC_PROG_CXX
AX_CXX_COMPILE_STDCXX_11([noext], [mandatory])
AC_PROG_RANLIB

dnl recursively configure smithlab_cpp
AC_CONFIG_SUBDIRS([src/smithlab_cpp])

dnl check for HTSLib if requested
hts_fail_msg="

Failed to locate HTSLib on your system. Please use the LDFLAGS and
CPPFLAGS variables to specify the directories where the HTSLib library
and headers can be found.
"
AC_ARG_ENABLE([hts],
  [AS_HELP_STRING([--enable-hts], [Enable HTSLib @<:@yes@:>@])],
  [enable_hts=yes], [enable_hts=no])
AS_IF([test "x$enable_hts" = "xyes"],
  [AC_CHECK_LIB([hts], [hts_version], [],
  [AC_MSG_FAILURE([$hts_fail_msg])])]
)
AM_CONDITIONAL([ENABLE_HTS], [test "x$enable_hts" = "xyes"])

AC_CONFIG_FILES([Makefile])
AC_OUTPUT
