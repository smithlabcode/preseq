# This file is part of preseq
#
# Copyright (C) 2018-2022: Andrew D. Smith
#
# Authors: Andrew D. Smith
#
# This is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This software is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.

EXTRA_DIST = README.md LICENSE preseqR data docs
ACLOCAL_AMFLAGS = -I m4

SUBDIRS := src/smithlab_cpp
install installdirs: SUBDIRS := $(filter-out src/smithlab_cpp, $(SUBDIRS))
AM_CPPFLAGS = -I $(top_srcdir)/src/smithlab_cpp

CXXFLAGS = -O2 # default has optimization on

LDADD = src/smithlab_cpp/libsmithlab_cpp.a

bin_PROGRAMS = preseq
if ENABLE_HTS
bin_PROGRAMS += to-mr
AM_CPPFLAGS += -DHAVE_HTSLIB
endif

preseq_SOURCES = \
        src/preseq.cpp \
        src/continued_fraction.hpp \
        src/continued_fraction.cpp \
        src/load_data_for_complexity.hpp \
        src/load_data_for_complexity.cpp \
        src/moment_sequence.hpp \
        src/moment_sequence.cpp

to_mr_SOURCES = \
        src/to-mr.cpp
