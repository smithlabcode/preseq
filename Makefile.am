# This file is part of preseq
#
# Copyright (C) 2018-2024: Andrew D. Smith
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

ACLOCAL_AMFLAGS = -I m4

SUBDIRS := src/smithlab_cpp
install installdirs: SUBDIRS := $(filter-out src/smithlab_cpp, $(SUBDIRS))
AM_CPPFLAGS = -I $(top_srcdir)/src/smithlab_cpp
AM_CPPFLAGS += -Wall -Wextra -Wpedantic -Wno-unknown-attributes

CXXFLAGS = -O2 # default has optimization on

EXTRA_DIST = \
	README.md \
	LICENSE \
	preseqR \
	data \
	docs \
	tests/md5sum.txt \
	tests/data/lc_extrap_input.vals \
	tests/data/gc_extrap_input.mr \
	tests/data/c_curve_input.hist \
	tests/scripts/test_c_curve.test \
	tests/scripts/test_gc_extrap.test \
	tests/scripts/test_lc_extrap.test

TESTS = \
	tests/scripts/test_c_curve.test \
	tests/scripts/test_lc_extrap.test \
	tests/scripts/test_gc_extrap.test

TEST_EXTENSIONS = .test

LDADD = src/smithlab_cpp/libsmithlab_cpp.a

bin_PROGRAMS = preseq
if ENABLE_HTS
bin_PROGRAMS += to-mr
AM_CPPFLAGS += -DHAVE_HTSLIB
endif

preseq_SOURCES = \
        src/preseq.cpp \
        src/common.hpp \
        src/common.cpp \
        src/c_curve.hpp \
        src/c_curve.cpp \
        src/gc_extrap.hpp \
        src/gc_extrap.cpp \
        src/lc_extrap.hpp \
        src/lc_extrap.cpp \
        src/bound_pop.hpp \
        src/bound_pop.cpp \
        src/pop_size.hpp \
        src/pop_size.cpp \
        src/continued_fraction.hpp \
        src/continued_fraction.cpp \
        src/load_data_for_complexity.hpp \
        src/load_data_for_complexity.cpp \
        src/moment_sequence.hpp \
        src/moment_sequence.cpp

to_mr_SOURCES = \
        src/to-mr.cpp
