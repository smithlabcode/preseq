# Copyright (C) 2011-2020 University of Southern California and
#                         Andrew D. Smith and Timothy Daley
#
# Authors: Timothy Daley and Andrew D. Smith
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.

PROGS = preseq
ifdef HAVE_HTSLIB
PROGS += to-mr
endif

INCLUDEDIRS = smithlab_cpp
INCLUDEARGS = $(addprefix -I, $(INCLUDEDIRS))
LIBS = -lz
CXX = g++
CXXFLAGS = -std=c++11 -Wall

ifdef DEBUG
CXXFLAGS += -g
else
CXXFLAGS += -O2
endif

ifdef HAVE_HTSLIB
CXXFLAGS += -DHAVE_HTSLIB
LIBS += -lhts
endif

all: $(PROGS)

$(PROGS): $(addprefix smithlab_cpp/, \
	smithlab_os.o smithlab_utils.o GenomicRegion.o \
	OptionParser.o MappedRead.o)

ifdef HAVE_HTSLIB
preseq to-mr: $(addprefix smithlab_cpp/, htslib_wrapper.o)
endif

preseq: continued_fraction.o load_data_for_complexity.o moment_sequence.o

%.o: %.cpp %.hpp
	$(CXX) $(CXXFLAGS) -c -o $@ $< $(INCLUDEARGS)

%: %.cpp
	$(CXX) $(CXXFLAGS) -o $@ $^ $(INCLUDEARGS) $(LIBS)

install: $(PROGS)
	@mkdir -p $(install_dir)/bin
	@install $(PROGS) $(install_dir)/bin

clean:
	@-make -C smithlab_cpp clean
	@-rm -f $(PROGS) *.o

.PHONY: install clean
