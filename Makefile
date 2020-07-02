#    Copyright (C) 2011-2104 University of Southern California and
#                            Andrew D. Smith and Timothy Daley
#
#    Authors: Timothy Daley and Andrew D. Smith
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
#


ifndef SRC_ROOT
SRC_ROOT = $(shell pwd)
endif

ifndef PREFIX
PREFIX = $(SRC_ROOT)
endif

ifndef SMITHLAB_CPP
SMITHLAB_CPP = $(SRC_ROOT)/smithlab_cpp
endif

SOURCES = $(wildcard *.cpp)
OBJECTS = $(patsubst %.cpp,%.o,$(SOURCES))
PROGS = preseq
ifdef HAVE_HTSLIB
PROGS += to-mr
endif
INCLUDEDIRS = $(SMITHLAB_CPP)
INCLUDEARGS = $(addprefix -I,$(INCLUDEDIRS))

LIBS += -lgsl -lgslcblas -lz

CXX = g++
CXXFLAGS = -std=c++11 -Wall -fPIC

OPTFLAGS = -O2
DEBUGFLAGS = -g

ifdef DEBUG
CXXFLAGS += $(DEBUGFLAGS)
endif

ifdef HAVE_HTSLIB
CXXFLAGS += -DHAVE_HTSLIB
LIBS += -lhts
endif

ifdef OPT
CXXFLAGS += $(OPTFLAGS)
endif

all: $(PROGS)

$(PROGS): $(addprefix $(SMITHLAB_CPP)/, \
	smithlab_os.o smithlab_utils.o GenomicRegion.o \
	OptionParser.o MappedRead.o)

preseq: continued_fraction.o load_data_for_complexity.o moment_sequence.o

ifdef HAVE_HTSLIB
preseq to-mr: $(addprefix $(SMITHLAB_CPP)/, htslib_wrapper.o)
endif

%.o: %.cpp %.hpp
	$(CXX) $(CXXFLAGS) -c -o $@ $< $(INCLUDEARGS)

%: %.cpp
	$(CXX) $(CXXFLAGS) -o $@ $^ $(INCLUDEARGS) $(LIBS)

install: $(PROGS)
	@mkdir -p $(PREFIX)/bin
	@install -m 755 $(PROGS) $(PREFIX)/bin

clean:
	@-rm -f $(PROGS) *.o *~
	@-rm -f $(SMITHLAB_CPP)*.o $(SMITHLAB_CPP)*~

.PHONY: clean
