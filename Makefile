#    Copyright (C) 2011 University of Southern California and
#                       Andrew D. Smith and Timothy Daley
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


ifndef ROOT
ROOT = $(shell pwd)
endif

ifndef SMITHLAB_CPP
SMITHLAB_CPP=$(ROOT)/smithlab_cpp/
endif


SOURCES = $(wildcard *.cpp)
OBJECTS = $(patsubst %.cpp,%.o,$(SOURCES))
PROGS =  gc_extrap lc_extrap c_curve mincount_c_curve count_extrap count_c_curve mincount_extrap saturation_extrap test_quad_replace
INCLUDEDIRS = $(SMITHLAB_CPP)
INCLUDEARGS = $(addprefix -I,$(INCLUDEDIRS))

LIBS += -lgsl -lgslcblas 

CXX = g++
CXXFLAGS = -Wall -fPIC -fmessage-length=50
OPTFLAGS = -O2
DEBUGFLAGS = -g -lefence -lpthread -L/home/cmb-01/as/andrewds/lib/

ifdef DEBUG
CXXFLAGS += $(DEBUGFLAGS)
endif

ifdef BAMTOOLS_ROOT
INCLUDEDIRS += $(BAMTOOLS_ROOT)/include
LIBS += -L$(BAMTOOLS_ROOT)/lib -lz -lbamtools
CXXFLAGS += -DHAVE_BAMTOOLS
endif

ifdef OPT
CXXFLAGS += $(OPTFLAGS)
endif

all: $(PROGS)

$(PROGS): $(addprefix $(SMITHLAB_CPP)/, GenomicRegion.o smithlab_os.o \
	smithlab_utils.o OptionParser.o MappedRead.o RNG.o)

gc_extrap saturation_extrap lc_extrap count_extrap mincount_extrap: continued_fraction.o

test_quad_replace: library_size_estimates.o newtons_method.o moment_sequence.o ZTNB.o

%.o: %.cpp %.hpp
	$(CXX) $(CXXFLAGS) -c -o $@ $< $(INCLUDEARGS)

%: %.cpp
	$(CXX) $(CXXFLAGS) -o $@ $^ $(INCLUDEARGS) $(LIBS)

install: $(PROGS)
	@mkdir -p $(ROOT)/bin
	@install -m 755 $(PROGS) $(ROOT)/bin

clean:
	@-rm -f $(PROGS) *.o *~

.PHONY: clean
