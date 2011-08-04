#    Copyright (C) 2010 University of Southern California and
#                       Andrew D. Smith
#                       Timothy Daley
#
#    Authors: Andrew D. Smith and Timothy Daley
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

ifndef SMITHLAB_CPP
$(error SMITHLAB_CPP variable undefined)
endif

PROGS = complexity_plot etcomplexity

LIBS = -lgsl -lgslcblas
LIBDIR = $(SMITHLAB_CPP)/lib
INCLUDEDIR = $(SMITHLAB_CPP)/
CXX = g++
CFLAGS = -Wall -fPIC -fmessage-length=50
OPTFLAGS = -O2
DEBUGFLAGS = -g
COMMON_DIR = $(SMITHLAB_CPP)/
TEST_DIR = $(SMITHLAB_CPP)/test

ifeq "$(shell uname)" "Darwin"
CFLAGS += -arch x86_64
endif

ifdef DEBUG
CFLAGS += $(DEBUGFLAGS)
endif

ifdef OPT
CFLAGS += $(OPTFLAGS)
endif

all:	$(PROGS)

$(PROGS): \
	$(addprefix $(COMMON_DIR), GenomicRegion.o rmap_os.o \
	rmap_utils.o OptionParser.o MappedRead.o)

%.o: %.cpp %.hpp
	$(CXX) $(CFLAGS) -c -o $@ $< -I$(COMMON_DIR)


install: all
	@mkdir -p $(SMITHLAB_CPP)/bin
	@install -m 755 $(PROGS) $(SMITHLAB_CPP)/bin

%: %.cpp
	$(CXX) $(CFLAGS) -o $@ $^ -I$(COMMON_DIR) $(LIBS)

# test_%:	%
# 	@$(TEST_DIR)/$@ $(TEST_DIR)

test:	$(addprefix test_, $(PROGS))

clean:
	@-rm -f $(PROGS) *.o *.so *.a *~

.PHONY: clean
