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
PROGS = library_complexity complexity_plot extrapolate_library_complexity poisson_generator poisson_mixture_estimation

LIBS = -lgsl -lgslcblas
LIBDIR = $(SMITHLAB_CPP)/
INCLUDEDIR = $(SMITHLAB_CPP)/
CXX = g++
CFLAGS = -Wall -fPIC -fmessage-length=50
OPTFLAGS = -O2
DEBUGFLAGS = -g
COMMON_DIR = $(SMITHLAB_CPP)/

ifeq "$(shell uname)" "Darwin"
CFLAGS += -arch x86_64
endif

ifdef DEBUG
CFLAGS += $(DEBUGFLAGS)
endif

ifdef OPT
CFLAGS += $(OPTFLAGS)
endif



library_complexity complexity_plot extrapolate_library_complexity poisson_generator poisson_mixture_estimation cutBEDfile poisson_mixture_estimation_clusterin2 poisson_estimation_hist poisson_estimation_hist_fake poisson_mix_estimation poisson_estimation_hist_AIC extrapolate_library_complexity_mixture_not_given_bootstrap: GenomicRegion.o rmap_utils.o OptionParser.o

GenomicRegion.o: GenomicRegion.cpp GenomicRegion.hpp
	$(CXX) $(CFLAGS) -c -o $@ $< 

rmap_utils.o: rmap_utils.cpp rmap_utils.hpp
	$(CXX) $(CFLAGS) -c -o $@ $< 

OptionParser.o: OptionParser.cpp OptionParser.hpp
	$(CXX) $(CFLAGS) -c -o $@ $< 


%: %.cpp
	$(CXX) $(CFLAGS) -o $@ $^ -I$(COMMON_DIR) -L$(LIBDIR) $(LIBS)

test_%:	%
	@$(TEST_DIR)/$@ $(TEST_DIR)

test:	$(addprefix test_, $(PROGS))

clean:
	@-rm -f $(PROGS) *.o *.so *.a *~

.PHONY: clean
