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
PROGS = library_complexity complexity_plot extrapolate_library_complexity poisson_generator poisson_mixture_estimation extrapolate_library_complexity_byexpectation test_euler

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

pade_test_fake2 pade_test2 coverage distinct: my_pade.o OptionParser.o rmap_utils.o GenomicRegion.o

etcomplexity: OptionParser.o rmap_utils.o GenomicRegion.o MappedRead.o

pade_test pade_test_fake pade_test_true_hist: gsl_pade.o OptionParser.o rmap_utils.o GenomicRegion.o

negbin_true_hist: OptionParser.o NBD_mixture.o rmap_utils.o

library_complexity complexity_plot extrapolate_library_complexity poisson_generator poisson_mixture_estimation cutBEDfile poisson_mixture_estimation_clusterin2 poisson_estimation_hist poisson_estimation_hist_fake poisson_mix_estimation poisson_estimation_hist_AIC extrapolate_library_complexity_mixture_not_given BEDsample combine_mixture extrapolate_library_complexity_byexpectation poisson_complexity fit_trunc_negbin negbin_generator poisson_estimation_hist_bootstrap negbin_estimation_AIC negbin_estimation_AIC_fake negbin_estimation_timer poisson_estimation_timer extrapolate_library_complexity_negbin: GenomicRegion.o rmap_utils.o OptionParser.o

poisson_estimation_extrapolation poisson_estimation_extrapolation_fake poiss_genome_size_simul poiss_extrap_simul poisson_generator: GenomicRegion.o rmap_utils.o OptionParser.o Poiss_mixture.o

negbin_estimation_extrapolation negbin_estimation_extrapolation_fake fit_ztnbd_output_val_loglike negbin_genome_size_simul: GenomicRegion.o rmap_utils.o OptionParser.o NBD_mixture.o

test_euler: OptionParser.o GenomicRegion.o rmap_utils.o NBD_mixture.o euler_series_transform.o

test_euler3: gsl_pade.o NBD_mixture.o OptionParser.o

NBD_mixture.o:NBD_mixture.cpp NBD_mixture.hpp
	$(CXX) $(CFLAGS) -c -o $@ $<

euler_series_transform.o: euler_series_transform.cpp euler_series_transform.hpp
	$(CXX) $(CFLAGS) -c -o $@ $<

Poiss_mixture.o:Poiss_mixture.cpp Poiss_mixture.hpp
	$(CXX) $(CFLAGS) -c -o $@ $<

GenomicRegion.o: GenomicRegion.cpp GenomicRegion.hpp
	$(CXX) $(CFLAGS) -c -o $@ $< 

rmap_utils.o: rmap_utils.cpp rmap_utils.hpp
	$(CXX) $(CFLAGS) -c -o $@ $< 

gsl_pade.o: gsl_pade.cpp gsl_pade.hpp
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
