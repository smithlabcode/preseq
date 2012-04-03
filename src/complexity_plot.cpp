/*    complexity_plot: 
 *
 *    Copyright (C) 2010 University of Southern California and
 *                       Andrew D. Smith and Timothy Dailey
 *
 *    Authors: Andrew D. Smith and Timothy Dailey
 *
 *    This program is free software: you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation, either version 3 of the License, or
 *    (at your option) any later version.
 *
 *    This program is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public License
 *    along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <fstream>

#include "OptionParser.hpp"
#include "smithlab_utils.hpp"
#include "GenomicRegion.hpp"
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <iomanip>
#include <numeric>

using std::string;
using std::vector;
using std::ostream;
using std::endl;
using std::cerr;
using std::sort;

static size_t
sample_and_unique(const gsl_rng *rng,
		  const size_t sample_size, 
		  const vector<size_t> &orig,
		  const vector<SimpleGenomicRegion> &regions) {
  vector<size_t> sample(sample_size);
  gsl_ran_choose(rng, (size_t *)&sample.front(), sample_size, 
		 (size_t *)&orig.front(), orig.size(), sizeof(size_t));
  size_t count = 1;
  size_t prev_idx = sample.front();
  for (size_t i = 1; i < sample.size(); ++i) {
    const size_t idx = sample[i];
    count += (regions[idx] != regions[prev_idx]);
    prev_idx = idx;
  }
  return count;
}


int main(int argc, const char **argv) {

  try {
    /* FILES */
    string outfile;
    size_t lower_limit = 1000000;
    size_t upper_limit = 0;
    size_t step_size = 1000000;

    bool VERBOSE = false;


    /****************** GET COMMAND LINE ARGUMENTS ***************************/
    OptionParser opt_parse("complexity_plot", "",
			   "<bed-format-file>");
    opt_parse.add_opt("output", 'o', "Name of output file (default: stdout)", 
		      false , outfile);
    opt_parse.add_opt("lower", 'l', "lower limit for samples", 
		      false , lower_limit);
    opt_parse.add_opt("upper", 'u', "upper limit for samples", 
		      false , upper_limit);
    opt_parse.add_opt("step", 's', "step size for samples", 
		      false , step_size);
    opt_parse.add_opt("verbose", 'v', "print more run information", 
		      false , VERBOSE);
    vector<string> leftover_args;
    opt_parse.parse(argc, argv, leftover_args);
    if (argc == 1 || opt_parse.help_requested()) {
      cerr << opt_parse.help_message() << endl;
      return EXIT_SUCCESS;
    }
    if (opt_parse.about_requested()) {
      cerr << opt_parse.about_message() << endl;
      return EXIT_SUCCESS;
    }
    if (opt_parse.option_missing()) {
      cerr << opt_parse.option_missing_message() << endl;
      return EXIT_SUCCESS;
    }
    if (leftover_args.empty()) {
      cerr << opt_parse.help_message() << endl;
      return EXIT_SUCCESS;
    }
    const string input_file_name = leftover_args.front();
    /**********************************************************************/
    
    // Setup the random number generator
    gsl_rng_env_setup();
    gsl_rng *rng = gsl_rng_alloc(gsl_rng_default);
    srand(time(0) + getpid());
    gsl_rng_set(rng, rand());
    
    if (VERBOSE)
      cerr << "loading mapped locations" << endl;
    vector<SimpleGenomicRegion> regions;
    ReadBEDFile(input_file_name, regions);
    if (!check_sorted(regions))
      throw SMITHLABException("regions not sorted");
    
    vector<size_t> orig(regions.size());
    for (size_t i = 0; i < orig.size(); ++i)
      orig[i] = i;

    if (upper_limit == 0)
      upper_limit = orig.size();

    std::ofstream of;
    if (!outfile.empty()) of.open(outfile.c_str());
    ostream out(outfile.empty() ? std::cout.rdbuf() : of.rdbuf());

    for (size_t i = lower_limit; i <= upper_limit; i += step_size) {
      if (VERBOSE)
	cerr << "sample size: " << i << endl;
      out << i << "\t" << sample_and_unique(rng, i, orig, regions) << endl;
    }
    
  }
  catch (SMITHLABException &e) {
    cerr << "ERROR:\t" << e.what() << endl;
    return EXIT_FAILURE;
  }
  catch (std::bad_alloc &ba) {
    cerr << "ERROR: could not allocate memory" << endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
