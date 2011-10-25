//test code


#ifndef My_PADE_HPP
#define MY_PADE_HPP

#include <fstream>
#include <iomanip>
#include <numeric>
#include <vector>

void compute_pade_coeffs(const std::vector<double> &coeffs, 
                         const size_t numer_size, const size_t denom_size,
                         std::vector<double> &num_coeffs, 
                         std::vector<double> &denom_coeffs);
double compute_pade_approx_numerator(const double time,
                                     const std::vector<double> &numerator);
double compute_pade_approx_denominator(const double time,
                                       const std::vector<double> &denominator);
double locate_polynomial_zero(const std::vector<double> &coeffs,
                              const double lower_lim, const double upper_lim,
                              const double tolerance);
void compute_pade_curve(const std::vector<double> &coeffs,
                        const double max_time,
                        const double time_step,
                        const double allowable_defect_error,
                        const double tolerance,
                        const size_t denom_size,
                        const bool VERBOSE,
                        std::vector<double> &numerator_approx,
                        std::vector<double> &denominator_approx,
                        bool &defect_flag);
void cont_frac_pd(const std::vector<double> &coeffs,
                  const size_t depth,
                  std::vector<double> &cf_coeffs); //implements pd algorithm to compute cont frac coeffs
double compute_cf_approx(const std::vector<double> &cf_coeffs,
                         const double time);
double compute_cf_approx_euler(const std::vector<double> &cf_coeffs,
                               const double time); //euler's recursion for cf's


#endif