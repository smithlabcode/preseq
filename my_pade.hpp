//test code


#ifndef My_PADE_HPP
#define MY_PADE_HPP

#include <fstream>
#include <iomanip>
#include <numeric>
#include <vector>
#include <complex>

class cont_frac {
public:
  cont_frac(const std::vector<double> in_cfs, const std::vector<double> in_offset,
            const size_t in_lower, const size_t in_upper):
    cf_coeffs(in_cfs), offset_coeffs(in_offset), lower_offset(in_lower), 
    upper_offset(in_upper) {;}
  void set_lower_offset(const size_t l_o) {lower_offset = l_o; upper_offset = 0;}
  void set_upper_offset(const size_t u_o) {upper_offset = u_o; lower_offset = 0;}
  void set_offset_coeffs(const std::vector<double> coeffs) {offset_coeffs = coeffs;}
  void set_cf_coeffs(const std::vector<double> coeffs) {cf_coeffs = coeffs;}
  void get_cf_coeffs(std::vector<double> &return_coeffs) {return_coeffs = cf_coeffs;}
  void get_offset_coeffs(std::vector<double> &return_coeffs) {return_coeffs = offset_coeffs;}
    
  void compute_cf_coeffs(const std::vector<double> ps_coeffs, const size_t depth);
  double cf_approx(const double time, const double tolerance);
  double cf_deriv_complex(const double val, const double dx,
                            const double tolerance);
  double locate_zero_cf_deriv(const double val, const double prev_val,
                                const double dx, const double tolerance);

private:
  std::vector<double> cf_coeffs;
  std::vector<double> offset_coeffs;
  size_t lower_offset;
  size_t upper_offset;
};
  

bool compute_pade_coeffs(const std::vector<double> &coeffs, 
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
double evaluate_polynomial(const std::vector<double> &coeffs,
                           const double val);
bool compute_pade_curve(const std::vector<double> &coeffs,
                        const double max_time,
                        const double time_step,
                        const double allowable_defect_error,
                        const double tolerance,
                        const size_t denom_size,
                        const size_t numer_size,
                        const bool VERBOSE,
                        std::vector<double> &numerator_approx,
                        std::vector<double> &denominator_approx,
                        bool &defect_flag);
  
double log_sum_log_vec(const std::vector<double> &vals, size_t limit);
double movement(const double a, const double b) ;



#endif