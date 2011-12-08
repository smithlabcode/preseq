
#include "my_pade.hpp"
#include "rmap_utils.hpp"


#include <gsl/gsl_vector_double.h>  
#include <gsl/gsl_matrix_double.h>  
#include <gsl/gsl_linalg.h>  
#include <math.h>  

#include <fstream>
#include <iomanip>
#include <numeric>
#include <vector>
#include <complex>

using std::string;
using std::vector;
using std::ostream;
using std::endl;
using std::cerr;
using std::pair;
using std::make_pair;
using std::sort;

using std::max;
using std::setw;
using std::fabs;
using std::ceil;
using std::greater;
using std::numeric_limits;
using std::complex;
using std::norm;
using std::real;
using std::imag;

void
solve_linear_system(const vector<vector<double> > &U, 
                    const vector<double> &v, vector<double> &b) {
  
  // compute Ub = v
  
  gsl_vector *v_gsl = gsl_vector_alloc(v.size());
  copy(v.begin(), v.end(), v_gsl->data);
  
  gsl_vector *b_gsl = gsl_vector_alloc(U.front().size());
  gsl_vector_set_zero(b_gsl);
  
  gsl_matrix *U_gsl = gsl_matrix_alloc(U.size(), U.front().size());
  for (size_t i = 0; i < U.size(); ++i)
    for (size_t j = 0; j < U[i].size(); ++j)
      gsl_matrix_set(U_gsl, i, j, U[i][j]);
  
  gsl_linalg_HH_solve(U_gsl, v_gsl, b_gsl);
  gsl_matrix_free(U_gsl);
  gsl_vector_free(v_gsl);
  
  copy(b_gsl->data, b_gsl->data + b_gsl->size, back_inserter(b));
  gsl_vector_free(b_gsl);
}

void
compute_denom_ceoffs(const vector<double> &coeffs, const size_t numer_size,
                     const size_t denom_size,
                     vector<double> &denom_coeffs) {
  
  //use the formula from Baker, Pade Approximants page 2 and 3
  vector<double> v(denom_size, 0.0);
  for (size_t j = 0; j < denom_size; j++){
    v[j] = -coeffs[j + numer_size];
  }
  
  size_t n = denom_size;
  size_t m = numer_size-1;
  vector<vector<double> > U(n, vector<double>(n, 0.0));
  for (size_t j = m + 1; j < m + n + 1; j++){
    for (size_t k = 1; k <= (j > n ? n : j); k++){
        U[j - m - 1][n-k] = coeffs[j - k];
    }
  }

  solve_linear_system(U, v, denom_coeffs);

}

void
test_coefficients(const vector<double> &coeffs, 
                  const vector<double> &num_coeffs, 
                  const vector<double> &denom_coeffs) {
  
  static const double COEFF_TEST_TOLERANCE = 1e-10;
  
  for (size_t i = 0; i < num_coeffs.size(); ++i) {
    double sum = coeffs[i];
    size_t upper_lim = std::min(num_coeffs.size(), std::min(i, denom_coeffs.size()));
    for (size_t j = 0; j < upper_lim; ++j)
      sum += denom_coeffs[denom_coeffs.size()-1-j]*coeffs[i - j - 1];
    assert(num_coeffs[i] == sum);
  }
  for (size_t i = 0; i < denom_coeffs.size(); ++i) {
    const size_t offset = num_coeffs.size() + i;
    double sum = coeffs[offset];
    for (size_t j = 0; j < denom_coeffs.size(); j++){
      if(offset-denom_coeffs.size()+j >= 0)
        sum += coeffs[offset-denom_coeffs.size()+j]*denom_coeffs[j];  // if the coeffs index < 0, coeff = 0
    }

    assert(fabs(sum) < COEFF_TEST_TOLERANCE);
  }
}

void
compute_pade_coeffs(const vector<double> &coeffs,
                    const size_t numer_size, const size_t denom_size, //numer_size = L+1
                    vector<double> &num_coeffs, 
                    vector<double> &denom_coeffs) {

  compute_denom_ceoffs(coeffs, numer_size, denom_size, denom_coeffs);  

  for (size_t i = 0; i < numer_size; ++i) {
    num_coeffs.push_back(coeffs[i]);
    size_t upper_lim = std::min(numer_size, denom_coeffs.size());
    for (size_t j = 0; j < std::min(i, upper_lim); ++j){
      num_coeffs[i] += denom_coeffs[denom_coeffs.size()-1-j]*coeffs[i - j-1];
    }
  }
  
  
  test_coefficients(coeffs, num_coeffs, denom_coeffs);
  //denom is backwards, rearrange it
  vector<double> denom_coeffs_copy(denom_coeffs);
  for(size_t j = 0; j < denom_coeffs.size(); j++)
    denom_coeffs[j] = denom_coeffs_copy[denom_coeffs_copy.size()-1-j];
  
}

double
compute_pade_approx_numerator(const double t,
                              const vector<double> &numerator) {

  double total = 0.0; 
  for (size_t j = 0; j < numerator.size(); ++j) {
    total += numerator[j]*pow(t, j);
  }
  return total;
}

double
compute_pade_approx_denominator(const double t,
                                const vector<double> &denominator) {
  //   for(j = 0; j <= b->size; j++)  
  //     den += (j == 0 ? 1.0 : gsl_vector_get(b, j - 1)) * pow(x, j);  
  double total = 1.0;
  for(size_t j = 0; j < denominator.size(); ++j)
    total += denominator[j]*pow(t, j + 1);
  return total;
}

double
movement(const double a, const double b) {
  return fabs((a - b)/max(a, b)); //delta
}

double
evaluate_polynomial(const vector<double> &coeffs,
                    const double val){
  double return_val = coeffs[0];
  for(size_t i = 1; i < coeffs.size(); i++)
    return_val += coeffs[i]*pow(val, i);
  return(return_val);
}

double
locate_polynomial_zero(const vector<double> &coeffs, //by bisection
                       const double lower_lim,
                       const double upper_lim,
                       const double tolerance){
  double z_low = lower_lim;
  double z_high = upper_lim;
  double z_mid = upper_lim;
  double mid_val;
  
  double lower_lim_sign = (evaluate_polynomial(coeffs, lower_lim) >= 0 ? 1.0 : -1.0);
  
  double diff = numeric_limits<double>::max();
  double prev_val = numeric_limits<double>::max();
  
  while(diff > tolerance && movement(z_high, z_low) > tolerance){
    z_mid = (z_high + z_low)/2.0;
    mid_val = evaluate_polynomial(coeffs, z_mid);
    if(mid_val*lower_lim_sign >= 0) z_low = z_mid;  //if f(z_mid) is same sign as f(z_low)
    else z_high = z_mid;
    diff = fabs((prev_val - mid_val)/prev_val);
    prev_val = mid_val;
  }
  return(z_mid);
}

void
compute_pade_curve(const std::vector<double> &coeffs,
                   const double max_time,
                   const double time_step,
                   const double allowable_defect_error,
                   const double tolerance,
                   const size_t denom_size,
                   const size_t numer_size,
                   const bool VERBOSE,
                   vector<double> &numerator_approx,
                   vector<double> &denominator_approx,
                   bool &defect_flag){  
  vector<double> denom_vec;
  vector<double> num_vec;
    
  compute_pade_coeffs(coeffs, numer_size, denom_size, num_vec, denom_vec); 
  
  vector<double> full_denom_vec(denom_vec);
  full_denom_vec.insert(full_denom_vec.begin(), 1.0);
  double t = 0.0;
  
  double prev_denom_val = 1.0;
  double current_denom_val = 1.0;
  double zero_location = 0.0; 
  
  numerator_approx.clear();
  denominator_approx.clear();

  while(t <= max_time){
    numerator_approx.push_back(compute_pade_approx_numerator(t, num_vec));
    denominator_approx.push_back(compute_pade_approx_denominator(t, denom_vec));
    current_denom_val = denominator_approx.back();
    if(current_denom_val*prev_denom_val < 0){
      double denom_zero = locate_polynomial_zero(full_denom_vec, t-time_step, t, tolerance);
      double numer_zero = locate_polynomial_zero(num_vec, t-time_step, t+time_step, tolerance);
      if(VERBOSE){
        cerr << "zero found, denom location = " << denom_zero << ", numerator location = "
        << numer_zero << "\n";
      }
      if(fabs(denom_zero - numer_zero) > allowable_defect_error)
        defect_flag = true;
    }
    prev_denom_val = current_denom_val;
    t += time_step;
  }
  if(defect_flag == true && VERBOSE)
    cerr << "defect found \n";
}


void 
cont_frac_pd(const vector<double> &coeffs, //product difference
             const size_t depth,
             vector<double> &cf_coeffs){
  vector< vector<double> > p_table(2*depth-2, vector<double>(2*depth-2, 0.0));
  p_table[0][0] = 1;
  for(size_t i = 0; i < p_table.size(); i++)
    p_table[i][1] = coeffs[i];
  for(size_t j = 2; j < p_table[0].size(); j++){
    for(size_t i = 0; i < p_table.size()-1; i++){
      p_table[i][j] = p_table[0][j-1]*p_table[i+1][j-2] - p_table[0][j-2]*p_table[i+1][j-1];
    }
  }
  cf_coeffs.push_back(coeffs[0]);
  for(size_t i = 1; i < depth; i++)
    cf_coeffs.push_back(p_table[0][i+1]/(p_table[0][i-1]*p_table[0][i]));
}

void
cont_frac_qd(const vector<double> &coeffs, //quotient-difference
             const size_t depth,  // odd depth gives M,M pade even depth gives M-1,M pade
             vector<double> &cf_coeffs){
  vector< vector<double> > q_table(depth, vector<double>(depth, 0.0));
  vector< vector<double> > e_table(depth, vector<double>(depth, 0.0));
  for(size_t i = 0; i < q_table[1].size(); i++)
    q_table[1][i] = coeffs[i+1]/coeffs[i];
  
  for(size_t j = 0; j < depth-1; j++)
    e_table[1][j] = q_table[1][j+1] - q_table[1][j] + e_table[0][j+1];
  
  for(size_t i = 2; i < depth; i++){
    for(size_t j = 0; j < depth; j++)
      q_table[i][j] = q_table[i-1][j+1]*e_table[i-1][j+1]/e_table[i-1][j];
      
    for(size_t j = 0; j < depth; j++)
      e_table[i][j] = q_table[i][j+1] - q_table[i][j] + e_table[i-1][j+1];
  }
  cf_coeffs.push_back(coeffs[0]);
  for(size_t i = 1; i < depth; i++){
    if(i % 2 == 0)
      cf_coeffs.push_back(-e_table[i/2][0]);
    else
      cf_coeffs.push_back(-q_table[(i+1)/2][0]);
  }
}

void
cont_frac_upper_offdiagonal(const vector<double> &coeffs,
                            const size_t depth,
                            const size_t offset,
                            vector<double> &offset_cf_coeffs,
                            vector<double> &cf_coeffs){ //first offset coefficients set to coeffs
  vector<double> holding_coeffs;
  for(size_t i = offset; i < depth; i++)
    holding_coeffs.push_back(coeffs[i]);
  cont_frac_qd(holding_coeffs, depth-offset, cf_coeffs);
  for(size_t i =0; i < offset; i++)
    offset_cf_coeffs.push_back(coeffs[i]);
}

double
compute_upper_offdiag_cf_approx(const vector<double> &cf_coeffs,
                                const vector<double> &offset_coeffs,
                                const double time,
                                const double tolerance){
  if(time == 0)
    return 0.0;
  else{
    double current_num = 0.0;
    double prev_num1 = cf_coeffs[0];
    double prev_num2 = 0.0;
    double current_denom = 0.0;
    double prev_denom1 = 1.0;
    double prev_denom2 = 1.0; 
    for(size_t i = 1; i < cf_coeffs.size(); i++){
      current_num = prev_num1 + cf_coeffs[i]*time*prev_num2;
      current_denom = prev_denom1 + cf_coeffs[i]*time*prev_denom2;
      prev_num2 = prev_num1;
      prev_num1 = current_num;
      prev_denom2= prev_denom1;
      prev_denom1 = current_denom;
      if(i % 10 == 0){ //rescale every 10th iter
        double rescale_val = fabs(current_num) + fabs(current_denom);
        if(rescale_val > 1/tolerance)
          rescale_val = 1/rescale_val;
        else if(rescale_val < tolerance)
          rescale_val = 1/rescale_val;
        else
          rescale_val = 1.0;
        current_num = current_num*rescale_val;
        current_denom = current_denom*rescale_val;
        prev_num1 = prev_num1*rescale_val;
        prev_num2 = prev_num2*rescale_val;
        prev_denom1 = prev_denom1*rescale_val;
        prev_denom2 = prev_denom2*rescale_val;
      }
    }
    double offset_terms = 0.0;
    for(size_t i = 0; i < offset_coeffs.size(); i++)
      offset_terms += offset_coeffs[i]*pow(time, i);
    return(time*(offset_terms + pow(time, offset_coeffs.size())*current_num/current_denom));
  }
}   

void
cont_frac_lower_offdiagonal(const vector<double> &coeffs,
                            const size_t depth,
                            const size_t offset,
                            vector<double> &offset_cf_coeffs,
                            vector<double> &cf_coeffs){  //remember to invert the formulas.
  //need to work with reciprocal series g = f^{-1}
  vector<double> reciprocal_coeffs;
  reciprocal_coeffs.push_back(1/coeffs[0]);
  for(size_t i = 1; i < depth; i++){
    double holding_val = 0.0;
    for(size_t j = 0; j < i; j++)
      holding_val += coeffs[i-j]*reciprocal_coeffs[j];

    reciprocal_coeffs.push_back(-holding_val/coeffs[0]);
  }
  vector<double> holding_coeffs;
  for(size_t i = offset; i < depth; i++)
    holding_coeffs.push_back(reciprocal_coeffs[i]);
  cont_frac_qd(holding_coeffs, depth-offset, cf_coeffs);
  for(size_t i = 0; i < offset; i++)
    offset_cf_coeffs.push_back(reciprocal_coeffs[i]);
}

double 
compute_lower_offdiag_cf_approx(const vector<double> &cf_coeffs,
                                const vector<double> &offset_coeffs,
                                const double time,
                                const double tolerance){
  if(time == 0)
    return 0.0;
  else{
    double current_num = 0.0;
    double prev_num1 = cf_coeffs[0];
    double prev_num2 = 0.0;
    double current_denom = 0.0;
    double prev_denom1 = 1.0;
    double prev_denom2 = 1.0; 
    for(size_t i = 1; i < cf_coeffs.size(); i++){
      current_num = prev_num1 + cf_coeffs[i]*time*prev_num2;
      current_denom = prev_denom1 + cf_coeffs[i]*time*prev_denom2;
      prev_num2 = prev_num1;
      prev_num1 = current_num;
      prev_denom2= prev_denom1;
      prev_denom1 = current_denom;
      if(i % 10 == 0){ //rescale every 10th iter
        double rescale_val = fabs(current_num) + fabs(current_denom);
        if(rescale_val > 1/tolerance)
          rescale_val = 1/rescale_val;
        else if(rescale_val < tolerance)
          rescale_val = 1/rescale_val;
        else
          rescale_val = 1.0;
        current_num = current_num*rescale_val;
        current_denom = current_denom*rescale_val;
        prev_num1 = prev_num1*rescale_val;
        prev_num2 = prev_num2*rescale_val;
        prev_denom1 = prev_denom1*rescale_val;
        prev_denom2 = prev_denom2*rescale_val;
      }
    }
    double offset_terms = 0.0;
    for(size_t i = 0; i < offset_coeffs.size(); i++)
      offset_terms += offset_coeffs[i]*pow(time, i);
    return(time/(offset_terms + pow(time, offset_coeffs.size())*current_num/current_denom));
  }
}
  
  
/*
double 
compute_cf_approx(const vector<double> &cf_coeffs,
                  const double time){
  if(time == 0.0){
    return 0.0;
  }
  else{
    double return_val = time*cf_coeffs.back();
    for(size_t i = cf_coeffs.size()-2; i > 0; i--){
      return_val = cf_coeffs[i]*time/(1.0+return_val);
    }
    return_val = cf_coeffs[0]*time/(1.0+return_val);
    return(return_val);
  }
}

double
compute_cf_approx_euler(const vector<double> &cf_coeffs, //uses euler's recursion
                        const double time){
  if(time == 0.0){
    return 0.0;
  }
  else{
    double current_num = 0.0;
    double prev_num1 = cf_coeffs[0];
    double prev_num2 = 0.0;
    double current_denom = 0.0;
    double prev_denom1 = 1.0;
    double prev_denom2 = 1.0; 
    for(size_t i = 1; i < cf_coeffs.size(); i++){
      current_num = prev_num1 + cf_coeffs[i]*time*prev_num2;
      current_denom = prev_denom1 + cf_coeffs[i]*time*prev_denom2;
      prev_num2 = prev_num1;
      prev_num1 = current_num;
      prev_denom2= prev_denom1;
      prev_denom1 = current_denom;
    }
    return(time*current_num/current_denom);
  }
}
 */

double
compute_cf_approx_euler(const vector<double> &cf_coeffs, //uses euler's recursion
                        const double time,
                        const double tolerance){
  if(time == 0.0){
    return 0.0;
  }
  else{
    double current_num = 0.0;
    double prev_num1 = cf_coeffs[0];
    double prev_num2 = 0.0;
    double current_denom = 0.0;
    double prev_denom1 = 1.0;
    double prev_denom2 = 1.0; 
    for(size_t i = 1; i < cf_coeffs.size(); i++){
      current_num = prev_num1 + cf_coeffs[i]*time*prev_num2;
      current_denom = prev_denom1 + cf_coeffs[i]*time*prev_denom2;
      prev_num2 = prev_num1;
      prev_num1 = current_num;
      prev_denom2= prev_denom1;
      prev_denom1 = current_denom;
      if(i % 10 == 0){ //rescale every 10th iter
        double rescale_val = fabs(current_num) + fabs(current_denom);
        if(rescale_val > 1/tolerance)
          rescale_val = 1/rescale_val;
        else if(rescale_val < tolerance)
          rescale_val = 1/rescale_val;
        else
          rescale_val = 1.0;
        current_num = current_num*rescale_val;
        current_denom = current_denom*rescale_val;
        prev_num1 = prev_num1*rescale_val;
        prev_num2 = prev_num2*rescale_val;
        prev_denom1 = prev_denom1*rescale_val;
        prev_denom2 = prev_denom2*rescale_val;
      }
    }
    return(time*current_num/current_denom);
  }
}

double 
log_sum_log_vec(const vector<double> &vals, size_t limit){
  const size_t max_idx = 
  max_element(vals.begin(), vals.begin() + limit) - 
  vals.begin();
  const double max_val = vals[max_idx];
  double sum = 1.0;
  for (size_t i = 0; i < limit; ++i) {
    if (i != max_idx) {
      sum += exp(vals[i] - max_val);
#ifdef DEBUG
      assert(finite(sum)); 
      // abort if the sum is infinte //
#endif
    }
  }
  return(max_val + log(sum));
}

void
cont_frac::compute_cf_coeffs(const vector<double> ps_coeffs,
                             const size_t depth){
  cf_coeffs.clear();
  offset_coeffs.clear();
  if(upper_offset != 0 && lower_offset != 0)
    cerr << "at least one offset must be zero, reset offset.\n";
  else if(upper_offset == 0 && lower_offset == 0){    
    vector<double> temp_cf_coeffs;
    cont_frac_qd(ps_coeffs, depth, temp_cf_coeffs);
    set_cf_coeffs(temp_cf_coeffs);
  }
  else if(upper_offset > 0){
    vector<double> temp_cf_coeffs;
    vector<double> temp_offset_coeffs;
    cont_frac_upper_offdiagonal(ps_coeffs, depth, upper_offset, 
                                temp_cf_coeffs, temp_offset_coeffs);
    set_offset_coeffs(temp_offset_coeffs);
    set_cf_coeffs(temp_cf_coeffs);
  }
  else if(lower_offset > 0){
    vector<double> temp_cf_coeffs;
    vector<double> temp_offset_coeffs;
    cont_frac_lower_offdiagonal(ps_coeffs, depth, lower_offset,
                                temp_offset_coeffs, temp_cf_coeffs);
    set_offset_coeffs(temp_offset_coeffs);
    set_cf_coeffs(temp_cf_coeffs);
  }
}

double
cont_frac::cf_approx(const double time, const double tolerance){
  if(upper_offset != 0 && lower_offset != 0){
    cerr << "at least one offset must be zero, reset offset.\n";
    return(0.0);
  }
  else if(upper_offset == 0 && lower_offset == 0)
    return(compute_cf_approx_euler(cf_coeffs, time, tolerance));
  else if(upper_offset > 0)
    return(compute_upper_offdiag_cf_approx(cf_coeffs, offset_coeffs,
                                           time, tolerance));
  else if(lower_offset > 0)
    return(compute_lower_offdiag_cf_approx(cf_coeffs, offset_coeffs,
                                           time, tolerance));
  else
    return 0.0;
}

static void
cf_approx_euler_complex(const vector<double> &cf_coeffs,
                        const complex<double> perturbed_val,
                        const double tolerance,
                        complex<double> &approx){
  const complex<double> i(0.0,1.0);
  if(norm(perturbed_val) == 0.0)
    approx = 0.0*i;
  else{
    complex<double> current_num(0.0, 0.0);
    complex<double> prev_num1(cf_coeffs[0], 0.0);
    complex<double> prev_num2(0.0, 0.0);
    complex<double> current_denom(0.0, 0.0);
    complex<double> prev_denom1(1.0, 0.0);
    complex<double> prev_denom2(1.0, 0.0);
    for(size_t j = 1; j < cf_coeffs.size(); j++){
      complex<double> coeff(cf_coeffs[j], 0.0);
      current_num = prev_num1 + coeff*perturbed_val*prev_num2;
      current_denom = prev_denom1 + coeff*perturbed_val*prev_denom2;
      prev_num2 = prev_num1;
      prev_num1 = current_num;
      prev_denom2 = prev_denom1;
      prev_denom1 = current_denom;
      if(j % 10 == 0){ //rescale every 10th iter
        double rescale_val = norm(current_num) + norm(current_denom);
        if(rescale_val > 1/tolerance)
          rescale_val = 1/rescale_val;
        else if(rescale_val < tolerance)
          rescale_val = 1/rescale_val;
        else
          rescale_val = 1.0;
        current_num = current_num*rescale_val;
        current_denom = current_denom*rescale_val;
        prev_num1 = prev_num1*rescale_val;
        prev_num2 = prev_num2*rescale_val;
        prev_denom1 = prev_denom1*rescale_val;
        prev_denom2 = prev_denom2*rescale_val;
      }
    }
    approx = perturbed_val*current_num/current_denom;
  }
}

double
compute_upper_offdiag_cf_approx_complex(const vector<double> &cf_coeffs,
                                        const vector<double> &offset_coeffs,
                                        const complex<double> perturbed_val,
                                        const double tolerance,
                                        complex<double> &approx){
  const complex<double> i(0.0,1.0);
  if(norm(perturbed_val) == 0.0)
    approx = 0.0*i;
  else{
    complex<double> current_num(0.0, 0.0);
    complex<double> prev_num1(cf_coeffs[0], 0.0);
    complex<double> prev_num2(0.0, 0.0);
    complex<double> current_denom(0.0, 0.0);
    complex<double> prev_denom1(1.0, 0.0);
    complex<double> prev_denom2(1.0, 0.0);
    for(size_t j = 1; j < cf_coeffs.size(); j++){
      complex<double> coeff(cf_coeffs[j], 0.0);
      current_num = prev_num1 + coeff*perturbed_val*prev_num2;
      current_denom = prev_denom1 + coeff*perturbed_val*prev_denom2;
      prev_num2 = prev_num1;
      prev_num1 = current_num;
      prev_denom2 = prev_denom1;
      prev_denom1 = current_denom;
      if(j % 10 == 0){ //rescale every 10th iter
        double rescale_val = norm(current_num) + norm(current_denom);
        if(rescale_val > 1/tolerance)
          rescale_val = 1/rescale_val;
        else if(rescale_val < tolerance)
          rescale_val = 1/rescale_val;
        else
          rescale_val = 1.0;
        current_num = current_num*rescale_val;
        current_denom = current_denom*rescale_val;
        prev_num1 = prev_num1*rescale_val;
        prev_num2 = prev_num2*rescale_val;
        prev_denom1 = prev_denom1*rescale_val;
        prev_denom2 = prev_denom2*rescale_val;
      }
    }
    complex<double> offset_terms(0.0, 0.0);
    for(size_t i = 0; i < offset_coeffs.size(); i++)
      offset_terms += offset_coeffs[i]*pow(perturbed_val, i);
    approx = perturbed_val*(offset_terms + pow(perturbed_val, offset_coeffs.size())*current_num/current_denom);
  }
}  

double
compute_lower_offdiag_cf_approx_complex(const vector<double> &cf_coeffs,
                                        const vector<double> &offset_coeffs,
                                        const complex<double> perturbed_val,
                                        const double tolerance,
                                        complex<double> &approx){
  const complex<double> i(0.0,1.0);
  if(norm(perturbed_val) == 0.0)
    approx = 0.0*i;
  else{
    complex<double> current_num(0.0, 0.0);
    complex<double> prev_num1(cf_coeffs[0], 0.0);
    complex<double> prev_num2(0.0, 0.0);
    complex<double> current_denom(0.0, 0.0);
    complex<double> prev_denom1(1.0, 0.0);
    complex<double> prev_denom2(1.0, 0.0);
    for(size_t j = 1; j < cf_coeffs.size(); j++){
      complex<double> coeff(cf_coeffs[j], 0.0);
      current_num = prev_num1 + coeff*perturbed_val*prev_num2;
      current_denom = prev_denom1 + coeff*perturbed_val*prev_denom2;
      prev_num2 = prev_num1;
      prev_num1 = current_num;
      prev_denom2 = prev_denom1;
      prev_denom1 = current_denom;
      if(j % 10 == 0){ //rescale every 10th iter
        double rescale_val = norm(current_num) + norm(current_denom);
        if(rescale_val > 1/tolerance)
          rescale_val = 1/rescale_val;
        else if(rescale_val < tolerance)
          rescale_val = 1/rescale_val;
        else
          rescale_val = 1.0;
        current_num = current_num*rescale_val;
        current_denom = current_denom*rescale_val;
        prev_num1 = prev_num1*rescale_val;
        prev_num2 = prev_num2*rescale_val;
        prev_denom1 = prev_denom1*rescale_val;
        prev_denom2 = prev_denom2*rescale_val;
      }
    }
    complex<double> offset_terms(0.0, 0.0);
    for(size_t i = 0; i < offset_coeffs.size(); i++)
      offset_terms += offset_coeffs[i]*pow(perturbed_val, i);
    approx = perturbed_val/(offset_terms + pow(perturbed_val, offset_coeffs.size())*current_num/current_denom);
  }
}  
    
double
cont_frac::cf_deriv_complex(const double val,
                            const double dx,
                            const double tolerance){
  const complex<double> i(0.0,1.0);
  complex<double> df(0.0, 0.0);
  complex<double> value(val, 0.0);
  if(upper_offset != 0 && lower_offset != 0)
    cerr << "at least one offset must be zero, reset offset.\n";
  if(upper_offset == 0 && lower_offset == 0)
    cf_approx_euler_complex(cf_coeffs, value + dx*i, tolerance, df);
  else if(upper_offset > 0)
    compute_upper_offdiag_cf_approx_complex(cf_coeffs, offset_coeffs,
                                            value + dx*i, tolerance, df);
  else if(lower_offset > 0)
    compute_lower_offdiag_cf_approx_complex(cf_coeffs, offset_coeffs,
                                            value + dx*i, tolerance, df);
  return(imag(df)/dx);
}
      
double
cont_frac::locate_zero_cf_deriv(const double val, const double prev_val,
                                const double dx, const double tolerance){
  double val_low = prev_val;
  double deriv_low = cf_deriv_complex(val_low, dx, tolerance);
  double val_high = val;
  double deriv_high = cf_deriv_complex(val_high, dx, tolerance);
  double val_mid = val;
  double mid_deriv = 0.0;
  double diff = numeric_limits<double>::max();
  double prev_deriv = numeric_limits<double>::max();
  while(diff > tolerance && movement(val_low, val_high) > tolerance){
    val_mid = (val_low + val_high)/2.0;
    mid_deriv = cf_deriv_complex(val_mid, dx, tolerance);
    if(mid_deriv*deriv_low  > 0)
      val_low = val_mid;
    else 
      val_high = val_mid;
    deriv_low = cf_deriv_complex(val_low, dx, tolerance);
    deriv_high = cf_deriv_complex(val_high, dx, tolerance);
    diff = fabs((prev_deriv - mid_deriv)/prev_deriv);
    prev_deriv = mid_deriv;
  }
  return(val_mid);
}




