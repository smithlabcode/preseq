
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

static inline double
movement(const double a, const double b) {
  return fabs((a - b)/max(a, b)); //delta
}

static inline double
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
                   const bool VERBOSE,
                   vector<double> &numerator_approx,
                   vector<double> &denominator_approx,
                   bool &defect_flag){
  size_t numer_size = coeffs.size()-denom_size;  //numer_size = L+1, denom_size = M
  
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
  //  numerator_approx.push_back(evaluate_polynomial(num_vec, t));
  //  denominator_approx.push_back(evaluate_polynomial(full_denom_vec, t));
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
cont_frac_pd(const vector<double> &coeffs,
             const size_t depth,
             vector<double> &cf_coeffs){
  vector< vector<double> > p_vec(2*depth-2, vector<double>(2*depth-2, 0.0));
  p_vec[0][0] = 1;
  for(size_t i = 0; i < p_vec.size(); i++)
    p_vec[i][1] = coeffs[i];
  for(size_t j = 2; j < p_vec[0].size(); j++){
    for(size_t i = 0; i < p_vec.size()-1; i++){
      p_vec[i][j] = p_vec[0][j-1]*p_vec[i+1][j-2] - p_vec[0][j-2]*p_vec[i+1][j-1];
    }
  }
  cerr << "p calculated \n";
  cf_coeffs.push_back(coeffs[0]);
  for(size_t i = 1; i < depth; i++)
    cf_coeffs.push_back(p_vec[0][i+1]/(p_vec[0][i-1]*p_vec[0][i]));
}

double 
compute_cf_approx(const vector<double> &cf_coeffs,
                  const double time){
  if(time == 0.0){
    return 0.0;
  }
  else{
    double return_val = 0.0;
    for(size_t i = cf_coeffs.size()-1; i > 0; i--){
      return_val = cf_coeffs[i]*time/(1.0+return_val);
    }
    return_val = cf_coeffs[0]*time/(1.0+return_val);
    return(return_val);
  }
}

double
compute_cf_approx_euler(const vector<double> &cf_coeffs,  //failure
                        const double time){
  if(time == 0.0){
    return 0.0;
  }
  else{
    double current_num = 0.0;
    double prev_num1 = 0.0;
    double prev_num2 = 1.0;
    double current_denom = 0.0;
    double prev_denom1 = 0.0;
    double prev_denom2 = 1.0; 
    for(size_t i = 0; i < cf_coeffs.size(); i++){
      current_num = prev_num1 + cf_coeffs[i]*time*prev_num2;
      current_denom = prev_denom1 + cf_coeffs[i]*time*prev_denom2;
      prev_num2 = prev_num1;
      prev_num1 = current_num;
      prev_denom2= prev_denom1;
      prev_denom1 = current_denom;
      cerr << current_num << ", " << current_denom << "\n";
    }
    return(current_num/current_denom);
  }
}
    
      
  
      





