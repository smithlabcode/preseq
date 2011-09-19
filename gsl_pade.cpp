
#include "gsl_pade.hpp"
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

double gsl_pade_approx(double x, const gsl_vector *a, const gsl_vector *b){  
  double num = 0.0, den = 1.0;  
  unsigned int j;  
  for(j = 0; j < a->size; j++)  
    num += gsl_vector_get(a, j) * pow(x, j);  
  for(j = 0; j < b->size; j++)  
    den += gsl_vector_get(b, j)*pow(x, j+1);  
  return num / den;  
}   

void gsl_pade_coeff(gsl_vector *a, gsl_vector *b, const gsl_vector *c){  
  if(a->size + b->size > c->size){  
    fprintf(stderr, "Insufficient data, for Pade's approximation!\n");  
    exit(2);  
  }  
  gsl_vector_set_zero(a);  
  gsl_vector_set_zero(b);  
  unsigned int m = a->size - 1, n = b->size;  
  unsigned int j, k;  
  //////////  
  //Ux = v  
  gsl_matrix *U  
  = gsl_matrix_alloc(n, n);  
  gsl_vector *v  
  = gsl_vector_alloc(n);  
  gsl_permutation *p  
  = gsl_permutation_alloc(n);  
  
  //Indices are not for us puny humans!  
  for(j = m + 1; j < m + n + 1; j++){  
    gsl_vector_set(v, j - m - 1, - gsl_vector_get(c, j));  
    for(k = 1; k <= (j > n ? n : j); k++){  
      gsl_matrix_set(U, j - m - 1, k - 1, gsl_vector_get(c, j - k));  
    }  
  } 
  
  
  gsl_linalg_HH_solve(U, v, b);//Solve Ux = v  
  gsl_matrix_free(U);  
  gsl_permutation_free(p);  
  gsl_vector_free(v);  
  ///////  
  //#define B(x) (x == 0 ? 1.0 : gsl_vector_get(b, x))  
  for(j = 0; j <= m; j++){  
    for(k = 0; k <= (j > n ? n : j); k++)  
      *gsl_vector_ptr(a, j) += (k == 0 ? 1.0 : gsl_vector_get(b, k - 1)) * gsl_vector_get(c, j - k);  
  }    
  
  return;  
}  