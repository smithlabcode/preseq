//test code


#ifndef GSL_PADE_HPP
#define GSL_PADE_HPP

#include <stdio.h>  
#include <stdint.h> 
#include <gsl/gsl_vector_double.h>  
#include <gsl/gsl_matrix_double.h>  
#include <gsl/gsl_linalg.h>  
#include <math.h>  

double gsl_pade_approx(double x, const gsl_vector *a, const gsl_vector *b);  
void gsl_pade_coeff(gsl_vector *a, gsl_vector *b, const gsl_vector *c);  

#endif