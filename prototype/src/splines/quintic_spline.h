//#ifndef QUINTIC_SPLINE_H_
//#define QUINTIC_SPLINE_H_

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "sym_penta_diagonal.c"

typedef struct quintic_spline

{

  int n;
  double xmin, xmax;
  double *p;

} quintic_spline_plan;

  quintic_spline_plan *new_quintic_spline(const int n, 
                const double xmin, const double xmax);
  double B_spline(const double x); 
  // Quintic B-spline 
  double *coeffs(const double *f);
  // Coefficients of the spline
  double spline(const double x, const double *f); 
  // The interpolator spline function
  quintic_spline_plan *delete_quintic_spline();

