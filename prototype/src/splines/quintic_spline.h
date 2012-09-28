//#ifndef QUINTIC_SPLINE_H_
//#define QUINTIC_SPLINE_H_

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "toep_penta_diagonal.h"

typedef struct quintic_spline

{

  int n;
  double xmin, xmax;
  double *coeffs;

} quintic_spline_plan;

  quintic_spline_plan *new_quintic_spline(const int n, 
               const double xmin, const double xmax, const double *f);
  double B(int j, int i, double x, quintic_spline_plan *plan_spline);
  void compute_coeffs(const double *f, quintic_spline_plan *plan_spline);
  // The interpolator spline function with boundary conditions
  double spline(const double x, quintic_spline_plan *plan);
  void delete_quintic_spline(quintic_spline_plan *plan);

