//#ifndef QUINTIC_SPLINE_H_
//#define QUINTIC_SPLINE_H_

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "sym_penta_diagonal.h"

typedef struct quintic_spline_hermite

{

  int n;
  double xmin, xmax;
  double *coeffs;

} quintic_spline_hermite_plan;

  quintic_spline_hermite_plan *new_quintic_spline_hermite(const int n, 
               const double xmin, const double xmax, const double *f);
  double B(int j, int i, double x, quintic_spline_hermite_plan *plan_spline);
  void compute_coeffs_hermite(const double *f, quintic_spline_hermite_plan *plan_spline);
  // The interpolator spline function with Hermite boundary conditions
  double spline_hermite(const double x, quintic_spline_hermite_plan *plan);
  void delete_quintic_spline_hermite(quintic_spline_hermite_plan *plan);

