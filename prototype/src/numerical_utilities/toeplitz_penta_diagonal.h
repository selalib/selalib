
// This is a Toeplitz penta-diagonal solver header

#ifndef toep_PENTA_DIAGONAL_H_
#define toep_PENTA_DIAGONAL_H_

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#ifdef __cplusplus
extern "C" {
#endif

  typedef struct toep_penta_diagonal_plan
  {
    int n;
    double *e1;
    double *e2;
    double *y;
    double *z;
    double *w;
    double *x;
    double *for_subsys;
  } toep_penta_diagonal_plan;


  double *alloc_double(int n);

  toep_penta_diagonal_plan *new_toep_penta_diagonal(int n);

  double *solve_subsystem(double l1, double l2, double *b, int n);

  double *solve_toep_penta_diagonal(double a, double b, double c, 
	        const double *f, toep_penta_diagonal_plan *plan);

  void dealloc_double(double *ptr);

  void delete_toep_penta_diagonal(toep_penta_diagonal_plan *plan);

#ifdef __cplusplus
}
#endif

#endif /* toep_PENTA_DIAGONAL_H_ */
