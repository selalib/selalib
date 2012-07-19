
// This is a symmetric penta-diagonal solver header

#ifndef SYM_PENTA_DIAGONAL_H_
#define SYM_PENTA_DIAGONAL_H_

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#ifdef __cplusplus
extern "C" {
#endif

  typedef struct sym_penta_diagonal_plan
  {
    int n;
    double *e1;
    double *e2;
    double *y;
    double *z;
    double *w;
    double *x;
    double *for_subsys;
  } sym_penta_diagonal_plan;


  double *alloc_double(int n);

  sym_penta_diagonal_plan *new_sym_penta_diagonal(int n);

  double *solve_subsystem(double l1, double l2, double *b, int n);

  double *solve_sym_penta_diagonal(double a, double b, double c, 
	        const double *f, sym_penta_diagonal_plan *plan);

  void dealloc_double(double *ptr);

  void delete_sym_penta_diagonal(sym_penta_diagonal_plan *plan);

#ifdef __cplusplus
}
#endif

#endif /* SYM_PENTA_DIAGONAL_H_ */
