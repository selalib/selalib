
// This is a symmetric penta-diagonal solver module

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "sym_penta_diagonal.h"


  double *alloc_double(int n)

  {

    double *ptr = malloc( n*sizeof(double) );

    if ( ptr == NULL )
    {
      fprintf(stderr,"Impossible to allocate pointer\n");
      exit(EXIT_FAILURE);
    }

    return ptr;

  }


  sym_penta_diagonal_plan *new_sym_penta_diagonal(int n)

  {
    // Declaration
    sym_penta_diagonal_plan *plan;

    if (n<3)
    {
      printf("Matrix size must be at least 3x3\n");
      printf("Exiting...\n");
      exit(0);
    }

    // Plan allocation
    plan = malloc(sizeof(sym_penta_diagonal_plan));

    if ( plan == NULL )
    {
      fprintf(stderr,"Impossible to allocate pointer\n");
      exit(EXIT_FAILURE);
    }

    // Plan components allocation
    plan->n  = n; 
    plan->e1 = alloc_double(n);
    plan->e2 = alloc_double(n);
    plan->y  = alloc_double(n);
    plan->z  = alloc_double(n);
    plan->w  = alloc_double(n);
    plan->x  = alloc_double(n);
    plan->for_subsys  = alloc_double(plan->n);//plan->x=malloc( n*sizeof(double) );

    return plan;

  }


  void solve_subsystem(double l1, double l2, double *b,
                       int n, double *tmp, double* y)

  {

    int i;
     
    tmp[0] = b[0];
    tmp[1] = b[1] - l1*tmp[0];
    for (i = 2; i < n; i++)
      tmp[i] = b[i] - ( l2*tmp[i-2] + l1*tmp[i-1] );

    y[n-1] = tmp[n-1];
    y[n-2] = tmp[n-2] - l1*y[n-1];
    for (i = n-3; i >= 0; i--)
      y[i] = tmp[i] - ( l1*y[i+1] + l2*y[i+2] );

  }


  double *solve_sym_penta_diagonal(double a, double b, double c, 
  		 const double *f, sym_penta_diagonal_plan *plan)

  {

    double s, t, p, l1, l2;
    double d, d1, d2;
    double *e1, *e2;
    double *y, *z, *w;
    double *x, *tmp;
    int n, i;
    int sign_of_a;

    if ( fabs(a) <= 2*(fabs(b)+fabs(c)) )
    {
      printf("a, b, and c must be such that: |a| > 2(|b|+|c|)\n");
      printf("Exiting...\n");
      exit(0);
    }

    e1  = plan->e1;
    e2  = plan->e2;
    y   = plan->y;
    z   = plan->z;
    w   = plan->w;
    x   = plan->x;
    tmp = plan->for_subsys;
    n   = plan->n;

    s = (a/2+c)*(a/2+c) - b*b;
    t = a*a/2 - b*b - 2*c*c;
    sign_of_a = (a!=0)? a/fabs(a):0;
    p = (a-2*c)/4 + sign_of_a*sqrt(sign_of_a*(a-2*c)*sqrt(s)+t)/2 
                                           + sign_of_a*sqrt(s)/2;
    l1 = b/(p+c);
    l2 = c/p;

    for (i = 0 ; i < n ; i++)
    {
      tmp[i] = f[i]/p;
      e1[i] = 0.;
      e2[i] = 0.;
    }
    e1[0] = 1.;
    e2[1] = 1.;

    solve_subsystem(l1, l2, tmp, n, tmp, y); 
    solve_subsystem(l1, l2, e1,  n, tmp, z);
    solve_subsystem(l1, l2, e2,  n, tmp, w);

    d = 1. + l1*l2*(z[1]+w[0]) + l2*l2*(z[0]+w[1]) + 
       l1*l1*z[0] + pow(l2,4)*(z[0]*w[1]-z[1]*w[0]);
    d1 = pow(l2,4)*(y[0]*w[1]-y[1]*w[0]) + 
          (l1*l1+l2*l2)*y[0] + l1*l2*y[1];
    d2 = pow(l2,4)*(y[1]*z[0]-y[0]*z[1]) + 
                  l1*l2*y[0] + l2*l2*y[1];

    for (i = 0 ; i < n ; i++)
      x[i] = y[i] - (d1*z[i]+d2*w[i])/d;	     

    return x;

  }


  void dealloc_double(double *ptr)

  {

    free(ptr);
    ptr = NULL;

  }


  void delete_sym_penta_diagonal(sym_penta_diagonal_plan *plan)

  {
 
    // Plan components deallocation 
    dealloc_double(plan->e1);
    dealloc_double(plan->e2);
    dealloc_double(plan->y);
    dealloc_double(plan->z);
    dealloc_double(plan->w);

    // Plan deallocation
    free(plan);
    plan = NULL;
 
  }

