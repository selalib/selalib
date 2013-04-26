
// This is a Toeplitz penta-diagonal solver module

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "toep_penta_diagonal.h"


  double *alloc_double(int n)
  {
    double *ptr = malloc( n*sizeof(double) );
    if ( ptr == NULL )    {
      fprintf(stderr,"Impossible to allocate pointer\n");
      exit(EXIT_FAILURE);
    }
    return ptr;
  }


  toep_penta_diagonal_plan *new_toep_penta_diagonal(int n)

  {
    // Declaration
    toep_penta_diagonal_plan *plan;

    if (n<3)
    {
      printf("Matrix size must be at least 3x3\n");
      printf("Exiting...\n");
      exit(0);
    }

    // Plan allocation
    plan = malloc(sizeof(toep_penta_diagonal_plan));

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
    plan->x  = alloc_double(n);

    return plan;

  }


  double *solve_subsystem(double l1, double l2, double *b, int n)

  {

    int i;
    double *x, *y;

    x = alloc_double(n);
    y = alloc_double(n);
     
    y[0] = b[0];
    y[1] = b[1] - l1*y[0];
    for (i = 2; i < n; i++)
      y[i] = b[i] - ( l2*y[i-2] + l1*y[i-1] );

    x[n-1] = y[n-1];
    x[n-2] = y[n-2] - l1*y[n-1];
    for (i = n-3; i >= 0; i--)
      x[i] = y[i] - ( l1*x[i+1] + l2*x[i+2] );

    return x;

    dealloc_double(y);

  }


  double *solve_toep_penta_diagonal(double a, double b, double c, 
  		 const double *f, toep_penta_diagonal_plan *plan)

  {

    double s, t, p, l1, l2;
    double d, d1, d2;
    double *e1, *e2;
    double *x, *y, *z, *w;
    int n, i;
    int sign_of_a;

    if ( fabs(a) <= 2*(fabs(b)+fabs(c)) )
    {
      printf("a, b, and c must be such that: |a| > 2(|b|+|c|)\n");
      printf("a=%f, b=%f, c=%f\n", a, b, c);
      printf("Exiting...\n");
      exit(0);
    }

    e1  = plan->e1;
    e2  = plan->e2;
    y   = plan->y;
    x   = plan->x;
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
      y[i] = f[i]/p;
      e1[i] = 0.;
      e2[i] = 0.;
    }
    e1[0] = 1.;
    e2[1] = 1.;

    y = solve_subsystem(l1, l2, y, n); 
    z = solve_subsystem(l1, l2, e1,  n);
    w = solve_subsystem(l1, l2, e2,  n);

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


  void delete_toep_penta_diagonal(toep_penta_diagonal_plan *plan)

  {
 
    // Plan components deallocation 
    dealloc_double(plan->e1);
    dealloc_double(plan->e2);
    dealloc_double(plan->y);

    // Plan deallocation
    free(plan);
    plan = NULL;
 
  }

