
// This is a draft symmetric penta-diagonal solver

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

typedef struct sym_penta_diagonal_plan
{
  unsigned n;
  double *e1;
  double *e2;
  double *y;
  double *z;
  double *w;
  double *x;
  double *for_subsys;
} plan_t;


double *alloc_double(double *ptr, unsigned n)

{

  ptr = malloc( n*sizeof(double) );

  if ( ptr == NULL )
    {
      fprintf(stderr,"Impossible to allocate pointer\n");
      exit(EXIT_FAILURE);
    }

}


plan_t *new_sym_penta_diagonal(unsigned n)

{
  // Declaration
  plan_t *plan;

  if (n<3)
     {
       printf("Matrix size must be at least 3x3\n");
       printf("Exiting...\n");
       exit(0);
     }

  // Plan allocation
  plan = malloc(sizeof(plan_t));

  if ( plan == NULL )
    {
      fprintf(stderr,"Impossible to allocate pointer\n");
      exit(EXIT_FAILURE);
    }

  // Plan components allocation
  plan->n  = n; 
  plan->e1 = alloc_double(plan->e1, n);
  plan->e2 = alloc_double(plan->e2, n);
  plan->y  = alloc_double(plan->y,  n);
  plan->z  = alloc_double(plan->z,  n);
  plan->w  = alloc_double(plan->w,  n);
  plan->x  = alloc_double(plan->x,  n);
  plan->for_subsys  = alloc_double(plan->for_subsys,  n);

  return plan;

}


void solve_subsystem(double l1, double l2, double *b,
                       int n, double *tmp, double* y)

{

  int i;
     
  tmp[0] = b[0];
  tmp[1] = b[1] - l1*tmp[0];
  for (i = 2; i < n; i++)
    {
      tmp[i] = b[i] - ( l2*tmp[i-2] + l1*tmp[i-1] );
    }

  y[n-1] = tmp[n-1];
  y[n-2] = tmp[n-2] - l1*y[n-1];
  for (i = n-3; i >= 0; i--)
    {
      y[i] = tmp[i] - ( l1*y[i+1] + l2*y[i+2] );
    }

}


double *solve_sym_penta_diagonal(double a, double b, double c, 
				     double *f, plan_t *plan)

{

  double s, t, p, l1, l2;
  double d, d1, d2;
  double *e1, *e2;
  double *y, *z, *w;
  double *x, *tmp;
  unsigned n, i;
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

  // Let define here sign(a)=1
  s = (a/2+c)*(a/2+c) - b*b;
  t = a*a/2 - b*b - 2*c*c;
  if ( a > 0 )
     {
       sign_of_a = 1;
     }
  else if ( a < 0 )
     {
        sign_of_a = -1;
     }
  else
     {
       sign_of_a = 0;
     }
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
    {
      x[i] = y[i] - (d1*z[i]+d2*w[i])/d;	     
    }

  return x;

}


void dealloc_double(double *ptr)

{

  free(ptr);
  ptr = NULL;

}


plan_t *delete_sym_penta_diagonal(plan_t *plan)

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

  return plan;
 
}



main(void)

{

  int n_max = 1000, nb_test = 100;
  int n, i_test, i;
  double a, b, c;
  double *f;
  double *x, *x_exact;
  double err, norm;
  plan_t *plan;

  for (i_test = 1; i_test <= nb_test; i_test++)
    {
      for (n = 3; n <= n_max; n++)
	{

          a = n*i_test; 
          b = a/4; 
          c = a/4 - 1;
	  f = alloc_double(f, n);
	  x = alloc_double(x, n);
	  x_exact = alloc_double(x_exact, n);

	  for (i = 0; i < n; i++)
	    {
	      x_exact[i] = i*n*i_test;//rand();
	    }

	  f[0] = a*x_exact[0] + b*x_exact[1] + c*x_exact[2]; 
	  f[1] = b*x_exact[0] + a*x_exact[1] + b*x_exact[2] + c*x_exact[3];

	  for (i = 2; i < n-2; i++)
	    {
	      f[i] = c*x_exact[i-2] + b*x_exact[i-1] + a*x_exact[i] + 
                                     b*x_exact[i+1] + c*x_exact[i+2];        
	    }

	  f[n-2] = c*x_exact[n-4] + b*x_exact[n-3] + a*x_exact[n-2] + b*x_exact[n-1];
	  f[n-1] = c*x_exact[n-3] + b*x_exact[n-2] + a*x_exact[n-1];

	  plan = new_sym_penta_diagonal(n);
	  x = solve_sym_penta_diagonal(a, b, c, f, plan);

	  err = 0.;
          norm = 0.;
	  for (i = 0; i < n; i++)
	    {
	      err = err + pow(x_exact[i]-x[i], 2);
              norm = norm + pow(x_exact[i], 2);
	    }
	  err = err/norm;

	  if (err > pow(10,-16))
	    {
	      printf("NOT PASS\n");
	      printf("Error = %f\n", err);
	      printf("i_test = %d\n", i_test);
	      printf("n = %d\n", n);
	      exit(0);
	    }

	  plan = delete_sym_penta_diagonal(plan);

	}
    }

  printf("PASS\n");

}

