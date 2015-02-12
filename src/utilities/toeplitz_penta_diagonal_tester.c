// This is a Toeplitz penta-diagonal solver tester

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "toep_penta_diagonal.c"
#include <time.h>

main(void)

{

  int n_max = 1000, nb_test = 100;
  int n, i_test, i;
  double a, b, c;
  double *f;
  double *x, *x_exact;
  double err, norm;
  toep_penta_diagonal_plan *plan;

  printf("Testing Toeplitz penta-diagonal solver...\n");
  /* initialize random seed: */
  srand ( time(NULL) );

  for (i_test = 1; i_test <= nb_test; i_test++)
  {
    for (n = 3; n <= n_max; n++)
    {

      a = n*i_test; 
      b = a/4; 
      c = a/4 - 1;

      f = alloc_double(n);
      x = alloc_double(n);
      x_exact = alloc_double(n);

      for (i = 0; i < n; i++)
      {
        x_exact[i] = log( fabs( rand() ) ); 
      }

      f[0] = a*x_exact[0] + b*x_exact[1] + c*x_exact[2]; 
      f[n-1] = c*x_exact[n-3] + b*x_exact[n-2] + a*x_exact[n-1];

      if ( n > 3 )
      {
         f[1] = b*x_exact[0] + a*x_exact[1] + b*x_exact[2] + c*x_exact[3];
         f[n-2] = c*x_exact[n-4] + b*x_exact[n-3] + a*x_exact[n-2] + b*x_exact[n-1];
      }
      else
      {
         f[1] = b*x_exact[0] + a*x_exact[1] + b*x_exact[2]
      }

      for (i = 2; i < n-2; i++)
      {
        f[i] = c*x_exact[i-2] + b*x_exact[i-1] + a*x_exact[i] + 
                               b*x_exact[i+1] + c*x_exact[i+2];        
      }

      plan = new_toep_penta_diagonal(n);
      x = solve_toep_penta_diagonal(a, b, c, f, plan);

      err = 0.;
      norm = 0.;
      for (i = 0; i < n; i++)
      {
        err = err + pow(x_exact[i]-x[i], 2);
        norm = norm + pow(x_exact[i], 2);
        //printf("%f\n", x[i]);
      }
      err = err/norm;

      if (err > pow(10,-15))
      {
        printf("Program stopped by Toeplitz penta-diagonal solver failure \n");
        printf("Error = %f\n", err);
        printf("i_test = %d\n", i_test);
        printf("n = %d\n", n);
        exit(0);
      }

      delete_toep_penta_diagonal(plan);

    }
  }

  printf("Toeplitz penta-diagonal solver unit test: PASS\n");

}

