// This is a symmetric penta-diagonal solver module

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "sym_penta_diagonal.c"

main(void)

{

  int n_max = 6, nb_test = 1;
  int n, i_test, i;
  double a, b, c;
  double *f;
  double *x, *x_exact;
  double err, norm;
  sym_penta_diagonal_plan *plan;




  for (i_test = 1; i_test <= nb_test; i_test++)
  {
    for (n = n_max; n <= n_max; n++)
    {

      /*a = n*i_test; 
      b = a/4; 
      c = a/4 - 1;*/
      a = -14;
      b = 4;
      c = -2;

      f = alloc_double(n);
      x = alloc_double(n);
      x_exact = alloc_double(n);

      for (i = 0; i < n; i++)
        x_exact[i] = i*n*i_test;//rand();

      /*f[0] = 1000.;//a*x_exact[0] + b*x_exact[1] + c*x_exact[2]; 
      f[1] = 1000.;//b*x_exact[0] + a*x_exact[1] + b*x_exact[2] + c*x_exact[3];

      for (i = 2; i < n-2; i++)
      {
        f[i] = 1000.;//c*x_exact[i-2] + b*x_exact[i-1] + a*x_exact[i] + 
                               b*x_exact[i+1] + c*x_exact[i+2];        
      }

      f[n-2] = 1000.;//c*x_exact[n-4] + b*x_exact[n-3] + a*x_exact[n-2] + b*x_exact[n-1];
      f[n-1] = 1000.;//c*x_exact[n-3] + b*x_exact[n-2] + a*x_exact[n-1];*/

      f[0] = -12; f[1] = -8; f[2] = -10; f[3] = -10; f[4] = -8; f[5] = -12;

      plan = new_sym_penta_diagonal(n);
      x = solve_sym_penta_diagonal(a, b, c, f, plan);

      err = 0.;
      norm = 0.;
      for (i = 0; i < n; i++)
      {
        err = err + pow(x_exact[i]-x[i], 2);
        norm = norm + pow(x_exact[i], 2);
        printf("%f\n", x[i]);
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

     // delete_sym_penta_diagonal(plan);

    }
  }

  printf(" PASS\n");

}

