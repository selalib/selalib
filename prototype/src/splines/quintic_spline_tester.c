

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "quintic_spline.h"
#include <time.h>

int main(void)
{
  
  #define pi 3.14159265
  quintic_spline_plan *plan1, *plan2;
  double *f1, *f2;
  double h, xmin, xmax, mu;
  double x, y, norm1, norm2;
  double err1, err2;
  int n, nmax, i, j, left;

  nmax = 1000;
  printf("Testing quintic_spline...\n");

  for(n=5; n<=nmax; n++)
  {

    /* initialize random seed: */
    srand ( time(NULL) );

    xmin=-10;
    xmax=10;
    h = (xmax-xmin)/n;
    mu = (xmin+xmax) / 2;

    f1 = malloc( (n+6) * sizeof(double));
    f2 = malloc( (n+6) * sizeof(double));
    f1[0] = 0.;f1[1] = 0.;f1[n+3] = 0.;f1[n+4] = 0.;f1[n+5] = 0.;
    f2[0] = 0.;f2[1] = 0.;f2[n+3] = 0.;f2[n+4] = 0.;f2[n+5] = 0.;

    for (i=0; i<=n; i++)
    {
      x = xmin + (double)i*h;
      f1[i+2] = exp( - .5*pow( x - mu, 2) ); 
      f2[i+2] = f2[i] = x*(x-xmin)*(xmax-x);
    }

    plan1 = new_quintic_spline(n, xmin, xmax, f1);
    plan2 = new_quintic_spline(n, xmin, xmax, f2);

    err1 = 0.; 
    norm1 = 0.;
    err2 = 0.; 
    norm2 = 0.;

    for (i=0; i<=n; i++)
    {

      x = xmin + (double)i*h;
      err1 += ( f1[i+2] - spline(x, plan1) ) * ( f1[i+2] - spline(x, plan1) );
      norm1 += f1[i+2] * f1[i+2];

      x = xmin + n*h*( (rand()%11) / 10. );
      left = (int)((x-xmin)/h);//Determine the leftmost support index 'i' of x

      y = x*(x-xmin)*(xmax-x);
    
      err2 += ( y - spline(x, plan2) ) * ( y - spline(x, plan2) );
      norm2 += y*y;

    }

    err1 = sqrt(err1/norm1);
    err2 = sqrt(err2/norm2);
    printf("Relative errors: %e, %e \n", err1, err2);

    if ( (err1 >= pow(10,-12)) || (err2 >= pow(10,-12)) )
    {
      printf("Program stopped by iteration number %d of quintic ", n);
      printf("splines unit test:\n");
      printf("Exiting...\n");
      exit(0);
    }

    free(f1);
    f1 = NULL;
    free(f2);
    f2 = NULL;
    delete_quintic_spline(plan1);
    delete_quintic_spline(plan2);

  }

  printf("Quitinc splines unit test: PASS\n");

}
