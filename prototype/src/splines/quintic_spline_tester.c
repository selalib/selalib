

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "quintic_spline.h"
#include <time.h>

int main(void)
{
  
  #define pi 3.14159265
  quintic_spline_plan *s;
  double *f;
  double h, xmin, xmax, mu;
  double x, y, norm1;
  double err1, err2;
  int n, nmax, i;

  nmax = 100;
  printf("Testing quintic_spline...\n");

  for(n=1; n<=nmax; n++)
  {

    /* initialize random seed: */
    srand ( time(NULL) );
    xmin=0;xmax=100;
    h = (xmax-xmin)/n;
    mu = ( (xmin-5*h) + (xmin+(n+6)*h) ) / 2;

    f = malloc( (n+6) * sizeof(double));

    for (i=-2; i<=n+3; i++)
    {
      x = xmin + (double)i*h;
      f[i+2] = exp( - .5*pow( x - mu, 2) ); 
    }

    s = new_quintic_spline(n, xmin, xmax, f);

    err1 = 0.; 
    norm1 = 0.;

    for (i=0; i<=n; i++)
    {

      x = xmin + (double)i*h;
      err1 += ( f[i+2] - spline(x, s) ) * ( f[i+2] - spline(x, s) );
      norm1 += f[i+2] * f[i+2];

      x = xmin + n*h*( (rand()%11) / 10. );

      y = exp( - .5*pow( x - mu, 2) ); 

    }

    err1 = sqrt(err1/norm1);
    printf("Average absolute errors: %f, %f \n", err1, err2);

    if ( (err1 > pow(10,-9)) )
    {
      printf("Program stopped by iteration number %d of quintic ", n);
      printf("splines unit test:\n");
    printf("%f, %f\n", y, spline(x, s) );
      printf("Exiting...\n");
      exit(0);
    }

    free(f);
    f = NULL;
    delete_quintic_spline(s);

  }

  printf("Quitinc splines unit test: PASS\n");

}
