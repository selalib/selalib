

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "quintic_spline.h"
#include <time.h>

int main(void)
{
  
  #define pi 3.14159265
  quintic_spline_hermite_plan *s;
  double *f;
  double h, xmin, xmax, x, x_rand, y, err1, err2, max, bound;
  int n, nmax, i;

  nmax = 1000;
  printf("Testing quintic_spline_hermite...\n");

  for(n=1; n<=nmax; n++)
  {

    /* initialize random seed: */
    srand ( time(NULL) );
    xmin = log( fabs(rand()) + 1 );
    xmax = xmin + fabs( log( fabs(rand()) + 1 ) );
    h = (xmax-xmin)/n;

    f = malloc( (n+6) * sizeof(double));

    for (i=-2; i<=n+3; i++)
    {
      x = xmin + (double)i*h;
      f[i+2] = log( (double)n ) * ( 2*pi*x/(n*h) - sin(2*pi*(x-xmin)/(n*h)) ); 
    }

    s = new_quintic_spline_hermite(n, xmin, xmax, f);

    err1 = 0.; 
    err2 = 0.;

    for (i=0; i<=n; i++)
    {

      x = xmin + (double)i*h;
      err1 += fabs( f[i+2] - spline_hermite(x, s) );

      x_rand = xmin + n*h*( (rand()%11) / 10. );

      y = log( (double)n ) * ( 2*pi*x_rand/(n*h) - sin( 2*pi*(x_rand-xmin)/(n*h) ) ); 
      max = fabs( y - spline_hermite(x_rand, s) );
      if (err2 < max )
      {
        err2 = max; 
      }

    }

    bound = log( (double)n ) * 8*pow(pi,4) / ( 2*pow(n,3)*pow(h,3/2) );
    printf("Average absolute errors: %f, %f, %f \n", err1/(n+1), err2, bound);

    if ( (err1/(n+1) > pow(10,-9)) || (err2 > bound) )
    {
      printf("Program stopped by iteration number %d of quintic ");
      printf("splines unit test:\n", n);
      printf("Exciting...\n");
      exit(0);
    }

    free(f);
    f = NULL;
    delete_quintic_spline_hermite(s);

  }

  printf("Quitinc splines unit test: PASS\n");

}
