

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

  for(n=1; n<=nmax; n++)
  {

    /* initialize random seed: */
    srand ( time(NULL) );
    xmin = (rand()%11) - 5;//Generating random number between -5 and 5
    xmax = xmin + (rand()%100) + 1;//Generating random number between 1 and 100
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

    for (i=-2; i<=n+3; i++)
    {

      x = xmin + (double)i*h;
      err1 += fabs( f[i+2] - spline_hermite(x, s) );

      x_rand = xmin-2*h + (n+5)*h*( (rand()%11) / 10. );

      y = log( (double)n ) * ( 2*pi*x_rand/(n*h) - sin( 2*pi*(x_rand-xmin)/(n*h) ) ); 
      max = fabs( y - spline_hermite(x_rand, s) );
      if (err2 < max )
      {
        err2 = max; 
      }

    }

    bound = log( (double)n ) * 2*pi/n;
    printf("Average absolute errors: %f, %f, %f \n", err1/(n+5), err2, bound);

    if ( (err1/(n+5) > pow(10,-9)) || (err2 > bound) )
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
