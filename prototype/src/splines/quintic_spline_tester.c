

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "quintic_spline.h"
#include <time.h>


  double B_test(int j, int i, double x, double xmin, double h)
  {

    int n;

    /*
                 x-t(i)                       t(i+j+1)-x
    B[j,i](x) = ----------- * B[j-1,i](x) + ----------------* B[j-1,i+1](x)
               t(i+j)-t(i)                  t(i+j+1)-t(i+1)
    And
  
    B[0,i] = 1 if t(i) <= x < t(i+1), and 0 otherwise.
    
    t(i) = xmin + i*h

    */

    if (j!=0)
    {                                         
      return ((x-xmin-(double)i*h)/((double)j*h)) * B_test(j-1, i, x, xmin, h) + 
        ((xmin+(double)(i+j+1)*h-x)/((double)j*h)) * B_test(j-1, i+1, x, xmin, h);
    }
    else
    {
      if ( ( xmin + (double)i*h <= x ) && ( x < xmin + (double)(i+1)*h ) )
      {
        return 1.;
      }
      else
      {
        return 0.;
      }

    }

  }



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

  for(n=1; n<=nmax; n++)
  {

    /* initialize random seed: */
    srand ( time(NULL) );

    xmin=0;
    xmax=100;
    h = (xmax-xmin)/n;
    mu = ( (xmin-5*h) + (xmin+(n+6)*h) ) / 2;

    f1 = malloc( (n+6) * sizeof(double));
    f2 = malloc( (n+6) * sizeof(double));

    for (i=-2; i<=n+3; i++)
    {

      x = xmin + (double)i*h;
      f1[i+2] = exp( - .5*pow( x - mu, 2) ); 

      f2[i+2] = 0.;
      for(j=i-5; j<=i; j++)
      {
         if ( (j>=-5) && (j<=n) )
         {
            f2[i+2] = f2[i+2] + B_test(5, j, x, xmin, h);
         }
      }

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

      y = 0.;
      for(j=left-5; j<=left; j++)
      {
         if ( (j>=-5) && (j<=n) )
         {
            y += B_test(5, j, x, xmin, h);
         }
      }
    
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
