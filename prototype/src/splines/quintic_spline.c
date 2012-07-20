
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "toep_penta_diagonal.h"
#include "quintic_spline.h"


  quintic_spline_hermite_plan *new_quintic_spline_hermite(const int n,
                const double xmin, const double xmax, const double *f)
  {

    quintic_spline_hermite_plan *plan;

    // Plan allocation
    plan = malloc(sizeof(quintic_spline_hermite_plan));
    if ( plan == NULL )
    {
      fprintf(stderr,"Impossible to allocate pointer\n");
      exit(EXIT_FAILURE);
    }

    // plan component allocation
    plan->coeffs = malloc( (n+6) * sizeof(double) );
    if ( plan->coeffs == NULL )    
    {
      fprintf(stderr,"Impossible to allocate pointer\n");
      exit(EXIT_FAILURE);
    }

    plan->n = n;
    plan->xmin = xmin;
    plan->xmax = xmax;
    compute_coeffs_hermite(f, plan);

    return plan;

  }


  void compute_coeffs_hermite(const double *f, 
        quintic_spline_hermite_plan *plan_spline)
  /* f is the vector of the values of the function 
     in the nodes of the mesh*/
 
  {

    double a, b, c, xmin, xmax, h;
    int n, i;
    toep_penta_diagonal_plan *plan_pent;

    n = plan_spline->n;
    xmin = plan_spline->xmin;
    xmax = plan_spline->xmax;
    h = (xmax-xmin)/n;

    a = B(5, n-3, xmax, plan_spline);
    b = B(5, n-4, xmax, plan_spline); 
    c = B(5, n-5, xmax, plan_spline); 

    plan_pent = new_toep_penta_diagonal(n+6);
    plan_spline->coeffs = solve_toep_penta_diagonal(a, b, c, f, plan_pent);
    delete_toep_penta_diagonal(plan_pent);    

  }

  double B(int j, int i, double x, quintic_spline_hermite_plan *plan_spline)
  {

    double xmin, xmax, h;
    int n;

    xmin = plan_spline->xmin;
    xmax = plan_spline->xmax;
    n = plan_spline->n;
    h = (xmax-xmin)/n;

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
      return ((x-xmin-(double)i*h)/((double)j*h)) * B(j-1, i, x, plan_spline) + 
        ((xmin+(double)(i+j+1)*h-x)/((double)j*h)) * B(j-1, i+1, x, plan_spline);
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


  double spline_hermite(const double x, quintic_spline_hermite_plan *plan_spline) 
  // The interpolator spline function

  {

    int n, i_x, i, j;
    double xmin, xmax, y, h, s = 0.;

    xmin = plan_spline->xmin;
    xmax = plan_spline->xmax;
    n = plan_spline->n;
    h = (xmax-xmin)/n;
    i = (int)((x-xmin)/h); // Determine the leftmost support index 'i' of x

    for (j=i-5; j<=i; j++)
    {
      if( (j>=-5) && (j<=n) )
      {
        s += plan_spline->coeffs[j+5] * B(5, j, x, plan_spline);
      }
    }
    return s;

  }


  void delete_quintic_spline_hermite(quintic_spline_hermite_plan *plan)

  {
 
    // Plan components deallocation 
    free(plan->coeffs);
    plan->coeffs = NULL;

    // Plan deallocation
    free(plan);
    plan = NULL;
 
  }


