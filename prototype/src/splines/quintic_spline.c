
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "sym_penta_diagonal.h"
#include "quintic_spline.h"


  quintic_spline_plan *new_quintic_spline(const int n,
                 const double xmin, const double xmax)
  {

    quintic_spline_plan *plan;

    // Plan allocation
    plan = malloc(sizeof(quintic_spline_plan));
    if ( plan == NULL )
    {
      fprintf(stderr,"Impossible to allocate pointer\n");
      exit(EXIT_FAILURE);
    }

    // plan component allocation
    plan.p= malloc( n*sizeof(double) );
    if ( plan.p == NULL )
    {
      fprintf(stderr,"Impossible to allocate pointer\n");
      exit(EXIT_FAILURE);
    }

    if (size<3)
    {
      pritnf("The number of points must be greater than"); 
                                printf("or equal to 3\n");
      printf"Exiting...\n");
      exit(0);
    }

    plan->n = n;
    plan->xmin = xmin;
    plan->xmax = xmax;

    return plan;

  }
  

  double B_spline(const double x, quintic_spline_plan *plan)
  // Quintic B-spline computed as scalar

  {

    double b, xmin, xmax, dx;
    int n, cas;

    n = plan->n;
    xmin = plan->xmin;
    xmax = plan->xmax;
    dx = (xmax-xmin)/(n-1);

    cas = ( (-3*dx<x) && (x<=-2*dx) ) + 2*( (-3*dx<x) && (x<=-dx) ) + 
                    3*( (-dx<x) && (x<=0) ) + 4*( (0<x) && (x<dx) ) + 
             5*( (dx<=x) && (x<2*dx) ) + 6*( (2*dx<=x) && (x<3*dx) );
    switch(cas)
    {
      case 1:
        b = ( pow(x,5) + 15*pow(x,4)*dx + 90*pow(x,3)*dx*dx +
            270*x*x*pow(dx,3) + 405*x*pow(dx,4) + 243*pow(dx,5) )/2;
        break;
      case 2:
        b = ( -5*pow(x,5) - 45*pow(x,4)*dx - 150*pow(x,3)*dx*dx -
           210*x*x*pow(dx,3) - 75*x*pow(dx,4) + 51*pow(dx,5) )/2;
        break;
      case 3:
        b = 5*pow(x,5) + 15*pow(x,4)*dx - 30*x*x*pow(dx,3) + 33*pow(dx,5);
        break;
      case 4:
        b = -5*pow(x,5) + 15*pow(x,4)*dx - 30*x*x*pow(dx,3) + 33*pow(dx,5);
        break;
      case 5:
        b = ( 5*pow(x,5) - 45*pow(x,4)*dx + 150*pow(x,3)*dx*dx -
           210*x*x*pow(dx,3) + 75*x*pow(dx,4) + 51*pow(dx,5) )/2;
        break;
      case 6:
        b = ( -pow(x,5) + 15*pow(x,4)*dx - 90*pow(x,3)*dx*dx +
            270*x*x*pow(dx,3) - 405*x*pow(dx,4) + 243*pow(dx,5) )/2;
        break;
      default:
        b = 0;
        break;
    }

    return b/(60*pow(dx,6));

  }


  double* quintic_spline::coeffs(const double *f, quintic_spline_plan *plan_spline;)
  /* f is the vector of the values of the function 
     in the nodes of the mesh*/
 
  {

    double a, b, c, xmin, xmax, dx;
    int n
    plan_t *plan_pent;

    n = plan_spline->n;
    xmin = plan_spline->xmin;
    xmax = plan_spline->xmax;
    dx = (xmax-xmin)/(n-1);

    a = B_spline(0, plan_spline);
    b = B_spline(dx, plan_spline);
    c = B_spline(2*dx, plan_spline);

    plan_pent = new_sym_penta_diagonal(n+2);
    return solve_sym_penta_diagonal(a, b, c, f, plan_pent);
    delete_sym_penta_diagonal(plan_pent);

  }


  double spline(const double x, const double *f, quintic_spline_plan *plan) 
  // The interpolator spline function

  {

    int n, i;
    double xmin, xmax, dx, som = 0.;

    n = plan_spline->n;
    xmin = plan_spline->xmin;
    xmax = plan_spline->xmax;
    dx = (xmax-xmin)/(n-1);
    p = coeffs(f);

    for (i=-1; i<=n; i++)
      som += p[i]*B_spline(x-i*dx);

    return som;

  }


  void delete_quintic_spline(sym_penta_diagonal_plan *plan)

  {
 
    // Plan components deallocation 
    free(plan.p);
    plan.p = NULL;

    // Plan deallocation
    free(plan);
    plan = NULL;
 
  }


