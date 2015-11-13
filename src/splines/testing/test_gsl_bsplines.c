#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_bspline.h>
#include <gsl/gsl_spline.h>

/* number of data points to fit */
#define N        5

/* number of fit coefficients */
#define NCOEFFS  5

/* nbreak = ncoeffs + 2 - k = ncoeffs - 2 since k = 4 */
#define NBREAK   (NCOEFFS - 2)

int
main (void)
{
  const size_t n = N;
  const size_t ncoeffs = NCOEFFS;
  const size_t nbreak = NBREAK;
  size_t i, j;
  gsl_bspline_workspace *bw;
  gsl_bspline_deriv_workspace *dw;
  gsl_vector *B;
  double xmin=0.0, xmax=2.0*M_PI;
  const size_t nderiv=1;
  gsl_matrix *dB;
  const size_t k = 4;
  double xi;
  gsl_vector *x;
  gsl_matrix *X;

  x = gsl_vector_alloc(n);
  

  /* allocate a cubic bspline workspace */
  bw = gsl_bspline_alloc(k, nbreak);
  dw = gsl_bspline_deriv_alloc(k);
  B  = gsl_vector_alloc(ncoeffs);
  dB = gsl_matrix_alloc(nbreak+k-2, nderiv+1);
  X = gsl_matrix_alloc(n, ncoeffs);

  /* use uniform breakpoints on [xmin, xmax] */
  gsl_bspline_knots_uniform(xmin, xmax, bw);

  /* construct the fit matrix X */
  for (i = 0; i < n; ++i)
    {
      double xi = ((xmax-xmin) / (N - 1)) * i;
      double ti = gsl_vector_get(bw->knots, i);

      gsl_vector_set(x, i, xi);

      /* compute B_j(xi) for all j */
      gsl_bspline_eval(xi, B, bw);
      int out = gsl_bspline_deriv_eval (xi, nderiv, dB, bw, dw);

      printf("%f %f\n", xi, ti);

      /* fill in row i of X */
      for (j = 0; j < ncoeffs; ++j)
        {
          double Bj = gsl_vector_get(B, j);
          gsl_matrix_set(X, i, j, Bj);
        }

    }

  gsl_bspline_free(bw);
  gsl_bspline_deriv_free(dw);
  gsl_vector_free(B);
  return 0;

} 


