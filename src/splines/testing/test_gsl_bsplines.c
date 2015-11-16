#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_bspline.h>
#include <gsl/gsl_spline.h>

/* number of data points to fit */
#define N        5

/* number of fit coefficients */
#define NCOEFFS  5

/* nbreak = ncoeffs + 2 - k */
#define NBREAK   (NCOEFFS - 2)

int
main (void)
{
  const size_t n = N;
  const size_t ncoeffs = NCOEFFS;
  const size_t nbreak = NBREAK;
  size_t i, j;
  gsl_bspline_workspace *bw;
  gsl_vector *B;
  double xmin=0.0, xmax=2.0;
  const size_t k = 4;
  double xi;
  gsl_vector *x;
  gsl_matrix *A;

  x = gsl_vector_alloc(n);
  for (i = 0; i < n; ++i)
    {
       
    }
  

  /* allocate a cubic bspline workspace */
  bw = gsl_bspline_alloc(k, nbreak);
  B  = gsl_vector_alloc(ncoeffs);
  A  = gsl_matrix_alloc(n, ncoeffs);

  /* gsl_bspline_knots(breakpts, bw) 
  Compute the knots from the given breakpoints:
   knots(1:k) = breakpts(1)
   knots(k+1:k+l-1) = breakpts(i), i = 2 .. l
   knots(n+1:n+k) = breakpts(l + 1) */

  /* use uniform breakpoints on [xmin, xmax] */
  gsl_bspline_knots_uniform(xmin, xmax, bw);
  //knots(1:k) = xmin
  //knots(k+1:k+l-1) = xmin + i*delta, i = 1 .. l - 1
  //knots(n+1:n+k) = xmax

  /* construct the fit matrix X */
  for (i = 0; i < n; ++i)
    {
      double xi = ((xmax-xmin) / (N - 1)) * i;
      double ti = gsl_vector_get(bw->knots, i);

      gsl_vector_set(x, i, xi);

      /* compute B_j(xi) for all j */
      gsl_bspline_eval(xi, B, bw);

      printf("%f %f\n", xi, ti);

      /* fill in row i of X */
      for (j = 0; j < ncoeffs; ++j)
        {
          double Bj = gsl_vector_get(B, j);
          gsl_matrix_set(A, i, j, Bj);
        }

    }
  for (i = 0; i < n; i++)  
    { 
      for (j = 0; j < ncoeffs; j++)
        printf ("%6.3f", gsl_matrix_get (A, i, j));
      printf("\n");
    }

/*
  gsl_permutation * p = gsl_permutation_alloc (ncoeffs);

  gsl_linalg_LU_decomp (A, p, &s);

  gsl_linalg_LU_solve (A, p, b, x);

  printf ("x = \n");
  gsl_vector_fprintf (stdout, x, "%g");

  gsl_permutation_free (p);
  gsl_vector_free (x);
*/

  //gsl_vector *d; d = gsl_vector_alloc(ncoeffs);
  //gsl_vector *e; e = gsl_vector_alloc(ncoeffs);
  //gsl_vector *f; f = gsl_vector_alloc(ncoeffs);
  //gsl_vector *b; b = gsl_vector_alloc(ncoeffs);
  //gsl_vector *x; x = gsl_vector_alloc(ncoeffs);

  // A = ( d_0 e_0  0  f_3 )
  //     ( f_0 d_1 e_1  0  )
  //     (  0  f_1 d_2 e_2 )
  //     ( e_3  0  f_2 d_3 )
 
  // gsl_linalg_solve_cyc_tridiag (d, e, f, b, x)

  gsl_bspline_free(bw);
  gsl_vector_free(B);
  return 0;

} 

