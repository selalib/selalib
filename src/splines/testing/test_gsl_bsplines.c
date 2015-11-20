#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_bspline.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_errno.h>

/* number of data points to fit */
#define N        5

/* number of fit coefficients */
#define NCOEFFS  5


/* Spline order */
#define K  3

/* nbreak = ncoeffs + 2 - k */
#define NBREAK   (NCOEFFS + 2 - K)

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
  const size_t k = K;
  double xi;
  gsl_vector *x;
  gsl_matrix *A;

  double *a;
  double *b;
  double *c;
  double *r;
  double *u;

  a = (double *) malloc(n*sizeof(double));
  b = (double *) malloc(n*sizeof(double));
  c = (double *) malloc(n*sizeof(double));
  r = (double *) malloc(n*sizeof(double));
  u = (double *) malloc(n*sizeof(double));

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

/*
void tridag(double a[], 
            double b[], 
            double c[], 
            double r[], 
            double u[], 
            unsigned long n)
{
  unsigned long j; 
  float bet;
  double *gam;

  gam = gsl_vector_alloc(n);

  if (b[1] == 0.0) 
    printf("Error 1 in tridag");

  u[1]=r[1]/(bet=b[1]);

  for (j=2;j<=n;j++) { 
   gsl_vector_set(gam, j, c[j-1]/bet);
   bet=b[j]-a[j]*gam[j];
   if (bet == 0.0) 
     printf("Error 2 in tridag");
   u[j]=(r[j]-a[j]*u[j-1])/bet;

  }

  for (j=(n-1);j>=1;j--)
    u[j] -= gsl_vector_get (gam,j+1)*u[j+1]; 

  gsl_vector_free(gam);

}

void cyclic(double a[], 
            double b[], 
            double c[], 
            double alpha, 
            double beta, 
            double r[], 
            double x[], unsigned long n)

{
  unsigned long i;
  double fact,gamma;
  double *bb;
  double *u;
  double *z;
  
  if (n <= 2) printf("n too small in cyclic"); 
  bb = (double *) malloc(n*sizeof(double));
  u  = (double *) malloc(n*sizeof(double));
  z  = (double *) malloc(n*sizeof(double));
  
  gamma = -b[1]; 
  bb[1]=b[1]-gamma; 
  bb[n]=b[n]-alpha*beta/gamma; 
  for (i=2;i<n;i++) 
    bb[i]=b[i]; 
  
  tridag(a,bb,c,r,x,n); 
  u[1]=gamma;
  u[n]=alpha;
  for (i=2;i<n;i++) 
    u[i]=0.0; 

  tridag(a,bb,c,u,z,n); 

  fact=(x[1]+beta*x[n]/gamma)/(1.0+z[1]+beta*z[n]/gamma);

  for (i=1;i<=n;i++) 
    x[i] -= fact*z[i]; 
  
  gsl_vector_free(z); 
  gsl_vector_free(u); 
  gsl_vector_free(bb);
}
*/
