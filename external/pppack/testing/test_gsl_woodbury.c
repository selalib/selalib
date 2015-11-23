#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_bspline.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_linalg.h>
#include <assert.h>

/* number of data points to fit */
#define N       8 

void tensor_product( gsl_vector * u, 
                     gsl_vector * v,
                     gsl_matrix * w) 
{

 assert( u->size == v->size); 
 size_t i, j, n;
 n = u->size;
 w  = gsl_matrix_alloc(n,n);
 for (i = 0; i < u->size; i++)
   for (j = 0; j < v->size; j++)
      gsl_matrix_set(w,i,j,gsl_vector_get(u,i)*gsl_vector_get(v,j));
  
}


int
main (void)
{
  const size_t n = N;
  size_t i, j;

  gsl_matrix *A;
  gsl_vector *b;

  printf("n       = %d \n", n);

  A  = gsl_matrix_alloc(n, n);

  for (i = 0; i < n; ++i)
    for (j = 0; j < n; ++j)
    {
      gsl_matrix_set(A, i, (i-2)%n, +1);
      gsl_matrix_set(A, i, (i-1)%n, -1);
      gsl_matrix_set(A, i, i, 4);
      gsl_matrix_set(A, i, (i+1)%n, -1);
      gsl_matrix_set(A, i, (i+2)%n, +1);
    }

  for ( i = 0; i < n; i++)  
    { for ( j = 0; j < n; j++)
        {
          printf ("%6.3f", gsl_matrix_get (A, i, j));
        }
      printf("\n");
    }

  b = gsl_vector_alloc(n);
  for ( i=0; i<n;i++)
      gsl_vector_set(b,i,1.);

  int s;

  gsl_permutation * p = gsl_permutation_alloc (n);

  gsl_linalg_LU_decomp (A, p, &s);

  gsl_vector *x; x = gsl_vector_alloc(n);

  gsl_linalg_LU_solve (A, p, b, x);

  printf ("x = \n");
  for ( i=0; i<n; i++)
    printf ("x[%d] = %f \n", i, gsl_vector_get(x,i));

  gsl_vector *d; d = gsl_vector_alloc(n);
  gsl_vector *e; e = gsl_vector_alloc(n-1);
  gsl_vector *f; f = gsl_vector_alloc(n-1);

  for (i = 0; i < n; i++)  
    { 
      gsl_vector_set(d, i,  4);
      if (i>0) gsl_vector_set(f, i-1, -1.);
      if (i<n-1) gsl_vector_set(e, i, -1.);
    }

  gsl_vector * bb; bb = gsl_vector_alloc(n);
  gsl_vector * u;  u  = gsl_vector_alloc(n);
  gsl_vector * z;  z  = gsl_vector_alloc(n);
  
  double alpha = -1.;
  double beta  = -1.;
  double gamma   = - gsl_vector_get(d,0); 
  gsl_vector_set(bb,0,   gsl_vector_get(d,0)-gamma); 
  gsl_vector_set(bb,n-1, gsl_vector_get(d,n-1)-alpha*beta/gamma); 

  for (int i=1;i<n-1;i++) 
    gsl_vector_set(bb,i,gsl_vector_get(b,i)); 
  
  gsl_linalg_solve_tridiag ( d, e, f, b, x);

  gsl_vector_set(u,0,gamma);
  gsl_vector_set(u,n-1,alpha);

  for (int i=1;i<n-1;i++) 
    gsl_vector_set(u,i,0.0); 

  gsl_linalg_solve_tridiag ( d, e, f, u, z);

  double fact;
  fact=( x->data[0]+beta*x->data[n-1]/gamma)
       /(1.0+z->data[0]+beta*z->data[n-1]/gamma);

  for (int i=0;i<n;i++) 
    x->data[i] -= fact*z->data[i]; 

  printf ("x = \n");
  for ( i=0; i<n; i++)
    printf ("x[%d] = %f \n", i, gsl_vector_get(x,i));

  /*
  gsl_vector *u;
  gsl_vector *v;
  gsl_matrix *w;
  u  = gsl_vector_alloc(n);
  v  = gsl_vector_alloc(n);
  w  = gsl_matrix_alloc(n,n);

  for ( i=0; i<n; i++) gsl_vector_set(u, i, 0.0);

  gsl_vector_set(u, 0, 1.0);
  gsl_vector_set(u, n-1, 1.0);
  for ( j=0; j<n; j++) gsl_vector_set(v, j, 0.0);
  gsl_vector_set(v, 0, 1.0);
  gsl_vector_set(v, n-1, 1.0);

  tensor_product(u,v,w);
  
  for (i = 0; i < n; i++)
    {
      for (j = 0; j < n; j++)
        printf ("%7.3f", gsl_matrix_get (w, i, j));
      printf("\n");
    }
  */

  gsl_matrix_free(A);
  return 0;

} 
