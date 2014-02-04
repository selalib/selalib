#include <cusp/csr_matrix.h>
#include <cusp/print.h>
#include <cusp/multiply.h>
#include <iostream>

//int main(void)
// function called from main fortran program
extern "C" void cusp_csr_mv_(int *api_ia, int *api_ja, float *apr_a, int *iparam, float *apr_x, float *apr_y)
{
    int nR  = iparam[0];
    int nC  = iparam[1];
    int nnz = iparam[2];

    // allocate storage for (nR,nC) matrix with nnz nonzeros
    cusp::csr_matrix<int,float,cusp::host_memory> A(nR,nC,nnz);

    // initialize matrix entries on host
    int i;
    for (i=0; i<nR+1; i++){
    	A.row_offsets[i] = api_ia[i];
    }
    for (i=0; i<nnz; i++){
    	A.column_indices[i] = api_ja[i]; A.values[i] = apr_a[i]; 
    }
    
    // get the vector x
    cusp::array1d<float, cusp::host_memory> x(nC);
    cusp::array1d<float, cusp::host_memory> y(nR);
    for (i=0; i<nC; i++){
    	x[i] = apr_x[i];
    }

    // compute y = A* x 
    cusp::multiply(A, x, y);

    // print matrix entries
    cusp::print(A);
    cusp::print(x);
    cusp::print(y);

    for (i=0; i<nR; i++){
    	apr_y[i]=y[i];
    }

    return ;
}

