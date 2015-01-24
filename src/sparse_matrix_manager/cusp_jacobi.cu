#include <cusp/csr_matrix.h>
#include <cusp/print.h>
#include <cusp/relaxation/jacobi.h>

// where to perform the computation
typedef cusp::device_memory MemorySpace;

// which floating point type to use
typedef float ValueType;

//int main(void)
// function called from main fortran program
extern "C" void cusp_jacobi_(int *api_ia, int *api_ja, float *apr_a, int *iparam, float *rparam, float *apr_x, float *apr_y)
{
    int nR  = iparam[0];
    int nC  = iparam[1];
    int nnz = iparam[2];

    float omega = rparam[0];

    typedef typename Matrix::memory_space Space;

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
    for (i=0; i<nR; i++){
    	y[i] = apr_y[i];
    }

    Matrix J(A);
    cusp::relaxation::jacobi<float, Space> relax(J, omega);
    relax(J, y, x);

    for (i=0; i<nR; i++){
    	apr_x[i]=x[i];
    }

    return ;
}

