#include <cusp/csr_matrix.h>
#include <cusp/print.h>
#include <cusp/krylov/cg.h>

// where to perform the computation
//typedef cusp::device_memory MemorySpace;
typedef cusp::host_memory MemorySpace;

// which floating point type to use
typedef float ValueType;
typedef int   IndexType;

//int main(void)
// function called from main fortran program
extern "C" void cusp_cg_(int *api_ia, int *api_ja, float *apr_a, int *iparam, float *rparam, float *apr_x, float *apr_y)
{
    int nR  = iparam[0];
    int nC  = iparam[1];
    int nnz = iparam[2];
    int niter = iparam[3];

    float tol = rparam[0];

    // allocate storage for (nR,nC) matrix with nnz nonzeros
    cusp::csr_matrix<IndexType,ValueType,MemorySpace> A(nR,nC,nnz);

    // initialize matrix entries on host
    int i;
    for (i=0; i<nR+1; i++){
    	A.row_offsets[i] = api_ia[i];
    }
    for (i=0; i<nnz; i++){
    	A.column_indices[i] = api_ja[i]; A.values[i] = apr_a[i]; 
    }

    // get the vector x
    cusp::array1d<ValueType, MemorySpace> x(nC);
    cusp::array1d<ValueType, MemorySpace> y(nR);
    for (i=0; i<nR; i++){
    	y[i] = apr_y[i];
    }

    // set stopping criteria:
//    cusp::verbose_monitor<ValueType> monitor(y, niter, tol);
    cusp::default_monitor<ValueType> monitor(y, niter, tol);

    // solve the linear system A * x = b with the Conjugate Gradient method
    // set preconditioner (identity)
    cusp::identity_operator<ValueType, MemorySpace> M(A.num_rows, A.num_rows);
    cusp::krylov::cg(A, x, y, monitor, M);

    for (i=0; i<nR; i++){
    	apr_x[i]=x[i];
    }

    return ;
}

