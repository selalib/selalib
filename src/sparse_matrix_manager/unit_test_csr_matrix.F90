program test_csr_matrix
#include "sll_working_precision.h"
#include "sll_memory.h"
use sll_sparse_matrix_module, only: sll_csr_matrix, uni2d
!
!  Solves the discrete Laplacian on the unit square by simple call to agmg.
!  The right-hand-side is such that the exact solution is the vector of all 1.
!
implicit none
real (kind(0d0)),allocatable :: f(:),x(:)
integer :: n,iter,iprint,nhinv
real (kind(0d0)) :: tol
type(sll_csr_matrix), pointer :: mat

!
!       set inverse of the mesh size (feel free to change)
nhinv=500
!
!       maximal number of iterations
iter=50
!
!       tolerance on relative residual norm
tol=1.e-6
!
!       unit number for output messages: 6 => standard output
iprint=6
!
!       generate the matrix in required format (CSR)
!
!         first allocate the vectors with correct size
N=(nhinv-1)**2
allocate(mat)
allocate (mat%val(5*N),mat%col_ind(5*N),mat%row_ptr(N+1),f(N),x(N))
!         next call subroutine to set entries


mat%num_rows = nhinv-1
mat%num_cols = nhinv-1

!call uni2d(nhinv-1,f,mat%val,mat%col_ind,mat%row_ptr)
call uni2d(mat, f)

!
!       call agmg
!         argument 5 (ijob)  is 0 because we want a complete solve
!         argument 7 (nrest) is 1 because we want to use flexible CG
!                            (the matrix is symmetric positive definite)
!
call dagmg(N,mat%val,mat%col_ind,mat%row_ptr,f,x,0,iprint,1,iter,tol)
!
!      uncomment the following lines to write solution on disk for checking
!
!       open(10,file='sol.out',form='formatted')
!       write(10,'(e22.15)') f(1:n)
!       close(10)
end program test_csr_matrix
!----------------------------------------------------------------------

