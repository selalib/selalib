program test_csr_matrix
#include "sll_working_precision.h"
#include "sll_memory.h"
use sll_m_sparse_matrix
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
tol=1.d-6
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

call uni_laplace_2d(nhinv-1,f,mat%val,mat%col_ind,mat%row_ptr)

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
subroutine uni_laplace_2d(m,f,a,ja,ia)
!
! Fill a matrix in CSR format corresponding to a constant coefficient
! five-point stencil on a square grid
!
implicit none
real (kind(0d0)) :: f(*),a(*)
integer :: m,ia(*),ja(*)
integer :: k,l,i,j
real (kind(0d0)), parameter :: zero=0.0d0,cx=-1.0d0,cy=-1.0d0, cd=4.0d0
!
k=0
l=0
ia(1)=1
do i=1,m
  do j=1,m
    k=k+1
    l=l+1
    a(l)=cd
    ja(l)=k
    f(k)=zero
    if(j < m) then
       l=l+1
       a(l)=cx
       ja(l)=k+1
      else
       f(k)=f(k)-cx
    end if
    if(i < m) then
       l=l+1
       a(l)=cy
       ja(l)=k+m
      else
       f(k)=f(k)-cy
    end if
    if(j > 1) then
       l=l+1
       a(l)=cx
       ja(l)=k-1
      else
       f(k)=f(k)-cx
    end if
    if(i >  1) then
       l=l+1
       a(l)=cy
       ja(l)=k-m
      else
       f(k)=f(k)-cy
    end if
    ia(k+1)=l+1
  end do
end do

return
end subroutine uni_laplace_2D


