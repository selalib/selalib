program test_paralution_solver

#include "sll_working_precision.h"

use, intrinsic :: ISO_C_BINDING, only : C_INT,       &
                                        C_PTR,       &
                                        C_DOUBLE,    &
                                        C_CHAR,      &
                                        C_NULL_CHAR, &
                                        C_LOC

use sll_m_paralution

implicit none


integer, parameter    :: infile = 10
integer(kind=C_INT)   :: n, m, nnz
integer, dimension(8) :: tbegin, tend
real(kind=8)          :: tsolver
integer               :: nhinv
integer               :: n2

real(kind=C_DOUBLE), allocatable, target :: x(:)
real(kind=C_DOUBLE), allocatable, target :: rhs(:)

type(paralution_solver) :: mat


nhinv = 400
n2    = (nhinv-1)**2
n     = n2
m     = n2
nnz   = 5*n2

allocate(rhs(n))
allocate(x(n))

rhs = 1.0_f64
x   = 0.0_f64

! Print L2 norm of solution vector
write(*,fmt='(A,F0.2)') '(Fortran) Initial L2 Norm(x) = ', sqrt( sum( x**2 ) )

call initialize(mat,n,nnz)

call uni2d(nhinv-1, mat, rhs)

call date_and_time(values = tbegin)

call solve(mat, rhs, x)

call date_and_time(values = tend)

tbegin = tend - tbegin
tsolver = 0.001_f64 * tbegin(8) + tbegin(7) + 60 * tbegin(6) + 3600 * tbegin(5)
write(*,fmt='(A,F0.2,A)') '(Fortran) Solver ended after ', tsolver,'sec.'


write(*,fmt='(A,F0.2)') '(Fortran) Final L2 Norm(x)   = ', sqrt( sum( x**2 ) )
write(*,*) "error =", maxval(abs(x-1))

deallocate( rhs, x )

call delete(mat)

contains

!> @brief
!> Test function to initialize a CSR matrix
!> @details
!> Fill a matrix in CSR format corresponding to a constant coefficient
!> five-point stencil on a square grid
!> RHS is set to get solution = 1
subroutine uni2d(m, mat,f)
type(paralution_solver) :: mat
real(kind=C_DOUBLE)     :: f(:)
integer                 :: i, j, k, l, m

real(kind(0d0)), parameter :: zero =  0.0d0
real(kind(0d0)), parameter :: cx   = -1.0d0
real(kind(0d0)), parameter :: cy   = -1.0d0
real(kind(0d0)), parameter :: cd   =  4.0d0


k=0
l=0
mat%row_ptr(1)=1
do i=1,m
  do j=1,m
    k=k+1
    l=l+1
    mat%val(l)=cd
    mat%col_ind(l)=k
    f(k)=zero
    if(j < m) then
       l=l+1
       mat%val(l)=cx
       mat%col_ind(l)=k+1
      else
       f(k)=f(k)-cx
    end if
    if(i < m) then
       l=l+1
       mat%val(l)=cy
       mat%col_ind(l)=k+m
      else
       f(k)=f(k)-cy
    end if
    if(j > 1) then
       l=l+1
       mat%val(l)=cx
       mat%col_ind(l)=k-1
      else
       f(k)=f(k)-cx
    end if
    if(i >  1) then
       l=l+1
       mat%val(l)=cy
       mat%col_ind(l)=k-m
      else
       f(k)=f(k)-cy
    end if
    mat%row_ptr(k+1)=l+1
  end do
end do

mat%num_nz = l

return
end subroutine uni2D


end program test_paralution_solver
