module sll_m_paralution
use, intrinsic :: ISO_C_BINDING

implicit none

type sll_paralution_solver

  integer(kind=C_INT)          :: num_rows
  integer(kind=C_INT)          :: num_cols
  integer(kind=C_INT)          :: num_nz
  integer(kind=C_INT), pointer :: row_ptr(:)
  integer(kind=C_INT), pointer :: col_ind(:)
  real(kind=C_DOUBLE), pointer :: val(:)

end type sll_paralution_solver

interface

  subroutine paralution_fortran_solve_csr( n, m, nnz, solver,       &
                                           mformat, preconditioner, &
                                           pformat, rows, cols, &
                                           rval, rhs, atol, rtol,   &
                                           div, maxiter, basis,     &
                                           p, q, x, iter, resnorm,  &
                                           ierr ) BIND(C)

  use, intrinsic :: ISO_C_BINDING, only : C_INT, C_PTR, C_DOUBLE, C_CHAR

  integer(kind=C_INT), value, intent(in)  :: n, m, nnz, maxiter, basis, p, q
  real(kind=C_DOUBLE), value, intent(in)  :: atol, rtol, div
  integer(kind=C_INT),        intent(out) :: iter, ierr
  real(kind=C_DOUBLE),        intent(out) :: resnorm
  type(C_PTR),         value, intent(in)  :: rows, cols, rval, rhs
  type(C_PTR),         value              :: x
  character(kind=C_CHAR)                  :: solver
  character(kind=C_CHAR)                  :: mformat
  character(kind=C_CHAR)                  :: preconditioner
  character(kind=C_CHAR)                  :: pformat

  end subroutine paralution_fortran_solve_csr

end interface

contains

!> @brief
!> Test function to initialize a CSR matrix
!> @details
!> Fill a matrix in CSR format corresponding to a constant coefficient
!> five-point stencil on a square grid
!> RHS is set to get solution = 1
subroutine uni2d(mat,f)
type(sll_paralution_solver) :: mat
real(kind=C_DOUBLE)         :: f(:)
integer                     :: i, j, k, l, m

real(kind(0d0)), parameter :: zero =  0.0d0
real(kind(0d0)), parameter :: cx   = -1.0d0
real(kind(0d0)), parameter :: cy   = -1.0d0
real(kind(0d0)), parameter :: cd   =  4.0d0

m = mat%num_rows

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

end module sll_m_paralution

