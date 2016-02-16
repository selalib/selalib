module sll_m_paralution
#include "sll_working_precision.h"
#include "sll_memory.h"
use, intrinsic :: ISO_C_BINDING

implicit none

type paralution_solver

  integer(kind=C_INT)          :: num_rows
  integer(kind=C_INT)          :: num_cols
  integer(kind=C_INT)          :: num_nz
  integer(kind=C_INT), pointer :: row_ptr(:)
  integer(kind=C_INT), pointer :: col_ind(:)
  real(kind=C_DOUBLE), pointer :: val(:)

end type paralution_solver

interface

  subroutine paralution_init() BIND(C)
  end subroutine paralution_init

  subroutine paralution_stop() BIND(C)
  end subroutine paralution_stop

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

interface initialize
   module procedure init_paralution
end interface initialize

interface solve
   module procedure solve_paralution_with_rhs
end interface solve

interface factorize
   module procedure factorize_paralution
end interface factorize

interface delete
   module procedure free_paralution
end interface delete

contains

subroutine init_paralution(self,n,nnz)

  type(paralution_solver)          :: self
  sll_int32,    intent(in)         :: n
  sll_int32,    intent(in)         :: nnz
  sll_int32                        :: error

  ! Allocate memory for CSR format specific arrays
  self%num_rows = n
  self%num_cols = n
  self%num_nz   = nnz
  SLL_ALLOCATE(self%val(nnz), error)
  SLL_ALLOCATE(self%col_ind(nnz), error)
  SLL_ALLOCATE(self%row_ptr(n+1), error)

  ! Initialize PARALUTION backend
  call paralution_init

end subroutine init_paralution

subroutine factorize_paralution(self)

  type(paralution_solver) :: self

end subroutine factorize_paralution

subroutine solve_paralution_with_rhs(self, rhs, sol)

  type(paralution_solver)  :: self
  real(kind=C_DOUBLE), target :: sol(:)
  real(kind=C_DOUBLE), target :: rhs(:)

  integer(kind=C_INT)   :: iter
  integer(kind=C_INT)   :: ierr
  real(kind=C_DOUBLE)   :: resnorm

  ! Run paralution C function for CSR matrices
  ! Doing a GMRES with MultiColored ILU(1,2) preconditioner
  ! Check paralution documentation for a detailed argument explanation
  call paralution_fortran_solve_csr(       &
    self%num_rows,                         &
    self%num_cols,                         &
    self%num_nz,                           &
    'CG' // C_NULL_CHAR,                   &
    'CSR' // C_NULL_CHAR,                  & 
    'MultiColoredILU' // C_NULL_CHAR,      &
    'CSR' // C_NULL_CHAR,                  &
    C_LOC(self%row_ptr(1)),                &
    C_LOC(self%col_ind(1)),                &
    C_LOC(self%val(1)),                    &
    C_LOC(rhs),                            &
    1e-15_C_DOUBLE,                        &
    1e-8_C_DOUBLE,                         &
    1e+8_C_DOUBLE,                         &
    5000,                                  &
    30,                                    &
    0,                                     &
    1,                                     &
    C_LOC(sol),                            &
    iter,                                  &
    resnorm,                               &
    ierr )

  ! Print solver details
  if ( ierr .eq. 0 ) then
    write(*,fmt='(A,I0,A,E12.5,A)') '(Fortran) Solver took ', iter, ' iterations with residual norm ', resnorm, '.'
  else
    write(*,fmt='(A,I0)') '(Fortran) Solver returned status code ', ierr
  end if

end subroutine solve_paralution_with_rhs

subroutine free_paralution(self)

  type(paralution_solver) :: self

  ! Stop PARALUTION backend
  call paralution_stop

end subroutine free_paralution

end module sll_m_paralution
