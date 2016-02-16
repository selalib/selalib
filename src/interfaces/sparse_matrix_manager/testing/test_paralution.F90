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
integer(kind=C_INT)   :: n, m, nnz, iter, ierr
real(kind=C_DOUBLE)   :: resnorm
integer, dimension(8) :: tbegin, tend
real(kind=8)          :: tsolver
integer               :: nhinv
integer               :: n2

real(kind=C_DOUBLE), allocatable, target :: x(:)
real(kind=C_DOUBLE), allocatable, target :: rhs(:)

type(sll_paralution_solver) :: mat


nhinv = 400
n2    = (nhinv-1)**2
n     = n2
m     = n2
nnz   = 5*n2

allocate(mat%val(nnz))
allocate(mat%col_ind(nnz))
allocate(mat%row_ptr(n+1))
allocate(rhs(n))
allocate(x(n))

! Allocate memory for CSR format specific arrays
mat%num_rows = nhinv-1
mat%num_cols = nhinv-1
mat%num_nz   = nnz

rhs = 1.0_f64
x   = 0.0_f64

call uni2d(mat, rhs)

! Print L2 norm of solution vector
write(*,fmt='(A,F0.2)') '(Fortran) Initial L2 Norm(x) = ', sqrt( sum( x**2 ) )

call date_and_time(values = tbegin)

! Run paralution C function for CSR matrices
! Doing a GMRES with MultiColored ILU(1,2) preconditioner
! Check paralution documentation for a detailed argument explanation
call paralution_fortran_solve_csr(       &
  n,                                     &
  m,                                     &
  nnz,                                   &
  'CG' // C_NULL_CHAR,                   &
  'CSR' // C_NULL_CHAR,                  & 
  'MultiColoredILU' // C_NULL_CHAR,      &
  'CSR' // C_NULL_CHAR,                  &
  C_LOC(mat%row_ptr(1)),                 &
  C_LOC(mat%col_ind(1)),                 &
  C_LOC(mat%val(1)),                     &
  C_LOC(rhs),                            &
  1e-15_C_DOUBLE,                        &
  1e-8_C_DOUBLE,                         &
  1e+8_C_DOUBLE,                         &
  5000,                                  &
  30,                                    &
  0,                                     &
  1,                                     &
  C_LOC(x),                              &
  iter,                                  &
  resnorm,                               &
  ierr )

call date_and_time(values = tend)

tbegin = tend - tbegin
tsolver = 0.001_f64 * tbegin(8) + tbegin(7) + 60 * tbegin(6) + 3600 * tbegin(5)
write(*,fmt='(A,F0.2,A)') '(Fortran) Solver ended after ', tsolver,'sec.'

! Print solver details
if ( ierr .eq. 0 ) then
  write(*,fmt='(A,I0,A,E12.5,A)') '(Fortran) Solver took ', iter, ' iterations with residual norm ', resnorm, '.'
  write(*,fmt='(A,F0.2)') '(Fortran) Final L2 Norm(x)   = ', sqrt( sum( x**2 ) )
else
  write(*,fmt='(A,I0)') '(Fortran) Solver returned status code ', ierr
end if

write(*,*) "error =", maxval(abs(x-1))

deallocate( rhs, x )

end program test_paralution_solver
