program test_paralution_solver

  use, intrinsic :: ISO_C_BINDING, only : C_INT,       &
                                          C_PTR,       &
                                          C_DOUBLE,    &
                                          C_CHAR,      &
                                          C_NULL_CHAR, &
                                          C_LOC

  use sll_sparse_matrix_module

  implicit none

  interface

    subroutine paralution_fortran_solve_csr( n, m, nnz, solver,       &
                                             mformat, preconditioner, &
                                             pformat, rows_ptr, cols, &
                                             rval, rhs, atol, rtol,   &
                                             div, maxiter, basis,     &
                                             p, q, x, iter, resnorm,  &
                                             ierr ) BIND(C)

    use, intrinsic :: ISO_C_BINDING, only : C_INT, C_PTR, C_DOUBLE, C_CHAR

    integer(kind=C_INT), value, intent(in)  :: n, m, nnz, maxiter, basis, p, q
    real(kind=C_DOUBLE), value, intent(in)  :: atol, rtol, div
    integer(kind=C_INT),        intent(out) :: iter, ierr
    real(kind=C_DOUBLE),        intent(out) :: resnorm
    type(C_PTR),         value, intent(in)  :: rows_ptr, cols, rval, rhs
    type(C_PTR),         value              :: x
    character(kind=C_CHAR)                  :: solver
    character(kind=C_CHAR)                  :: mformat
    character(kind=C_CHAR)                  :: preconditioner
    character(kind=C_CHAR)                  :: pformat

    end subroutine paralution_fortran_solve_csr

  end interface

  integer, parameter    :: infile = 10
  integer(kind=C_INT)   :: n, m, nnz, i, j, iter, ierr
  real(kind=C_DOUBLE)   :: resnorm
  integer, dimension(8) :: tbegin, tend
  real(kind=8)          :: tsolver
  integer               :: nhinv

  logical               :: sym = .false.

  real(kind=C_DOUBLE), allocatable, target :: x(:)
  real(kind=C_DOUBLE), allocatable, target :: rhs(:)

  character(len=10)  :: rep
  character(len=7)   :: field
  character(len=19)  :: symm
  character(len=128) :: arg

  type(sll_csr_matrix), pointer :: mat

  nhinv = 500
  n=(nhinv-1)**2
  allocate(mat)
  allocate (mat%val(5*N),mat%col_ind(5*N),mat%row_ptr(N+1),rhs(N),x(N))

  ! Allocate memory for CSR format specific arrays
  mat%num_rows = n
  mat%num_cols = m
  mat%num_nz   = nnz

  ! Allocate and initialize rhs and solution vector
  !do i = 1, n
  !  rhs(i) = 1._C_DOUBLE
  !  x(i)   = 0._C_DOUBLE
  !end do
  rhs = 1.
  x   = 0.

  call uni2d(mat, rhs)

  ! Print L2 norm of solution vector
  write(*,fmt='(A,F0.2)') '(Fortran) Initial L2 Norm(x) = ', sqrt( sum( x**2 ) )

  call date_and_time(values = tbegin)

  ! Run paralution C function for CSR matrices
  ! Doing a GMRES with MultiColored ILU(1,2) preconditioner
  ! Check paralution documentation for a detailed argument explanation
!  call paralution_fortran_solve_csr(       &
!    n,                                     &
!    m,                                     &
!    nnz,                                   &
!    'CG' // C_NULL_CHAR,                   &
!    'CSR' // C_NULL_CHAR,                  & 
!    'MultiColoredILU' // C_NULL_CHAR,      &
!    'CSR' // C_NULL_CHAR,                  &
!    C_LOC(rows),                           &
!    C_LOC(cols),                           &
!    C_LOC(rval),                           &
!    C_LOC(rhs),                            &
!    1e-15_C_DOUBLE,                        &
!    1e-8_C_DOUBLE,                         &
!    1e+8_C_DOUBLE,                         &
!    5000,                                  &
!    30,                                    &
!    0,                                     &
!    1,                                     &
!    C_LOC(x),                              &
!    iter,                                  &
!    resnorm,                               &
!    ierr )
!
!  call date_and_time(values = tend)
!
!  tbegin = tend - tbegin
!  tsolver = 0.001 * tbegin(8) + tbegin(7) + 60 * tbegin(6) + 3600 * tbegin(5)
!  write(*,fmt='(A,F0.2,A)') '(Fortran) Solver ended after ', tsolver,'sec.'
!
!  ! Print solver details
!  if ( ierr .eq. 0 ) then
!    write(*,fmt='(A,I0,A,E11.5,A)') '(Fortran) Solver took ', iter, ' iterations with residual norm ', resnorm, '.'
!    write(*,fmt='(A,F0.2)') '(Fortran) Final L2 Norm(x)   = ', sqrt( sum( x**2 ) )
!  else
!    write(*,fmt='(A,I0)') '(Fortran) Solver returned status code ', ierr
!  end if
!
!  do i=1,n
!    write(10,*) x(i)
!  end do
!
!  deallocate( rows, cols, rval, rhs, x )

end program test_paralution_solver

