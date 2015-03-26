program test_paralution_solver

  use, intrinsic :: ISO_C_BINDING, only : C_INT,       &
                                          C_PTR,       &
                                          C_DOUBLE,    &
                                          C_CHAR,      &
                                          C_NULL_CHAR, &
                                          C_LOC

  implicit none

  interface

    subroutine paralution_fortran_solve_coo( n, m, nnz, solver,       &
                                             mformat, preconditioner, &
                                             pformat, rows, cols,     &
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

    end subroutine paralution_fortran_solve_coo

  end interface

  !> @brief type for CSR format
  type sll_coo_matrix
    integer(kind=C_INT)              :: num_rows    !< rows, public
    integer(kind=C_INT)              :: num_cols    !< columns
    integer(kind=C_INT)              :: num_nz      !< non zeros
    integer(kind=C_INT), allocatable :: row_ind(:)
    integer(kind=C_INT), allocatable :: col_ind(:)
    real(kind=C_DOUBLE), allocatable :: val(:)
  end type sll_coo_matrix

  integer, parameter    :: infile = 10
  integer(kind=C_INT)   :: n, m, nnz, fnz, i, j, iter, ierr
  real(kind=C_DOUBLE)   :: resnorm
  integer, dimension(8) :: tbegin, tend
  real(kind=8)          :: tsolver

  logical               :: sym = .false.

  real(kind=C_DOUBLE), allocatable, target :: rhs(:), x(:)

  character(len=10)  :: rep
  character(len=7)   :: field
  character(len=19)  :: symm
  character(len=128) :: arg

  type(sll_coo_matrix) :: mat


  nnz = fnz
  if ( symm .eq. 'symmetric' .or. symm .eq. 'hermitian' ) then
    nnz = 2 * ( fnz - n ) + n
    sym = .true.
  end if

  ! Allocate memory for COO format specific arrays
  allocate( mat%row_ind(nnz), mat%col_ind(nnz), mat%val(nnz))

  ! Allocate and initialize rhs and solution vector
  allocate( rhs(n), x(n) )
  do i = 1, n
    rhs(i) = 1._C_DOUBLE
    x(i)   = 0._C_DOUBLE
  end do

  ! Print L2 norm of solution vector
  write(*,fmt='(A,F0.2)') '(Fortran) Initial L2 Norm(x) = ', sqrt( sum( x**2 ) )

  call date_and_time(values = tbegin)

!  ! Run paralution C function for COO matrices
!  ! Doing a GMRES with MultiColored ILU(1,2) preconditioner
!  ! Check paralution documentation for a detailed argument explanation
!  call paralution_fortran_solve_coo(       &
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

subroutine uni2d(m,f,a,ja,ia)
!
! Fill a matrix in COO format corresponding to a constant coefficient
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
    ia(l)=l+1
  end do
end do

return
end subroutine uni2D
