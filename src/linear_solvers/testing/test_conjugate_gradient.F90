program test_conjugate_gradient
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  use sll_m_working_precision, only: f64

  use sll_m_vector_space_real_arrays, only: &
    sll_t_vector_space_real_1d, &
    sll_t_vector_space_real_2d

  use sll_m_linear_operator_matrix_dense, only: sll_t_linear_operator_matrix_dense

  use sll_m_conjugate_gradient, only: sll_t_conjugate_gradient

  implicit none
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  ! Working precision
  integer, parameter :: wp = f64

  integer :: i, j, k, n

  ! 2D and 1D real arrays
  real(wp), allocatable :: A(:,:)
  real(wp), allocatable :: x(:)
  real(wp), allocatable :: b(:)
  real(wp), allocatable :: z(:) ! 1D array of zeros

  ! 2D and 1D equivalent vector space real arrays
  type(sll_t_vector_space_real_2d) :: AA
  type(sll_t_vector_space_real_1d) :: xx
  type(sll_t_vector_space_real_1d) :: bb
  type(sll_t_vector_space_real_1d) :: x0

  ! Linear operator constructed from AA
  type(sll_t_linear_operator_matrix_dense) :: AA_linear_operator

  ! Conjugate gradient solver
  type(sll_t_conjugate_gradient) :: conjugate_gradient

  real(wp) :: error
  real(wp), parameter :: tol = 1.0e-14_wp
  logical , parameter :: verbose = .true.

  ! For CTest
  logical :: passed
  logical :: success

  passed  = .true.

  !-----------------------------------------------------------------------------
  ! Test #1: diagonal matrices of increasing size
  !-----------------------------------------------------------------------------

  write(*,*)
  write(*,'(a)') " **************************"
  write(*,'(a)') " Test #1: diagonal matrices"
  write(*,'(a)') " **************************"

  do k = 1, 8

    n = 2**k

    allocate( A(n,n) )
    allocate( x(n) )
    allocate( b(n) )
    allocate( z(n) )

    ! Output
    write(*,*)
    write(*,'(a,i0,a,i0)') " Matrix A: ", size(A,1), " x ", size(A,2)

    ! Initialize matrix A diagonal
    A(:,:) = 0.0_wp
    do i = 1, n
      call random_number( A(i,i) )
      A(i,i) = 0.5_wp + A(i,i) ! 0.5 <= A(i,i) < 1.5
    end do

    ! Initialize vector x
    do i = 1, n
      call random_number( x(i) )
      x(i) = -0.5_wp + x(i) ! -0.5 <= x < 0.5
    end do

    ! Compute b=Ax
    b(:) = 0.0_wp
    do i = 1, n
      do j = 1, n
        b(i) = b(i) + A(i,j)*x(j)
      end do
    end do

    ! Array of zeros
    z(:) = 0.0_wp

    ! Construct linear operator from matrix A
    call AA % attach( A )
    call AA_linear_operator % init( AA )

    ! Construct vector space from vector b
    call bb % attach( b )

    ! Construct vector space object for initial guess
    call x0 % attach( z )

    ! Construct vector space object for solution
    call xx % attach( z )

    ! Solve linear system Ax=b for x using conjugate gradient method
    call conjugate_gradient % solve( &
      A       = AA_linear_operator, &
      b       = bb                , &
      x0      = x0                , &
      tol     = tol               , &
      verbose = verbose           , &
      x       = xx )

    ! Check error and write to output
    error = maxval( abs( xx % array(:) - x(:) ) )
    write(*,*)
    write(*,'(a,es8.2/)') " Maximum absolute error (numerical vs. exact solution): ", error

    if ( error <= tol ) then
      success = .true.
    else
      success = .false.
      write(*,'(a/)') " Test FAILED"
    end if

    passed = passed .and. success

    ! Deallocate allocatables and free objects
    call AA_linear_operator % free()
    call AA % delete()
    call xx % delete()
    call bb % delete()
    call x0 % delete()
    deallocate( A, x, b, z )

  end do

  !-----------------------------------------------------------------------------
  ! Test #2: banded matrices
  !-----------------------------------------------------------------------------

  write(*,*)
  write(*,'(a)') " ************************"
  write(*,'(a)') " Test #2: banded matrices"
  write(*,'(a)') " ************************"

  do k = 1, 8

    n = 2**k

    allocate( A(n,n) )
    allocate( x(n) )
    allocate( b(n) )
    allocate( z(n) )

    ! Output
    write(*,*)
    write(*,'(a,i0,a,i0)') " Matrix A: ", size(A,1), " x ", size(A,2)

    ! Initialize matrix A
    A(:,:) = 0.0_wp
    ! diagonal
    do i = 1, n
      call random_number( A(i,i) )
      A(i,i) = 2.0_wp + A(i,i) ! 2.0 <= A(i,i) < 3.0
    end do
    ! upper-diagonal
    do i = 1, n-1
      A(i,i+1) = -1.0_wp
    end do
    ! lower-diagonal
    do i = 2, n
      A(i,i-1) = -1.0_wp
    end do

    ! Initialize vector x
    do i = 1, n
      call random_number( x(i) )
      x(i) = -0.5_wp + x(i) ! -0.5 <= x < 0.5
    end do

    ! Compute b=Ax
    b(:) = 0.0_wp
    do i = 1, n
      do j = 1, n
        b(i) = b(i) + A(i,j)*x(j)
      end do
    end do

    ! Array of zeros
    z(:) = 0.0_wp

    ! Construct linear operator from matrix A
    call AA % attach( A )
    call AA_linear_operator % init( AA )

    ! Construct vector space from vector b
    call bb % attach( b )

    ! Construct vector space object for initial guess
    call x0 % attach( z )

    ! Construct vector space object for solution
    call xx % attach( z )

    ! Solve linear system Ax=b for x using conjugate gradient method
    call conjugate_gradient % solve( &
      A       = AA_linear_operator, &
      b       = bb                , &
      x0      = x0                , &
      tol     = tol               , &
      verbose = verbose           , &
      x       = xx )

    ! Check error and write to output
    error = maxval( abs( xx % array(:) - x(:) ) )
    write(*,*)
    write(*,'(a,es8.2/)') " Maximum absolute error (numerical vs. exact solution): ", error

    if ( error <= tol ) then
      success = .true.
    else
      success = .false.
      write(*,'(a/)') " !!! Test FAILED !!!"
    end if

    passed = passed .and. success

    ! Deallocate allocatables and free objects
    call AA_linear_operator % free()
    call AA % delete()
    call xx % delete()
    call bb % delete()
    call x0 % delete()
    deallocate( A, x, b, z )

  end do

  !-----------------------------------------------------------------------------
  ! Check if all tests passed
  !-----------------------------------------------------------------------------
  if ( passed ) then
    write(*,'(/a/)') "PASSED"
  else
    write(*,*)
  end if

end program test_conjugate_gradient
