program test_conjugate_gradient
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  use sll_m_working_precision, only: f64

  use sll_m_vector_space_real_array_1d, only: sll_t_vector_space_real_array_1d

  use sll_m_linear_operator_matrix_dense, only: sll_t_linear_operator_matrix_dense

  use sll_m_conjugate_gradient, only: sll_t_conjugate_gradient

  implicit none
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  ! Working precision
  integer, parameter :: wp = f64

  integer  :: i, j, n, p, m, ki, kj, li, lj
  real(wp) :: Amin, Amax

  ! 2D and 1D real arrays
  real(wp), allocatable :: A(:,:)
  real(wp), allocatable :: La(:,:)
  real(wp), allocatable :: Id(:,:)
  real(wp), allocatable :: x(:)
  real(wp), allocatable :: b(:)
  real(wp), allocatable :: z(:) ! 1D array of zeros

  ! 2D and 1D equivalent vector space real arrays
  type(sll_t_vector_space_real_array_1d) :: xx
  type(sll_t_vector_space_real_array_1d) :: bb

  ! Linear operator constructed from A
  type(sll_t_linear_operator_matrix_dense) :: A_linop_dense

  ! Conjugate gradient solver
  type(sll_t_conjugate_gradient) :: conjugate_gradient

  real(wp) :: error
  real(wp), parameter :: tol = 1.0e-14_wp
  logical , parameter :: verbose = .true.

  ! For CTest
  logical :: passed
  logical :: success

  passed = .true.

  !-----------------------------------------------------------------------------
  ! Test #1: diagonal matrices of increasing size
  !-----------------------------------------------------------------------------

  write(*,*)
  write(*,'(a)') " **************************"
  write(*,'(a)') " Test #1: diagonal matrices"
  write(*,'(a)') " **************************"

  do m = 1, 10

    n = 2**m

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
    b = matmul( A, x )

    ! Array of zeros
    z(:) = 0.0_wp

    ! Construct linear operator from matrix A
    call A_linop_dense % init( A )

    ! Construct vector space from vector b
    allocate( bb % array( size( b ) ), source=b )

    ! Construct vector space object for solution
    allocate( xx % array( size( z ) ), source=z )

    ! Initialize conjugate gradient solver
    call conjugate_gradient % init( tol=tol, verbose=verbose, template_vector=xx )

    ! Solve linear system Ax=b for x using conjugate gradient method
    call conjugate_gradient % solve( &
      A       = A_linop_dense, &
      b       = bb           , &
      x       = xx )

    ! Free conjugate gradient solver
    call conjugate_gradient % free()

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
    call A_linop_dense % free()
    deallocate( bb % array )
    deallocate( xx % array )
    deallocate( A, x, b, z )

  end do

  !-----------------------------------------------------------------------------
  ! Test #2: Laplacian matrices
  !-----------------------------------------------------------------------------

  write(*,*)
  write(*,'(a)') " ***************************"
  write(*,'(a)') " Test #2: Laplacian matrices"
  write(*,'(a)') " ***************************"

  do m = 3, 12

    ! Auxiliary matrices
    allocate( La(m,m) )
    allocate( Id(m,m) )

    n = m**2

    allocate( A(n,n) )
    allocate( x(n) )
    allocate( b(n) )
    allocate( z(n) )

    ! Output
    write(*,*)
    write(*,'(a,i0,a,i0)') " Matrix A: ", size(A,1), " x ", size(A,2)

    ! Initialize Laplacian (m x m) matrix La
    La(:,:) = 0.0_wp
    ! diagonal
    do i = 1, m
      La(i,i) = 2.0_wp
    end do
    ! upper-diagonal
    do i = 1, m-1
      La(i,i+1) = -1.0_wp
    end do
    ! lower-diagonal
    do i = 2, m
      La(i,i-1) = -1.0_wp
    end do

    ! Initialize identity (m x m) matrix Id
    Id(:,:) = 0.0_wp
    do i = 1, m
      Id(i,i) = 1.0_wp
    end do

    ! Initialize (m^2 x m^2) matrix A = La x Id + Id x La (Kronecker products)
    do j = 1, m
      do i = 1, m
        lj = 1
        do kj = (j-1)*m+1, j*m
          li = 1
          do ki = (i-1)*m+1, i*m
            A(ki,kj) = La(i,j)*Id(li,lj) + Id(i,j)*La(li,lj)
            li = li + 1
          end do
          lj = lj + 1
        end do
      end do
    end do

    ! Initialize vector x
    do i = 1, n
      call random_number( x(i) )
      x(i) = -0.5_wp + x(i) ! -0.5 <= x < 0.5
    end do

    ! Compute b=Ax
    b = matmul( A, x )

    ! Array of zeros
    z(:) = 0.0_wp

    ! Construct linear operator from matrix A
    call A_linop_dense % init( A )

    ! Construct vector space from vector b
    allocate( bb % array( size( b ) ), source=b )

    ! Construct vector space object for solution
    allocate( xx % array( size( z ) ), source=z )

    ! Initialize conjugate gradient solver
    call conjugate_gradient % init( tol=tol, verbose=verbose, template_vector=xx )

    ! Solve linear system Ax=b for x using conjugate gradient method
    call conjugate_gradient % solve( &
      A       = A_linop_dense, &
      b       = bb           , &
      x       = xx )

    ! Free conjugate gradient solver
    call conjugate_gradient % free()

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
    call A_linop_dense % free()
    deallocate( bb % array )
    deallocate( xx % array )
    deallocate( A, La, Id, x, b, z )

  end do

  !-----------------------------------------------------------------------------
  ! Test #3: multi-diagonal matrices
  !-----------------------------------------------------------------------------

  write(*,*)
  write(*,'(a)') " ********************************"
  write(*,'(a)') " Test #3: multi-diagonal matrices"
  write(*,'(a)') " ********************************"

  ! Following Gershgorin circle theorem, we construct A as follows:
  ! - random diagonal elements in a certain positive interval (does not need to be positive maybe)
  ! - p upper and lower diagonal elements: sum of absolute values of non-diagonal elements on each
  !   row kept smaller than value of corresponding diagonal element (=> eigenvalues all positive)
  ! - lower diagonal elements constructed by symmetry

  do n = 3, 12

    allocate( A(n,n) )
    allocate( x(n) )
    allocate( b(n) )
    allocate( z(n) )

    do p = 1, n-2

      ! Output
      write(*,*)
      write(*,'(a,i0,a,i0,a,i0,a,i0,a)') " Matrix A: ", size(A,1), " x ", size(A,2), &
                                         " with 2*", p, "+1=", 2*p+1, " diagonals"

      ! Initialize matrix A diagonal
      A(:,:) = 0.0_wp
      ! diagonal
      do i = 1, n
        call random_number( A(i,i) )
        A(i,i) = 0.5_wp + A(i,i) ! 0.5 <= A(i,i) < 1.5
      end do
      ! upper diagonals
      do j = 1, p
        do i = 1, n-j
          Amin = - A(i,i) / ( 2.0_wp*p )
          Amax =   A(i,i) / ( 2.0_wp*p )
          call random_number( A(i,i+j) )
          A(i,i+j) = Amin + A(i,i+j) * (Amax-Amin)
        end do
      end do
      ! lower diagonals (apply symmetry)
      do j = 1, p
        do i = 1+j, n
          A(i,i-j) = A(i-j,i)
        end do
      end do

      ! Initialize vector x
      do i = 1, n
        call random_number( x(i) )
        x(i) = -0.5_wp + x(i) ! -0.5 <= x < 0.5
      end do

      ! Compute b=Ax
      b = matmul( A, x )

      ! Array of zeros
      z(:) = 0.0_wp

      ! Construct linear operator from matrix A
      call A_linop_dense % init( A )

      ! Construct vector space from vector b
      allocate( bb % array( size( b ) ), source=b )

      ! Construct vector space object for solution
      allocate( xx % array( size( z ) ), source=z )

      ! Initialize conjugate gradient solver
      call conjugate_gradient % init( tol=tol, verbose=verbose, template_vector=xx )

      ! Solve linear system Ax=b for x using conjugate gradient method
      call conjugate_gradient % solve( &
        A       = A_linop_dense, &
        b       = bb           , &
        x       = xx )

      ! Free conjugate gradient solver
      call conjugate_gradient % free()

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
      call A_linop_dense % free()
      deallocate( bb % array )
      deallocate( xx % array )

    end do

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
