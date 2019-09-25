program test_dot_product_stencil
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  use sll_m_working_precision, only: f64

  use sll_m_vector_space_base, only: sll_c_vector_space

  use sll_m_vector_space_real_array_1d, only: sll_t_vector_space_real_array_1d

  use sll_m_vector_space_real_array_2d, only: sll_t_vector_space_real_array_2d

  use sll_m_linear_operator_matrix_stencil_to_stencil, only: sll_t_linear_operator_matrix_stencil_to_stencil

  implicit none
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  ! Working precision
  integer, parameter :: wp = f64

  integer  :: n1, n2, p1, p2
  integer  :: j, i1, i2, j1, j2, k1, k2, l1, l2

  ! Matrices
  real(wp), allocatable :: Ad(:,:)     ! dense
  real(wp), allocatable :: As(:,:,:,:) ! stencil

  ! Vectors
  real(wp), allocatable :: Vd(:)   ! dense
  real(wp), allocatable :: Wd(:)   ! dense
  real(wp), allocatable :: Rd(:)   ! dense
  real(wp), allocatable :: Vs(:,:) ! stencil

  ! Linear operators
  type(sll_t_linear_operator_matrix_stencil_to_stencil) :: As_linop

  ! Vector spaces
  type(sll_t_vector_space_real_array_2d) :: Vs_vecsp

  ! Vector spaces for results of dot product
  class(sll_c_vector_space), allocatable :: Ws_vecsp

  real(wp) :: error
  real(wp), parameter :: tol = 1.0e-13_wp

  ! For CTest
  logical :: passed
  logical :: success

  passed = .true.

  n1 = 10
  n2 = 20

  p1 = 3
  p2 = 5
 
  ! Allocations of dense matrices and vectors
  allocate( Ad( n1*n2, n1*n2 ) )
  allocate( Vd( n1*n2 ) )
  allocate( Wd( n1*n2 ) )
  allocate( Rd( n1*n2 ) )
  ! Allocations of stencil matrices and vectors
  allocate( As( -p1:p1, -p2:p2, n1, n2 ) )
  allocate( Vs( 1-p1:n1+p1, 1-p2:n2+p2 ) )

  ! Initialize stencil matrix As
  As = 0.0_wp
  do i2 = 1, n2
    do i1 = 1, n1
      do k2 = -p2, p2
        do k1 = -p1, p1
          call random_number( As(k1,k2,i1,i2) )
        end do
      end do
    end do
  end do

  ! Construct linear operator from stencil matrix As
  call As_linop % init( n1, n2, p1, p2 )
  As_linop % A = As

  ! Convert stencil matrix to dense matrix
  Ad = 0.0_wp
  call As_linop % to_array( Ad )

  ! Initialize stencil vector Vs
  do l2 = 1-p2, n2+p2
    do l1 = 1-p1, n1+p1
      call random_number( Vs(l1,l2) )
    end do
  end do

  ! Construct vector space from stencil vector Vs
  associate( l1 => lbound( Vs, 1 ), &
             l2 => lbound( Vs, 2 ), &
             u1 => ubound( Vs, 1 ), &
             u2 => ubound( Vs, 2 ) )
    allocate( Vs_vecsp % array(l1:u1,l2:u2), source=Vs ) 
  end associate

  ! Convert stencil vector to dense vector
  do j2 = 1, n2
    do j1 = 1, n1
      j = (j1-1) * n2 + j2
      Vd(j) = Vs(j1,j2)
    end do
  end do

  ! Compute dot product in dense format
  Wd = matmul( Ad, Vd )

  ! Construct vector space for stencil result
  allocate( sll_t_vector_space_real_array_2d :: Ws_vecsp )
  call Vs_vecsp % source( Ws_vecsp )

  ! Compute dot product in stencil format
  call As_linop % dot( Vs_vecsp, Ws_vecsp )

  ! Copy result from stencil to dense format
  select type( Ws_vecsp )

  type is( sll_t_vector_space_real_array_2d )

    do j2 = 1, n2
      do j1 = 1, n1
        j = (j1-1) * n2 + j2
        Rd(j) = Ws_vecsp % array(j1,j2)
      end do
    end do

  end select

  ! Compare results in dense format
  error = 0.0_wp
  do j = 1, n1*n2
    error = max( error, abs( Wd(j) - Rd(j) ) )
  end do

  write(*,'(/a)') " *********************************************"
  write(*,'(a)' ) " Test dot product in dense and stencil formats"
  write(*,'(a)' ) " *********************************************"
  write(*,'(/a,es8.2/)') " L_inf norm of error = ", error

  if ( error <= tol ) then
    success = .true.
  else
    success = .false.
    write(*,'(a/)') " Test FAILED"
  end if

  passed = passed .and. success

  ! Check if all tests passed
  if ( passed ) write(*,'(a/)') "PASSED"

  ! Deallocations
  deallocate( Ad )
  deallocate( As )
  deallocate( Vd )
  deallocate( Wd )
  deallocate( Rd )
  deallocate( Vs )

  ! Free objects
  call As_linop % free()

end program test_dot_product_stencil
