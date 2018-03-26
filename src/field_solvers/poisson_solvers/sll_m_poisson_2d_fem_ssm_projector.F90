module sll_m_poisson_2d_fem_ssm_projector
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_assert.h"

  use sll_m_working_precision, only: f64

  implicit none

  public :: sll_t_poisson_2d_fem_ssm_projector

  private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  ! Working precision
  integer, parameter :: wp = f64

  type :: sll_t_poisson_2d_fem_ssm_projector

    integer :: n1
    integer :: n2

    ! Temporary storage needed for projection
    real(wp), allocatable :: Qp_temp(:,:)
    real(wp), allocatable :: L(:,:)
    real(wp), allocatable :: Lt(:,:) ! transpose

  contains

    procedure :: init         => s_poisson_2d_fem_ssm_projector__init
    procedure :: change_basis => s_poisson_2d_fem_ssm_projector__change_basis
    procedure :: free         => s_poisson_2d_fem_ssm_projector__free

  end type sll_t_poisson_2d_fem_ssm_projector

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
contains
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  ! Initializer
  subroutine s_poisson_2d_fem_ssm_projector__init( self, n1, n2, L )
    class(sll_t_poisson_2d_fem_ssm_projector), intent(inout) :: self
    integer                                  , intent(in   ) :: n1
    integer                                  , intent(in   ) :: n2
    real(wp)                                 , intent(in   ) :: L(:,:) ! matrix of barycentric coordinates

    integer :: nn

    self % n1 = n1
    self % n2 = n2

    nn = (n1-2)*n2

    ! Allocate temporary storage
    allocate( self % Qp_temp( 2*n2, 3 ) )
    allocate( self % L ( size(L,1), size(L,2) ) )
    allocate( self % Lt( size(L,2), size(L,1) ) )

    SLL_ASSERT( size(L,1) == 2*n2  )
    SLL_ASSERT( size(L,2) == 3     )

    self % L  = L
    self % Lt = transpose( L )

  end subroutine s_poisson_2d_fem_ssm_projector__init

  ! Change basis: C1 projection of stiffness and mass matrices
  ! NOTE: 'self' has intent(inout) because temporary storage has to be assigned
  subroutine s_poisson_2d_fem_ssm_projector__change_basis( self, Q, Qp )
    class(sll_t_poisson_2d_fem_ssm_projector), intent(inout) :: self
    real(wp)                                 , intent(in   ) :: Q(:,:)
    real(wp)                                 , intent(inout) :: Qp(:,:)

    integer :: nn

    associate( n1 => self % n1, n2 => self % n2 )

      nn = (n1-2)*n2

      ! Checks
      SLL_ASSERT( size(Q ,1) == n1*n2 )
      SLL_ASSERT( size(Q ,2) == n1*n2 )
      SLL_ASSERT( size(Qp,1) == 3+nn  )
      SLL_ASSERT( size(Qp,2) == 3+nn  )

      ! Fill block 1: 3 x 3
      self % Qp_temp = matmul( Q(1:2*n2,1:2*n2), self % L )
      Qp(1:3,1:3)    = matmul( self % Lt, self % Qp_temp )

      ! Fill block 2: 3 x (n1-2)*n2
      Qp(1:3,4:3+nn) = matmul( self % Lt, Q(1:2*n2,2*n2+1:n1*n2) )

      ! Fill block 3: (n1-2)*n2 x 3
      Qp(4:3+nn,1:3) = matmul( Q(2*n2+1:n1*n2,1:2*n2), self % L )

      ! Fill block 4: (n1-2)*n2 x (n1-2)*n2
      Qp(4:3+nn,4:3+nn) = Q(2*n2+1:n1*n2,2*n2+1:n1*n2)

    end associate

  end subroutine s_poisson_2d_fem_ssm_projector__change_basis

  ! Deallocate allocatables
  subroutine s_poisson_2d_fem_ssm_projector__free( self )
    class(sll_t_poisson_2d_fem_ssm_projector), intent(inout) :: self

    deallocate( self % Qp_temp )
    deallocate( self % L       )
    deallocate( self % Lt      )

  end subroutine s_poisson_2d_fem_ssm_projector__free

end module sll_m_poisson_2d_fem_ssm_projector
