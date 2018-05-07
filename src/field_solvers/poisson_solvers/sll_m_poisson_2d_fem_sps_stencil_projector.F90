module sll_m_poisson_2d_fem_sps_stencil_projector
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_assert.h"

  use sll_m_working_precision, only: f64

  use sll_m_vector_space_c1_block, only: sll_t_vector_space_c1_block

  use sll_m_linear_operator_matrix_stencil_to_stencil, only: sll_t_linear_operator_matrix_stencil_to_stencil

  use sll_m_linear_operator_matrix_c1_block, only: sll_t_linear_operator_matrix_c1_block

  implicit none

  public :: sll_t_poisson_2d_fem_sps_stencil_projector

  private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  ! Working precision
  integer, parameter :: wp = f64

  type :: sll_t_poisson_2d_fem_sps_stencil_projector

    integer :: n1
    integer :: n2
    integer :: p1
    integer :: p2

    ! Temporary storage needed for projection
!    real(wp), allocatable :: Q_temp(:,:)
    real(wp), allocatable :: Qp_temp(:,:)
    real(wp), allocatable :: L(:,:)
    real(wp), allocatable :: Lt(:,:) ! transpose

  contains

    procedure :: init                    => s_poisson_2d_fem_sps_stencil_projector__init
    procedure :: change_basis_matrix     => s_poisson_2d_fem_sps_stencil_projector__change_basis_matrix
    procedure :: change_basis_vector     => s_poisson_2d_fem_sps_stencil_projector__change_basis_vector
    procedure :: change_basis_vector_inv => s_poisson_2d_fem_sps_stencil_projector__change_basis_vector_inv
    procedure :: free                    => s_poisson_2d_fem_sps_stencil_projector__free

  end type sll_t_poisson_2d_fem_sps_stencil_projector

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
contains
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  ! Initializer
  subroutine s_poisson_2d_fem_sps_stencil_projector__init( self, n1, n2, p1, p2, L )
    class(sll_t_poisson_2d_fem_sps_stencil_projector), intent(inout) :: self
    integer                                          , intent(in   ) :: n1
    integer                                          , intent(in   ) :: n2
    integer                                          , intent(in   ) :: p1
    integer                                          , intent(in   ) :: p2
    real(wp)                                         , intent(in   ) :: L(:,:) ! matrix of barycentric coordinates

    integer :: nn

    self % n1 = n1
    self % n2 = n2
    self % p1 = p1
    self % p2 = p2

    nn = (n1-2)*n2

    ! Allocate temporary storage
!    allocate( self % Q_temp( n1*n2, n1*n2 ) )
    allocate( self % Qp_temp( 2*n2, 3 ) )
    allocate( self % L ( size(L,1), size(L,2) ) )
    allocate( self % Lt( size(L,2), size(L,1) ) )

    SLL_ASSERT( size(L,1) == 2*n2  )
    SLL_ASSERT( size(L,2) == 3     )

    self % L  = L
    self % Lt = transpose( L )

  end subroutine s_poisson_2d_fem_sps_stencil_projector__init

  ! Change basis: C1 projection of stiffness and mass matrices
  ! NOTE: 'self' has intent(inout) because temporary storage has to be assigned
  subroutine s_poisson_2d_fem_sps_stencil_projector__change_basis_matrix( self, Ql, Qp )
    class(sll_t_poisson_2d_fem_sps_stencil_projector)    , intent(inout) :: self
    type(sll_t_linear_operator_matrix_stencil_to_stencil), intent(inout) :: Ql
    type(sll_t_linear_operator_matrix_c1_block)          , intent(inout) :: Qp

    integer :: i1, i2, j1, j2, k1, k2, i, j, ip, jp, ll

    associate( n1 => self % n1, &
               n2 => self % n2, &
               p1 => self % p1, &
               p2 => self % p2 )

!      self % Q_temp = 0.0_wp
!      call Ql % to_array( self % Q_temp )
!
!      ! Fill block 1: 3 x 3
!      self % Qp_temp = matmul( self % Q_temp(1:2*n2,1:2*n2), self % L )
!      Qp % block1 % A(:,:) = matmul( self % Lt, self % Qp_temp )
!
!      ! Fill block 2: 3 x (n1-2)*n2
!      Qp % block2 % A(:,:) = matmul( self % Lt, self % Q_temp(1:2*n2,2*n2+1:n1*n2) )
!
!      ! Fill block 3: (n1-2)*n2 x 3
!      Qp % block3 % A(:,:) = matmul( self % Q_temp(2*n2+1:n1*n2,1:2*n2), self % L )
!
!      ! Fill block 4: (n1-2)*n2 x (n1-2)*n2
!      do i2 = 1, n2
!        do i1 = 1, n1-2
!          do k2 = -p2, p2
!            do k1 = -p1, p1
!              j1 = modulo( i1 - 1 + k1, n1-2 ) + 1
!              j2 = modulo( i2 - 1 + k2, n2   ) + 1
!              i  = (i1-1) * n2 + i2
!              j  = (j1-1) * n2 + j2
!              Qp % block4 % A(k1,k2,i1,i2) = self % Q_temp(2*n2+i,2*n2+j)
!            end do
!          end do
!        end do
!      end do

      ! Fill blocks
      self % Qp_temp (:,:) = 0.0_wp
      Qp % block1 % A(:,:) = 0.0_wp
      Qp % block2 % A(:,:) = 0.0_wp
      Qp % block3 % A(:,:) = 0.0_wp
      do ll = 1, 3
        do i2 = 1, n2
          do i1 = 1, n1
            do k2 = -p2, p2
              do k1 = -p1, p1
                j1 = modulo( i1 - 1 + k1, n1 ) + 1
                j2 = modulo( i2 - 1 + k2, n2 ) + 1
                i  = (i1-1) * n2 + i2
                j  = (j1-1) * n2 + j2
                ! block 1: 3 x 3
                if ( 1 <= i .and. i <= 2*n2 .and. 1 <= j .and. j <= 2*n2 ) then
                  self % Qp_temp(i,ll) = self % Qp_temp(i,ll) + Ql % A(k1,k2,i1,i2) * self % L(j,ll)
                ! block 2: 3 x (n1-2)*n2
                else if ( 1 <= i .and. i <= 2*n2 .and. 2*n2+1 <= j .and. j <= n1*n2 ) then
                  Qp % block2 % A(ll,j-2*n2) = Qp % block2 % A(ll,j-2*n2) + self % Lt(ll,i) * Ql % A(k1,k2,i1,i2)
                ! block 3: (n1-2)*n2 x 3
                else if ( 2*n2+1 <= i .and. i <= n1*n2 .and. 1 <= j .and. j <= 2*n2 ) then
                  Qp % block3 % A(i-2*n2,ll) = Qp % block3 % A(i-2*n2,ll) + Ql % A(k1,k2,i1,i2) * self % L(j,ll)
                end if
              end do
            end do
          end do
        end do
      end do
      ! complete calculation of block 1
      Qp % block1 % A(:,:) = matmul( self % Lt, self % Qp_temp )

      ! block 4: (n1-2)*n2 x (n1-2)*n2
      Qp % block4 % A(:,:,:,:) = 0.0_wp
      do ll = 1, n1-2
        do i2 = 1, n2
          do i1 = 1, n1
            do k2 = -p2, p2
              do k1 = -p1, p1
                j1 = modulo( ll - 1 + k1, n1-2 ) + 1
                j2 = modulo( i2 - 1 + k2, n2   ) + 1
                i  = (ll-1) * n2 + i2
                j  = (j1-1) * n2 + j2
                j1 = modulo( i1 - 1 + k1, n1 ) + 1
                ip = (i1-1) * n2 + i2
                jp = (j1-1) * n2 + j2
                if ( ip == 2*n2+i .and. jp == 2*n2+j ) &
                  Qp % block4 % A(k1,k2,ll,i2) = Ql % A(k1,k2,i1,i2)
              end do
            end do
          end do
        end do
      end do

    end associate

  end subroutine s_poisson_2d_fem_sps_stencil_projector__change_basis_matrix

  ! Change basis: C1 projection of vectors
  subroutine s_poisson_2d_fem_sps_stencil_projector__change_basis_vector( self, V, Vp )
    class(sll_t_poisson_2d_fem_sps_stencil_projector), intent(in   ) :: self
    real(wp)                                         , intent(in   ) :: V (:)
    type(sll_t_vector_space_c1_block)                , intent(inout) :: Vp

    integer :: j, j1, j2

    associate( n1 => self % n1, &
               n2 => self % n2, &
               p1 => self % p1, &
               p2 => self % p2 )

      ! Checks
      SLL_ASSERT( size( V ) == n1*n2 )
      SLL_ASSERT( size( Vp % vd % array    ) == 3 )
      SLL_ASSERT( size( Vp % vs % array, 1 ) == n1-2+2*p1 )
      SLL_ASSERT( size( Vp % vs % array, 2 ) == n2  +2*p2 )

      Vp % vd % array(1:3) = matmul( self % Lt, V(1:2*n2) )

      do j2 = 1, n2
        do j1 = 1, n1-2
          j = (j1-1) * n2 + j2
          Vp % vs % array(j1,j2) = V(2*n2+j)
        end do
      end do

    end associate

  end subroutine s_poisson_2d_fem_sps_stencil_projector__change_basis_vector

  ! Change basis: C1 projection of vectors
  subroutine s_poisson_2d_fem_sps_stencil_projector__change_basis_vector_inv( self, Vp, V )
    class(sll_t_poisson_2d_fem_sps_stencil_projector), intent(in   ) :: self
    type(sll_t_vector_space_c1_block)                , intent(inout) :: Vp
    real(wp)                                         , intent(inout) :: V (:)

    integer :: j, j1, j2

    associate( n1 => self % n1, &
               n2 => self % n2, &
               p1 => self % p1, &
               p2 => self % p2 )

      ! Checks
      SLL_ASSERT( size( Vp % vd % array    ) == 3 )
      SLL_ASSERT( size( Vp % vs % array, 1 ) == n1-2+2*p1 )
      SLL_ASSERT( size( Vp % vs % array, 2 ) == n2  +2*p2 )
      SLL_ASSERT( size( V ) == n1*n2 )

      V(1:2*n2) = matmul( self % L, Vp % vd % array(1:3) )

      do j2 = 1, n2
        do j1 = 1, n1-2
          j = (j1-1) * n2 + j2
          V(2*n2+j) = Vp % vs % array(j1,j2)
        end do
      end do

    end associate

  end subroutine s_poisson_2d_fem_sps_stencil_projector__change_basis_vector_inv

  ! Deallocate allocatables
  subroutine s_poisson_2d_fem_sps_stencil_projector__free( self )
    class(sll_t_poisson_2d_fem_sps_stencil_projector), intent(inout) :: self

!    deallocate( self % Q_temp  )
    deallocate( self % Qp_temp )
    deallocate( self % L       )
    deallocate( self % Lt      )

  end subroutine s_poisson_2d_fem_sps_stencil_projector__free

end module sll_m_poisson_2d_fem_sps_stencil_projector
