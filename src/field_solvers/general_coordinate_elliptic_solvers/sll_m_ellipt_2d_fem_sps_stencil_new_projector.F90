module sll_m_ellipt_2d_fem_sps_stencil_new_projector
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_assert.h"

   use sll_m_working_precision, only: f64

   use sll_m_vector_space_c1_block, only: sll_t_vector_space_c1_block

   use sll_m_linear_operator_matrix_stencil_to_stencil, only: sll_t_linear_operator_matrix_stencil_to_stencil

   use sll_m_linear_operator_matrix_c1_block_new, only: sll_t_linear_operator_matrix_c1_block_new

   implicit none

   public :: sll_t_ellipt_2d_fem_sps_stencil_new_projector

   private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   ! Working precision
   integer, parameter :: wp = f64

   type :: sll_t_ellipt_2d_fem_sps_stencil_new_projector

      integer :: n1
      integer :: n2
      integer :: p1
      integer :: p2

      ! Temporary storage needed for projection
      real(wp), allocatable :: Qp_temp(:, :, :)
      real(wp), allocatable :: L(:, :, :)

   contains

      procedure :: init => s_ellipt_2d_fem_sps_stencil_new_projector__init
      procedure :: change_basis_matrix => s_ellipt_2d_fem_sps_stencil_new_projector__change_basis_matrix
      procedure :: change_basis_vector => s_ellipt_2d_fem_sps_stencil_new_projector__change_basis_vector
      procedure :: change_basis_vector_inv => s_ellipt_2d_fem_sps_stencil_new_projector__change_basis_vecinv
      procedure :: free => s_ellipt_2d_fem_sps_stencil_new_projector__free

   end type sll_t_ellipt_2d_fem_sps_stencil_new_projector

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
contains
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   ! Initializer
   subroutine s_ellipt_2d_fem_sps_stencil_new_projector__init(self, n1, n2, p1, p2, L)
      class(sll_t_ellipt_2d_fem_sps_stencil_new_projector), intent(inout) :: self
      integer, intent(in) :: n1
      integer, intent(in) :: n2
      integer, intent(in) :: p1
      integer, intent(in) :: p2
      real(wp), intent(in) :: L(:, :, :) ! matrix of barycentric coordinates

      self%n1 = n1
      self%n2 = n2
      self%p1 = p1
      self%p2 = p2

      SLL_ASSERT(size(L, 1) == 2)
      SLL_ASSERT(size(L, 2) == n2)
      SLL_ASSERT(size(L, 3) == 3)

      ! Allocate temporary storage
      allocate (self%L(size(L, 1), size(L, 2), size(L, 3)), source=L)
      allocate (self%Qp_temp(size(L, 1), size(L, 2), size(L, 3)))

   end subroutine s_ellipt_2d_fem_sps_stencil_new_projector__init

   ! Change basis: C1 projection of stiffness and mass matrices
   ! NOTE: 'self' has intent(inout) because temporary storage has to be assigned
   subroutine s_ellipt_2d_fem_sps_stencil_new_projector__change_basis_matrix(self, Ql, Qp)
      class(sll_t_ellipt_2d_fem_sps_stencil_new_projector), intent(inout) :: self
      type(sll_t_linear_operator_matrix_stencil_to_stencil), intent(inout) :: Ql
      type(sll_t_linear_operator_matrix_c1_block_new), intent(inout) :: Qp

      integer :: i1, i2, j1, j2, k1, k2, ip, jp, ll

      associate (n1 => self%n1, &
                 n2 => self%n2, &
                 p1 => self%p1, &
                 p2 => self%p2)

         ! Fill blocks
         self%Qp_temp(:, :, :) = 0.0_wp
         Qp%block1%A(:, :) = 0.0_wp
         Qp%block2%A(:, :, :) = 0.0_wp
         Qp%block3%A(:, :, :) = 0.0_wp
         do ll = 1, 3
            do i2 = 1, n2
               do i1 = 1, n1
                  do k2 = -p2, p2
                     do k1 = -p1, p1

                        j1 = i1 + k1
                        j2 = i2 + k2

                        ! Hand-made modulo operation
                        if (j1 < 1) then
                           j1 = j1 + n1
                        else if (j1 > n1) then
                           j1 = j1 - n1
                        end if

                        ! Hand-made modulo operation
                        if (j2 < 1) then
                           j2 = j2 + n2
                        else if (j2 > n2) then
                           j2 = j2 - n2
                        end if

                        ! block 1: 3 x 3
                        if (i1 <= 2 .and. j1 <= 2) then
                           self%Qp_temp(i1, i2, ll) = self%Qp_temp(i1, i2, ll) + Ql%A(k1, k2, i1, i2)*self%L(j1, j2, ll)
                           ! block 2: 3 x (n1-2)*n2
                        else if (i1 <= 2 .and. j1 <= 2 + p1) then
                           Qp%block2%A(ll, j1 - 2, j2) = Qp%block2%A(ll, j1 - 2, j2) + self%L(i1, i2, ll)*Ql%A(k1, k2, i1, i2)
                           ! block 3: (n1-2)*n2 x 3
                        else if (i1 <= 2 + p1 .and. j1 <= 2) then
                           Qp%block3%A(i1 - 2, i2, ll) = Qp%block3%A(i1 - 2, i2, ll) + Ql%A(k1, k2, i1, i2)*self%L(j1, j2, ll)
                        end if

                     end do
                  end do
               end do
            end do
         end do
         ! complete calculation of block 1
         do ip = 1, 3
            do jp = 1, 3
               do i2 = 1, n2
                  do i1 = 1, 2
                     Qp%block1%A(ip, jp) = Qp%block1%A(ip, jp) + self%L(i1, i2, ip)*self%Qp_temp(i1, i2, jp)
                  end do
               end do
            end do
         end do

         ! block 4: (n1-2)*n2 x (n1-2)*n2
         Qp%block4%A(:, :, :, :) = 0.0_wp
         do i2 = 1, n2
            do i1 = 1, n1 - 2
               do k2 = -p2, p2
                  do k1 = -p1, p1
                     j1 = modulo(i1 - 1 + k1, n1 - 2)
                     ll = modulo(i1 + 1 + k1, n1)
                     if (ll == 2 + j1) Qp%block4%A(k1, k2, i1, i2) = Ql%A(k1, k2, i1 + 2, i2)
                  end do
               end do
            end do
         end do

      end associate

   end subroutine s_ellipt_2d_fem_sps_stencil_new_projector__change_basis_matrix

   ! Change basis: C1 projection of vectors
   subroutine s_ellipt_2d_fem_sps_stencil_new_projector__change_basis_vector(self, V, Vp)
      class(sll_t_ellipt_2d_fem_sps_stencil_new_projector), intent(in) :: self
      real(wp), intent(in) :: V(:, :)
      type(sll_t_vector_space_c1_block), intent(inout) :: Vp

      integer :: i1, i2, ll

      associate (n1 => self%n1, &
                 n2 => self%n2, &
                 p1 => self%p1, &
                 p2 => self%p2)

         ! Checks
         SLL_ASSERT(size(V, 1) == n1)
         SLL_ASSERT(size(V, 2) == n2)
         SLL_ASSERT(size(Vp%vd%array) == 3)
         SLL_ASSERT(size(Vp%vs%array, 1) == n1 - 2 + 2*p1)
         SLL_ASSERT(size(Vp%vs%array, 2) == n2 + 2*p2)

         Vp%vd%array(:) = 0.0_wp
         do ll = 1, 3
            do i2 = 1, n2
               do i1 = 1, 2
                  Vp%vd%array(ll) = Vp%vd%array(ll) + self%L(i1, i2, ll)*V(i1, i2)
               end do
            end do
         end do

         do i2 = 1, n2
            do i1 = 3, n1
               Vp%vs%array(i1 - 2, i2) = V(i1, i2)
            end do
         end do

      end associate

   end subroutine s_ellipt_2d_fem_sps_stencil_new_projector__change_basis_vector

   ! Change basis: C1 projection of vectors
   subroutine s_ellipt_2d_fem_sps_stencil_new_projector__change_basis_vecinv(self, Vp, V)
      class(sll_t_ellipt_2d_fem_sps_stencil_new_projector), intent(in) :: self
      type(sll_t_vector_space_c1_block), intent(inout) :: Vp
      real(wp), intent(inout) :: V(:)

      integer :: i, i1, i2, j, j1, j2, ll

      associate (n1 => self%n1, &
                 n2 => self%n2, &
                 p1 => self%p1, &
                 p2 => self%p2)

         ! Checks
         SLL_ASSERT(size(Vp%vd%array) == 3)
         SLL_ASSERT(size(Vp%vs%array, 1) == n1 - 2 + 2*p1)
         SLL_ASSERT(size(Vp%vs%array, 2) == n2 + 2*p2)
         SLL_ASSERT(size(V) == n1*n2)

         V(:) = 0.0_wp
         do ll = 1, 3
            do i2 = 1, n2
               do i1 = 1, 2
                  i = (i1 - 1)*n2 + i2
                  V(i) = V(i) + self%L(i1, i2, ll)*Vp%vd%array(ll)
               end do
            end do
         end do

         do j2 = 1, n2
            do j1 = 1, n1 - 2
               j = (j1 - 1)*n2 + j2
               V(2*n2 + j) = Vp%vs%array(j1, j2)
            end do
         end do

      end associate

   end subroutine s_ellipt_2d_fem_sps_stencil_new_projector__change_basis_vecinv

   ! Deallocate allocatables
   subroutine s_ellipt_2d_fem_sps_stencil_new_projector__free(self)
      class(sll_t_ellipt_2d_fem_sps_stencil_new_projector), intent(inout) :: self

      deallocate (self%Qp_temp)
      deallocate (self%L)

   end subroutine s_ellipt_2d_fem_sps_stencil_new_projector__free

end module sll_m_ellipt_2d_fem_sps_stencil_new_projector
