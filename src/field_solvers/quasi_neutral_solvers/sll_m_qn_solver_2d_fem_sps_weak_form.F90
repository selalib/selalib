module sll_m_qn_solver_2d_fem_sps_weak_form
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   use sll_m_working_precision, only: f64

   use sll_m_ellipt_2d_fem_sps_weak_form, only: sll_c_ellipt_2d_fem_sps_weak_form

   implicit none

   public :: sll_t_qn_solver_2d_fem_sps_weak_form

   private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   ! Working precision
   integer, parameter :: wp = f64

   type, extends(sll_c_ellipt_2d_fem_sps_weak_form) :: sll_t_qn_solver_2d_fem_sps_weak_form

   contains

      procedure :: element_mat => s_qn_solver_2d_fem_sps_weak_form__element_mat
      procedure :: element_rhs => s_qn_solver_2d_fem_sps_weak_form__element_rhs

   end type sll_t_qn_solver_2d_fem_sps_weak_form

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
contains
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   ! Compute elementary contribution to stiffness and mass matrices
   subroutine s_qn_solver_2d_fem_sps_weak_form__element_mat( &
      self, &
      test_values_and_derivs_eta1, &
      test_values_and_derivs_eta2, &
      trial_values_and_derivs_eta1, &
      trial_values_and_derivs_eta2, &
      int_volume, &
      inv_metric, &
      coeffs1, &
      coeffs2, &
      Aij, &
      Mij)
      class(sll_t_qn_solver_2d_fem_sps_weak_form), intent(in) :: self
      real(wp), intent(in) :: test_values_and_derivs_eta1(:, :)
      real(wp), intent(in) :: test_values_and_derivs_eta2(:, :)
      real(wp), intent(in) :: trial_values_and_derivs_eta1(:, :)
      real(wp), intent(in) :: trial_values_and_derivs_eta2(:, :)
      real(wp), intent(in) :: int_volume(:, :)
      real(wp), intent(in) :: inv_metric(:, :, :, :)
      real(wp), intent(in) :: coeffs1(:, :)
      real(wp), intent(in) :: coeffs2(:, :)
      real(wp), intent(out) :: Aij
      real(wp), intent(out) :: Mij

      integer :: q1, q2, Nq1, Nq2

      real(wp) :: v(2), w(2), t(2)

      ! Extract number of quadrature points
      Nq1 = size(test_values_and_derivs_eta1, 1)
      Nq2 = size(test_values_and_derivs_eta2, 1)

      Aij = 0.0_wp
      Mij = 0.0_wp

      ! Quadrature
      do q2 = 1, Nq2
         do q1 = 1, Nq1

            ! 2D gradient in logical test space
            v = [test_values_and_derivs_eta1(q1, 2)*test_values_and_derivs_eta2(q2, 1), &
                 test_values_and_derivs_eta1(q1, 1)*test_values_and_derivs_eta2(q2, 2)]

            ! 2D gradient in logical trial space
            w = [trial_values_and_derivs_eta1(q1, 2)*trial_values_and_derivs_eta2(q2, 1), &
                 trial_values_and_derivs_eta1(q1, 1)*trial_values_and_derivs_eta2(q2, 2)]

            t = matmul(inv_metric(q1, q2, :, :), w)

            ! Elementary contribution to stiffness matrix
            Aij = Aij + (coeffs1(q1, q2)*dot_product(v, t) + coeffs2(q1, q2)* &
                         test_values_and_derivs_eta1(q1, 1)*test_values_and_derivs_eta2(q2, 1)* &
                         trial_values_and_derivs_eta1(q1, 1)*trial_values_and_derivs_eta2(q2, 1))* &
                  int_volume(q1, q2)

            ! Elementary contribution to mass matrix
            Mij = Mij + &
                  test_values_and_derivs_eta1(q1, 1)*test_values_and_derivs_eta2(q2, 1)* &
                  trial_values_and_derivs_eta1(q1, 1)*trial_values_and_derivs_eta2(q2, 1)* &
                  int_volume(q1, q2)

         end do
      end do

   end subroutine s_qn_solver_2d_fem_sps_weak_form__element_mat

   ! Compute elementary contribution to right hand side
   subroutine s_qn_solver_2d_fem_sps_weak_form__element_rhs( &
      self, &
      test_values_and_derivs_eta1, &
      test_values_and_derivs_eta2, &
      data_2d_rhs, &
      int_volume, &
      bi)
      class(sll_t_qn_solver_2d_fem_sps_weak_form), intent(in) :: self
      real(wp), intent(in) :: test_values_and_derivs_eta1(:, :)
      real(wp), intent(in) :: test_values_and_derivs_eta2(:, :)
      real(wp), intent(in) :: data_2d_rhs(:, :)
      real(wp), intent(in) :: int_volume(:, :)
      real(wp), intent(out) :: bi

      integer :: q1, q2, Nq1, Nq2

      ! Extract number of quadrature points
      Nq1 = size(test_values_and_derivs_eta1, 1)
      Nq2 = size(test_values_and_derivs_eta2, 1)

      bi = 0.0_wp

      ! Quadrature
      do q2 = 1, Nq2
         do q1 = 1, Nq1

            ! Elementary contribution to mass matrix
            bi = bi + &
                 test_values_and_derivs_eta1(q1, 1)*test_values_and_derivs_eta2(q2, 1)* &
                 data_2d_rhs(q1, q2)*int_volume(q1, q2)

         end do
      end do

   end subroutine s_qn_solver_2d_fem_sps_weak_form__element_rhs

end module sll_m_qn_solver_2d_fem_sps_weak_form
