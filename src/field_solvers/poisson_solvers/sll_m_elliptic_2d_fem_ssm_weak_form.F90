module sll_m_elliptic_2d_fem_ssm_weak_form
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  use sll_m_working_precision, only: f64

  implicit none

  public :: sll_c_elliptic_2d_fem_ssm_weak_form

  private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  ! Working precision
  integer, parameter :: wp = f64

  type, abstract :: sll_c_elliptic_2d_fem_ssm_weak_form

  contains

    ! Deferred procedures
    procedure(i_sub_element), deferred :: element

  end type sll_c_elliptic_2d_fem_ssm_weak_form

  ! Interfaces for deferred procedures
  abstract interface

    subroutine i_sub_element( &
      self        , &
      test_values_and_derivs_eta1 , &
      test_values_and_derivs_eta2 , &
      trial_values_and_derivs_eta1, &
      trial_values_and_derivs_eta2, &
      int_volume                  , &
      inv_metric                  , &
      Aij, &
      Mij )
      import sll_c_elliptic_2d_fem_ssm_weak_form, wp
      class(sll_c_elliptic_2d_fem_ssm_weak_form), intent(in   ) :: self
      real(wp)                                  , intent(in   ) :: test_values_and_derivs_eta1 (:,:)
      real(wp)                                  , intent(in   ) :: test_values_and_derivs_eta2 (:,:)
      real(wp)                                  , intent(in   ) :: trial_values_and_derivs_eta1(:,:)
      real(wp)                                  , intent(in   ) :: trial_values_and_derivs_eta2(:,:)
      real(wp)                                  , intent(in   ) :: int_volume(:,:)
      real(wp)                                  , intent(in   ) :: inv_metric(:,:,:,:)
      real(wp)                                  , intent(  out) :: Aij
      real(wp)                                  , intent(  out) :: Mij
    end subroutine i_sub_element

  end interface

end module sll_m_elliptic_2d_fem_ssm_weak_form
