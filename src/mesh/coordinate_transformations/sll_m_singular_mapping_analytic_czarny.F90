! This module implements the analytical mapping given by
! equation (3) of https://doi.org/10.1016/j.jcp.2019.108889:
!
! x(s,theta)=(1-sqrt(1+epsilon*(epsilon+2*s*cos(theta))))/epsilon
! y(s,theta)=y0+(e*xi*s*sin(theta))/(1+epsilon*x(s,theta))

module sll_m_singular_mapping_analytic_czarny
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_assert.h"

   use sll_m_working_precision, only: f64

   use sll_m_constants, only: sll_p_twopi

   use sll_m_singular_mapping_analytic, only: sll_c_singular_mapping_analytic

   implicit none

   public :: sll_t_singular_mapping_analytic_czarny

   private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   !> Working precision
   integer, parameter :: wp = f64

   !> Concrete type, analytical singular mapping
   type, extends(sll_c_singular_mapping_analytic) :: sll_t_singular_mapping_analytic_czarny

      ! Default values (no specific meaning), can be overwritten from 'init' method
      real(wp) :: x0(2) = [0.0_wp, 0.0_wp]
      real(wp) :: e = 1.0_wp
      real(wp) :: eps = 1.0_wp

   contains

      procedure :: init => s_singular_mapping_analytic_czarny__init
      procedure :: eval => f_singular_mapping_analytic_czarny__eval
      procedure :: jmat => f_singular_mapping_analytic_czarny__jmat ! Jacobian matrix
      procedure :: jmat_comp => f_singular_mapping_analytic_czarny__jmat_comp

   end type sll_t_singular_mapping_analytic_czarny

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
contains
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   subroutine s_singular_mapping_analytic_czarny__init(self, x0, e, eps)
      class(sll_t_singular_mapping_analytic_czarny), intent(inout) :: self
      real(wp), optional, intent(in) :: x0(2)
      real(wp), optional, intent(in) :: e
      real(wp), optional, intent(in) :: eps

      ! Overwrite parameters
      if (present(x0)) self%x0 = x0
      if (present(e)) self%e = e
      if (present(eps)) self%eps = eps

   end subroutine s_singular_mapping_analytic_czarny__init

   !-----------------------------------------------------------------------------
   SLL_PURE function f_singular_mapping_analytic_czarny__eval(self, eta) result(x)
      class(sll_t_singular_mapping_analytic_czarny), intent(in) :: self
      real(wp), intent(in) :: eta(2)
      real(wp) :: x(2)

      associate (s => eta(1), t => eta(2), x0 => self%x0, e => self%e, eps => self%eps)

         ! TODO: check if x0(1) can be added to x(1) as well
         x(1) = (1.0_wp - sqrt(1.0_wp + eps*(eps + 2.0_wp*s*cos(t))))/eps
         x(2) = x0(2) + e*s*sin(t)/((1.0_wp + eps*x(1))*sqrt(1.0_wp - eps**2*0.25_wp))

      end associate

   end function f_singular_mapping_analytic_czarny__eval

   !-----------------------------------------------------------------------------
   SLL_PURE function f_singular_mapping_analytic_czarny__jmat(self, eta) result(jmat)
      class(sll_t_singular_mapping_analytic_czarny), intent(in) :: self
      real(wp), intent(in) :: eta(2)
      real(wp) :: jmat(2, 2)

      real(wp) :: tmp1, tmp2

      associate (s => eta(1), t => eta(2), e => self%e, eps => self%eps)

         tmp1 = sqrt(1.0_wp + eps*(eps + 2.0_wp*s*cos(t)))
         tmp2 = sqrt(1.0_wp - eps**2*0.25_wp)

         ! J_11 = d(x1)/d(eta1)
         ! J_12 = d(x1)/d(eta2)
         ! J_21 = d(x2)/d(eta1)
         ! J_22 = d(x2)/d(eta2)
         jmat(1, 1) = -cos(t)/tmp1
         jmat(1, 2) = s*sin(t)/tmp1
         jmat(2, 1) = e*sin(t)/((2.0_wp - tmp1)*tmp2) + eps*e*s*sin(t)*cos(t)/(tmp1*tmp2*(2.0_wp - tmp1)**2)
         jmat(2, 2) = e*s*cos(t)/((2.0_wp - tmp1)*tmp2) - eps*e*s**2*sin(t)**2/(tmp1*tmp2*(2.0_wp - tmp1)**2)

      end associate

   end function f_singular_mapping_analytic_czarny__jmat

   !-----------------------------------------------------------------------------
   SLL_PURE function f_singular_mapping_analytic_czarny__jmat_comp(self, eta) result(jmat_comp)
      class(sll_t_singular_mapping_analytic_czarny), intent(in) :: self
      real(wp), intent(in) :: eta(2)
      real(wp) :: jmat_comp(2, 2)

      real(wp) :: tmp1, tmp2

      associate (s => eta(1), t => eta(2), e => self%e, eps => self%eps)

         tmp1 = sqrt(1.0_wp + eps*(eps + 2.0_wp*s*cos(t)))
         tmp2 = sqrt(1.0_wp - eps**2*0.25_wp)

         jmat_comp(1, 1) = e/(tmp2*(2.0_wp - tmp1))
         jmat_comp(1, 2) = 0.0_wp
         jmat_comp(2, 1) = -e*eps*s*sin(t)/(tmp1*tmp2*(2.0_wp - tmp1)**2)
         jmat_comp(2, 2) = -1.0_wp/tmp1

         jmat_comp = jmat_comp/(-e/(tmp1*tmp2*(2.0_wp - tmp1)))

      end associate

   end function f_singular_mapping_analytic_czarny__jmat_comp

end module sll_m_singular_mapping_analytic_czarny
