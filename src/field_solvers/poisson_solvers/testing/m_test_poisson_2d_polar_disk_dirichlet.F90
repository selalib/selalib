!> @authors Yaman Güçlü, IPP Garching
!> @authors Edoardo Zoni, IPP Garching
!>
!> @brief
!> Test-cases with Dirichlet BCs for Poisson's equation in 2D polar coordinates.
!>
!> @details
!> This module contains:
!> - Dirichlet base class with default parameters (all overwritable except for BCs)
!> - Test-case with parabolic radial profile:
!>   phi(r,th) = a * (1-(r/rmax)^2) + b * 4(r/rmax)(1-r/rmax)cos(k(th-th0))
!>
!> For a solver employing a 2nd-order scheme in the radial direction r and a
!> spectral Fourier method in the angular direction theta, the numerical error
!> expected on the first test-case (parabolic radial profile) is zero (machine
!> precision), if k (Fourier mode) is <= ntheta/2.

module m_test_poisson_2d_polar_disk_dirichlet
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_working_precision.h"

   use m_test_poisson_2d_polar_base, only: &
      c_test_poisson_2d_polar_base

   use sll_m_boundary_condition_descriptors, only: &
      sll_p_polar_origin, &
      sll_p_dirichlet

   implicit none

   public :: &
      t_test_poisson_2d_polar_disk_dirichlet_quadratic

   private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   !-----------------------------------------------------------------------------
   ! Dirichlet base class
   !-----------------------------------------------------------------------------
   type, extends(c_test_poisson_2d_polar_base), abstract :: c_test_disk_dirichlet

      ! Boundary conditions (not overwritable)
      sll_real64, private :: rmin = 0.0_f64
      sll_int32, private :: bc_rmin = sll_p_polar_origin
      sll_int32, private :: bc_rmax = sll_p_dirichlet

      ! Default parameters (may be overwritten by user)
      sll_real64 :: rmax = 2.0_f64
      sll_real64 :: a = 1.0_f64
      sll_real64 :: b = 0.7_f64
      sll_real64 :: th0 = 0.8_f64
      sll_int32  :: k = 3

   contains

      ! Get domain limits and boundary conditions
      procedure :: get_rlim => dirichlet_get_rlim
      procedure :: get_bcs => dirichlet_get_bcs

   end type c_test_disk_dirichlet

   !-----------------------------------------------------------------------------
   ! Test-case with expected zero numerical error (parabolic radial profile).
   ! Solution is linear combination of two terms with unitary amplitude:
   ! \phi(r,th) = a * (1-(r/rmax)^2) + b * 4(r/rmax)(1-r/rmax)cos(k(th-th0))
   !-----------------------------------------------------------------------------
   type, extends(c_test_disk_dirichlet) :: t_test_poisson_2d_polar_disk_dirichlet_quadratic

   contains
      ! 2D manufactured solution
      procedure :: phi_ex => dirichlet_zero_error_phi_ex
      procedure :: phi_ex_diff1_r => dirichlet_zero_error_phi_ex_d1r
      procedure :: phi_ex_diff2_r => dirichlet_zero_error_phi_ex_d2r
      procedure :: phi_ex_diff2_th => dirichlet_zero_error_phi_ex_d2th

   end type t_test_poisson_2d_polar_disk_dirichlet_quadratic

!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
contains
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

   !=============================================================================
   ! Dirichlet base class
   !=============================================================================

   pure function dirichlet_get_rlim(self) result(rlim)
      class(c_test_disk_dirichlet), intent(in) :: self
      sll_real64 :: rlim(2)

      rlim(1) = self%rmin
      rlim(2) = self%rmax

   end function dirichlet_get_rlim

   !-----------------------------------------------------------------------------
   pure function dirichlet_get_bcs(self) result(bcs)
      class(c_test_disk_dirichlet), intent(in) :: self
      sll_int32 :: bcs(2)

      bcs(1) = self%bc_rmin
      bcs(2) = self%bc_rmax

   end function dirichlet_get_bcs

   !=============================================================================
   ! Test-case with expected zero numerical error (parabolic radial profile).
   ! Solution is linear combination of two terms with unitary amplitude:
   ! \phi(r,th) = a * (1-(r/rmax)^2) + b * 4(r/rmax)(1-r/rmax)cos(k(th-th0))
   !=============================================================================

   pure function dirichlet_zero_error_phi_ex(self, r, th) result(val)
      class(t_test_poisson_2d_polar_disk_dirichlet_quadratic), intent(in) :: self
      sll_real64, intent(in) :: r
      sll_real64, intent(in) :: th
      sll_real64 :: val

      associate (rbar => r/self%rmax)
         val = self%a*(1.0_f64 - rbar**2) + &
               self%b*4.0_f64*rbar*(1.0_f64 - rbar)*cos(self%k*(th - self%th0))
      end associate

   end function dirichlet_zero_error_phi_ex

   !-----------------------------------------------------------------------------
   pure function dirichlet_zero_error_phi_ex_d1r(self, r, th) result(val)
      class(t_test_poisson_2d_polar_disk_dirichlet_quadratic), intent(in) :: self
      sll_real64, intent(in) :: r
      sll_real64, intent(in) :: th
      sll_real64 :: val

      associate (rbar => r/self%rmax)
         val = self%a*(-2.0_f64/self%rmax)*rbar + &
               self%b*(4.0_f64/self%rmax)*(1.0_f64 - 2.0_f64*rbar)*cos(self%k*(th - self%th0))
      end associate

   end function dirichlet_zero_error_phi_ex_d1r

   !-----------------------------------------------------------------------------
   pure function dirichlet_zero_error_phi_ex_d2r(self, r, th) result(val)
      class(t_test_poisson_2d_polar_disk_dirichlet_quadratic), intent(in) :: self
      sll_real64, intent(in) :: r
      sll_real64, intent(in) :: th
      sll_real64 :: val

      val = self%a*(-2.0_f64/self%rmax**2) + &
            self%b*(-8.0_f64/self%rmax**2)*cos(self%k*(th - self%th0))

   end function dirichlet_zero_error_phi_ex_d2r

   !-----------------------------------------------------------------------------
   pure function dirichlet_zero_error_phi_ex_d2th(self, r, th) result(val)
      class(t_test_poisson_2d_polar_disk_dirichlet_quadratic), intent(in) :: self
      sll_real64, intent(in) :: r
      sll_real64, intent(in) :: th
      sll_real64 :: val

      associate (rbar => r/self%rmax)
         val = self%b*(-self%k**2)*4.0_f64*rbar*(1.0_f64 - rbar)*cos(self%k*(th - self%th0))
      end associate

   end function dirichlet_zero_error_phi_ex_d2th

end module m_test_poisson_2d_polar_disk_dirichlet
