!> @authors Yaman Güçlü, IPP Garching
!> @authors Edoardo Zoni, IPP Garching
!>
!> @brief
!> Test-cases with Dirichlet BCs for quasi-neutrality equation in 2D polar coordinates.
!>
!> @details
!> This module contains:
!> - Dirichlet base class with default parameters (all overwritable except for BCs)
!> - Test-case with parabolic radial profile:
!>   phi(r,th) = a * (1-(r/rmax)^2) + b * 4(r/rmax)(1-r/rmax)cos(k(th-th0))
!>
!> For a solver employing a 2nd-order scheme in the radial direction r and a
!> spectral Fourier method in the angular direction theta, the numerical error
!> expected on this test-case is zero (machine precision), if k (Fourier mode)
!> is <= ntheta/2.

module m_test_qn_solver_2d_polar_disk_dirichlet
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_working_precision.h"

   use m_test_qn_solver_2d_polar_base, only: &
      c_test_qn_solver_2d_polar_base

   use sll_m_boundary_condition_descriptors, only: &
      sll_p_polar_origin, &
      sll_p_dirichlet

   implicit none

   public :: &
      t_test_qn_solver_2d_polar_disk_dirichlet_quadratic

   private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   !-----------------------------------------------------------------------------
   ! Dirichlet base class
   !-----------------------------------------------------------------------------
   type, extends(c_test_qn_solver_2d_polar_base), abstract :: c_test_disk_dirichlet

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
      logical    :: adiabatic_electrons = .true.
      logical    :: use_zonal_flow = .true.
      sll_real64 :: epsilon_0 = 1.0_f64

   contains

      ! Get domain limits and boundary conditions
      procedure :: get_rlim => dirichlet_get_rlim
      procedure :: get_bcs => dirichlet_get_bcs
      procedure :: get_parameters => dirichlet_get_parameters

   end type c_test_disk_dirichlet

   !-----------------------------------------------------------------------------
   ! Test-case with expected zero numerical error (parabolic radial profile)
   ! Solution is linear combination of two terms with unitary amplitude:
   ! \phi(r,th) = a * (1-(r/rmax)^2) + b * 4(r/rmax)(1-r/rmax)cos(k(th-th0))
   !-----------------------------------------------------------------------------
   type, extends(c_test_disk_dirichlet) :: t_test_qn_solver_2d_polar_disk_dirichlet_quadratic

   contains
      ! 1D input profiles
      procedure :: rho_m0 => f_test__rho_m0
      procedure :: rho_m0_diff1_r => f_test__rho_m0_diff1_r
      procedure :: b_magn => f_test__b_magn
      procedure :: b_magn_diff1_r => f_test__b_magn_diff1_r
      procedure :: lambda => f_test__lambda
      ! 2D manufactured solution
      procedure :: phi_ex => f_test__phi_ex
      procedure :: phi_ex_avg_th => f_test__phi_ex_avg_th
      procedure :: phi_ex_diff1_r => f_test__phi_ex_diff1_r
      procedure :: phi_ex_diff2_r => f_test__phi_ex_diff2_r
      procedure :: phi_ex_diff2_th => f_test__phi_ex_diff2_th

   end type t_test_qn_solver_2d_polar_disk_dirichlet_quadratic

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

   !-----------------------------------------------------------------------------
   pure subroutine dirichlet_get_parameters(self, adiabatic_electrons, &
                                            use_zonal_flow, epsilon_0)
      class(c_test_disk_dirichlet), intent(in) :: self
      logical, intent(out) :: adiabatic_electrons
      logical, intent(out) :: use_zonal_flow
      sll_real64, intent(out) :: epsilon_0

      adiabatic_electrons = self%adiabatic_electrons
      use_zonal_flow = self%use_zonal_flow
      epsilon_0 = self%epsilon_0

   end subroutine dirichlet_get_parameters

   !=============================================================================
   ! Test-case with expected zero numerical error (parabolic radial profile)
   ! Solution is linear combination of two terms with unitary amplitude:
   ! \phi(r,th) = a * (1-(r/rmax)^2) + b * 4(r/rmax)(1-r/rmax)cos(k(th-th0))
   !=============================================================================

   !-----------------------------------------------------------------------------
   ! 1D PROFILES
   !-----------------------------------------------------------------------------

   pure function f_test__rho_m0(self, r) result(val)
      class(t_test_qn_solver_2d_polar_disk_dirichlet_quadratic), intent(in) :: self
      sll_real64, intent(in) :: r
      sll_real64 :: val
      val = 1.0_f64
   end function f_test__rho_m0

   !-----------------------------------------------------------------------------
   pure function f_test__rho_m0_diff1_r(self, r) result(val)
      class(t_test_qn_solver_2d_polar_disk_dirichlet_quadratic), intent(in) :: self
      sll_real64, intent(in) :: r
      sll_real64 :: val
      val = 0.0_f64
   end function f_test__rho_m0_diff1_r

   !-----------------------------------------------------------------------------
   pure function f_test__b_magn(self, r) result(val)
      class(t_test_qn_solver_2d_polar_disk_dirichlet_quadratic), intent(in) :: self
      sll_real64, intent(in) :: r
      sll_real64 :: val
      val = 1.0_f64
   end function f_test__b_magn

   !-----------------------------------------------------------------------------
   pure function f_test__b_magn_diff1_r(self, r) result(val)
      class(t_test_qn_solver_2d_polar_disk_dirichlet_quadratic), intent(in) :: self
      sll_real64, intent(in) :: r
      sll_real64 :: val
      val = 0.0_f64
   end function f_test__b_magn_diff1_r

   !-----------------------------------------------------------------------------
   pure function f_test__lambda(self, r) result(val)
      class(t_test_qn_solver_2d_polar_disk_dirichlet_quadratic), intent(in) :: self
      sll_real64, intent(in) :: r
      sll_real64 :: val

      sll_real64 :: temp_e
      associate (minv => 0.01_f64, &
                 maxv => 1.0_f64, &
                 width => 0.1_f64*self%rmax, &
                 dens_e => self%rho_m0(r))
         temp_e = minv + (maxv - minv)*exp(-(r/2.0_f64*width)**2)  ! in [0.01, 1]
         val = 0.01_f64*sqrt(temp_e/dens_e)
      end associate
   end function f_test__lambda

   !-----------------------------------------------------------------------------
   ! 2D MANUFACTURED SOLUTION
   !-----------------------------------------------------------------------------

   pure function f_test__phi_ex(self, r, th) result(val)
      class(t_test_qn_solver_2d_polar_disk_dirichlet_quadratic), intent(in) :: self
      sll_real64, intent(in) :: r
      sll_real64, intent(in) :: th
      sll_real64 :: val

      associate (rbar => r/self%rmax)
         val = self%a*(1.0_f64 - rbar**2) + &
               self%b*4.0_f64*rbar*(1.0_f64 - rbar)*cos(self%k*(th - self%th0))
      end associate

   end function f_test__phi_ex

   !-----------------------------------------------------------------------------
   pure function f_test__phi_ex_avg_th(self, r, th) result(val)
      class(t_test_qn_solver_2d_polar_disk_dirichlet_quadratic), intent(in) :: self
      sll_real64, intent(in) :: r
      sll_real64, intent(in) :: th
      sll_real64 :: val

      associate (rbar => r/self%rmax)
         if (self%k == 0) then
            val = self%a*(1.0_f64 - rbar**2) + self%b*4.0_f64*rbar*(1.0_f64 - rbar)
         else
            val = self%a*(1.0_f64 - rbar**2)
         end if
      end associate

   end function f_test__phi_ex_avg_th

   !-----------------------------------------------------------------------------
   pure function f_test__phi_ex_diff1_r(self, r, th) result(val)
      class(t_test_qn_solver_2d_polar_disk_dirichlet_quadratic), intent(in) :: self
      sll_real64, intent(in) :: r
      sll_real64, intent(in) :: th
      sll_real64 :: val

      associate (rbar => r/self%rmax)
         val = self%a*(-2.0_f64/self%rmax)*rbar + &
               self%b*(4.0_f64/self%rmax)*(1.0_f64 - 2.0_f64*rbar)*cos(self%k*(th - self%th0))
      end associate

   end function f_test__phi_ex_diff1_r

   !-----------------------------------------------------------------------------
   pure function f_test__phi_ex_diff2_r(self, r, th) result(val)
      class(t_test_qn_solver_2d_polar_disk_dirichlet_quadratic), intent(in) :: self
      sll_real64, intent(in) :: r
      sll_real64, intent(in) :: th
      sll_real64 :: val

      val = self%a*(-2.0_f64/self%rmax**2) + &
            self%b*(-8.0_f64/self%rmax**2)*cos(self%k*(th - self%th0))

   end function f_test__phi_ex_diff2_r

   !-----------------------------------------------------------------------------
   pure function f_test__phi_ex_diff2_th(self, r, th) result(val)
      class(t_test_qn_solver_2d_polar_disk_dirichlet_quadratic), intent(in) :: self
      sll_real64, intent(in) :: r
      sll_real64, intent(in) :: th
      sll_real64 :: val

      associate (rbar => r/self%rmax)
         val = self%b*(-self%k**2)*4.0_f64*rbar*(1.0_f64 - rbar)*cos(self%k*(th - self%th0))
      end associate

   end function f_test__phi_ex_diff2_th

end module m_test_qn_solver_2d_polar_disk_dirichlet
