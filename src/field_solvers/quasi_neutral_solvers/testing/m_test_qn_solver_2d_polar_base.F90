!> @authors Yaman Güçlü, IPP Garching
!> @authors Edoardo Zoni, IPP Garching
!>
!> @brief
!> Method of manufactured solutions for quasi-neutrality equation in 2D polar coordinates.
!>
!> @details
!> This module defines an abstract interface which requires subclasses to implement
!> an analytical function phi(r,theta), its derivatives and its average over theta.
!> Additional analytical inputs are: rho_m0(r) and b_magn(r) with the respective
!> derivatives, and lambda(r).
!> The analytical rho(r,theta) is calculated here according to the quasi-neutrality
!> equation, and it is not overridable.
!> In the numerical tests rho(r,theta) will be given to the quasi-neutral solver,
!> and the resulting numerical phi will be compared to the exact solution.

module m_test_qn_solver_2d_polar_base
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_working_precision.h"

   implicit none

   public :: &
      c_test_qn_solver_2d_polar_base

   private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   type, abstract :: c_test_qn_solver_2d_polar_base

   contains
      ! Get domain limits, boundary conditions and parameters
      procedure(i_func_get_rlim), deferred :: get_rlim
      procedure(i_func_get_bcs), deferred :: get_bcs
      procedure(i_subr_get_parameters), deferred :: get_parameters
      ! 1D input profiles
      procedure(i_func_1d_real), deferred :: rho_m0
      procedure(i_func_1d_real), deferred :: rho_m0_diff1_r
      procedure(i_func_1d_real), deferred :: b_magn
      procedure(i_func_1d_real), deferred :: b_magn_diff1_r
      procedure(i_func_1d_real), deferred :: lambda
      ! 2D manufactured solution
      procedure(i_func_2d_real), deferred :: phi_ex
      procedure(i_func_2d_real), deferred :: phi_ex_avg_th
      procedure(i_func_2d_real), deferred :: phi_ex_diff1_r
      procedure(i_func_2d_real), deferred :: phi_ex_diff2_r
      procedure(i_func_2d_real), deferred :: phi_ex_diff2_th
      ! 2D right-hand side to solver
      procedure, non_overridable :: rho => f_test_case__rho

   end type c_test_qn_solver_2d_polar_base

   !-----------------------------------------------------------------------------
   abstract interface

      ! Get domain limits
      pure function i_func_get_rlim(self) result(rlim)
         use sll_m_working_precision
         import c_test_qn_solver_2d_polar_base
         class(c_test_qn_solver_2d_polar_base), intent(in) :: self
         sll_real64 :: rlim(2)
      end function i_func_get_rlim

      ! Get boundary conditions
      pure function i_func_get_bcs(self) result(bcs)
         use sll_m_working_precision
         import c_test_qn_solver_2d_polar_base
         class(c_test_qn_solver_2d_polar_base), intent(in) :: self
         sll_int32 :: bcs(2)
      end function i_func_get_bcs

      ! Get parameters
      pure subroutine i_subr_get_parameters(self, adiabatic_electrons, &
                                            use_zonal_flow, epsilon_0)
         use sll_m_working_precision
         import c_test_qn_solver_2d_polar_base
         class(c_test_qn_solver_2d_polar_base), intent(in) :: self
         logical, intent(out) :: adiabatic_electrons
         logical, intent(out) :: use_zonal_flow
         sll_real64, intent(out) :: epsilon_0
      end subroutine i_subr_get_parameters

      ! 1D radial profile, scalar real function
      pure function i_func_1d_real(self, r) result(val)
         use sll_m_working_precision
         import c_test_qn_solver_2d_polar_base
         class(c_test_qn_solver_2d_polar_base), intent(in) :: self
         sll_real64, intent(in) :: r
         sll_real64 :: val
      end function i_func_1d_real

      ! 2D polar profile, scalar real function
      pure function i_func_2d_real(self, r, th) result(val)
         use sll_m_working_precision
         import c_test_qn_solver_2d_polar_base
         class(c_test_qn_solver_2d_polar_base), intent(in) :: self
         sll_real64, intent(in) :: r
         sll_real64, intent(in) :: th
         sll_real64 :: val
      end function i_func_2d_real

   end interface

!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
contains
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

   pure function f_test_case__rho(self, r, th) result(rho_c1)
      class(c_test_qn_solver_2d_polar_base), intent(in) :: self
      sll_real64, intent(in) :: r
      sll_real64, intent(in) :: th
      sll_real64 :: rho_c1

      sll_real64 :: g
      sll_real64 :: dg_dr
      sll_real64 :: dphi

      logical    :: adiabatic_electrons
      logical    :: use_zonal_flow
      sll_real64 :: epsilon_0

      call self%get_parameters(adiabatic_electrons, use_zonal_flow, epsilon_0)

      associate (rho_m0 => self%rho_m0(r), &
                 rho_m0_diff1_r => self%rho_m0_diff1_r(r), &
                 b_magn => self%b_magn(r), &
                 b_magn_diff1_r => self%b_magn_diff1_r(r))

         ! Utility function: g(r) = rho_m0 / (B^2 epsilon_0)
         g = rho_m0/(b_magn**2*epsilon_0)
         dg_dr = (rho_m0_diff1_r*b_magn**2 &
                  - rho_m0*2.0_f64*b_magn*b_magn_diff1_r)/(b_magn**4*epsilon_0)
      end associate

      rho_c1 = g*self%phi_ex_diff2_r(r, th) &
               + (g/r + dg_dr)*self%phi_ex_diff1_r(r, th) &
               + g/r**2*self%phi_ex_diff2_th(r, th)

      if (adiabatic_electrons) then
         if (use_zonal_flow) then
            dphi = self%phi_ex(r, th) - self%phi_ex_avg_th(r, th)
         else
            dphi = self%phi_ex(r, th)
         end if
         rho_c1 = rho_c1 - dphi/self%lambda(r)**2
      end if

      rho_c1 = -epsilon_0*rho_c1

   end function f_test_case__rho

end module m_test_qn_solver_2d_polar_base
