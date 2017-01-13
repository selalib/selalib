module m_test_case_qn_par_2d_base
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_working_precision.h"

implicit none

public :: &
  c_test_case_qn_solver_2d_polar

private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  type, abstract :: c_test_case_qn_solver_2d_polar

    sll_real64 :: rmin
    sll_real64 :: rmax
    logical    :: adiabatic_electrons
    logical    :: use_zonal_flow
    sll_real64 :: epsilon_0
    sll_int32  :: bc_rmin
    sll_int32  :: bc_rmax

  contains
    ! 1D input profiles
    procedure( i_func_1d_real ), deferred :: rho_m0
    procedure( i_func_1d_real ), deferred :: rho_m0_diff1_r
    procedure( i_func_1d_real ), deferred :: b_magn
    procedure( i_func_1d_real ), deferred :: b_magn_diff1_r
    procedure( i_func_1d_real ), deferred :: lambda
    ! 2D manufactured solution
    procedure( i_func_2d_real ), deferred :: phi_ex
    procedure( i_func_2d_real ), deferred :: phi_ex_avg_th
    procedure( i_func_2d_real ), deferred :: phi_ex_diff1_r
    procedure( i_func_2d_real ), deferred :: phi_ex_diff2_r
    procedure( i_func_2d_real ), deferred :: phi_ex_diff2_th
    ! 2D right-hand side to solver
    procedure, non_overridable :: rhs => f_test_case__rhs

  end type c_test_case_qn_solver_2d_polar

  !-----------------------------------------------------------------------------
  abstract interface

    ! 1D radial profile, scalar real function
    pure function i_func_1d_real( self, r ) result( val )
      use sll_m_working_precision
      import c_test_case_qn_solver_2d_polar
      class( c_test_case_qn_solver_2d_polar ), intent(in) :: self
      sll_real64                             , intent(in) :: r
      sll_real64 :: val
    end function i_func_1d_real

    ! 2D polar profile, scalar real function
    pure function i_func_2d_real( self, r, th ) result( val )
      use sll_m_working_precision
      import c_test_case_qn_solver_2d_polar
      class( c_test_case_qn_solver_2d_polar ), intent(in) :: self
      sll_real64                             , intent(in) :: r
      sll_real64                             , intent(in) :: th
      sll_real64 :: val
    end function i_func_2d_real

  end interface

!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
contains
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

  pure function f_test_case__rhs( self, r, th ) result( rho_c1 )
    class(c_test_case_qn_solver_2d_polar), intent(in) :: self
    sll_real64                           , intent(in) :: r
    sll_real64                           , intent(in) :: th
    sll_real64 :: rho_c1

    sll_real64 :: g
    sll_real64 :: dg_dr
    sll_real64 :: dphi

    associate( rho_m0         => self%rho_m0        ( r ), &
               rho_m0_diff1_r => self%rho_m0_diff1_r( r ), &
               b_magn         => self%b_magn        ( r ), &
               b_magn_diff1_r => self%b_magn_diff1_r( r ) )

      ! Utility function: g(r) = rho_m0 / (B^2 epsilon_0)
      g     =  rho_m0 / (b_magn**2 * self%epsilon_0)
      dg_dr = (rho_m0_diff1_r * b_magn**2 &
        - rho_m0 * 2.0_f64 * b_magn * b_magn_diff1_r) &
        / (b_magn**4 * self%epsilon_0)
    end associate

    rho_c1 = g            * self%phi_ex_diff2_r ( r, th ) &
           +(g/r + dg_dr) * self%phi_ex_diff1_r ( r, th ) &
           + g/r**2       * self%phi_ex_diff2_th( r, th )

    if (self%adiabatic_electrons) then
      if (self%use_zonal_flow) then
        dphi = self%phi_ex( r, th ) - self%phi_ex_avg_th( r, th )
      else
        dphi = self%phi_ex( r, th )
      end if
      rho_c1 = rho_c1 - dphi / self%lambda( r )**2
    end if

    rho_c1 = -self%epsilon_0 * rho_c1

  end function f_test_case__rhs

  
end module m_test_case_qn_par_2d_base
