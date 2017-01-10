module m_test_case_2d_neumann_mode0
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_working_precision.h"

  use m_test_case_2d_base, only: c_test_case_qn_solver_2d_polar

  implicit none

  public :: t_test_neumann_mode0_zero_error

  private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  type, extends(c_test_case_qn_solver_2d_polar) :: t_test_neumann_mode0_zero_error

    !---------------------------------------------------------------------------
    ! \phi(r,th) = a*Q(r) + b*cos(m*(th-th0))*R(r)
    !
    ! Q(r) = (rmax-r)(r+rmax-2rmin)/(rmax-rmin)^2
    ! R(r) = 4(rmax-r)(r-rmin)/(rmax-rmin)^2
    !---------------------------------------------------------------------------
    sll_real64, private :: a = 0.0001_f64
    sll_real64, private :: b = 1000.0_f64
    sll_real64, private :: th0 = 0.3_f64
    sll_int32,  private :: m = 3

  contains
    ! 1D input profiles
    procedure :: rho_m0          => f_test__rho_m0
    procedure :: rho_m0_diff1_r  => f_test__rho_m0_diff1_r
    procedure :: b_magn          => f_test__b_magn
    procedure :: b_magn_diff1_r  => f_test__b_magn_diff1_r
    procedure :: lambda          => f_test__lambda
    ! 2D manufactured solution
    procedure :: phi_ex          => f_test__phi_ex
    procedure :: phi_ex_avg_th   => f_test__phi_ex_avg_th
    procedure :: phi_ex_diff1_r  => f_test__phi_ex_diff1_r
    procedure :: phi_ex_diff2_r  => f_test__phi_ex_diff2_r
    procedure :: phi_ex_diff2_th => f_test__phi_ex_diff2_th

  end type t_test_neumann_mode0_zero_error

!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
contains
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

  !=============================================================================
  ! 1D PROFILES
  !=============================================================================

  pure function f_test__rho_m0( self, r ) result( val )
    class(t_test_neumann_mode0_zero_error), intent(in) :: self
    sll_real64, intent(in) :: r
    sll_real64 :: val

    associate( rmin => self%rmin, rmax => self%rmax )
      val = 10.0_f64*(rmax-r)/(rmax-rmin) + 2.0_f64*(r-rmin)/(rmax-rmin)
    end associate

  end function f_test__rho_m0

  pure function f_test__rho_m0_diff1_r( self, r ) result( val )
    class(t_test_neumann_mode0_zero_error), intent(in) :: self
    sll_real64, intent(in) :: r
    sll_real64 :: val

    associate( rmin => self%rmin, rmax => self%rmax )
      val = -8.0_f64/(rmax-rmin)
    end associate

  end function f_test__rho_m0_diff1_r 

  pure function f_test__b_magn( self, r ) result( val )
    class(t_test_neumann_mode0_zero_error), intent(in) :: self
    sll_real64, intent(in) :: r
    sll_real64 :: val

    val = 1.0_f64

  end function f_test__b_magn

  pure function f_test__b_magn_diff1_r( self, r ) result( val )
    class(t_test_neumann_mode0_zero_error), intent(in) :: self
    sll_real64, intent(in) :: r
    sll_real64 :: val

    val = 0.0_f64

  end function f_test__b_magn_diff1_r

  pure function f_test__lambda( self, r ) result( val )
    class(t_test_neumann_mode0_zero_error), intent(in) :: self
    sll_real64, intent(in) :: r
    sll_real64 :: val

    associate( rmean => 0.5_f64*(self%rmax+self%rmin), &
               width => 0.1_f64*(self%rmax-self%rmin) )
      val = 0.01_f64 * (1.0_f64 + exp( -(r-rmean)**2/width**2 ))
    end associate

  end function f_test__lambda

  !=============================================================================
  ! 2D MANUFACTURED SOLUTION
  !=============================================================================

  pure function f_test__phi_ex( self, r, th ) result( val )
    class(t_test_neumann_mode0_zero_error), intent(in) :: self
    sll_real64, intent(in) :: r
    sll_real64, intent(in) :: th
    sll_real64 :: val

    associate( rmin => self%rmin, rmax => self%rmax, &
               a => self%a, b => self%b, m => self%m, th0 => self%th0 )

      val = a*(rmax-r)*(r+rmax-2.0_f64*rmin)/(rmax-rmin)**2 &
            +b*cos(m*(th-th0))*4.0_f64*(rmax-r)*(r-rmin)/(rmax-rmin)**2

    end associate

  end function f_test__phi_ex

  pure function f_test__phi_ex_avg_th( self, r, th ) result( val )
    class(t_test_neumann_mode0_zero_error), intent(in) :: self
    sll_real64, intent(in) :: r
    sll_real64, intent(in) :: th
    sll_real64 :: val

    associate( rmin => self%rmin, rmax => self%rmax, a => self%a )

      val = a*(rmax-r)*(r+rmax-2.0_f64*rmin)/(rmax-rmin)**2

    end associate

  end function f_test__phi_ex_avg_th

  pure function f_test__phi_ex_diff1_r( self, r, th ) result( val )
    class(t_test_neumann_mode0_zero_error), intent(in) :: self
    sll_real64, intent(in) :: r
    sll_real64, intent(in) :: th
    sll_real64 :: val

    associate( rmin => self%rmin, rmax => self%rmax, &
               a => self%a, b => self%b, m => self%m, th0 => self%th0 )

      val = a*2.0_f64*(rmin-r)/(rmax-rmin)**2 &
            +b*cos(m*(th-th0))*4.0_f64*(-2.0_f64*r+rmin+rmax)/(rmax-rmin)**2

    end associate

  end function f_test__phi_ex_diff1_r

  pure function f_test__phi_ex_diff2_r( self, r, th ) result( val )
    class(t_test_neumann_mode0_zero_error), intent(in) :: self
    sll_real64, intent(in) :: r
    sll_real64, intent(in) :: th
    sll_real64 :: val

    associate( rmin => self%rmin, rmax => self%rmax, &
               a => self%a, b => self%b, m => self%m, th0 => self%th0 )

      val = -a*2.0_f64/(rmax-rmin)**2-b*cos(m*(th-th0))*8.0_f64/(rmax-rmin)**2

    end associate

  end function f_test__phi_ex_diff2_r

  pure function f_test__phi_ex_diff2_th( self, r, th ) result( val )
    class(t_test_neumann_mode0_zero_error), intent(in) :: self
    sll_real64, intent(in) :: r
    sll_real64, intent(in) :: th
    sll_real64 :: val

    associate( rmin => self%rmin, rmax => self%rmax, &
               b => self%b, m => self%m, th0 => self%th0 )

      val = -b*real(m**2,f64)*cos(m*(th-th0))*4.0_f64*(rmax-r)*(r-rmin)/(rmax-rmin)**2

    end associate

  end function f_test__phi_ex_diff2_th

end module m_test_case_2d_neumann_mode0
