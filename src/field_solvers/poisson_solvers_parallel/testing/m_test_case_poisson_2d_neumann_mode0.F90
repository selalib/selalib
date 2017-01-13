module m_test_case_poisson_2d_neumann_mode0
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_working_precision.h"

  use m_test_case_poisson_2d_base, only: c_test_case_poisson_2d_polar

  implicit none

  public :: t_test_neumann_mode0_zero_error

  private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  type, extends(c_test_case_poisson_2d_polar) :: t_test_neumann_mode0_zero_error

    !---------------------------------------------------------------------------
    ! \phi(r,th) = a(r-rmax)(r+rmax-2rmin) + b(r-rmin)(r-rmax)sin(k(th-th0))
    !---------------------------------------------------------------------------
    sll_real64, private :: a   = 0.1_f64
    sll_real64, private :: b   = 1.0_f64
    sll_real64, private :: th0 = 0.3_f64
    sll_int32,  private :: k   = 3

  contains
    ! 2D manufactured solution
    procedure :: phi_ex          => neumann_mode0_zero_error_phi_ex
    procedure :: phi_ex_diff1_r  => neumann_mode0_zero_error_phi_ex_d1r
    procedure :: phi_ex_diff2_r  => neumann_mode0_zero_error_phi_ex_d2r
    procedure :: phi_ex_diff2_th => neumann_mode0_zero_error_phi_ex_d2th

  end type t_test_neumann_mode0_zero_error

!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
contains
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

  pure function neumann_mode0_zero_error_phi_ex( self, r, th ) result( val )
    class(t_test_neumann_mode0_zero_error), intent(in) :: self
    sll_real64                        , intent(in) :: r
    sll_real64                        , intent(in) :: th
    sll_real64 :: val

    val = self%a*(r-self%rmax)*(r+self%rmax-2.0_f64*self%rmin) &
          + self%b*(r-self%rmin)*(r-self%rmax)*sin( self%k*(th-self%th0) )

  end function neumann_mode0_zero_error_phi_ex


  pure function neumann_mode0_zero_error_phi_ex_d1r( self, r, th ) result( val )
    class(t_test_neumann_mode0_zero_error), intent(in) :: self
    sll_real64                        , intent(in) :: r
    sll_real64                        , intent(in) :: th
    sll_real64 :: val

    associate( rmin => self%rmin, rmax => self%rmax, &
               a => self%a, b => self%b, k => self%k, th0 => self%th0 )

      val = a*2.0_f64*(r-rmin) + b*(2.0_f64*r-rmin-rmax)*sin( k*(th-th0) )

    end associate

  end function neumann_mode0_zero_error_phi_ex_d1r


  pure function neumann_mode0_zero_error_phi_ex_d2r( self, r, th ) result( val )
    class(t_test_neumann_mode0_zero_error), intent(in) :: self
    sll_real64                        , intent(in) :: r
    sll_real64                        , intent(in) :: th
    sll_real64 :: val

    val = self%a*2.0_f64 + self%b*2.0_f64*sin( self%k*(th-self%th0) )

  end function neumann_mode0_zero_error_phi_ex_d2r


  pure function neumann_mode0_zero_error_phi_ex_d2th( self, r, th ) result( val )
    class(t_test_neumann_mode0_zero_error), intent(in) :: self
    sll_real64                        , intent(in) :: r
    sll_real64                        , intent(in) :: th
    sll_real64 :: val

    associate( rmin => self%rmin, rmax => self%rmax, &
               k => self%k, th0 => self%th0 )

      val = -self%b*(r-rmin)*(r-rmax)*k**2*sin( k*(th-th0) )

    end associate

  end function neumann_mode0_zero_error_phi_ex_d2th

end module m_test_case_poisson_2d_neumann_mode0
