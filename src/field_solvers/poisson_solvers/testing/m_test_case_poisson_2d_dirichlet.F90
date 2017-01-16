module m_test_case_poisson_2d_dirichlet
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_working_precision.h"

  use m_test_case_poisson_2d_base, only: c_test_case_poisson_2d_polar

  implicit none

  public :: t_test_dirichlet_zero_error
  public :: t_test_dirichlet

  private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  type, extends(c_test_case_poisson_2d_polar) :: t_test_dirichlet_zero_error

    !---------------------------------------------------------------------------
    ! \phi(r,th) = (r-rmin)(r-rmax)(a+sin(k(th-th0))
    !---------------------------------------------------------------------------
    sll_real64, private :: a   = 0.1_f64
    sll_int32 , private :: k   = 3
    sll_real64, private :: th0 = 0.3_f64

  contains
    ! 2D manufactured solution
    procedure :: phi_ex          => dirichlet_zero_error_phi_ex
    procedure :: phi_ex_diff1_r  => dirichlet_zero_error_phi_ex_d1r
    procedure :: phi_ex_diff2_r  => dirichlet_zero_error_phi_ex_d2r
    procedure :: phi_ex_diff2_th => dirichlet_zero_error_phi_ex_d2th

  end type t_test_dirichlet_zero_error

  type, extends(c_test_case_poisson_2d_polar) :: t_test_dirichlet

    !---------------------------------------------------------------------------
    ! \phi(r,th) = r(r-rmin)(r-rmax)(a+sin(k(th-th0))
    !---------------------------------------------------------------------------
    sll_real64, private :: a   = 0.1_f64
    sll_int32 , private :: k   = 3
    sll_real64, private :: th0 = 0.3_f64

  contains
    ! 2D manufactured solution
    procedure :: phi_ex          => f_test__phi_ex
    procedure :: phi_ex_diff1_r  => f_test__phi_ex_diff1_r
    procedure :: phi_ex_diff2_r  => f_test__phi_ex_diff2_r
    procedure :: phi_ex_diff2_th => f_test__phi_ex_diff2_th

  end type t_test_dirichlet

!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
contains
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

  pure function dirichlet_zero_error_phi_ex( self, r, th ) result( val )
    class(t_test_dirichlet_zero_error), intent(in) :: self
    sll_real64                        , intent(in) :: r
    sll_real64                        , intent(in) :: th
    sll_real64 :: val

    val = (r-self%rmin)*(r-self%rmax)*(self%a+sin( self%k*(th-self%th0) ))

  end function dirichlet_zero_error_phi_ex


  pure function dirichlet_zero_error_phi_ex_d1r( self, r, th ) result( val )
    class(t_test_dirichlet_zero_error), intent(in) :: self
    sll_real64                        , intent(in) :: r
    sll_real64                        , intent(in) :: th
    sll_real64 :: val

    associate( rmin => self%rmin, rmax => self%rmax, &
               a => self%a, k => self%k, th0 => self%th0 )

      val = (2.0_f64*r-rmin-rmax)*(a+sin( k*(th-th0) ))

    end associate

  end function dirichlet_zero_error_phi_ex_d1r


  pure function dirichlet_zero_error_phi_ex_d2r( self, r, th ) result( val )
    class(t_test_dirichlet_zero_error), intent(in) :: self
    sll_real64                        , intent(in) :: r
    sll_real64                        , intent(in) :: th
    sll_real64 :: val

    val = 2.0_f64*(self%a+sin( self%k*(th-self%th0) ))

  end function dirichlet_zero_error_phi_ex_d2r


  pure function dirichlet_zero_error_phi_ex_d2th( self, r, th ) result( val )
    class(t_test_dirichlet_zero_error), intent(in) :: self
    sll_real64                        , intent(in) :: r
    sll_real64                        , intent(in) :: th
    sll_real64 :: val

    associate( rmin => self%rmin, rmax => self%rmax, &
               k => self%k, th0 => self%th0 )

      val = -(r-rmin)*(r-rmax)*k**2*sin( k*(th-th0) )

    end associate

  end function dirichlet_zero_error_phi_ex_d2th

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  pure function f_test__phi_ex( self, r, th ) result( val )
    class(t_test_dirichlet), intent(in) :: self
    sll_real64                        , intent(in) :: r
    sll_real64                        , intent(in) :: th
    sll_real64 :: val

    val = r*(r-self%rmin)*(r-self%rmax)*(self%a+sin( self%k*(th-self%th0) ))

  end function f_test__phi_ex


  pure function f_test__phi_ex_diff1_r( self, r, th ) result( val )
    class(t_test_dirichlet), intent(in) :: self
    sll_real64                        , intent(in) :: r
    sll_real64                        , intent(in) :: th
    sll_real64 :: val

    associate( rmin => self%rmin, rmax => self%rmax, &
               a => self%a, k => self%k, th0 => self%th0 )

      val = (3.0_f64*r**2-2.0_f64*(rmin+rmax)*r+rmin*rmax)*(a+sin( k*(th-th0) ))

    end associate

  end function f_test__phi_ex_diff1_r


  pure function f_test__phi_ex_diff2_r( self, r, th ) result( val )
    class(t_test_dirichlet), intent(in) :: self
    sll_real64                        , intent(in) :: r
    sll_real64                        , intent(in) :: th
    sll_real64 :: val

    associate( rmin => self%rmin, rmax => self%rmax, &
               a => self%a, k => self%k, th0 => self%th0 )

      val = (6.0_f64*r-2.0_f64*(rmin+rmax))*(a+sin( k*(th-th0) ))

    end associate

  end function f_test__phi_ex_diff2_r


  pure function f_test__phi_ex_diff2_th( self, r, th ) result( val )
    class(t_test_dirichlet), intent(in) :: self
    sll_real64                        , intent(in) :: r
    sll_real64                        , intent(in) :: th
    sll_real64 :: val

    associate( rmin => self%rmin, rmax => self%rmax, &
               k => self%k, th0 => self%th0 )

      val = -r*(r-rmin)*(r-rmax)*k**2*sin( k*(th-th0) )

    end associate

  end function f_test__phi_ex_diff2_th

end module m_test_case_poisson_2d_dirichlet
