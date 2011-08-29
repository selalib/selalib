module geometry_functions
#include "sll_working_precision.h"

use numeric_constants
implicit none

sll_real64, parameter :: c1_test = 0.1_f64
sll_real64, parameter :: c2_test = 0.0_f64
contains
!
! geometry functions should provide direct mapping, inverse mapping and jacobian
! all of these should be implement following the examples using a common name to identity
! on specific mapping

  ! identity function
  !-------------------
  ! direct mapping
  function identity_x1 ( eta1, eta2 )
    sll_real64  :: identity_x1
    sll_real64, intent(in)   :: eta1
    sll_real64, intent(in)   :: eta2
    identity_x1 = eta1
  end function identity_x1

  function identity_x2 ( eta1, eta2 )
    sll_real64  :: identity_x2
    sll_real64, intent(in)   :: eta1
    sll_real64, intent(in)   :: eta2
    identity_x2 = eta2
  end function identity_x2

  ! inverse mapping
  function identity_eta1 ( x1, x2 )
    sll_real64  :: identity_eta1
    sll_real64, intent(in)   :: x1
    sll_real64, intent(in)   :: x2
    identity_eta1 = x1
  end function identity_eta1

  function identity_eta2 ( x1, x2 )
    sll_real64  :: identity_eta2
    sll_real64, intent(in)   :: x1
    sll_real64, intent(in)   :: x2
    identity_eta2 = x2
  end function identity_eta2

  ! jacobian maxtrix
  function identity_jac11 ( eta1, eta2 )
    sll_real64  :: identity_jac11
    sll_real64, intent(in)   :: eta1
    sll_real64, intent(in)   :: eta2
    identity_jac11 = 1.0_f64
  end function identity_jac11

    function identity_jac12 ( eta1, eta2 )
    sll_real64  :: identity_jac12
    sll_real64, intent(in)   :: eta1
    sll_real64, intent(in)   :: eta2
    identity_jac12 = 0.0_f64
  end function identity_jac12

  function identity_jac21 ( eta1, eta2 )
    sll_real64  :: identity_jac21
    sll_real64, intent(in)   :: eta1
    sll_real64, intent(in)   :: eta2
    identity_jac21 = 0.0_f64
  end function identity_jac21

  function identity_jac22 ( eta1, eta2 )
    sll_real64  :: identity_jac22
    sll_real64, intent(in)   :: eta1
    sll_real64, intent(in)   :: eta2
    identity_jac22 = 1.0_f64
  end function identity_jac22

  ! jacobian ie determinant of jacobian matrix
  function identity_jac ( eta1, eta2 )
    sll_real64  :: identity_jac
    sll_real64, intent(in)   :: eta1
    sll_real64, intent(in)   :: eta2
    identity_jac = 1.0_f64
  end function identity_jac

  ! polar coordinates (r = eta1, theta = eta2)
  ! x1 = eta1 * cos (eta2)
  ! x2 = eta1 * sin (eta2)
  !-------------------
  ! direct mapping
  function polar_x1 ( eta1, eta2 )
    sll_real64  :: polar_x1
    sll_real64, intent(in)   :: eta1
    sll_real64, intent(in)   :: eta2
    polar_x1 = eta1 * cos( eta2 )
  end function polar_x1

  function polar_x2 ( eta1, eta2 )
    sll_real64  :: polar_x2
    sll_real64, intent(in)   :: eta1
    sll_real64, intent(in)   :: eta2
    polar_x2 = eta1 * sin( eta2 )
  end function polar_x2

  ! inverse mapping
  function polar_eta1 ( x1, x2 )
    sll_real64  :: polar_eta1
    sll_real64, intent(in)   :: x1
    sll_real64, intent(in)   :: x2
    polar_eta1 = sqrt( x1*x1 + x2*x2 )
  end function polar_eta1

  function polar_eta2 ( x1, x2 )
    sll_real64  :: polar_eta2
    sll_real64, intent(in)   :: x1
    sll_real64, intent(in)   :: x2
    polar_eta2 = atan( x2 / x1 ) 
  end function polar_eta2

  ! jacobian matrix
  function polar_jac11 ( eta1, eta2 )
    sll_real64  :: polar_jac11
    sll_real64, intent(in)   :: eta1
    sll_real64, intent(in)   :: eta2
    polar_jac11 = cos ( eta2) 
  end function polar_jac11

    function polar_jac12 ( eta1, eta2 )
    sll_real64  :: polar_jac12
    sll_real64, intent(in)   :: eta1
    sll_real64, intent(in)   :: eta2
    polar_jac12 = - eta1 * sin( eta2)
  end function polar_jac12

  function polar_jac21 ( eta1, eta2 )
    sll_real64  :: polar_jac21
    sll_real64, intent(in)   :: eta1
    sll_real64, intent(in)   :: eta2
    polar_jac21 = sin ( eta2 )
  end function polar_jac21

  function polar_jac22 ( eta1, eta2 )
    sll_real64  :: polar_jac22
    sll_real64, intent(in)   :: eta1
    sll_real64, intent(in)   :: eta2
    polar_jac22 = eta1 * cos ( eta2 )
  end function polar_jac22

 ! jacobian ie determinant of jacobian matrix
  function polar_jac ( eta1, eta2 )
    sll_real64  :: polar_jac
    sll_real64, intent(in)   :: eta1
    sll_real64, intent(in)   :: eta2
    polar_jac = eta1
  end function polar_jac

  ! sinusoidal product (see P. Colella et al. JCP 230 (2011) formula (102) p 2968)
  ! x1 = eta1 + 0.1 * sin(2*pi*eta1) * sin(2*pi*eta2)
  ! x2 = eta2 + 0.1 * sin(2*pi*eta1) * sin(2*pi*eta2)
  !-------------------
  ! direct mapping
  function sinprod_x1 ( eta1, eta2 )
    sll_real64  :: sinprod_x1
    sll_real64, intent(in)   :: eta1
    sll_real64, intent(in)   :: eta2
    sinprod_x1 = eta1 + 0.1_f64 * sin(2*sll_pi*eta1) * sin(2*sll_pi*eta2)
  end function sinprod_x1

  function sinprod_x2 ( eta1, eta2 )
    sll_real64  :: sinprod_x2
    sll_real64, intent(in)   :: eta1
    sll_real64, intent(in)   :: eta2
    sinprod_x2 = eta2 + 0.1_f64 * sin(2*sll_pi*eta1) * sin(2*sll_pi*eta2)
  end function sinprod_x2

  ! inverse mapping 
  ! cannot be computed analytically in this case. Use fixed point iterations.
  function sinprod_eta1 ( x1, x2 )
    sll_real64  :: sinprod_eta1
    sll_real64, intent(in)   :: x1
    sll_real64, intent(in)   :: x2
    ! NEEDS TO BE IMPLEMENTED
    STOP 'function not implemented'
    sinprod_eta1 = x1
  end function sinprod_eta1

  function sinprod_eta2 ( x1, x2 )
    sll_real64  :: sinprod_eta2
    sll_real64, intent(in)   :: x1
    sll_real64, intent(in)   :: x2
    ! NEEDS TO BE IMPLEMENTED
    STOP 'function not implemented'
    sinprod_eta2 = x2
  end function sinprod_eta2

  ! jacobian matrix
  function sinprod_jac11 ( eta1, eta2 )
    sll_real64  :: sinprod_jac11
    sll_real64, intent(in)   :: eta1
    sll_real64, intent(in)   :: eta2
    sinprod_jac11 = 1.0_f64 + 0.2_f64 *sll_pi * cos (2*sll_pi*eta1) * sin (2*sll_pi*eta2)
  end function sinprod_jac11

    function sinprod_jac12 ( eta1, eta2 )
    sll_real64  :: sinprod_jac12
    sll_real64, intent(in)   :: eta1
    sll_real64, intent(in)   :: eta2
    sinprod_jac12 = 0.2_f64 *sll_pi * sin (2*sll_pi*eta1) * cos (2*sll_pi*eta2)
  end function sinprod_jac12

  function sinprod_jac21 ( eta1, eta2 )
    sll_real64  :: sinprod_jac21
    sll_real64, intent(in)   :: eta1
    sll_real64, intent(in)   :: eta2
    sinprod_jac21 = 0.2_f64 * sll_pi * cos (2*sll_pi*eta1) * sin (2*sll_pi*eta2)
  end function sinprod_jac21

  function sinprod_jac22 ( eta1, eta2 )
    sll_real64  :: sinprod_jac22
    sll_real64, intent(in)   :: eta1
    sll_real64, intent(in)   :: eta2
    sinprod_jac22 = 1.0_f64 + 0.2_f64 * sll_pi * sin (2*sll_pi*eta1) * cos (2*sll_pi*eta2)
  end function sinprod_jac22

   ! jacobian ie determinant of jacobian matrix
  function sinprod_jac ( eta1, eta2 )
    sll_real64  :: sinprod_jac
    sll_real64, intent(in)   :: eta1
    sll_real64, intent(in)   :: eta2
    !sinprod_jac = 1.0_f64 + 0.2_f64 *sll_pi * sin (2*sll_pi**(eta1+eta2)) 
    sinprod_jac = (1.0_f64 + 0.2_f64 *sll_pi * cos (2*sll_pi*eta1) * sin (2*sll_pi*eta2)) * &
         (1.0_f64 + 0.2_f64 * sll_pi * sin (2*sll_pi*eta1) * cos (2*sll_pi*eta2)) - &
         0.2_f64 *sll_pi * sin (2*sll_pi*eta1) * cos (2*sll_pi*eta2) * &
         0.2_f64 * sll_pi * cos (2*sll_pi*eta1) * sin (2*sll_pi*eta2)
    
  end function sinprod_jac

  ! test function
  !-------------------
  ! direct mapping
  function test_x1 ( eta1, eta2 )
    sll_real64  :: test_x1
    sll_real64, intent(in)   :: eta1
    sll_real64, intent(in)   :: eta2
    test_x1 = eta1 + c1_test * sin( 2.0_f64* sll_pi * eta1 )
    !test_x1 = eta1**2
  end function test_x1

  function test_x2 ( eta1, eta2 )
    sll_real64  :: test_x2
    sll_real64, intent(in)   :: eta1
    sll_real64, intent(in)   :: eta2
    test_x2 = eta2 + c2_test * sin( 2.0_f64* sll_pi * eta2 )
  end function test_x2

  ! inverse mapping
  function test_eta1 ( x1, x2 )
    sll_real64  :: test_eta1
    sll_real64, intent(in)   :: x1
    sll_real64, intent(in)   :: x2
    test_eta1 = x1 / c1_test
  end function test_eta1

  function test_eta2 ( x1, x2 )
    sll_real64  :: test_eta2
    sll_real64, intent(in)   :: x1
    sll_real64, intent(in)   :: x2
    test_eta2 = x2 / c2_test
  end function test_eta2

  ! inverse jacobian matrix
  function test_jac11 ( eta1, eta2 )
    sll_real64  :: test_jac11
    sll_real64, intent(in)   :: eta1
    sll_real64, intent(in)   :: eta2
    test_jac11 = 1.0_f64 / (1.0_f64 + 2.0_f64 * sll_pi* c1_test * cos( 2.0_f64* sll_pi * eta1))
  end function test_jac11

    function test_jac12 ( eta1, eta2 )
    sll_real64  :: test_jac12
    sll_real64, intent(in)   :: eta1
    sll_real64, intent(in)   :: eta2
    test_jac12 = 0.0_f64
  end function test_jac12

  function test_jac21 ( eta1, eta2 )
    sll_real64  :: test_jac21
    sll_real64, intent(in)   :: eta1
    sll_real64, intent(in)   :: eta2
    test_jac21 = 0.0_f64
  end function test_jac21

  function test_jac22 ( eta1, eta2 )
    sll_real64  :: test_jac22
    sll_real64, intent(in)   :: eta1
    sll_real64, intent(in)   :: eta2
    test_jac22 = 1.0_f64 / (1.0_f64 + 2.0_f64 * sll_pi* c2_test * cos( 2.0_f64* sll_pi * eta2))
  end function test_jac22

  ! jacobian ie determinant of jacobian matrix
  function test_jac ( eta1, eta2 )
    sll_real64  :: test_jac
    sll_real64, intent(in)   :: eta1
    sll_real64, intent(in)   :: eta2
    test_jac =  (1.0_f64 + 2.0_f64 * sll_pi* c1_test * cos( 2.0_f64* sll_pi * eta1)) * &
         (1.0_f64 + 2.0_f64 * sll_pi* c2_test * cos( 2.0_f64* sll_pi * eta2))
    !test_jac =  2 * eta1!
  end function test_jac
end module geometry_functions

