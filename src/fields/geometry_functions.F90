module geometry_functions
#include "sll_working_precision.h"
  
#include "sll_assert.h"

!  use sll_splines
  use sll_constants
  implicit none
  
  sll_real64, parameter :: c1_test = 0.1_f64
  sll_real64, parameter :: c2_test = 0.1_f64

!  type(sll_spline_2D), pointer :: spl2D_x1, spl2D_x2 

contains
  
  ! geometry functions should provide direct mapping, inverse mapping and 
  ! jacobian all of these should be implement following the examples using 
  ! a common name to identity on specific mapping

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

  ! affine function (logical mesh is [0,1]x[0,1])
  ! x1 = (b1-a1)*eta1 + a1; x2=(b2-a2)*eta2 + a2
  !-------------------
#define A1 (-1.0_f64)
#define B1  1.0_f64
#define A2 (-1.0_f64)
#define B2  1.0_f64
  ! direct mapping
  function affine_x1 ( eta1, eta2 )
    sll_real64  :: affine_x1
    sll_real64, intent(in)   :: eta1
    sll_real64, intent(in)   :: eta2
    affine_x1 = (B1-A1)*eta1 + A1
  end function affine_x1

  function affine_x2 ( eta1, eta2 )
    sll_real64  :: affine_x2
    sll_real64, intent(in)   :: eta1
    sll_real64, intent(in)   :: eta2
    affine_x2 = (B2-A2)*eta2 + A2
  end function affine_x2

  ! jacobian maxtrix
  function affine_jac11 ( eta1, eta2 )
    sll_real64  :: affine_jac11
    sll_real64, intent(in)   :: eta1
    sll_real64, intent(in)   :: eta2
    affine_jac11 = B1-A1
  end function affine_jac11

    function affine_jac12 ( eta1, eta2 )
    sll_real64  :: affine_jac12
    sll_real64, intent(in)   :: eta1
    sll_real64, intent(in)   :: eta2
    affine_jac12 = 0.0_f64
  end function affine_jac12

  function affine_jac21 ( eta1, eta2 )
    sll_real64  :: affine_jac21
    sll_real64, intent(in)   :: eta1
    sll_real64, intent(in)   :: eta2
    affine_jac21 = 0.0_f64
  end function affine_jac21

  function affine_jac22 ( eta1, eta2 )
    sll_real64  :: affine_jac22
    sll_real64, intent(in)   :: eta1
    sll_real64, intent(in)   :: eta2
    affine_jac22 = B2-A2
  end function affine_jac22

  ! jacobian ie determinant of jacobian matrix
  function affine_jac ( eta1, eta2 )
    sll_real64  :: affine_jac
    sll_real64, intent(in)   :: eta1
    sll_real64, intent(in)   :: eta2
    affine_jac = (B1-A1) * (B2-A2)
  end function affine_jac

#undef A1
#undef B1
#undef A2
#undef B2
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

  ! sinusoidal product (see P. Colella et al. JCP 230 (2011) formula 
  ! (102) p 2968)
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

  ! Alternative formulation for the polar coordinate transformation:
  !
  ! X1 = (Rmin + (Rmax-Rmin)*eta1)*cos(2*pi*eta2)
  ! X2 = (Rmin + (Rmax-Rmin)*eta1)*sin(2*pi*eta2)
  !
  ! Where eta1 and eta2 are defined in the interval [0,1]. Following are
  ! the functions that embody this transformation. This is used for testing
  ! purposes, hence the parameters R1 and R2 do not form part of the functions'
  ! interfaces. This may become a limitation and should be discussed further.
#define R1 0.1_f64
#define R2 1.0_f64
  function x1_polar_f( eta1, eta2 )
    sll_real64 :: x1_polar_f
    sll_real64, intent(in) :: eta1, eta2
    x1_polar_f = (R1 + (R2-R1)*eta1)*cos(2.0_f64*sll_pi*eta2)
  end function x1_polar_f

  function x2_polar_f( eta1, eta2 )
    sll_real64 :: x2_polar_f
    sll_real64, intent(in) :: eta1, eta2
    x2_polar_f = (R1 + (R2-R1)*eta1)*sin(2.0_f64*sll_pi*eta2)
  end function x2_polar_f

  function deriv_x1_polar_f_eta1( eta1, eta2 )
    sll_real64 :: deriv_x1_polar_f_eta1
    sll_real64, intent(in) :: eta1, eta2
    deriv_x1_polar_f_eta1 = (R2-R1)*cos(2.0_f64*sll_pi*eta2)
  end function deriv_x1_polar_f_eta1

  function deriv_x1_polar_f_eta2( eta1, eta2 )
    sll_real64 :: deriv_x1_polar_f_eta2
    sll_real64, intent(in) :: eta1, eta2
    sll_real64 :: k
    k = 2.0_f64*sll_pi
    deriv_x1_polar_f_eta2 = -(R1+(R2-R1)*eta1)*sin(k*eta2)*k
  end function deriv_x1_polar_f_eta2

  function deriv_x2_polar_f_eta1( eta1, eta2 )
    sll_real64 :: deriv_x2_polar_f_eta1
    sll_real64, intent(in) :: eta1, eta2
    deriv_x2_polar_f_eta1 = (R2-R1)*sin(2.0_f64*sll_pi*eta2)
  end function deriv_x2_polar_f_eta1

  function deriv_x2_polar_f_eta2( eta1, eta2 )
    sll_real64 :: deriv_x2_polar_f_eta2
    sll_real64, intent(in) :: eta1, eta2
    sll_real64 :: k
    k = 2.0_f64*sll_pi
    deriv_x2_polar_f_eta2 = (R1+(R2-R1)*eta1)*cos(k*eta2)*k
  end function deriv_x2_polar_f_eta2

  function jacobian_polar_f( eta1, eta2 ) result(jac)
    sll_real64             :: jac
    sll_real64, intent(in) :: eta1, eta2
    jac = 2.0_f64*sll_pi*(R1+(R2-R1)*eta1)*(R2-R1)
  end function jacobian_polar_f

  function deriv1_jacobian_polar_f( ) result(deriv)
    sll_real64             :: deriv
    deriv = 2.0_f64*sll_pi*(R2-R1)**2
  end function deriv1_jacobian_polar_f


#undef R1
#undef R2

  function zero_function(eta1, eta2)
    sll_real64, intent(in) :: eta1, eta2
    sll_real64             :: zero_function
    zero_function = 0.0_f64
  end function zero_function
    
  !************************************************************************
  ! 1D maps
  !************************************************************************

#define A (-1.0_f64)
#define B  1.0_f64

  function linear_map_f( eta ) result(val)
    sll_real64 :: val
    sll_real64, intent(in) :: eta
    val = (B-A)*eta + A
  end function linear_map_f

  function linear_map_jac_f( eta ) result(val)
    sll_real64 :: val
    sll_real64, intent(in) :: eta
    val = (B-A)
  end function linear_map_jac_f

#undef A
#undef B

#define A 0.0_f64
#define B 6.2831853071795862_f64

  function linear_map_poisson_f( eta ) result(val)
    sll_real64 :: val
    sll_real64, intent(in) :: eta
    val = (B-A)*eta + A
  end function linear_map_poisson_f

  function linear_map_poisson_jac_f( eta ) result(val)
    sll_real64 :: val
    sll_real64, intent(in) :: eta
    val = (B-A)
  end function linear_map_poisson_jac_f

#undef A
#undef B
  

end module geometry_functions

