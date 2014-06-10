module sll_common_coordinate_transformations
#include "sll_working_precision.h"
#include "sll_assert.h"
  use sll_constants
  implicit none
  

contains
  
  ! This module provides some common coordinate transformations in terms of the
  ! direct mapping, inverse mapping and jacobian.  All of these should be 
  ! implement following similar naming conventions.

  ! **************************************************************************
  !
  !                         Identity transformation
  !
  ! **************************************************************************

  ! direct mapping
  function identity_x1 ( eta1, eta2, params )
    sll_real64  :: identity_x1
    sll_real64, intent(in)   :: eta1
    sll_real64, intent(in)   :: eta2
    sll_real64, dimension(:), intent(in) :: params
    identity_x1 = eta1
  end function identity_x1

  function identity_x2 ( eta1, eta2, params )
    sll_real64  :: identity_x2
    sll_real64, intent(in)   :: eta1
    sll_real64, intent(in)   :: eta2
    sll_real64, dimension(:), intent(in) :: params
    identity_x2 = eta2
  end function identity_x2

  ! inverse mapping
  function identity_eta1 ( x1, x2, params )
    sll_real64  :: identity_eta1
    sll_real64, intent(in)   :: x1
    sll_real64, intent(in)   :: x2
    sll_real64, dimension(:), intent(in) :: params
    identity_eta1 = x1
  end function identity_eta1

  function identity_eta2 ( x1, x2, params )
    sll_real64  :: identity_eta2
    sll_real64, intent(in)   :: x1
    sll_real64, intent(in)   :: x2
    sll_real64, dimension(:), intent(in) :: params
    identity_eta2 = x2
  end function identity_eta2

  ! jacobian maxtrix
  function identity_jac11 ( eta1, eta2, params )
    sll_real64  :: identity_jac11
    sll_real64, intent(in)   :: eta1
    sll_real64, intent(in)   :: eta2
    sll_real64, dimension(:), intent(in) :: params
    identity_jac11 = 1.0_f64
  end function identity_jac11

  function identity_jac12 ( eta1, eta2, params )
    sll_real64  :: identity_jac12
    sll_real64, intent(in)   :: eta1
    sll_real64, intent(in)   :: eta2
    sll_real64, dimension(:), intent(in) :: params
    identity_jac12 = 0.0_f64
  end function identity_jac12

  function identity_jac21 ( eta1, eta2, params )
    sll_real64  :: identity_jac21
    sll_real64, intent(in)   :: eta1
    sll_real64, intent(in)   :: eta2
    sll_real64, dimension(:), intent(in) :: params
    identity_jac21 = 0.0_f64
  end function identity_jac21

  function identity_jac22 ( eta1, eta2, params )
    sll_real64  :: identity_jac22
    sll_real64, intent(in)   :: eta1
    sll_real64, intent(in)   :: eta2
    sll_real64, dimension(:), intent(in) :: params
    identity_jac22 = 1.0_f64
  end function identity_jac22

  ! jacobian ie determinant of jacobian matrix
  function identity_jac ( eta1, eta2, params )
    sll_real64  :: identity_jac
    sll_real64, intent(in)   :: eta1
    sll_real64, intent(in)   :: eta2
    sll_real64, dimension(:), intent(in) :: params
    identity_jac = 1.0_f64
  end function identity_jac

  ! **************************************************************************
  !
  !        affine transformation (logical mesh is [0,1]x[0,1]):
  !
  !        x1 = (b1-a1)*(cos(alpha)*eta1-sin(alpha)*eta2) + a1
  !        x2 = (b2-a2)*(sin(alpha)*eta1+cos(alpha)*eta2) + a2
  !
  ! **************************************************************************

  ! developer's note: made the choice of the params array as the full 
  ! sequence (A1 B1 A2 B2) as the same params array may thus be passed as
  ! argument any call related with this transformation. 
  ! Delete these formerly default values:
  ! A1 = -1.0
  ! B1 =  1.0
  ! A2 = -1.0
  ! B2 =  1.0
  ! While convenient, there is a risk associated with this: 


  ! direct mapping
  function affine_x1 ( eta1, eta2, params )
    sll_real64  :: affine_x1
    sll_real64, intent(in)   :: eta1
    sll_real64, intent(in)   :: eta2
    sll_real64, dimension(:), intent(in) :: params
    sll_real64 :: A1
    sll_real64 :: B1
    sll_real64 :: alpha

    SLL_ASSERT(size(params) >= 5)
    A1 = params(1)
    B1 = params(2)
    alpha = params(5)
    affine_x1 = (B1-A1)*(cos(alpha)*eta1-sin(alpha)*eta2) + A1
  end function affine_x1

  function affine_x2 ( eta1, eta2, params )
    sll_real64  :: affine_x2
    sll_real64, intent(in)   :: eta1
    sll_real64, intent(in)   :: eta2
    sll_real64, dimension(:), intent(in) :: params
    sll_real64 :: A2
    sll_real64 :: B2
    sll_real64 :: alpha

    SLL_ASSERT(size(params) >= 5)
    A2 = params(3)
    B2 = params(4)
    alpha = params(5)
    affine_x2 = (B2-A2)*(sin(alpha)*eta1+cos(alpha)*eta2) + A2
  end function affine_x2

  ! jacobian maxtrix
  function affine_jac11 ( eta1, eta2, params )
    sll_real64  :: affine_jac11
    sll_real64, intent(in)   :: eta1
    sll_real64, intent(in)   :: eta2
    sll_real64, dimension(:), intent(in) :: params
    sll_real64 :: A1
    sll_real64 :: B1
    sll_real64 :: alpha

    SLL_ASSERT(size(params) >= 5)
    A1 = params(1)
    B1 = params(2)
    alpha = params(5)
    affine_jac11 = (B1-A1)*cos(alpha)
  end function affine_jac11

  function affine_jac12 ( eta1, eta2, params )
    sll_real64  :: affine_jac12
    sll_real64, intent(in)   :: eta1
    sll_real64, intent(in)   :: eta2
    sll_real64, dimension(:), intent(in) :: params
    sll_real64 :: alpha, a1, b1
    SLL_ASSERT(size(params) >= 5)
    a1 = params(1)
    b1 = params(2)
    alpha = params(5)
    affine_jac12 = -(B1-A1)*sin(alpha)
  end function affine_jac12

  function affine_jac21 ( eta1, eta2, params )
    sll_real64  :: affine_jac21
    sll_real64, intent(in)   :: eta1
    sll_real64, intent(in)   :: eta2
    sll_real64 :: alpha, a2, b2
    sll_real64, dimension(:), intent(in) :: params

    SLL_ASSERT(size(params) >= 5)
    a2 = params(3)
    b2 = params(4)
    alpha = params(5)
    affine_jac21 = (b2-a2)*sin(alpha)
  end function affine_jac21

  function affine_jac22 ( eta1, eta2, params )
    sll_real64  :: affine_jac22
    sll_real64, intent(in)   :: eta1
    sll_real64, intent(in)   :: eta2
    sll_real64, dimension(:), intent(in) :: params
    sll_real64 :: A2
    sll_real64 :: B2
    sll_real64 :: alpha

    SLL_ASSERT(size(params) >= 5)
    A2 = params(3)
    B2 = params(4)
    alpha = params(5)

    affine_jac22 = (B2-A2)*cos(alpha)

  end function affine_jac22

  ! jacobian ie determinant of jacobian matrix
  function affine_jac ( eta1, eta2, params )
    sll_real64  :: affine_jac
    sll_real64, intent(in)   :: eta1
    sll_real64, intent(in)   :: eta2
    sll_real64, dimension(:), intent(in) :: params
    sll_real64 :: A1
    sll_real64 :: B1
    sll_real64 :: A2
    sll_real64 :: B2

    SLL_ASSERT(size(params) >= 4)
    A1    = params(1)
    B1    = params(2)
    A2    = params(3)
    B2    = params(4)

    affine_jac = (B1-A1)*(B2-A2)

  end function affine_jac


  ! **************************************************************************
  !
  !        homography transformation (logical mesh is [0,1]x[0,1]):
  !
  !  the homography is used to rectify a perspective image, for example to 
  !  generate a "plan" view of a building from a "perspective" photo. 
  !  We can also convert a square to a trapeze.
  !
  !        x1 = (a*eta1*b*eta2+c)/(g*eta1+h*eta2+1) 
  !        x2 = (d*eta1*e*eta2+f)/(g*eta1+h*eta2+1) 
  !
  !  a = fixed scale factor in x1 direction with scale x2 unchanged.
  !  b = scale factor in x1 direction proportional to x2 distance from origin.
  !  c = origin translation in x1 direction.
  !  d = scale factor in x2 direction proportional to x1 distance from origin.
  !  e = fixed scale factor in x2 direction with scale x1 unchanged.
  !  f = origin translation in x2 direction.
  !  g = proportional scale factors x1 and x2 in function of x1.
  !  h = proportional scale factors x1 and x2 in function of x2.
  !   
  ! **************************************************************************

  ! direct mapping
  function homography_x1 ( eta1, eta2, params )
    sll_real64  :: homography_x1
    sll_real64, intent(in)   :: eta1
    sll_real64, intent(in)   :: eta2
    sll_real64, dimension(:), intent(in) :: params
    sll_real64 :: a, b, c, g, h

    SLL_ASSERT(size(params) >= 8)

    a = params(1)
    b = params(2)
    c = params(3)
    g = params(7)
    h = params(8)

    homography_x1 = (a*eta1+b*eta2+c)/(g*eta1+h*eta2+1)

  end function homography_x1

  function homography_x2 ( eta1, eta2, params )
    sll_real64  :: homography_x2
    sll_real64, intent(in)   :: eta1
    sll_real64, intent(in)   :: eta2
    sll_real64, dimension(:), intent(in) :: params
    sll_real64 :: d, e, f, g, h

    SLL_ASSERT(size(params) >= 8)

    d = params(4)
    e = params(5)
    f = params(6)
    g = params(7)
    h = params(8)

    homography_x2 = (d*eta1+e*eta2+f)/(g*eta1+h*eta2+1)

  end function homography_x2

  ! jacobian maxtrix
  function homography_jac11 ( eta1, eta2, params )
    sll_real64  :: homography_jac11
    sll_real64, intent(in)   :: eta1
    sll_real64, intent(in)   :: eta2
    sll_real64 :: a, b, c, g, h
    sll_real64, dimension(:), intent(in) :: params

    SLL_ASSERT(size(params) >= 8)

    a = params(1)
    b = params(2)
    c = params(3)
    g = params(7)
    h = params(8)


    homography_jac11 = a/(eta1*g+eta2*h+1)-(a*eta1+b*eta2+c)*g/(eta1*g+eta2*h+1)**2

  end function homography_jac11

  function homography_jac12 ( eta1, eta2, params )
    sll_real64  :: homography_jac12
    sll_real64, intent(in)   :: eta1
    sll_real64, intent(in)   :: eta2
    sll_real64, dimension(:), intent(in) :: params
    sll_real64 :: a, b, c, g, h

    SLL_ASSERT(size(params) >= 8)

    a = params(1)
    b = params(2)
    c = params(3)
    g = params(7)
    h = params(8)

    homography_jac12 = b/(eta1*g+eta2*h+1)-(a*eta1+b*eta2+c)*h/(eta1*g+eta2*h+1)**2

  end function homography_jac12

  function homography_jac21 ( eta1, eta2, params )
    sll_real64  :: homography_jac21
    sll_real64, intent(in)   :: eta1
    sll_real64, intent(in)   :: eta2
    sll_real64, dimension(:), intent(in) :: params
    sll_real64 :: d, e, f, g, h

    SLL_ASSERT(size(params) >= 8)

    d = params(4)
    e = params(5)
    f = params(6)
    g = params(7)
    h = params(8)

    homography_jac21 = d/(eta1*g+eta2*h+1)-(d*eta1+e*eta2+f)*g/(eta1*g+eta2*h+1)**2

  end function homography_jac21

  function homography_jac22 ( eta1, eta2, params )
    sll_real64  :: homography_jac22
    sll_real64, intent(in)   :: eta1
    sll_real64, intent(in)   :: eta2
    sll_real64, dimension(:), intent(in) :: params
    sll_real64 :: d, e, f, g, h

    SLL_ASSERT(size(params) >= 8)

    d = params(4)
    e = params(5)
    f = params(6)
    g = params(7)
    h = params(8)

    homography_jac22 = e/(eta1*g+eta2*h+1)-(d*eta1+e*eta2+f)*h/(eta1*g+eta2*h+1)**2

  end function homography_jac22

  ! jacobian ie determinant of jacobian matrix
  function homography_jac ( eta1, eta2, params )
    sll_real64  :: homography_jac
    sll_real64, intent(in)   :: eta1
    sll_real64, intent(in)   :: eta2
    sll_real64, dimension(:), intent(in) :: params
    sll_real64 :: a, b, c, d, e, f, g, h

    SLL_ASSERT(size(params) >= 8)

    a = params(1)
    b = params(2)
    c = params(3)
    d = params(4)
    e = params(5)
    f = params(6)
    g = params(7)
    h = params(8)


    homography_jac = -(b*d-a*e+(c*e-b*f)*g-(c*d-a*f)*h) &
                     /(eta1**3*g**3+eta2**3*h**3+3*eta1**2*g**2 &
                     +3*(eta1*eta2**2*g+eta2**2)*h**2+3*eta1*g &
                     +3*(eta1**2*eta2*g**2+2*eta1*eta2*g+eta2)*h+1)

  end function homography_jac

  ! **************************************************************************
  !
  !        rubber_sheeting transformation 
  !
  !  the rubber-sheeting is similar to homography. The difference is that the 
  !  Homography gives priority to alignments, whereas the Rubber-Sheeting gives 
  !  priority to linear proportions.
  !
  !        x1 = a*eta1*eta2+b*eta1+c*eta2+d
  !        x2 = e*eta1*eta2+f*eta1+g*eta2+h
  !
  !  a = scale factor in x1 direction proportional to the multiplication x1 * x2.
  !  b = fixed scale factor in x1 direction with scale x2 unchanged.
  !  c = scale factor in x1 direction proportional to x2 distance from origin.
  !  d = origin translation in x1 direction.
  !  e = scale factor in x2 direction proportional to the multiplication x1 * x2.
  !  f = fixed scale factor in x2 direction with scale x1 unchanged.
  !  g = scale factor in x2 direction proportional to x1 distance from origin.
  !  h = origin translation in x2 direction.
  !   
  ! **************************************************************************

  ! direct mapping
  function rubber_sheeting_x1 ( eta1, eta2, params )
    sll_real64  :: rubber_sheeting_x1
    sll_real64, intent(in)   :: eta1
    sll_real64, intent(in)   :: eta2
    sll_real64, dimension(:), intent(in) :: params
    sll_real64 :: a, b, c, d

    SLL_ASSERT(size(params) >= 8)

    a = params(1)
    b = params(2)
    c = params(3)
    d = params(4)

    rubber_sheeting_x1 = a*eta1*eta2+b*eta1+c*eta2+d

  end function rubber_sheeting_x1

  function rubber_sheeting_x2 ( eta1, eta2, params )
    sll_real64  :: rubber_sheeting_x2
    sll_real64, intent(in)   :: eta1
    sll_real64, intent(in)   :: eta2
    sll_real64, dimension(:), intent(in) :: params
    sll_real64 :: e, f, g, h

    SLL_ASSERT(size(params) >= 8)

    e = params(5)
    f = params(6)
    g = params(7)
    h = params(8)

    rubber_sheeting_x2 = e*eta1*eta2+f*eta1+g*eta2+h

  end function rubber_sheeting_x2

  ! jacobian maxtrix
  function rubber_sheeting_jac11 ( eta1, eta2, params )
    sll_real64  :: rubber_sheeting_jac11
    sll_real64, intent(in)   :: eta1
    sll_real64, intent(in)   :: eta2
    sll_real64 :: a, b
    sll_real64, dimension(:), intent(in) :: params

    SLL_ASSERT(size(params) >= 8)

    a = params(1)
    b = params(2)

    rubber_sheeting_jac11 = a*eta2+b

  end function rubber_sheeting_jac11

  function rubber_sheeting_jac12 ( eta1, eta2, params )
    sll_real64  :: rubber_sheeting_jac12
    sll_real64, intent(in)   :: eta1
    sll_real64, intent(in)   :: eta2
    sll_real64, dimension(:), intent(in) :: params
    sll_real64 :: a, c

    SLL_ASSERT(size(params) >= 8)

    a = params(1)
    c = params(3)

    rubber_sheeting_jac12 = a*eta1+c

  end function rubber_sheeting_jac12

  function rubber_sheeting_jac21 ( eta1, eta2, params )
    sll_real64  :: rubber_sheeting_jac21
    sll_real64, intent(in)   :: eta1
    sll_real64, intent(in)   :: eta2
    sll_real64, dimension(:), intent(in) :: params
    sll_real64 :: e, f

    SLL_ASSERT(size(params) >= 8)

    e = params(5)
    f = params(6)

    rubber_sheeting_jac21 = e*eta2+f

  end function rubber_sheeting_jac21

  function rubber_sheeting_jac22 ( eta1, eta2, params )
    sll_real64  :: rubber_sheeting_jac22
    sll_real64, intent(in)   :: eta1
    sll_real64, intent(in)   :: eta2
    sll_real64, dimension(:), intent(in) :: params
    sll_real64 :: e, g

    SLL_ASSERT(size(params) >= 8)

    e = params(5)
    g = params(7)

    rubber_sheeting_jac22 = e*eta1+g

  end function rubber_sheeting_jac22

  ! jacobian ie determinant of jacobian matrix
  function rubber_sheeting_jac ( eta1, eta2, params )
    sll_real64  :: rubber_sheeting_jac
    sll_real64, intent(in)   :: eta1
    sll_real64, intent(in)   :: eta2
    sll_real64, dimension(:), intent(in) :: params
    sll_real64 :: a, b, c, e, f, g

    SLL_ASSERT(size(params) >= 8)

    a = params(1)
    b = params(2)
    c = params(3)
    e = params(5)
    f = params(6)
    g = params(7)

    rubber_sheeting_jac = (e*eta1+g)*(a*eta2+b) - (a*eta1+c)*(e*eta2+f)

  end function rubber_sheeting_jac



  ! **************************************************************************
  !
  !       polar coordinate transformation (r = eta1, theta = eta2):
  !
  !        x1 = eta1 * cos (eta2)
  !        x2 = eta1 * sin (eta2)
  !
  ! **************************************************************************

  ! direct mapping
  function polar_x1 ( eta1, eta2, params )
    sll_real64  :: polar_x1
    sll_real64, intent(in)   :: eta1
    sll_real64, intent(in)   :: eta2
    sll_real64, dimension(:), intent(in) :: params
    polar_x1 = eta1 * cos( eta2 )
  end function polar_x1

  function polar_x2 ( eta1, eta2, params )
    sll_real64  :: polar_x2
    sll_real64, intent(in)   :: eta1
    sll_real64, intent(in)   :: eta2
    sll_real64, dimension(:), intent(in) :: params
    polar_x2 = eta1 * sin( eta2 )
  end function polar_x2

  ! inverse mapping
  function polar_eta1 ( x1, x2, params )
    sll_real64  :: polar_eta1
    sll_real64, intent(in)   :: x1
    sll_real64, intent(in)   :: x2
    sll_real64, dimension(:), intent(in) :: params
    polar_eta1 = sqrt( x1*x1 + x2*x2 )
  end function polar_eta1

  function polar_eta2 ( x1, x2, params )
    sll_real64  :: polar_eta2
    sll_real64, intent(in)   :: x1
    sll_real64, intent(in)   :: x2
    sll_real64, dimension(:), intent(in) :: params
    polar_eta2 = atan( x2 / x1 ) 
  end function polar_eta2

  ! jacobian matrix
  function polar_jac11 ( eta1, eta2, params )
    sll_real64  :: polar_jac11
    sll_real64, intent(in)   :: eta1
    sll_real64, intent(in)   :: eta2
    sll_real64, dimension(:), intent(in) :: params
    polar_jac11 = cos ( eta2 ) 
  end function polar_jac11

    function polar_jac12 ( eta1, eta2, params )
    sll_real64  :: polar_jac12
    sll_real64, intent(in)   :: eta1
    sll_real64, intent(in)   :: eta2
    sll_real64, dimension(:), intent(in) :: params
    polar_jac12 = - eta1 * sin( eta2 )
  end function polar_jac12

  function polar_jac21 ( eta1, eta2, params )
    sll_real64  :: polar_jac21
    sll_real64, intent(in)   :: eta1
    sll_real64, intent(in)   :: eta2
    sll_real64, dimension(:), intent(in) :: params
    polar_jac21 = sin ( eta2 )
  end function polar_jac21

  function polar_jac22 ( eta1, eta2, params )
    sll_real64  :: polar_jac22
    sll_real64, intent(in)   :: eta1
    sll_real64, intent(in)   :: eta2
    sll_real64, dimension(:), intent(in) :: params
    polar_jac22 = eta1 * cos ( eta2 )
  end function polar_jac22

 ! jacobian ie determinant of jacobian matrix
  function polar_jac ( eta1, eta2, params )
    sll_real64  :: polar_jac
    sll_real64, intent(in)   :: eta1
    sll_real64, intent(in)   :: eta2
    sll_real64, dimension(:), intent(in) :: params
    polar_jac = eta1
  end function polar_jac


  ! **************************************************************************
  !
  ! "Colella transformation";
  ! sinusoidal product (see P. Colella et al. JCP 230 (2011) formula 
  ! (102) p 2968):
  !
  !        x1 = eta1 + alpha1 * sin(2*pi/L1 * eta1) * sin(2*pi/L2 * eta2)
  !        x2 = eta2 + alpha2 * sin(2*pi/L1 * eta1) * sin(2*pi/L2 * eta2)
  !
  ! Domain: [0,L1] X [0,L2]
  ! But we generalize it by taking a transformation 
  ! Domain: [a,b] X [c,d]  -----> Domain: [a',b'] X [c',d']
  ! so the transformation becomes 
  ! 
  !        x1 = (b' - a')/(b- a) * (  eta1 + alpha1 * sin(2*pi*( eta1-a)/( b-a)) * sin(2*pi*( eta2-c)/( d-c)) ) 
  !             + ( a' b - b' a)/(b- a)
  !        x2 = (d' - c')/(d- c) * (  eta2 + alpha2 * sin(2*pi*( eta1-a)/( b-a)) * sin(2*pi*( eta2-c)/( d-c)) ) 
  !             + ( c' d - d' c)/(d- c)
  ! 
  ! The parameters are:
  !      alpha1 = params(1)
  !      alpha2 = params(2)
  !      a      = params(3)
  !      b      = params(4)
  !      a'     = params(5)
  !      b'     = params(6)
  !      c      = params(7)
  !      d      = params(8)
  !      c'     = params(9)
  !      d'     = params(10)
  !
  !
  ! **************************************************************************

  ! direct mapping
  function sinprod_gen_x1 ( eta1, eta2, params )
    sll_real64  :: sinprod_gen_x1
    sll_real64, intent(in)   :: eta1
    sll_real64, intent(in)   :: eta2
    sll_real64, dimension(:), intent(in) :: params
    sll_real64  :: alpha1
    sll_real64  :: alpha2
    sll_real64  :: a_1 
    sll_real64  :: b_1
    sll_real64  :: a_2
    sll_real64  :: b_2
    sll_real64  :: c_1 
    sll_real64  :: d_1
    sll_real64  :: c_2
    sll_real64  :: d_2
    sll_real64  :: rl1 ! reciprocal of the length of the domain
    sll_real64  :: rl2
    sll_real64  :: pi2
    sll_real64  :: coef1, coef2

    SLL_ASSERT(size(params) >= 10)
    alpha1 = params(1)
    alpha2 = params(2)
    a_1    = params(3)
    b_1    = params(4)
    a_2    = params(5)
    b_2    = params(6)
    c_1    = params(7)
    d_1    = params(8)
    c_2    = params(9)
    d_2    = params(10)


    rl1    = ( eta1 - a_1 )/ (b_1 - a_1 ) 
    rl2    = ( eta2 - c_1 )/ (d_1 - c_1 )

    coef1 = ( b_2 - a_2 )/ (b_1 - a_1 )
    coef2 = ( a_2*b_1 - b_2*a_1 )/ (b_1 - a_1 )
 
    pi2 = 2.0_f64*sll_pi
    sinprod_gen_x1 = eta1 + alpha1 * sin(pi2*rl1)*sin(pi2*rl2)
    sinprod_gen_x1 = coef1 * sinprod_gen_x1 + coef2
  end function sinprod_gen_x1

  function sinprod_gen_x2 ( eta1, eta2, params )
    sll_real64  :: sinprod_gen_x2
    sll_real64, intent(in)   :: eta1
    sll_real64, intent(in)   :: eta2
    sll_real64, dimension(:), intent(in) :: params
    sll_real64  :: alpha1
    sll_real64  :: alpha2
    sll_real64  :: a_1 
    sll_real64  :: b_1
    sll_real64  :: a_2
    sll_real64  :: b_2
    sll_real64  :: c_1 
    sll_real64  :: d_1
    sll_real64  :: c_2
    sll_real64  :: d_2
    sll_real64  :: coef1, coef2
    sll_real64  :: rl1 ! reciprocal of the length of the domain
    sll_real64  :: rl2
    sll_real64  :: pi2
    

    SLL_ASSERT(size(params) >= 10)
    alpha1 = params(1)
    alpha2 = params(2)
    a_1    = params(3)
    b_1    = params(4)
    a_2    = params(5)
    b_2    = params(6)
    c_1    = params(7)
    d_1    = params(8)
    c_2    = params(9)
    d_2    = params(10)


    rl1    = ( eta1 - a_1 )/ (b_1 - a_1 ) 
    rl2    = ( eta2 - c_1 )/ (d_1 - c_1 )

    coef1 = ( d_2 - c_2 )/ (d_1 - c_1 )
    coef2 = ( c_2*d_1 - d_2*c_1 )/ (d_1 - c_1 )
 
    pi2 = 2.0_f64*sll_pi
    sinprod_gen_x2 = eta2 + alpha2 * sin(pi2*rl1)*sin(pi2*rl2)
    sinprod_gen_x2 = coef1 * sinprod_gen_x2 + coef2
  end function sinprod_gen_x2

  ! inverse mapping 
  ! cannot be computed analytically in this case. Use fixed point iterations.
  function sinprod_gen_eta1 ( x1, x2, params )
    sll_real64  :: sinprod_gen_eta1
    sll_real64, intent(in)   :: x1
    sll_real64, intent(in)   :: x2
    sll_real64, dimension(:), optional, intent(in) :: params
    ! NEEDS TO BE IMPLEMENTED
    STOP 'function not implemented'
    sinprod_gen_eta1 = x1
  end function sinprod_gen_eta1

  function sinprod_gen_eta2 ( x1, x2, params )
    sll_real64  :: sinprod_gen_eta2
    sll_real64, intent(in)   :: x1
    sll_real64, intent(in)   :: x2
    sll_real64, dimension(:), optional, intent(in) :: params
    ! NEEDS TO BE IMPLEMENTED
    STOP 'function not implemented'
    sinprod_gen_eta2 = x2
  end function sinprod_gen_eta2

  ! jacobian matrix
  function sinprod_gen_jac11 ( eta1, eta2, params )
    sll_real64  :: sinprod_gen_jac11
    sll_real64, intent(in)   :: eta1
    sll_real64, intent(in)   :: eta2
    sll_real64, dimension(:), intent(in) :: params
    sll_real64  :: alpha1
    sll_real64  :: alpha2
    sll_real64  :: a_1 
    sll_real64  :: b_1
    sll_real64  :: a_2
    sll_real64  :: b_2
    sll_real64  :: c_1 
    sll_real64  :: d_1
    sll_real64  :: c_2
    sll_real64  :: d_2
    sll_real64  :: coef1, coef2
    sll_real64  :: rl1 ! reciprocal of the length of the domain
    sll_real64  :: rl2
    sll_real64  :: pi2
    
    
    SLL_ASSERT(size(params) >= 10)
    alpha1 = params(1)
    alpha2 = params(2)
    a_1    = params(3)
    b_1    = params(4)
    a_2    = params(5)
    b_2    = params(6)
    c_1    = params(7)
    d_1    = params(8)
    c_2    = params(9)
    d_2    = params(10)
    
    rl1    = ( eta1 - a_1 )/ (b_1 - a_1 ) 
    rl2    = ( eta2 - c_1 )/ (d_1 - c_1 )

    coef1 = ( b_2 - a_2 )/ (b_1 - a_1 )
    coef2 = ( a_2*b_1 - b_2*a_1 )/ (b_1 - a_1 )
 
    pi2 = 2.0_f64*sll_pi

    sinprod_gen_jac11 = 1.0_f64 + alpha1*pi2/ (b_1 - a_1 ) *cos(pi2*rl1)*sin(pi2*rl2)
    sinprod_gen_jac11 = coef1 * sinprod_gen_jac11
  end function sinprod_gen_jac11

  function sinprod_gen_jac12 ( eta1, eta2, params )
    sll_real64  :: sinprod_gen_jac12
    sll_real64, intent(in)   :: eta1
    sll_real64, intent(in)   :: eta2
    sll_real64, dimension(:), intent(in) :: params
    sll_real64  :: alpha1
    sll_real64  :: alpha2
    sll_real64  :: a_1 
    sll_real64  :: b_1
    sll_real64  :: a_2
    sll_real64  :: b_2
    sll_real64  :: c_1 
    sll_real64  :: d_1
    sll_real64  :: c_2
    sll_real64  :: d_2
    sll_real64  :: coef1, coef2
    sll_real64  :: rl1 ! reciprocal of the length of the domain
    sll_real64  :: rl2
    sll_real64  :: pi2
    
    SLL_ASSERT(size(params) >= 10)
    alpha1 = params(1)
    alpha2 = params(2)
    a_1    = params(3)
    b_1    = params(4)
    a_2    = params(5)
    b_2    = params(6)
    c_1    = params(7)
    d_1    = params(8)
    c_2    = params(9)
    d_2    = params(10)
    
    rl1    = ( eta1 - a_1 )/ (b_1 - a_1 ) 
    rl2    = ( eta2 - c_1 )/ (d_1 - c_1 )

    coef1 = ( b_2 - a_2 )/ (b_1 - a_1 )
    coef2 = ( a_2*b_1 - b_2*a_1 )/ (b_1 - a_1 )
 
    pi2 = 2.0_f64*sll_pi

    sinprod_gen_jac12 = alpha1*pi2/ (d_1 - c_1 ) *sin(pi2*rl1)*cos(pi2*rl2)
    sinprod_gen_jac12 = coef1 * sinprod_gen_jac12
  end function sinprod_gen_jac12
  
  function sinprod_gen_jac21 ( eta1, eta2, params )
    sll_real64  :: sinprod_gen_jac21
    sll_real64, intent(in)   :: eta1
    sll_real64, intent(in)   :: eta2
    sll_real64, dimension(:), intent(in) :: params
    sll_real64  :: alpha1
    sll_real64  :: alpha2
    sll_real64  :: a_1 
    sll_real64  :: b_1
    sll_real64  :: a_2
    sll_real64  :: b_2
    sll_real64  :: c_1 
    sll_real64  :: d_1
    sll_real64  :: c_2
    sll_real64  :: d_2
    sll_real64  :: coef1, coef2
    sll_real64  :: rl1 ! reciprocal of the length of the domain
    sll_real64  :: rl2
    sll_real64  :: pi2
    
    SLL_ASSERT(size(params) >= 10)
    alpha1 = params(1)
    alpha2 = params(2)
    a_1    = params(3)
    b_1    = params(4)
    a_2    = params(5)
    b_2    = params(6)
    c_1    = params(7)
    d_1    = params(8)
    c_2    = params(9)
    d_2    = params(10)

    
    rl1    = ( eta1 - a_1 )/ (b_1 - a_1 ) 
    rl2    = ( eta2 - c_1 )/ (d_1 - c_1 )
    
    coef1 = ( d_2 - c_2 )/ (d_1 - c_1 )
    coef2 = ( c_2*d_1 - d_2*c_1 )/ (d_1 - c_1 )
 
    pi2 = 2.0_f64*sll_pi
    sinprod_gen_jac21 = alpha2*pi2/ (b_1 - a_1 ) *cos(pi2*rl1)*sin(pi2*rl2)
    sinprod_gen_jac21 = coef1 * sinprod_gen_jac21
  end function sinprod_gen_jac21

  function sinprod_gen_jac22 ( eta1, eta2, params )
    sll_real64  :: sinprod_gen_jac22
    sll_real64, intent(in)   :: eta1
    sll_real64, intent(in)   :: eta2
    sll_real64, dimension(:), intent(in) :: params
    sll_real64  :: alpha1
    sll_real64  :: alpha2
    sll_real64  :: a_1 
    sll_real64  :: b_1
    sll_real64  :: a_2
    sll_real64  :: b_2
    sll_real64  :: c_1 
    sll_real64  :: d_1
    sll_real64  :: c_2
    sll_real64  :: d_2
    sll_real64  :: coef1, coef2
    sll_real64  :: rl1 ! reciprocal of the length of the domain
    sll_real64  :: rl2
    sll_real64  :: pi2
    
    SLL_ASSERT(size(params) >= 10)
    alpha1 = params(1)
    alpha2 = params(2)
    a_1    = params(3)
    b_1    = params(4)
    a_2    = params(5)
    b_2    = params(6)
    c_1    = params(7)
    d_1    = params(8)
    c_2    = params(9)
    d_2    = params(10)

    
    rl1    = ( eta1 - a_1 )/ (b_1 - a_1 ) 
    rl2    = ( eta2 - c_1 )/ (d_1 - c_1 )
    
    coef1 = ( d_2 - c_2 )/ (d_1 - c_1 )
    coef2 = ( c_2*d_1 - d_2*c_1 )/ (d_1 - c_1 )
 
    pi2 = 2.0_f64*sll_pi

    sinprod_gen_jac22 = 1.0_f64 + alpha2*pi2/ (d_1 - c_1 )*sin(pi2*rl1)*cos(pi2*rl2)
    sinprod_gen_jac22 = coef1 * sinprod_gen_jac22
    
  end function sinprod_gen_jac22

   ! jacobian ie determinant of jacobian matrix
  function sinprod_gen_jac ( eta1, eta2, params )
    sll_real64  :: sinprod_gen_jac
    sll_real64, intent(in)   :: eta1
    sll_real64, intent(in)   :: eta2
    sll_real64, dimension(:), intent(in) :: params
    sll_real64  :: alpha1
    sll_real64  :: alpha2
    sll_real64  :: rl1 ! reciprocal of the length of the domain
    sll_real64  :: rl2
    sll_real64  :: pi2


    SLL_ASSERT(size(params) >= 10)
    !alpha1 = params(1)
    !alpha2 = params(2)
    !rl1    = 1.0_f64/params(3)
    !rl2    = 1.0_f64/params(4)
    !pi2 = 2.0_f64*sll_pi
    !sinprod_jac = 1.0_f64 + 0.2_f64 *sll_pi * sin (2*sll_pi**(eta1+eta2)) 
    !sinprod_gen_jac = 1.0_f64 + alpha2*pi2*rl2*sin(pi2*rl1*eta1)*cos(pi2*rl2*eta2) + &
    !                        alpha1*pi2*rl1*cos(pi2*rl1*eta1)*sin(pi2*rl2*eta2)

    sinprod_gen_jac=sinprod_gen_jac22( eta1, eta2, params )*sinprod_gen_jac11( eta1, eta2, params )&
                  - sinprod_gen_jac12( eta1, eta2, params )*sinprod_gen_jac21( eta1, eta2, params )
  end function sinprod_gen_jac

  ! **************************************************************************
  !
  ! "Colella transformation";
  ! sinusoidal product (see P. Colella et al. JCP 230 (2011) formula 
  ! (102) p 2968):
  !
  !        x1 = eta1 + alpha1 * sin(2*pi*eta1) * sin(2*pi*eta2)
  !        x2 = eta2 + alpha2 * sin(2*pi*eta1) * sin(2*pi*eta2)
  !
  ! Domain: [0,L1] X [0,L2]
  ! By default the values of the alpha parameters are:
  !      alpha1 = 0.1
  !      alpha2 = 0.1
  !      L1     = 1.0
  !      L2     = 1.0
  !
  ! These parameters are stored in the params array as: 
  !     ( alpha1, alpha2, L1, L2 )
  !
  ! **************************************************************************

  ! direct mapping
  function sinprod_x1 ( eta1, eta2, params )
    sll_real64  :: sinprod_x1
    sll_real64, intent(in)   :: eta1
    sll_real64, intent(in)   :: eta2
    sll_real64, dimension(:), intent(in) :: params
    sll_real64  :: alpha1
    sll_real64  :: rl1 ! reciprocal of the length of the domain
    sll_real64  :: rl2
    sll_real64  :: pi2

    SLL_ASSERT(size(params) >= 4)
    alpha1 = params(1)
    rl1    = 1.0_f64/params(3)
    rl2    = 1.0_f64/params(4)
    pi2 = 2.0_f64*sll_pi
    sinprod_x1 = eta1 + alpha1 * sin(pi2*rl1*eta1)*sin(pi2*rl2*eta2)
  end function sinprod_x1

  function sinprod_x2 ( eta1, eta2, params )
    sll_real64  :: sinprod_x2
    sll_real64, intent(in)   :: eta1
    sll_real64, intent(in)   :: eta2
    sll_real64, dimension(:),intent(in) :: params
    sll_real64  :: alpha2
    sll_real64  :: rl1 ! reciprocal of the length of the domain
    sll_real64  :: rl2
    sll_real64  :: pi2

    SLL_ASSERT(size(params) >= 4)
    alpha2 = params(2)
    rl1    = 1.0_f64/params(3)
    rl2    = 1.0_f64/params(4)
    pi2 = 2.0_f64*sll_pi
    sinprod_x2 = eta2 + alpha2*sin(pi2*rl1*eta1)*sin(pi2*rl2*eta2)
  end function sinprod_x2

  ! inverse mapping 
  ! cannot be computed analytically in this case. Use fixed point iterations.
  function sinprod_eta1 ( x1, x2, params )
    sll_real64  :: sinprod_eta1
    sll_real64, intent(in)   :: x1
    sll_real64, intent(in)   :: x2
    sll_real64, dimension(:), optional, intent(in) :: params
    ! NEEDS TO BE IMPLEMENTED
    STOP 'function not implemented'
    sinprod_eta1 = x1
  end function sinprod_eta1

  function sinprod_eta2 ( x1, x2, params )
    sll_real64  :: sinprod_eta2
    sll_real64, intent(in)   :: x1
    sll_real64, intent(in)   :: x2
    sll_real64, dimension(:), optional, intent(in) :: params
    ! NEEDS TO BE IMPLEMENTED
    STOP 'function not implemented'
    sinprod_eta2 = x2
  end function sinprod_eta2

  ! jacobian matrix
  function sinprod_jac11 ( eta1, eta2, params )
    sll_real64  :: sinprod_jac11
    sll_real64, intent(in)   :: eta1
    sll_real64, intent(in)   :: eta2
    sll_real64, dimension(:), intent(in) :: params
    sll_real64  :: alpha1
    sll_real64  :: rl1 ! reciprocal of the length of the domain
    sll_real64  :: rl2
    sll_real64  :: pi2

    SLL_ASSERT(size(params) >= 4)
    alpha1 = params(1)
    rl1    = 1.0_f64/params(3)
    rl2    = 1.0_f64/params(4)
    pi2 = 2.0_f64*sll_pi
    sinprod_jac11 = 1.0_f64 + alpha1*pi2*rl1*cos(pi2*rl1*eta1)*sin(pi2*rl2*eta2)
  end function sinprod_jac11

  function sinprod_jac12 ( eta1, eta2, params )
    sll_real64  :: sinprod_jac12
    sll_real64, intent(in)   :: eta1
    sll_real64, intent(in)   :: eta2
    sll_real64, dimension(:), intent(in) :: params
    sll_real64  :: alpha1
    sll_real64  :: rl1 ! reciprocal of the length of the domain
    sll_real64  :: rl2
    sll_real64  :: pi2

    SLL_ASSERT(size(params) >= 4)
    alpha1 = params(1)
    rl1    = 1.0_f64/params(3)
    rl2    = 1.0_f64/params(4)
    pi2 = 2.0_f64*sll_pi
    sinprod_jac12 = alpha1*pi2*rl2*sin(pi2*rl1*eta1)*cos(pi2*rl2*eta2)
  end function sinprod_jac12

  function sinprod_jac21 ( eta1, eta2, params )
    sll_real64  :: sinprod_jac21
    sll_real64, intent(in)   :: eta1
    sll_real64, intent(in)   :: eta2
    sll_real64, dimension(:), intent(in) :: params
    sll_real64  :: alpha2
    sll_real64  :: rl1 ! reciprocal of the length of the domain
    sll_real64  :: rl2
    sll_real64  :: pi2
    sll_real64  :: l1,l2

    SLL_ASSERT(size(params) >= 4)
    alpha2 = params(2)
    rl1    = 1.0_f64/params(3)
    rl2    = 1.0_f64/params(4)
    pi2 = 2.0_f64*sll_pi
    sinprod_jac21 = alpha2*pi2*rl1*cos(pi2*rl1*eta1)*sin(pi2*rl2*eta2)
  end function sinprod_jac21

  function sinprod_jac22 ( eta1, eta2, params )
    sll_real64  :: sinprod_jac22
    sll_real64, intent(in)   :: eta1
    sll_real64, intent(in)   :: eta2
    sll_real64, dimension(:), intent(in) :: params
    sll_real64  :: alpha2
    sll_real64  :: rl1 ! reciprocal of the length of the domain
    sll_real64  :: rl2
    sll_real64  :: pi2
    sll_real64  :: l1,l2

    SLL_ASSERT(size(params) >= 4)
    alpha2 = params(2)
    rl1    = 1.0_f64/params(3)
    rl2    = 1.0_f64/params(4)
    pi2 = 2.0_f64*sll_pi
    sinprod_jac22 = 1.0_f64 + alpha2*pi2*rl2*sin(pi2*rl1*eta1)*cos(pi2*rl2*eta2)
  end function sinprod_jac22

   ! jacobian ie determinant of jacobian matrix
  function sinprod_jac ( eta1, eta2, params )
    sll_real64  :: sinprod_jac
    sll_real64, intent(in)   :: eta1
    sll_real64, intent(in)   :: eta2
    sll_real64, dimension(:), intent(in) :: params
    sll_real64  :: alpha1
    sll_real64  :: alpha2
    sll_real64  :: rl1 ! reciprocal of the length of the domain
    sll_real64  :: rl2
    sll_real64  :: pi2


    SLL_ASSERT(size(params) >= 4)
    alpha1 = params(1)
    alpha2 = params(2)
    rl1    = 1.0_f64/params(3)
    rl2    = 1.0_f64/params(4)
    pi2 = 2.0_f64*sll_pi
    !sinprod_jac = 1.0_f64 + 0.2_f64 *sll_pi * sin (2*sll_pi**(eta1+eta2)) 
    sinprod_jac = 1.0_f64 + alpha2*pi2*rl2*sin(pi2*rl1*eta1)*cos(pi2*rl2*eta2) + &
                            alpha1*pi2*rl1*cos(pi2*rl1*eta1)*sin(pi2*rl2*eta2)
  end function sinprod_jac

#if 0
! Only one Colella transformation should survive, the parametrized one above...
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! direct mapping collela on (0,4*pi)x (0,1)
  ! Same story as above, we separate the alpha coefficients in each direction
  ! and give default values of 0.1. These and the previous transformation
  ! should be merged and parametrized with params.
  function sinprod_x1_rect ( eta1, eta2, params )
    real(8)  :: sinprod_x1_rect
    real(8), intent(in)   :: eta1
    real(8), intent(in)   :: eta2
    sll_real64, dimension(:), optional, intent(in) :: params
    sll_real64  :: alpha1
    sll_real64  :: pi2
    if(present(params)) then
       SLL_ASSERT(size(params) >= 2)
       alpha1 = params(1)
    else
       alpha1 = 0.1_f64
    end if 
    pi2 = 2.0_f64*sll_pi
    sinprod_x1_rect = eta1 + alpha1*sin(0.5*eta1)*sin(pi2*eta2)
  end function sinprod_x1_rect
  
  function sinprod_x2_rect ( eta1, eta2, params )
    real(8)  :: sinprod_x2_rect
    real(8), intent(in)   :: eta1
    real(8), intent(in)   :: eta2
    sll_real64, dimension(:), optional, intent(in) :: params
    sll_real64  :: alpha2
    sll_real64  :: pi2
    if(present(params)) then
       SLL_ASSERT(size(params) >= 2)
       alpha2 = params(2)
    else
       alpha2 = 0.1_f64
    end if
    pi2 = 2.0_f64*sll_pi
    sinprod_x2_rect = eta2 + alpha2*sin(0.5*eta1)*sin(pi2*eta2)
  end function sinprod_x2_rect
  
  
  ! jacobian matrix
  function sinprod_jac11_rect ( eta1, eta2, params )
    real(8)  :: sinprod_jac11_rect
    real(8), intent(in)   :: eta1
    real(8), intent(in)   :: eta2
    sll_real64, dimension(:), optional, intent(in) :: params
    if(present(params)) then
       SLL_ASSERT(size(params) >= 2)
       alpha1 = params(1)
    else
       alpha1 = 0.1_f64
    end if 
    sinprod_jac11_rect = 1.0_8 + 0.5_8*alpha1 * cos (0.5*eta1) * sin (2*sll_pi*eta2)
  end function sinprod_jac11_rect
  
  function sinprod_jac12_rect ( eta1, eta2, params )
    real(8)  :: sinprod_jac12_rect
    real(8), intent(in)   :: eta1
    real(8), intent(in)   :: eta2
    sll_real64, dimension(:), optional, intent(in) :: params
    if(present(params)) then
       SLL_ASSERT(size(params) >= 2)
       alpha1 = params(1)
    else
       alpha1 = 0.1_f64
    end if
    sinprod_jac12_rect = 2._8 *sll_pi*alpha1* sin (0.5*eta1) * cos (2*sll_pi*eta2)
  end function sinprod_jac12_rect
  
  function sinprod_jac21_rect ( eta1, eta2, params )
    real(8)  :: sinprod_jac21_rect
    real(8), intent(in)   :: eta1
    real(8), intent(in)   :: eta2
    sll_real64, dimension(:), optional, intent(in) :: params
    if(present(params)) then
       SLL_ASSERT(size(params) >= 2)
       alpha2 = params(2)
    else
       alpha2 = 0.1_f64
    end if
    sinprod_jac21_rect = 0.5_8*alpha2* cos (0.5*eta1) * sin (2*sll_pi*eta2)
  end function sinprod_jac21_rect

  function sinprod_jac22_rect ( eta1, eta2, params )
    real(8)  :: sinprod_jac22_rect
    real(8), intent(in)   :: eta1
    real(8), intent(in)   :: eta2
    sll_real64, dimension(:), optional, intent(in) :: params
    if(present(params)) then
       SLL_ASSERT(size(params) >= 2)
       alpha2 = params(2)
    else
       alpha2 = 0.1_f64
    end if
    sinprod_jac22_rect = 1.0_8 + &
         2.0_8*sll_pi*alpha2*sin(0.5*eta1)*cos(2*sll_pi*eta2)
  end function sinprod_jac22_rect
  
  ! jacobian ie determinant of jacobian matrix
  function sinprod_jac_rect ( eta1, eta2, params )
    real(8)  :: sinprod_jac_rect
    real(8), intent(in)   :: eta1
    real(8), intent(in)   :: eta2
    sll_real64, dimension(:), optional, intent(in) :: params
    !sinprod_jac = 1.0_f64 + 0.2_f64 *sll_pi * sin (2*sll_pi**(eta1+eta2)) 
    if(present(params)) then
       SLL_ASSERT(size(params) >= 2)
       alpha2 = params(2)
       alpha1 = params(1)
    else
       alpha2 = 0.1_f64
       alpha1 = 0.1_f64
    end if
    sinprod_jac_rect = &
         (1.0_8 + 0.5_8*alpha1 * cos (0.5*eta1) * sin (2*sll_pi*eta2))* &
         (1.0_8 + 2.0_8* sll_pi*alpha2 * sin(0.5*eta1) * cos(2*sll_pi*eta2) ) - &
         2.0_8 *sll_pi*alpha2 * sin (0.5*eta1) * cos (2*sll_pi*eta2) * &
         0.5_8*alpha1 * cos (0.5*eta1) * sin (2*sll_pi*eta2)
  end function sinprod_jac_rect
  


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! direct mapping collela on (0,4*pi)x (0,4*pi)
  function sinprod_x1_square ( eta1, eta2, params )
    real(8)  :: sinprod_x1_square
    real(8), intent(in)   :: eta1
    real(8), intent(in)   :: eta2
    sll_real64, dimension(:), optional, intent(in) :: params
    sinprod_x1_square = eta1 + 0.1_8 * sin(0.5*eta1) * sin(0.5*eta2)
  end function sinprod_x1_square
  
  function sinprod_x2_square ( eta1, eta2, params )
    real(8)  :: sinprod_x2_square
    real(8), intent(in)   :: eta1
    real(8), intent(in)   :: eta2
    sll_real64, dimension(:), optional, intent(in) :: params
    sinprod_x2_square = eta2 + 0.1_8 * sin(0.5*eta1) * sin(0.5*eta2)
  end function sinprod_x2_square
  
  
  ! jacobian matrix
  function sinprod_jac11_square ( eta1, eta2, params )
    real(8)  :: sinprod_jac11_square
    real(8), intent(in)   :: eta1
    real(8), intent(in)   :: eta2
    sll_real64, dimension(:), optional, intent(in) :: params
    sinprod_jac11_square = 1.0_8 + &
         0.5_8*0.1_8 * cos (0.5*eta1) * sin (0.5*eta2)
  end function sinprod_jac11_square
  
  function sinprod_jac12_square ( eta1, eta2, params )
    real(8)  :: sinprod_jac12_square
    real(8), intent(in)   :: eta1
    real(8), intent(in)   :: eta2
    sll_real64, dimension(:), optional, intent(in) :: params
    sinprod_jac12_square = 0.5_8*0.1_8 * sin (0.5*eta1) * cos (0.5*eta2)
  end function sinprod_jac12_square
  
  function sinprod_jac21_square ( eta1, eta2, params )
    real(8)  :: sinprod_jac21_square
    real(8), intent(in)   :: eta1
    real(8), intent(in)   :: eta2
    sll_real64, dimension(:), optional, intent(in) :: params
    sinprod_jac21_square = 0.5_8*0.1_8 * cos (0.5*eta1) * sin (0.5*eta2)
  end function sinprod_jac21_square
  
  function sinprod_jac22_square ( eta1, eta2, params )
    real(8)  :: sinprod_jac22_square
    real(8), intent(in)   :: eta1
    real(8), intent(in)   :: eta2
    sll_real64, dimension(:), optional, intent(in) :: params
    sinprod_jac22_square = 1.0_8 + &
         0.5_8*0.1_8*sin(0.5*eta1)*cos(0.5*eta2)
  end function sinprod_jac22_square
  
  ! jacobian ie determinant of jacobian matrix
  function sinprod_jac_square ( eta1, eta2, params )
    real(8)  :: sinprod_jac_square
    real(8), intent(in)   :: eta1
    real(8), intent(in)   :: eta2
    sll_real64, dimension(:), optional, intent(in) :: params
    !sinprod_jac = 1.0_f64 + 0.2_f64 *sll_pi * sin (2*sll_pi**(eta1+eta2)) 
    sinprod_jac_square = &
         (1.0_8 + 0.5_8*0.1_8 * cos (0.5*eta1) * sin (0.5*eta2))* &
         (1.0_8 + 0.5_8*0.1_8 * sin(0.5*eta1) * cos(0.5*eta2) ) - &
         0.5_8*0.1_8 * sin (0.5*eta1) * cos (0.5*eta2) * &
         0.5_8*0.1_8 * cos (0.5*eta1) * sin (0.5*eta2)
  end function sinprod_jac_square

#endif

#if 0
  ! What is this???
  ! test function
  !-------------------
  ! direct mapping
  function test_x1 ( eta1, eta2, params )
    sll_real64  :: test_x1
    sll_real64, intent(in)   :: eta1
    sll_real64, intent(in)   :: eta2
    sll_real64, dimension(:), optional, intent(in) :: params
    test_x1 = eta1 + 0.1_f64 * sin( 2.0_f64* sll_pi * eta1 )
    !test_x1 = eta1**2
  end function test_x1

  function test_x2 ( eta1, eta2, params )
    sll_real64  :: test_x2
    sll_real64, intent(in)   :: eta1
    sll_real64, intent(in)   :: eta2
    sll_real64, dimension(:), optional, intent(in) :: params
    test_x2 = eta2 + 0.1_f64 * sin( 2.0_f64* sll_pi * eta2 )
  end function test_x2

  ! inverse mapping
  function test_eta1 ( x1, x2, params )
    sll_real64  :: test_eta1
    sll_real64, intent(in)   :: x1
    sll_real64, intent(in)   :: x2
    sll_real64, dimension(:), optional, intent(in) :: params
    test_eta1 = x1 / 0.1_f64
  end function test_eta1

  function test_eta2 ( x1, x2, params )
    sll_real64  :: test_eta2
    sll_real64, intent(in)   :: x1
    sll_real64, intent(in)   :: x2
    sll_real64, dimension(:), optional, intent(in) :: params
    test_eta2 = x2 / 0.1_f64
  end function test_eta2

  ! inverse jacobian matrix
  function test_jac11 ( eta1, eta2, params )
    sll_real64  :: test_jac11
    sll_real64, intent(in)   :: eta1
    sll_real64, intent(in)   :: eta2
    sll_real64, dimension(:), optional, intent(in) :: params
    test_jac11 = 1.0_f64 / (1.0_f64 + 2.0_f64 * sll_pi* 0.1_f64 * cos( 2.0_f64* sll_pi * eta1))
  end function test_jac11

    function test_jac12 ( eta1, eta2, params )
    sll_real64  :: test_jac12
    sll_real64, intent(in)   :: eta1
    sll_real64, intent(in)   :: eta2
    sll_real64, dimension(:), optional, intent(in) :: params
    test_jac12 = 0.0_f64
  end function test_jac12

  function test_jac21 ( eta1, eta2, params )
    sll_real64  :: test_jac21
    sll_real64, intent(in)   :: eta1
    sll_real64, intent(in)   :: eta2
    sll_real64, dimension(:), optional, intent(in) :: params
    test_jac21 = 0.0_f64
  end function test_jac21

  function test_jac22 ( eta1, eta2, params )
    sll_real64  :: test_jac22
    sll_real64, intent(in)   :: eta1
    sll_real64, intent(in)   :: eta2
    sll_real64, dimension(:), optional, intent(in) :: params
    test_jac22 = 1.0_f64 / (1.0_f64 + 2.0_f64 * sll_pi* 0.1_f64 * cos( 2.0_f64* sll_pi * eta2))
  end function test_jac22

  ! jacobian ie determinant of jacobian matrix
  function test_jac ( eta1, eta2, params )
    sll_real64  :: test_jac
    sll_real64, intent(in)   :: eta1
    sll_real64, intent(in)   :: eta2
    sll_real64, dimension(:), optional, intent(in) :: params
    test_jac =  (1.0_f64 + 2.0_f64 * sll_pi* 0.1_f64 * cos( 2.0_f64* sll_pi * eta1)) * &
         (1.0_f64 + 2.0_f64 * sll_pi* 0.1_f64 * cos( 2.0_f64* sll_pi * eta2))
    !test_jac =  2 * eta1!
  end function test_jac
#endif

  ! ***************************************************************************
  !
  ! Alternative formulation for the polar coordinate transformation:
  !
  ! X1 = (Rmin + (Rmax-Rmin)*eta1)*cos(2*pi*eta2)
  ! X2 = (Rmin + (Rmax-Rmin)*eta1)*sin(2*pi*eta2)
  !
  ! Where eta1 and eta2 are defined in the interval [0,1]. The 'params' array
  ! contains the information (R1, R2). Typically:
  !
  ! R1 = 0.1
  ! R2 = 1.0
  !
  ! ***************************************************************************


  function x1_polar_f( eta1, eta2, params )
    sll_real64 :: x1_polar_f
    sll_real64, intent(in) :: eta1, eta2
    sll_real64, dimension(:), intent(in) :: params
    sll_real64 :: r1
    sll_real64 :: r2

    SLL_ASSERT(size(params) >= 2)
    r1 = params(1)
    r2 = params(2)
    x1_polar_f = (r1 + (r2-r1)*eta1)*cos(2.0_f64*sll_pi*eta2)
  end function x1_polar_f

  function x2_polar_f( eta1, eta2, params )
    sll_real64 :: x2_polar_f
    sll_real64, intent(in) :: eta1, eta2
    sll_real64, dimension(:), intent(in) :: params
    sll_real64 :: r1
    sll_real64 :: r2

    SLL_ASSERT(size(params) >= 2)
    r1 = params(1)
    r2 = params(2)
    x2_polar_f = (r1 + (r2-r1)*eta1)*sin(2.0_f64*sll_pi*eta2)
  end function x2_polar_f

  function deriv_x1_polar_f_eta1( eta1, eta2, params )
    sll_real64 :: deriv_x1_polar_f_eta1
    sll_real64, intent(in) :: eta1, eta2
    sll_real64, dimension(:), intent(in) :: params
    sll_real64 :: r1
    sll_real64 :: r2

    SLL_ASSERT(size(params) >= 2)
    r1 = params(1)
    r2 = params(2)
    deriv_x1_polar_f_eta1 = (r2-r1)*cos(2.0_f64*sll_pi*eta2)
  end function deriv_x1_polar_f_eta1

  function deriv_x1_polar_f_eta2( eta1, eta2, params )
    sll_real64 :: deriv_x1_polar_f_eta2
    sll_real64, intent(in) :: eta1, eta2
    sll_real64, dimension(:), intent(in) :: params
    sll_real64 :: k
    sll_real64 :: r1
    sll_real64 :: r2

    SLL_ASSERT(size(params) >= 2)
    r1 = params(1)
    r2 = params(2)
    k = 2.0_f64*sll_pi
    deriv_x1_polar_f_eta2 = -(r1+(r2-r1)*eta1)*sin(k*eta2)*k
  end function deriv_x1_polar_f_eta2

  function deriv_x2_polar_f_eta1( eta1, eta2, params )
    sll_real64 :: deriv_x2_polar_f_eta1
    sll_real64, intent(in) :: eta1, eta2
    sll_real64, dimension(:), intent(in) :: params
    sll_real64 :: r1
    sll_real64 :: r2

    SLL_ASSERT(size(params) >= 2)
    r1 = params(1)
    r2 = params(2)
    deriv_x2_polar_f_eta1 = (r2-r1)*sin(2.0_f64*sll_pi*eta2)
  end function deriv_x2_polar_f_eta1

  function deriv_x2_polar_f_eta2( eta1, eta2, params )
    sll_real64 :: deriv_x2_polar_f_eta2
    sll_real64, intent(in) :: eta1, eta2
    sll_real64, dimension(:), intent(in) :: params
    sll_real64 :: k
    sll_real64 :: r1
    sll_real64 :: r2

    SLL_ASSERT(size(params) >= 2)
    r1 = params(1)
    r2 = params(2)
    k = 2.0_f64*sll_pi
    deriv_x2_polar_f_eta2 = (r1+(r2-r1)*eta1)*cos(k*eta2)*k
  end function deriv_x2_polar_f_eta2

  function jacobian_polar_f( eta1, eta2, params ) result(jac)
    sll_real64             :: jac
    sll_real64, intent(in) :: eta1, eta2
    sll_real64, dimension(:), intent(in) :: params
    sll_real64 :: r1
    sll_real64 :: r2

    SLL_ASSERT(size(params) >= 2)
    r1 = params(1)
    r2 = params(2)
    jac = 2.0_f64*sll_pi*(r1+(r2-r1)*eta1)*(r2-r1)
  end function jacobian_polar_f

  ! what is the following used for? It is not meant or used for the 
  ! coordinate transformation class... this one is used in the unit_test_2d.F90
  ! file but, what else?
  function deriv1_jacobian_polar_f(eta1, eta2, params ) result(deriv)
    sll_real64             :: deriv
    sll_real64, intent(in) :: eta1, eta2
    sll_real64, dimension(:), intent(in) :: params
    sll_real64 :: r1
    sll_real64 :: r2

    SLL_ASSERT(size(params) >= 2)
    r1 = params(1)
    r2 = params(2)
    deriv = 2.0_f64*sll_pi*(r2-r1)**2
  end function deriv1_jacobian_polar_f

! why is this here?
!!$  function zero_function(eta1, eta2)
!!$    sll_real64, intent(in) :: eta1, eta2
!!$    sll_real64             :: zero_function
!!$    zero_function = 0.0_f64
!!$  end function zero_function
    

  !************************************************************************
  !
  ! 1D maps
  !
  !************************************************************************

  ! Linear map
  !
  ! x1 = A + (B-A)*eta1
  !
  ! The params array is organized as (A,B) with default values:
  !
  ! A = -1.0
  ! B =  1.0

  function linear_map_f( eta, params ) result(val)
    sll_real64 :: val
    sll_real64, intent(in) :: eta
    sll_real64, dimension(:), intent(in) :: params
    sll_real64 :: a
    sll_real64 :: b

    SLL_ASSERT(size(params) >= 2)
    a = params(1)
    b = params(2)
    val = (b-a)*eta + a
  end function linear_map_f

  function linear_map_jac_f( eta, params ) result(val)
    sll_real64 :: val
    sll_real64, intent(in) :: eta
    sll_real64, dimension(:), intent(in) :: params
    sll_real64 :: a
    sll_real64 :: b

    SLL_ASSERT(size(params) >= 2)
    a = params(1)
    b = params(2)
        val = (b-a)
  end function linear_map_jac_f


!!$#define A 0.0_f64
!!$#define B 6.2831853071795862_f64
!!$
!!$  function linear_map_poisson_f( eta, params ) result(val)
!!$    sll_real64 :: val
!!$    sll_real64, intent(in) :: eta
!!$    val = (B-A)*eta + A
!!$  end function linear_map_poisson_f
!!$
!!$  function linear_map_poisson_jac_f( eta, params ) result(val)
!!$    sll_real64 :: val
!!$    sll_real64, intent(in) :: eta
!!$    val = (B-A)
!!$  end function linear_map_poisson_jac_f
!!$
!!$#undef A
!!$#undef B
  

end module sll_common_coordinate_transformations

