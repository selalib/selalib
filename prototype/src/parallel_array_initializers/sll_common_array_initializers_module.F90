module sll_common_array_initializers_module

#include "sll_working_precision.h"
#include "sll_constants.h"

  implicit none

  ! The functions specified here are meant to have the specific signature
  ! described in the sll_parallel_array_initializer_module. Else, they could
  ! not be used with this module.

contains

  ! -------------------------------------------------------------------------
  !
  !             Landau damping 4d initialization function
  !
  ! -------------------------------------------------------------------------


  ! This is a simplistic initializer aimed at a 4d cartesian distribution
  ! function, periodic in x and y, and compact-ish in vx and vy.
  !
  ! Basically:
  !                 1                           -(vx^2 + vy^2)
  ! f(x,y,vx,vy) = ----(1+epsilon*cos(kx*x))*exp(------------- )
  !                 2*pi                              2
  !
  ! It is meant to be used in the intervals:
  ! x:  [ 0,2*pi/kx]
  ! y:  [ 0,2*pi/ky]
  ! vx: [-6,6]
  ! vy: [-6,6]

  ! convention for the params array:
  ! params(1) = eta1_min
  ! params(2) = eta1_max
  ! params(3) = eta2_min
  ! params(4) = eta2_max
  ! params(5) = epsilon
  !
  ! The params array is declared optional to conform with the expected 
  ! function signature of the initializer subroutines, but in the particular
  ! case of the landau initializer, the params array must be passed.

  function sll_gaussian_initializer_4d( x, y, vx, vy, params ) 
    sll_real64 :: sll_gaussian_initializer_4d
    sll_real64, intent(in) :: x
    sll_real64, intent(in) :: y
    sll_real64, intent(in) :: vx
    sll_real64, intent(in) :: vy

    sll_real64, dimension(:), intent(in), optional :: params
    sll_real64 :: xc
    sll_real64 :: yc
    sll_real64 :: vxc
    sll_real64 :: vyc

    sll_real64 :: alpha = 1.
    sll_real64 :: beta  = 0.0

    if( .not. present(params) ) then
       print *, 'sll_gaussian_initializer_4d, error: the params array must ', &
            'be passed. params(1) = xc, params(2) = yc, params(3) = vxc...'
       stop
    end if

    xc  = params(1)
    yc  = params(2)
    vxc = params(3)
    vyc = params(4)

    sll_gaussian_initializer_4d = alpha*exp(-0.5_f64*((x-xc)**2+(y-yc)**2)) + &
                                  beta *exp(-0.5_f64*((vx-vxc)**2+(vy-vyc)**2))

  end function sll_gaussian_initializer_4d

  function sll_landau_initializer_4d( x, y, vx, vy, params ) 
    sll_real64 :: sll_landau_initializer_4d
    sll_real64, intent(in) :: x
    sll_real64, intent(in) :: y
    sll_real64, intent(in) :: vx
    sll_real64, intent(in) :: vy

    sll_real64, dimension(:), intent(in), optional :: params
    sll_real64 :: eta1_min
    sll_real64 :: eta1_max
    sll_real64 :: eta2_min
    sll_real64 :: eta2_max

    sll_real64 :: eps
    sll_real64 :: kx
    sll_real64 :: factor1

    if( .not. present(params) ) then
       print *, 'sll_landau_initializer_4d, error: the params array must ', &
            'be passed. params(1) = epsilon, params(2) = kx, params(3) = ky.'
       stop
    end if

    eta1_min = params(1)
    eta1_max = params(2)
    eta2_min = params(3)
    eta2_max = params(4)

    eps = params(5)
    kx  = 2. * sll_pi / (eta1_max - eta1_min)

    !Normalization
    !sagemath command
    !sage : var('u v epsilon a b c d x y')
    !sage : f(a,b,c,d,epsilon) =integral(integral(integral(integral((1+epsilon*cos(2*pi/(b-a)*x))*exp(-(u*u+v*v)/2),u,-oo,oo),v,-oo,oo),x,a,b),y,c,d)
    
    factor1 =  1./( (eta2_min - eta2_max) &
               *(((eta1_min - eta1_max)* &
               sin(2*sll_pi*eta1_min/(eta1_min - eta1_max)) &
                - (eta1_min - eta1_max)* &
               sin(2*sll_pi*eta1_max/(eta1_min - eta1_max)))*eps  &
               + 2*sll_pi*eta1_min - 2*sll_pi*eta1_max))
    
    sll_landau_initializer_4d = factor1 * &
         (1.0_f64+eps*cos(kx*x))*exp(-0.5_f64*(vx**2+vy**2))

  end function sll_landau_initializer_4d

  ! this function is a 1D landau initializer used for debugging
  ! 4D drift kinetic simulations in variables x1,x2,x3 ,v1
  ! the function is constant with respect to x2 and x3

  function sll_landau_initializer_dk_test_4d(v1,x1,x2,x3,params ) 
    sll_real64 :: sll_landau_initializer_dk_test_4d
    sll_real64, intent(in) :: x1
    sll_real64, intent(in) :: x2
    sll_real64, intent(in) :: x3
    sll_real64, intent(in) :: v1
    sll_real64, dimension(:), intent(in), optional :: params

    sll_real64 :: epsilon
    sll_real64 :: kx
    sll_real64 :: factor1

    if( .not. present(params) ) then
       print *, 'sll_landau_initializer_dk_test_4d, error: the params array must ', &
            'be passed. params(1) = epsilon, params(2) = kx'
       stop
    end if

    epsilon = params(1)
    kx      = params(2)
    factor1 = 0.5_f64/sll_pi

    !write(*,*) 'x1 v1',x1,v1

    sll_landau_initializer_dk_test_4d = factor1*&
         (1.0_f64+epsilon*cos(kx*x1))*exp(-0.5_f64*(v1**2))
  end function sll_landau_initializer_dk_test_4d

end module sll_common_array_initializers_module
