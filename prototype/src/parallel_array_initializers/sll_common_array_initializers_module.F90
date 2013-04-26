module sll_common_array_initializers_module
#include "sll_working_precision.h"
  use sll_constants
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
  !                 1                                      -(vx^2 + vy^2)
  ! f(x,y,vx,vy) = ----(1+epsilon*cos(kx*x)*cos(ky*y))*exp(------------- )
  !                 2*pi                                         2
  !
  ! It is meant to be used in the intervals:
  ! x:  [ 0,2*pi/kx]
  ! y:  [ 0,2*pi/ky]
  ! vx: [-6,6]
  ! vy: [-6,6]

  ! convention for the params array:
  ! params(1) = epsilon
  ! params(2) = kx
  ! params(3) = ky
  ! The params array is declared optional to conform with the expected 
  ! function signature of the initializer subroutines, but in the particular
  ! case of the landau initializer, the params array must be passed.

  function sll_landau_initializer_4d( x, y, vx, vy, params ) 
    sll_real64 :: sll_landau_initializer_4d
    sll_real64, intent(in) :: x
    sll_real64, intent(in) :: y
    sll_real64, intent(in) :: vx
    sll_real64, intent(in) :: vy
    sll_real64, dimension(:), intent(in), optional :: params

    sll_real64 :: epsilon
    sll_real64 :: kx
    sll_real64 :: ky
    sll_real64 :: factor1

    if( .not. present(params) ) then
       print *, 'sll_landau_initializer_4d, error: the params array must ', &
            'be passed. params(1) = epsilon, params(2) = kx, params(3) = ky.'
       stop
    end if

    epsilon = params(1)
    kx      = params(2)
    ky      = params(3)
    factor1 = 0.5_f64/sll_pi

    sll_landau_initializer_4d = factor1*&
         (1.0_f64+cos(kx*x)*cos(ky*y)*exp(-0.5_f64*(vx**2+vy**2)))
  end function sll_landau_initializer_4d


end module sll_common_array_initializers_module
