module initial_distribution_functions
#include "sll_working_precision.h"
use sll_constants
  implicit none
  
contains 
  elemental function landau(x, v)
    sll_real64, intent(in) :: x, v
    sll_real64 :: landau
    sll_real64, parameter :: eps = 0.01_f64
    
    ! sll_kx is defined in the sll_constants module. 
    ! It is set at the mesh initialization 
    landau = ( 1 + eps * cos(sll_kx*x) ) / sqrt(2*sll_pi) * exp(-0.5_f64*v*v)
  end function landau
  
  elemental function two_stream(x,v) result(fval)
    sll_real64, intent(in) :: x, v
    sll_real64 :: fval
    ! local variables
    sll_real64, parameter :: eps = 0.01_f64
    sll_real64, parameter :: xi = 0.90_f64
    sll_real64, parameter :: v0 = 2.4_f64
    sll_real64 :: vv
    ! sll_kx is defined in the sll_constants module. 
    ! It is set at the mesh initialization 
    vv = v*v
    fval=(1+eps*((cos(2*sll_kx*x)+cos(3*sll_kx*x))/1.2_f64+cos(sll_kx*x)))* &
         (1/sqrt(2*sll_pi))*((2-2*xi)/(3-2*xi))*  &
         (1+.5_f64*vv/(1-xi))*exp(-.5_f64*vv)
    !fval=(1+eps*eps*cos(sll_kx*x)+eps*cos(2*sll_kx*x))*0.5_f64/sqrt(2*pi) &
    !     *(exp(-.5_f64*(vx-v0)**2)+ exp(-.5_f64*(vx+v0)**2))
    !fval=(1+eps*cos(sll_kx*x))*0.5_f64/sqrt(2*pi)*(exp(-.5_f64*(vx-v0)**2) &
    !     + exp(-.5_f64*(vx+v0)**2))
    !fval= 1/sqrt (2 * sll_pi) * ( 0.9_f64 * exp (-.5_f64 * vx * vx) &
    !     + 0.2_f64 * exp(-0.5_f64 * (vx - 4.5_f64)*(vx - 4.5_f64)   &
    !     /(0.5_f64 * 0.5_f64))) * (1. + 0.03_f64 * cos (0.3_f64 * x))
    ! fval=(1+eps*cos(sll_kx*x))*1/sqrt(2*pi)*exp(-.5_f64*vv)
    ! fval=exp(-.5_f64*(xx+vv))
  end function two_stream

  ! Samuel : I've change elemental by pure because there is an error with 
  ! gcc 4.7.0 in the unit_test because he does'nt like pointer to elemental 
  ! function.

  pure function gaussian(x,v) result(fval)
    sll_real64, intent(in) :: x, v
    sll_real64 :: fval
    ! local variables
    sll_real64, parameter :: xoffset = 0.5_f64
    sll_real64, parameter :: voffset = 0.5_f64
    sll_real64 :: xx, vv
  
    xx = (x - xoffset )**2
    vv = (v - voffset )**2
    fval =  exp(-0.5_f64*(xx+vv)*40.)
  end function gaussian

end module
