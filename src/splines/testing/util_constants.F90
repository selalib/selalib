#ifndef DOXYGEN_SHOULD_SKIP_THIS
!*********************************************************
!
! Selalib      
! Module: util_constants.F90
!
! @brief 
! Module of utilitity constants for sll_splines unit test
!  
! @authors                    
! Aliou DIOUF (aliou.l.diouf@inria.fr)
!                                  
!*********************************************************

module util_constants
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_working_precision.h"

  use sll_m_constants, only: &
    sll_p_pi

  implicit none

  public :: &
    np, &
    npx1, &
    npx2, &
    r1, &
    r2, &
    x1max, &
    x1min, &
    x2max, &
    x2min, &
    xmax, &
    xmin

  private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  sll_int32,  parameter :: np    =  65
  sll_int32,  parameter :: npx1  =  129
  sll_int32,  parameter :: npx2  =  65
  sll_real64, parameter :: x1min = -2.0_f64*sll_p_pi
  sll_real64, parameter :: x1max =  2.0_f64*sll_p_pi
  sll_real64, parameter :: x2min = -sll_p_pi
  sll_real64, parameter :: x2max =  sll_p_pi  
  sll_real64, parameter :: xmin  = -sll_p_pi
  sll_real64, parameter :: xmax  =  sll_p_pi
  ! parameters for 2d transformation with fine tuned boundary conditions.
  sll_real64, parameter :: r1 = 0.1_f64
  sll_real64, parameter :: r2 = 1.0_f64
  sll_real64, parameter :: theta_min = 0.0_f64
  sll_real64, parameter :: theta_max = 2.0_f64*sll_p_pi

end module util_constants

#endif /* DOXYGEN_SHOULD_SKIP_THIS */
