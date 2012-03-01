
!*********************************************************
!
! Selalib      
! Module: util_constants.F90
!
!> @brief 
!> Module of utilitity constants for sll_splines unit test
!   
!> @authors                    
!> Aliou DIOUF (aliou.l.diouf@inria.fr)
!                                  
!*********************************************************

module util_constants
#include "sll_working_precision.h"
  use numeric_constants
  implicit none

  sll_int32,  parameter :: NP    = 512
  sll_int32,  parameter :: NPX1  =  65
  sll_int32,  parameter :: NPX2  = 33
  sll_real64, parameter :: X1MIN = -2.0_f64*sll_pi
  sll_real64, parameter :: X1MAX = 2.0_f64*sll_pi
  sll_real64, parameter :: X2MIN = -sll_pi
  sll_real64, parameter :: X2MAX = sll_pi  
  sll_real64, parameter :: XMIN  = -sll_pi
  sll_real64, parameter :: XMAX  = sll_pi

end module util_constants
