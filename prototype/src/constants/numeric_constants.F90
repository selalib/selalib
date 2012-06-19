!------------------------------------------------------------------------------
! Selalib
!------------------------------------------------------------------------------
!
! MODULE: numeric_constants
!
!> @author
!> Module Author Name and Affiliation
!
! DESCRIPTION: 
!> numeric constants.
!>
!> In this file are defined the numeric constants.
!>
!> <b> How to use it : </b> \n
!> ****************
!>
!> Just add the line
!> \code use numeric_constants \endcode
!> and you have access to the constants defined in the module.
!>
! REVISION HISTORY:
! DD Mmm YYYY - Initial Version
! TODO_dd_mmm_yyyy - TODO_describe_appropriate_changes - TODO_name
!------------------------------------------------------------------------------
module numeric_constants
#include "sll_working_precision.h"
  implicit none


  sll_real64, parameter :: sll_pi = 3.1415926535897932384626433_f64
  ! sll_kx is the fundamental mode in the x-direction. 
  ! It should be set at mesh initialization
  sll_real64            :: sll_kx = 2*sll_pi 

end module numeric_constants
