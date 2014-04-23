!------------------------------------------------------------------------------
! Selalib
!------------------------------------------------------------------------------
!
! MODULE: pic_1d_particle_loading
!
! DESCRIPTION:
!> @author Jakob Ameres
!> @brief Loading utility for particle-in-cell method in 1d
!> @details This module provides tools as distribution sampling for the PIC loading mechanism.
!
! REVISION HISTORY:
! DD Mmm YYYY - Initial Version
! TODO_dd_mmm_yyyy - TODO_describe_appropriate_changes - TODO_name
!------------------------------------------------------------------------------
module pic_1d_particle_loading
#include "sll_working_precision.h"
#include "sll_assert.h"
  implicit none


contains

sll_real function gaussianrnd( mu , sigma  ) RESULT(ans)
#include "sll_working_precision.h"
use sll_constants
IMPLICIT NONE
sll_real :: mu    !< mean
sll_real :: sigma !< standard deviation
! Box Muller Wiener
sll_real :: R1, R2, X ,Y


CALL random_number(R1)
CALL random_number(R2)

R1= 1.0-R1
R1 = -ALOG(real(R1))
R1 = SQRT(2.0*R1)
R2 = 2.0*sll_pi*R2

X  = R1*COS(R2)- mu
Y  = R1*SIN(R2) -mu

ans=X
endfunction



end module
