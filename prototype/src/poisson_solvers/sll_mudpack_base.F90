!> @brief
!> Base module to provide interface to mudpack library
!> @details
!> This library contains multigrid solvers for PDE equations
module sll_mudpack_base

#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"
use sll_cubic_spline_interpolator_2d

implicit none
sll_int32, private :: i, j, k

!> Fishpack solver cartesian 2d
type, public :: mudpack_2d

   sll_real64, dimension(:), allocatable :: work !< array for tmp data
   sll_int32  :: mgopt(4) !< Option to control multigrid
   sll_int32  :: iprm(16) !< Indices to control grid sizes
   sll_real64 :: fprm(6)  !< Real to set boundary conditions
   sll_int32  :: iguess   !< Initial solution or loop over time
   !class(sll_interpolator_2d_base), pointer   :: cxx_interp
   sll_int32, dimension(:,:), allocatable :: iwork

end type mudpack_2d

integer, parameter :: CARTESIAN_2D = 2    !< geometry parameter
integer, parameter :: CARTESIAN_3D = 3    !< geometry parameter
integer, parameter :: POLAR        = 11   !< geometry parameter
integer, parameter :: CYLINDRICAL  = 12   !< geometry parameter

integer, parameter :: SLL_SEPARABLE  = 1    !< type of equation
integer, parameter :: SLL_NON_SEPARABLE_WITHOUT_CROSS_TERMS = 2    !< type of equation
integer, parameter :: SLL_NON_SEPARABLE_WITH_CROSS_TERMS = 3    !< type of equation


end module sll_mudpack_base
