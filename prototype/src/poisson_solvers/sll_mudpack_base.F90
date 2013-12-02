!> Base module to provide interface to mudpack library
!> This library contains multigrid solvers for PDE equations
module sll_mudpack_base

#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"

implicit none
sll_int32, private :: i, j, k

!> Fishpack solver cartesian 2d
type, public :: mudpack_2d

   sll_real64, dimension(:), allocatable :: work !< array for tmp data
   sll_int32  :: mgopt(4) !< Option to control multigrid
   sll_int32  :: iprm(16) !< Indices to control grid sizes
   sll_real64 :: fprm(6)  !< Real to set boundary conditions
   sll_int32  :: iguess   !< Initial solution or loop over time

end type mudpack_2d

integer, parameter :: CARTESIAN_2D = 2    !< geometry parameter
integer, parameter :: CARTESIAN_3D = 3    !< geometry parameter
integer, parameter :: POLAR        = 11   !< geometry parameter
integer, parameter :: CYLINDRICAL  = 12   !< geometry parameter
integer, parameter :: PERIODIC     = 0    !< boundary condition parameter
integer, parameter :: DIRICHLET    = 1    !< boundary condition parameter
integer, parameter :: NEUMANN      = 2    !< boundary condition parameter

end module sll_mudpack_base
