module sll_mudpack_base

#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"

implicit none
sll_int32, private :: i, j, k

!> Fishpack solver cartesian 2d
type, public :: mudpack_2d

   sll_real64, dimension(:), allocatable :: work
   sll_int32  :: mgopt(4)
   sll_int32  :: iprm(16)
   sll_real64 :: fprm(6)
   sll_int32  :: iguess

end type mudpack_2d

integer, parameter :: CARTESIAN_2D = 2
integer, parameter :: CARTESIAN_3D = 3
integer, parameter :: POLAR        = 11
integer, parameter :: CYLINDRICAL  = 12
integer, parameter :: PERIODIC     = 0
integer, parameter :: DIRICHLET    = 1
integer, parameter :: NEUMANN      = 2

end module sll_mudpack_base
