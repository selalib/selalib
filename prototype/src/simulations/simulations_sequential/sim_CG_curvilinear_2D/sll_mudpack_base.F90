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

enum, bind(C)
   enumerator :: CARTESIAN_2D = 2
   enumerator :: CARTESIAN_3D = 3
   enumerator :: POLAR        = 11
   enumerator :: CYLINDRICAL  = 12
   enumerator :: PERIODIC     = 0
   enumerator :: DIRICHLET    = 1
   enumerator :: NEUMANN      = 2
end enum

end module sll_mudpack_base
