!> @ingroup poisson_solvers
!> @brief
!> Base module to provide interface to mudpack library
!> @details
!> This library contains multigrid solvers for PDE equations
!> You have to download and install mudpack 
!> http://www2.cisl.ucar.edu/resources/legacy/mudpack
module sll_mudpack_base
#include "sll_working_precision.h"

implicit none
private

!> Mudpack solver cartesian 2d
type, public :: mudpack_2d

   sll_real64, dimension(:), allocatable :: work !< array for tmp data
   sll_int32  :: mgopt(4)           !< Option to control multigrid
   sll_int32  :: iprm(16)           !< Indices to control grid sizes
   sll_real64 :: fprm(6)            !< Real to set boundary conditions
   sll_int32  :: iguess             !< Initial solution or loop over time
   sll_int32, pointer :: iwork(:,:) !< Internal work array for mudpack library

end type mudpack_2d

integer, parameter, public :: CARTESIAN_2D = 2    !< geometry parameter
integer, parameter, public :: CARTESIAN_3D = 3    !< geometry parameter
integer, parameter, public :: POLAR        = 11   !< geometry parameter
integer, parameter, public :: CYLINDRICAL  = 12   !< geometry parameter

integer, parameter, public :: SLL_SEPARABLE  = 1                        !< type of equation
integer, parameter, public :: SLL_NON_SEPARABLE_WITHOUT_CROSS_TERMS = 2 !< type of equation
integer, parameter, public :: SLL_NON_SEPARABLE_WITH_CROSS_TERMS = 3    !< type of equation


end module sll_mudpack_base
