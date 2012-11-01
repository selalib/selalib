module simulation_VP1D_cartesian_non_unif

#include "sll_working_precision.h"
#include "sll_assert.h"
#include "sll_memory.h"

  use sll_simulation_base
  implicit none

  type, extends(sll_simulation_base_class) :: &
    sll_simulation_VP1D_cartesian_non_unif
    ! Numerical parameters
    sll_real64 :: dt
  contains
    procedure, pass(sim) :: run => run_VP1D_cartesian_non_unif
  end type sll_simulation_VP1D_cartesian_non_unif

contains


  ! Note that the following function has no local variables, which is silly...
  ! This just happened since the guts of the unit test were transplanted here
  ! directly, but this should be cleaned up.
  subroutine run_VP1D_cartesian_non_unif(sim)
    class(sll_simulation_VP1D_cartesian_non_unif), intent(inout) :: sim
  end subroutine run_VP1D_cartesian_non_unif


end module simulation_VP1D_cartesian_non_unif