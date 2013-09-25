program comm_unit_test
  use sll_collective
  use sll_comm_module
#include "sll_memory.h"
#include "sll_working_precision.h"
#include "sll_utilities.h"
  implicit none

  type(sll_comm_real64), pointer :: comm

  call sll_boot_collective()

  comm => new_comm_real64( sll_world_collective, 2, 128 )
  print *, 'created new comm...'

  call sll_create_comm_real64_ring( comm )
  print *, 'configured the comm as a ring'

  print *, 'proceeding to delete comm...'
  call delete_comm_real64( comm )

  print *, 'PASSED'

end program comm_unit_test
