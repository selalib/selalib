program test_descriptors
#include "sll_working_precision.h"
use sll_descriptors

implicit none

sll_int32 :: idx
type(sll_vlasovpoisson_sim) :: testcase


call testcase%parse('  SLL_LANDAU_DIAG       ')
 print *, testcase%name()
 print *, SLL_LANDAU_DIAG%name(), SLL_LANDAU_DIAG%id

print *, "PASSED"

print *, trim('  SLL_LANDAU_DIAG       ')


end program