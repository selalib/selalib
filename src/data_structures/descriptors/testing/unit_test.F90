program test_descriptors
#include "sll_working_precision.h"
use sll_m_descriptors

implicit none

type(sll_vlasovpoisson_sim) :: testcase


call testcase%parse('  SLL_LANDAU_DIAG       ')
 print *, testcase%name()
 print *, SLL_LANDAU_DIAG%name(), SLL_LANDAU_DIAG%id

print *, "PASSED"

print *, trim('  SLL_LANDAU_DIAG       ')


end program
