program test_descriptors
#include "sll_working_precision.h"
use sll_descriptors

implicit none

sll_int32 :: idx

do idx=1, 5
print *, sll_vp_descriptor_key(idx)

end do

print *, "PASSED"

end program