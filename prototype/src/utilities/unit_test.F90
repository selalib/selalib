program utils_tester
#include "sll_utilities.h"
#include "sll_working_precision.h"
  implicit none

  sll_real64 :: re64
  sll_real32 :: re32
  sll_int64  :: in64
  sll_int32  :: in32
  sll_int64  :: i
  sll_int32  :: acc

  re64 = 3.2
  re32 = 1.0
  in64 =  transfer(z'7fffffff',in64) ! largest 32-bit int, 2**31-1
  in32 = 2

  print *, in64
  print *, 'The size in bytes of an sll_real64 variable is: ', BYTE_SIZEOF(re64)
  print *, 'The size in bytes of an sll_real32 variable is: ', BYTE_SIZEOF(re32)
  print *, 'The size in bytes of an sll_int64 variable is: ', BYTE_SIZEOF(in64)
  print *, 'The size in bytes of an sll_int32 variable is: ', BYTE_SIZEOF(in32)

  print *, 'The size in int32 of an sll_real64 variable is: ',INT32_SIZEOF(re64)
  print *, 'The size in int32 of an sll_real32 variable is: ',INT32_SIZEOF(re32)
  print *, 'The size in int32 of an sll_int64 variable is: ',INT32_SIZEOF(in64)
  print *, 'The size in int32 of an sll_int32 variable is: ',INT32_SIZEOF(in32)

  print *, 'Counting the amount of numbers that are a power of two between 1 '
  print *, 'and 2^31-1... (should be 31, powers 0 through 30, inclusive)'
  acc = 0
  do i=1,in64
     if( is_power_of_two(i) ) then
        acc = acc+1
     end if
  end do
  print *, 'acc = ', acc
  ! Improve this, split the tests between the different sizeofs and the
  ! power of 2 function, so as to accomodate other tests. All in one place
  ! means that if something fails, one has no information on which test
  ! failed...
  if( (BYTE_SIZEOF(re64) .eq. 8) .and. &
      (BYTE_SIZEOF(re32) .eq. 4) .and. &
      (BYTE_SIZEOF(in64) .eq. 8) .and. &
      (BYTE_SIZEOF(in32) .eq. 4) .and. &
      (INT32_SIZEOF(re64).eq. 2) .and. &
      (INT32_SIZEOF(re32).eq. 1) .and. &
      (INT32_SIZEOF(in64).eq. 2) .and. &
      (INT32_SIZEOF(in32).eq. 1) .and. &
      acc .eq. 31 ) then
     print *, 'PASSED TEST'
  else
     print *, 'FAILED TEST'
  end if

end program utils_tester
