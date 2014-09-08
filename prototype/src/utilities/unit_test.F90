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


  ! for factorial
  sll_int32  :: n1
  sll_int64  :: n2
  logical    :: test_passed
  sll_int64  :: largest_int64
  sll_int64  :: largest_int32
  sll_int64  :: acc64

  re64 = 3.2
  re32 = 1.0
  in64 =  transfer(z'7fffffff',in64) ! largest 32-bit int, 2**31-1
  !largest 64-bit integer 2**63-1
  largest_int64 = transfer(z'7fffffffffffffff',largest_int64)
  ! largest 32-bit int, 2**31-1
  largest_int32 = transfer(z'7fffffff',largest_int32)
  in32 = 2


  test_passed = .true.

  ! **********************************************************************
  !
  !          First test: is_power_of_two()
  !
  ! **********************************************************************

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

  if(acc == 31) then
     test_passed = test_passed .and. .true.
     print *, '... yes, 31.'
  else
     test_passed = test_passed .and. .false.
     print *, 'Error with is_power_of_two(), acc = ', acc
  end if

  ! **********************************************************************
  !
  !          Second test: is_power_of_two()
  !
  ! **********************************************************************

  if( (BYTE_SIZEOF(re64) .eq. 8) .and. &
      (BYTE_SIZEOF(re32) .eq. 4) .and. &
      (BYTE_SIZEOF(in64) .eq. 8) .and. &
      (BYTE_SIZEOF(in32) .eq. 4) .and. &
      (INT32_SIZEOF(re64).eq. 2) .and. &
      (INT32_SIZEOF(re32).eq. 1) .and. &
      (INT32_SIZEOF(in64).eq. 2) .and. &
      (INT32_SIZEOF(in32).eq. 1) .and. &
      acc .eq. 31 ) then
     test_passed = test_passed .and. .true.
  else
     test_passed = test_passed .and. .false.
     print *, 'Test of BYTE_SIZEOF and INT32_SIZEOF FAILED'
  end if

  ! **********************************************************************
  !
  !          Third test: factorial()
  !
  ! **********************************************************************
!!$  do i=0,22
!!$     acc64 = sll_factorial(int(i,i64))
!!$     print *, 'n = ', i, 'factorial = ', acc64
!!$  end do

  print *, 'largest int32 = ', in64
  print *, 'largest int64 = ', largest_int64
  acc64 = sll_factorial(20_i64)
  if( acc64 .ne. 2432902008176640000_i64 ) then
     test_passed = test_passed .and. .false.
     print *, 'test of factorial function failed'
     print *, 'factorial = ', sll_factorial(22_i32)
  else
     test_passed = test_passed .and. .true.
     print *, 'factorial is OK...'
  end if

  ! **********************************************************************
  !
  !          Final checking of test_passed flag
  !
  ! **********************************************************************

  if( test_passed .eqv. .true. ) then
     print *, "PASSED"
  else
     print *, "FAILED"
  end if

end program utils_tester
