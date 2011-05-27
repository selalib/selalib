program assert_test
#include "sll_memory.h"
#include "sll_assert.h"
#include "sll_working_precision.h"
  implicit none
  sll_int32, dimension(:), allocatable :: a
  sll_int32 :: ierr
  SLL_ALLOCATE( a(1000), ierr )

  print *, 'Array values: '
  write (*, '(a, i4)') 'The size of a is: ', size(a)
  write (*, '(a, i4)') 'a(1) = ', get_val(a, 1)
  write (*, '(a, i4)') 'a(117) = ', get_val(a, 117)
  write (*, '(a, i4)') 'a(1001) = ', get_val(a, 500001)




contains

  function get_val( a, i )
    sll_int32 :: get_val
    sll_int32, intent(in) :: i
    sll_int32, dimension(:), intent(in) :: a
    SLL_ASSERT( (i .ge. 1) .and. (i .le. size(a)) )
    get_val = a(i)
  end function get_val


end program assert_test
