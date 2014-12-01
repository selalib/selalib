program assert_test
!#include "sll_memory.h"
#include "sll_assert.h"
!#include "sll_working_precision.h"
  implicit none
  integer, dimension(:), allocatable :: a
  integer :: b
  allocate( a(1000) )

  print *, 'Array values: '
  write (*, '(a, i4)') 'The size of a is: ', size(a)
  b = get_val(a, 1)
  write (*, '(a, i4)') 'a(1) = ', b
  b = get_val(a, 117)
  write (*, '(a, i4)') 'a(117) = ', b
  b = get_val(a, 500001) 
  write (*, '(a, i4)') 'a(1001) = ', b




contains

  function get_val( a, i )
    integer :: get_val
    integer, intent(in) :: i
    integer, dimension(:), intent(in) :: a
    SLL_ASSERT( (i .ge. 1) .and. (i .le. size(a)) )
    get_val = a(i)
  end function get_val


end program assert_test
