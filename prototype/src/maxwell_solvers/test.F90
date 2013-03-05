program test_allocate
#include "sll_memory.h"
implicit none

integer  :: error
real(8), dimension(:,:), allocatable :: ex

SLL_CLEAR_ALLOCATE(ex(16,16), error)

write(*,*) ex(1:2,:)
end program test_allocate
