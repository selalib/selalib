! unit test for the proposed memory allocator

program memory_tester
!  use sll_strings
#include "sll_memory.h"
  implicit none
  integer :: err
  integer :: i, j, k
  real, dimension(:), pointer :: a=>null()
  real, dimension(:,:,:), pointer :: b=>null()
  real, dimension(:), allocatable :: c
!  integer, dimension(:), pointer :: ai=>null()

  print *, "Memory allocator testing program"

  print *, "allocating something small... "
  SLL_CLEAR_ALLOCATE(b(1:4,1:3,1:2),err)
  print *, 'allocation successful.'
  print *, 'array just after allocation: '
  print *, b(:,:,:)
  print *, 'change the values of the elements...'
  forall(i=1:4, j=1:3, k=1:2) b(i,j,k) = 7 ! initialize manually...
  print *, b(:,:,:)
  print *,'clear the array with SLL_INIT_ARRAY():'
  SLL_INIT_ARRAY(b,0)
  print *, 'array after clearing:'
  print *, b(:,:,:)
  print *, 'allocate and initialize to zero a large array'
  SLL_CLEAR_ALLOCATE(c(1:100000000),err)
  do i=1,size(c)
     if (c(i) .ne. 0.0) then
        print *, 'non zero value found in cleared array. ERROR!'
        stop
     end if
  end do
  print *, 'successful call to SLL_CLEAR_ALLOCATE'
! turn the following on to test the deallocator
#if 0 
  print *, 'Deallocate an unitialized pointer (i.e. force a crash)'
  SLL_DEALLOCATE(a,err)
  print *, 'deallocation successful'
  write(*, '(a, l5)') 'Is the previous pointer associated?', associated(a)
#endif
  print *, 'Proceeding to break the allocator by requesting a gazillion bytes: '
  SLL_ALLOCATE(b(1:10000,1:10000,1:10000), err)

  print *, selected_int_kind(12)
  print *, huge(0_8)

!  call TEST_MACRO(a, 1,90000000000_8 )
#ifndef STDF95
  SLL_ALLOCATE(a(90000000000_8), err)
#endif
 


end program memory_tester
