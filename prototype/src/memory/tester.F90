! What we require from the memory allocators is very little:
! To allocate the memory and in case of failure, to provide information
! about the call that failed.

program memory_tester
!  use sll_strings
#include "sll_memory.h"
  implicit none
  integer :: err
  integer :: i, j, k
  real, dimension(:), pointer :: a=>null()
  real, dimension(:,:,:), pointer :: b=>null()
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

! turn the following on to test the deallocator
#if 0 
  print *, b(:,:,:)
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
  SLL_ALLOCATE(a(90000000000_8), err)
 

#if 0  
contains

subroutine sll_allocate2(a_ptr, limits, descriptor)
  intrinsic :: modulo, reshape
  integer, parameter :: long_int = selected_int_kind(12)
  real, dimension(:), pointer :: a_ptr
  integer, dimension(:), intent(in) :: limits
  character(len=*), intent(in) :: descriptor
  integer :: err
  integer :: num_dims
  integer :: i
  integer(kind=long_int) :: sz = 1
  num_dims = size(limits)
  if (modulo(num_dims,2) .ne. 0) then
     stop 'inconsistent limits passed to sll_allocate'
  end if
  num_dims = num_dims/2
  do i=1,num_dims
     sz = sz*(limits(i+1)-limits(i))
  end do
  call sll_allocate(a_ptr, 1, sz, descriptor)
!  call reshape(a_ptr, limits)
end subroutine sll_allocate2

#endif

end program memory_tester
