program comm_unit_test
  use sll_collective
  use sll_comm_module
#include "sll_memory.h"
#include "sll_working_precision.h"
#include "sll_utilities.h"
  implicit none
#define PROBLEM_SIZE 4

  type(sll_comm_real64), pointer :: comm
  sll_real64, dimension(:), allocatable, target :: array1
  sll_real64, dimension(:), allocatable, target :: array2
  sll_real64, dimension(:), pointer :: buf1
  sll_real64, dimension(:), pointer :: buf2
  sll_int32 :: count
  sll_int32 :: rank
  sll_int32 :: size
  sll_int32 :: ierr
  sll_int32 :: i

  call sll_boot_collective()

  rank = sll_get_collective_rank(sll_world_collective)
  size = sll_get_collective_size(sll_world_collective)
  comm => new_comm_real64( sll_world_collective, 2, PROBLEM_SIZE )
  print *, 'created new comm...'

  call sll_create_comm_real64_ring( comm )
  print *, 'configured the comm as a ring'
  
  SLL_ALLOCATE(array1(PROBLEM_SIZE),ierr)
  SLL_ALLOCATE(array2(PROBLEM_SIZE),ierr)

  do i=1,PROBLEM_SIZE
     array1(i) = rank*(PROBLEM_SIZE)+i
  end do

  print *, 'rank: ', rank, 'array = ', array1(:)

  buf1 => get_buffer(comm,1)
  print *, 'buf1 is asociated: ', associated(buf1)

  buf1(1:PROBLEM_SIZE) = array1(1:PROBLEM_SIZE)

  call comm_send_real64( comm, 1, PROBLEM_SIZE)
  print *, 'rank: ', rank, ' sent buffer 1'
!!$  call sleep(2)
!!$  buf2 => get_buffer(comm,1)

!  print *, 'rank: ', rank, 'contents of buffer in port 2 : ', buf2(:)
!  print *, 'rank: ', rank, ' reading from buffer in port 2: ', 
!  buf2 => get_buffer(comm,2)
!  call comm_receive_real64( comm, 2, count )
  buf2 => get_buffer(comm,2)
  buf2(1:PROBLEM_SIZE) = array1(1:PROBLEM_SIZE)
  call comm_send_real64( comm, 2, PROBLEM_SIZE )
  print *, 'rank: ', rank, 'sent buffer 2'

  call comm_receive_real64( comm, 1, count )
  print *, 'rank: ', rank, 'received buffer 1'
  buf1 => get_buffer(comm,1)
!  print *, 'rank: ', rank, ' received buffer 2, count = ', count
  print *, 'rank: ', rank, 'received contents of buffer in port 1 : ', buf1(:)

  call comm_receive_real64( comm, 2, count )
  print *, 'rank: ', rank, 'received buffer 2'
  buf2 => get_buffer(comm,2)
!  print *, 'rank: ', rank, ' received buffer 2, count = ', count
  print *, 'rank: ', rank, 'received contents of buffer in port 2 : ', buf2(:)
 

  print *, 'proceeding to delete comm...'
  call delete_comm_real64( comm )
  call sll_halt_collective( )
  print *, 'PASSED'

end program comm_unit_test
