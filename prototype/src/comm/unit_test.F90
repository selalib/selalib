program comm_unit_test
  use sll_collective
  use sll_comm_module
#include "sll_memory.h"
#include "sll_working_precision.h"
#include "sll_utilities.h"
  implicit none
#define PROBLEM_SIZE 8

  type(sll_comm_real64), pointer :: comm
  sll_real64, dimension(:), pointer :: array1
  sll_real64, dimension(:), pointer :: array2
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

  SLL_ALLOCATE(array1(PROBLEM_SIZE),ierr)
  SLL_ALLOCATE(array2(PROBLEM_SIZE),ierr)

  do i=1,PROBLEM_SIZE
     array1(i) = rank*size+i
  end do

  print *, 'rank: ', rank, 'array = ', array1(:)

  buf1 => get_buffer(comm,1)
  buf1(1:8) = array1(1:8)

  call comm_send_real64( comm, 1, PROBLEM_SIZE)
  print *, 'rank: ', rank, ' sent buffer 1'

  call comm_receive_real64( comm, 2, count )
  buf2 => get_buffer(comm,2)
  print *, 'rank: ', rank, ' received buffer 2, count = ', count

  call sll_create_comm_real64_ring( comm )
  print *, 'configured the comm as a ring'

  print *, 'proceeding to delete comm...'
  call delete_comm_real64( comm )

  print *, 'PASSED'

end program comm_unit_test
