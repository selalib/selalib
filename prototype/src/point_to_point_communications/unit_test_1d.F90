program comm_unit_test
  use sll_collective
  use sll_point_to_point_comms_module
#include "sll_memory.h"
#include "sll_working_precision.h"
#include "sll_utilities.h"
  implicit none
#define PROBLEM_SIZE 4

  type(sll_p2p_comm_real64), pointer :: comm
  sll_real64, dimension(:), pointer :: array1
  sll_real64, dimension(:), pointer :: array_left
  sll_real64, dimension(:), pointer :: array_right
  sll_real64, dimension(:), pointer :: buf1
  sll_real64, dimension(:), pointer :: buf2
  sll_int32 :: count
  sll_int32 :: rank
  sll_int32 :: size
  sll_int32 :: ierr
  sll_int32 :: i
  logical, dimension(1)   :: local_pass
  logical, dimension(1)   :: general_pass

  call sll_boot_collective()

  rank = sll_get_collective_rank(sll_world_collective)
  size = sll_get_collective_size(sll_world_collective)
  comm => new_comm_real64( sll_world_collective, 2, PROBLEM_SIZE )
  if(rank == 0) then
     print *, 'created new comm, size = ', size
     call flush()
  end if

  ! In this test the processors in the communicator are organized as a ring,
  ! thus there are two ports which are linked with the left(1) and right(2).
  call sll_create_comm_real64_ring( comm )

  if(rank == 0) then
     print *, 'configured the comm as a 1D ring'
     call flush()
  end if

  SLL_ALLOCATE(array1(PROBLEM_SIZE),ierr)
  SLL_ALLOCATE(array_left(PROBLEM_SIZE),ierr) 
  SLL_ALLOCATE(array_right(PROBLEM_SIZE),ierr) 

  do i=1,PROBLEM_SIZE
     array1(i)      = rank*PROBLEM_SIZE+i
     array_left(i)  = mod(rank+size-1,size)*PROBLEM_SIZE+i
     array_right(i) = mod(rank+size+1,size)*PROBLEM_SIZE+i
  end do

!  print *, 'rank: ', rank, 'problem size: ', PROBLEM_SIZE, 'array = ', &
!  array1(:)

!  call sll_view_port(comm,1)
!  call sll_view_port(comm,2)

  ! Load the buffer on port 1 with the data and send
  buf1 => get_buffer(comm,1)
  buf1(1:PROBLEM_SIZE) = array1(1:PROBLEM_SIZE)
  call comm_send_real64( comm, 1, PROBLEM_SIZE)
  print *, 'rank: ', rank, ' sent buffer 1'
!!$  call sll_view_port(comm, 1)
!!$  call sll_view_port(comm, 2)
  ! Just check that the proper behavior is obtained. i.e.: the port becomes
  ! unavailable for a write after a send.
  buf1 => get_buffer(comm,1)

!!$  if(.not. associated(buf1)) then
!!$     print *, 'rank: ', rank, 'buffer 1 is not present for writing'
!!$  else
!!$     print *, 'rank: ', rank, 'buffer 1 IS present for writing'
!!$  end if
!!$

  buf2 => get_buffer(comm,2)

!!$  if(.not. associated(buf2)) then
!!$     print *, 'rank: ', rank, 'buffer 2 is not present for writing'
!!$  else
!!$     print *, 'rank: ', rank, 'buffer 2 IS present for writing'
!!$  end if

  ! Load the buffer on port 2 with the data and send.
  buf2(1:PROBLEM_SIZE) = array1(1:PROBLEM_SIZE)
!  print *, 'rank: ', rank, 'sending buffer on port 2:', buf2
  call comm_send_real64( comm, 2, PROBLEM_SIZE)

!!$  call sll_view_port(comm, 1)
!!$  call sll_view_port(comm, 2)

!  print *, 'rank: ', rank, ' sent buffer 2'

  ! And now receive the data.
  call comm_receive_real64( comm, 1, count )
!  print *, 'rank ', rank, ' received count on port 1', count
  buf1 => get_buffer(comm,1)

  if( 0.0_f64 == sum(array_left(:)-buf1(1:PROBLEM_SIZE)) ) then
     local_pass(1) = .true.
  else
     local_pass(1) = .false.
  end if

!  print *, 'rank ', rank, 'buffer received on port 1: ', buf1
  call comm_receive_real64( comm, 2, count )
  print *, 'rank ', rank, ' received count on port 2 ', count
  buf2 => get_buffer(comm,2)
!  print *, 'rank: ', rank, ' reading from buffer in port 2: ', buf2
  if( 0.0_f64 == sum(array_right(:)-buf2(1:PROBLEM_SIZE)) ) then
     local_pass(1) = local_pass(1) .and. .true.
  else
     local_pass(1) = local_pass(1) .and. .false.
  end if

  call sll_collective_reduce( comm%collective, local_pass, 1, MPI_LAND, 0, &
       general_pass )

  print *, 'proceeding to delete comm...'
  call delete_comm_real64( comm )

  print *, "after deletion, is comm associated?", associated(comm)

  if( rank == 0 ) then
     if( general_pass(1) .eqv. .true.) then  
        print *, 'PASSED'
     else
        print *, 'FAILED'
     end if
  end if
  call sll_halt_collective()
end program comm_unit_test
