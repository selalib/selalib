program comm_unit_test_2d
  use sll_collective
  use sll_point_to_point_comms_module
#include "sll_memory.h"
#include "sll_working_precision.h"
#include "sll_utilities.h"
  implicit none

#define SZ_X 4
#define SZ_Y 5
! BUF_SIZE MUST BE THE GREATEST BETWEEN THE SZ_X AND SZ_Y
#define BUF_SIZE SZ_Y

  type(sll_p2p_comm_real64), pointer :: comm
  sll_real64, dimension(:,:), pointer :: main_array
  ! A single 'recyclable' pointer could be used of course...
  sll_real64, dimension(:), pointer :: buf
  sll_int32 :: count
  sll_int32 :: rank
  sll_int32 :: col_size   ! collective size
  sll_int32 :: ierr
  sll_int32 :: i, j
  logical, dimension(1)   :: local_pass
  logical, dimension(1)   :: general_pass

  call sll_boot_collective()

  rank = sll_get_collective_rank(sll_world_collective)
  col_size = sll_get_collective_size(sll_world_collective)
  comm => new_comm_real64( sll_world_collective, 4, BUF_SIZE )
  if(rank == 0) then
     print *, 'created new comm, size = ', col_size
     call flush()
  end if

  ! In this test the processors in the communicator are organized as a 2D ring,
  ! thus there are four ports which are linked with the left(1) and right(2)
  ! and to the bottom(3) and top(4). At the borders of the domain, the ports
  ! are linked periodically.
  call sll_configure_comm_real64_ring_2d( comm,2,2 ) ! not general! change!! this will break if the program is run on a number of procs different than 4

  if(rank == 0) then
     print *, 'configured the comm as a 2D ring'
     call flush()
  end if

  SLL_ALLOCATE(main_array(SZ_X,SZ_Y),ierr)

  ! just give the same values to every process!
  do j=1,SZ_Y
     do i=1,SZ_X
        main_array(i,j) = real((i-1)+SZ_X*(j-1),f64)
     end do
  end do
  !print *, 'rank: ', comm%rank, 'main array = ', main_array(:,:)
  ! Load the buffer on port 1 with the data and send
  buf => get_buffer(comm,1)
  buf(1:SZ_Y) = main_array(1,:)
  call comm_send_real64( comm, 1, SZ_Y)

  buf => get_buffer(comm,2)
  buf(1:SZ_Y) = main_array(SZ_X,:)
  call comm_send_real64( comm, 2, SZ_Y)

  buf => get_buffer(comm,3)
  buf(1:SZ_X) = main_array(:,1)
  call comm_send_real64( comm, 3, SZ_X)

  buf => get_buffer(comm,4)
  buf(1:SZ_X) = main_array(:,SZ_Y)
  call comm_send_real64( comm, 4, SZ_X)

  ! here in principle the user would do some useful work with the main core
  ! of the main array while the communications complete...

  ! And now receive the data.
  call comm_receive_real64( comm, 1, count )
  buf => get_buffer(comm,1)
  if( 0.0_f64 == sum(main_array(SZ_X,:)-buf(1:SZ_Y)) ) then
     local_pass(1) = .true.
  else
     local_pass(1) = .false.
  end if

  call comm_receive_real64( comm, 2, count )
  buf => get_buffer(comm,2)
  if( 0.0_f64 == sum(main_array(1,:)-buf(1:SZ_Y)) ) then
     local_pass(1) = local_pass(1) .and. .true.
  else
     local_pass(1) = local_pass(1) .and. .false.
  end if

  call comm_receive_real64( comm, 3, count )
  buf => get_buffer(comm,3)
  if( 0.0_f64 == sum(main_array(:,SZ_Y)-buf(1:SZ_X)) ) then
     local_pass(1) = local_pass(1) .and. .true.
  else
     local_pass(1) = local_pass(1) .and. .false.
  end if

  call comm_receive_real64( comm, 4, count )
  buf => get_buffer(comm,4)
  if( 0.0_f64 == sum(main_array(:,1)-buf(1:SZ_X)) ) then
     local_pass(1) = local_pass(1) .and. .true.
  else
     local_pass(1) = local_pass(1) .and. .false.
  end if

  call sll_collective_reduce( comm%collective, local_pass, 1, MPI_LAND, 0, &
       general_pass )

  call delete_comm_real64( comm )

  if( rank == 0 ) then
     if( general_pass(1) .eqv. .true.) then  
        print *, 'PASSED'
     else
        print *, 'FAILED'
     end if
  end if
  call sll_halt_collective()

end program comm_unit_test_2d
