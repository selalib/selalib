program comm_unit_test_2d
  use sll_collective
  use sll_point_to_point_comms_module
#include "sll_memory.h"
#include "sll_working_precision.h"
#include "sll_utilities.h"
  implicit none

  ! This test takes a 2D array and sends the nodes at the borders to the
  ! neighboring processors, mimicking the type of operations that one would do
  ! on a ghost point scheme.

#define SZ_X 256
#define SZ_Y 512
! BUF_SIZE MUST BE THE GREATEST BETWEEN THE SZ_X AND SZ_Y
#define BUF_SIZE SZ_Y

  type(sll_p2p_comm_real64), pointer :: comm
  sll_real64, dimension(:,:), pointer :: main_array
  sll_real64, dimension(:), pointer :: buf ! used to refer to port buffers
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

  ! each comm has, in this 2D example, 4 ports, with an associated memory
  ! buffer of BUF_SIZE real64 elements.
  !
  !                +------------------------+
  !                |         port 4         |
  !                |                        |
  !                |                        |
  !                |                        |
  !                |port 1            port 2|
  !                |                        |
  !                |                        |
  !                |                        |
  !                |         port 3         |
  !                +------------------------+
  !
  ! All processes will contain a 'comm' configured like this after new_comm
  ! is called.

  comm => new_comm_real64( sll_world_collective, 4, BUF_SIZE )
  if(rank == 0) then
     print *, 'created new comm, size = ', col_size
     call flush(6)
  end if

  ! In this test the processors in the communicator are organized as a 2D ring,
  ! thus there are four ports which are linked with the left (port 1) and 
  ! right(port 2) and to the bottom(port 3) and top(port 4). At the borders of 
  ! the domain, the ports are linked periodically. The 2D ring means that
  ! the processors are arranged as a 2D mesh, with port 2 connected to the
  ! port 3 of the right neighbor, port 4 to the port 3 of the upper neighbor
  ! and so on. At the edges of the processor mesh, the processes are connected
  ! periodically. Thus the processors are connected as a toroidal surface.

  call sll_configure_comm_real64_torus_2d( comm,2,2 ) ! not general! change!! this will break if the program is run on a number of procs different than 4

  if(rank == 0) then
     print *, 'configured the comm as a toroidal surface'
     call flush(6)
  end if

  SLL_ALLOCATE(main_array(SZ_X,SZ_Y),ierr)

  ! just give the same values to every process!
  do j=1,SZ_Y
     do i=1,SZ_X
        main_array(i,j) = real((i-1)+SZ_X*(j-1),f64)
     end do
  end do

  ! Load the buffer on port 1 with the data and send
  ! one should always first ask comm for the buffer to write or read. The
  ! second argument is the port number one is working with.
  buf => get_buffer(comm,1)    
  buf(1:SZ_Y) = main_array(1,:)
  ! Non-blocking send. We indicate the port to send and how many elements.
  call comm_send_real64( comm, 1, SZ_Y)
  ! After this send operation, a call to get_buffer would yield a null pointer.
  ! There would thus be no writable buffer at this point.

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

  ! And now receive the data. Just indicate the number of the port to receive.
  ! The count variable is an 'out' variable that has the information on how
  ! many elements were received.
  call comm_receive_real64( comm, 1, count )
!  print *, 'rank ', rank, 'received ', count, 'elements on port 1'
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
