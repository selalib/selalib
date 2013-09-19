module sll_comm_module
#include "sll_assert.h"
#include "sll_memory.h"
  use sll_collective
  implicit none

  ! ************************************************************************
  ! The comm abstraction has been adapted from technology present in some
  ! of the fastest parallel codes today, notably, the Desmond molecular
  ! dynamics code. Comm is an abstraction that facilitates the use of
  ! point-to-point communications in a parallel network by hiding many
  ! implementation details and exposing a simplified interface that is 
  ! reminiscent of a shared memory abstraction.
  !
  ! The sll_comm module's mission is to manage a set of ports. A port
  ! is a bidirectional communication path.
  !
  ! The comm can be represented as a graph on a set of nodes. The endpoints
  ! of each edge are represented by (rank, port) pairs. The rank refers to 
  ! a process and the port to an abstraction that acts like a shared memory 
  ! buffer being exchanged between the two processes that share the edge.
  !
  ! Upon initialization of an edge, there is a memory buffer present at
  ! each end. A process can thus write onto this buffer right away. At 
  ! some point, one of the processes (let's call it 'local') can comm_send 
  ! its buffer along the edge to its counterpart (the 'remote'). At this 
  ! instant, the buffer becomes "not-present": there is no writable buffer 
  ! at the local port. If the local process then executes a comm_receive
  ! (which can only be done when there is no buffer present), the local
  ! process blocks until the remote buffer becomes present locally (which
  ! can only happen if the remote has executed a comm_send). Once 
  ! comm_receive returns, there is a new buffer present locally which can
  ! be read/written at will. 
  !
  ! In summary: comm_send makes the local buffer "not-present" (or 
  ! in-transit) while comm_receive makes the buffer most recently held by 
  ! the remote process present (blocking until it is). Apart from each 
  ! edge within the same comm needing simultaneous initialization with
  ! comm_init, the buffers are otherwise independent. 
  !
  ! Notes for threaded code:
  !
  ! If this facility is used in a program threaded with the "threader"
  ! library, one should be aware that the communicators should be created
  ! from the Master thread only. Buffers on each port are hooked up and 
  ! sized from the Master thread only. 
  !
  ! When a buffer is present, any thread may load data. comm_send is to be
  ! called from the Master thread only, and similarly for the comm_receive
  ! call.

  ! The message padding is present to avoid sending zero-length messages.

#define SLL_MAX_NUM_PORTS 4096
#define MESSAGE_PADDING 4
  
  private
  
  ! Since this object contains an MPI_Request, it would be desirable to move 
  ! this function to the collective module. For now, MPI is allowed to spill 
  ! over...

  type buffer_real64
     sll_real64, dimension(:), pointer :: data
     MPI_Request                       :: request
  end type buffer_real64

  type sll_remote
     sll_int32 :: rank
     sll_int32 :: port
  end type sll_remote

  type port_real64
     type(buffer_real64), dimension(:), pointer :: buffer  ! 2 buffers only.
     sll_int32 :: bit
     sll_int32 :: other_rank
     sll_int32 :: other_port
  end type port_real64


  public
  

  ! Consider the dynamic allocation of the 'remotes' array...
  type sll_comm_real64
     sll_int32 :: comm_size =  0
     sll_int32 :: rank      = -1
     sll_int32 :: num_ports
     sll_int64 :: buffer_size
     type(sll_collective_t), pointer :: collective
     type(port_real64), pointer      :: ports  ! array of ports
     type(sll_remote), pointer       :: remotes  ! may not be needed
  end type sll_comm_real64

contains

  function mpi_request( comm, port )
    sll_int32                      :: mpi_request
    type(sll_comm_real64), pointer :: comm
    sll_int32, intent(in)          :: port
    sll_int32 :: bit
    bit = comm%ports(port)%bit
    mpi_request = comm%ports(port)%buffer(bit)%request
  end function mpi_request

  subroutine flip_buffer( comm, port )
    type(sll_comm_real64), pointer :: comm
    sll_int32, intent(in)          :: port
    comm%ports(port)%bit = not(comm%ports(port)%bit)
  end subroutine flip_buffer

  function get_buffer( comm, port )
    sll_real64, dimension(:), pointer :: get_buffer
    type(sll_comm_real64), pointer :: comm
    sll_int32, intent(in)          :: port
    sll_int32 :: bit
    bit = comm%ports(port)%bit
    get_buffer => comm%ports(port)%buffer(bit)%data
  end function get_buffer

  ! Tag generators. 
  ! We need unique integer identifiers that we can use as tags for the
  ! nonblocking communications. These tags should be determined by the
  ! bit flag and the ports of the intervening processes. The simplest way
  ! is to construct an identifier with all this information while keeping 
  ! the distinction between send and receive tags.
  !
  ! The tag is built within a 32-bit integer.

  function receive_tag( bit, my_p, other_port )
    sll_int32 :: receive_tag
    sll_int32 :: bit
    sll_int32 :: my_port
    sll_int32 :: other_port
    sll_int32 :: higher
    sll_int32 :: lower

    ! for the lower part of the tag, pick the lowest 14 bits of the integer
    ! that represents the local process's port together with the value of the
    ! bit field, pushed to the 15th bit position.
    ! In decimal notation, z'3fff' = 16,384
    lower  = ior( ishft(    bit, 14), ior(my_port,   z'3fff'))  
    higher = ior( ishft(not(bit),14), ior(other_port,z'3fff'))

    ! shift the higher part of the tag to the upper bits, starting at bit 16 
    ! and leaving the lower 15 bits available for the lower part of the tag.
    receive_tag = ior(ishft(higher,15),lower)
  end function receive_tag

  ! The send tag is analogous to the receive tag with the difference that
  ! the lower and higher sections of the tag are switched.
  function send_tag( bit, my_port, other_port )
    sll_int64 :: receive_tag
    sll_int64 :: bit
    sll_int64 :: my_port
    sll_int64 :: other_port
    sll_int64 :: higher
    sll_int64 :: lower

    lower  = ior( ishft(not(bit),14), ior(other_port, z'3fff'))  
    higher = ior( ishft(    bit, 14), ior(my_port,    z'3fff'))
    receive_tag = ior(ishft(higher,15),lower)
  end function send_tag

  subroutine initialize_buffer_real64( buff, num_elems )
    type(buffer_real64), intent(out) :: buff
    sll_int32 :: ierr
    SLL_ALLOCATE(buff%data(num_elems),ierr)
    buff%request = MPI_REQUEST_NULL
  end subroutine initialize_buffer_real64

  subroutine initialize_port_real64( port, buf_num_elems )
    type(port_real_64), intent(out) :: port
    sll_int32, intent(in)           :: buf_num_elems
    sll_int32 :: ierr
    port%bit     =  0
    port%other   = -1
    port%other_p =  0
    SLL_ALLOCATE( port%buffer(2), ierr )
    call initialize_buffer_real64( port%buffer(1), buf_num_elems )
    call initialize_buffer_real64( port%buffer(2), buf_num_elems )
  end subroutine initialize_port_real64

  subroutine check_size( comm, size )
    type(sll_comm_real64), pointer :: comm

  end subroutine check_size

  subroutine check_port( comm, port )
    type(sll_comm_real64), pointer :: comm

  end subroutine check_port

  function port_is_busy( comm, port )
    logical                        :: port_is_busy
    type(sll_comm_real64), pointer :: comm
    sll_int32, intent(in)          :: port

    if( mpi_request(comm,port) .ne. MPI_REQUEST_NULL) then
       port_is_busy = .true.
    else
       port_is_busy = .false.
    end if
  end function port_is_busy

  subroutine check_other_rank( comm, other_rank )
    type(sll_comm_real64), pointer :: comm

  end subroutine check_other_rank

  function comm_num_ports( comm )
    type(sll_comm_real64), pointer :: comm

  end function comm_num_ports

  function comm_buffer_size( comm )
    type(sll_comm_real64), pointer :: comm

  end function comm_buffer_size

  function new_comm_real64( collective, num_ports, buffer_size ) result(comm)
    type(sll_comm_real64), pointer :: comm
    sll_int32                      :: ierr
    sll_int32                      :: i
    SLL_ALLOCATE(comm, ierr)

    ! Check if the requested number of ports is greater than some maximum.

    SLL_ALLOCATE(comm%ports(num_ports), ierr)

    do i=1,num_ports
       call initialize_port_real64( comm%ports(i), buffer_size )
    end do

    ! it is a good idea probably to store a 'duplicate' collective,
    ! with an unerlying duplicate communicator generated by MPI_Comm_dup
    ! and use this as the base collective for the comm... this is pending.
  end function new_comm_real64

  subroutine connect_ports( comm, port, remote, remote_port )
    type(sll_comm_real64), pointer :: comm
    sll_int32, intent(in)          :: port
    sll_int32, intent(in)          :: remote
    sll_int32, intent(in)          :: remote_port
    sll_int64 :: tag
    sll_int32 :: ierr
    ! do diverse checks here ... pending

    ! this error must be greatly improved...
    if(comm%ports(port) >= 0) then
       print *, 'comm connect error; port already in use'
       stop
    end if

    comm%ports(port)%other_rank = remote
    comm%ports(port)%other_port = remote_port
    tag = receive_tag( comm%ports(port)%bit, port, comm%ports(port)%other_port)

    ! post a 'receive' on first buffer and then flip it. For now, we allow
    ! the mpi functions to be called directly, but it is desirable to send
    ! this back to the collective module, so a wrapper routine is necessary.
    call MPI_IRecv( &
         comm%ports(port)%buffer(bit)%data, &
         comm%buffer_size, &
         MPI_DOUBLE_PRECISION, &
         comm%ports(port)%other_rank, &
         tag, &
         comm%collective%comm, &
         mpi_request(comm,port), &
         ierr )
  end subroutine connect_ports

  subroutine comm_send_real64( comm, port, count )
    type(sll_comm_real64), pointer :: comm
    sll_int64, intent(in)  :: port
    sll_int64, intent(out) :: count
    sll_int64 :: tag
    sll_int32 :: ierr

    ! arguments tests here
    if(comm%ports(port)%other_rank < 0) then
       print *, 'comm_send_real64(), error, port not connected.'
       stop
    end if

    tag = send_tag( comm%ports(port)%bit, port, comm%ports(port)%other_port)

    call MPI_Isend( &
         comm%ports(port)%buffer(bit)%data, &
         size, &
         MPI_DOUBLE_PRECISION, &
         comm%ports(port)%other_rank, &
         tag, &
         comm%collective%comm, &
         mpi_request(comm,port), &
         ierr )
    if( ierr .ne. MPI_SUCCESS ) then
       print *, 'comm_send_real64() error in mpi call'
       stop
    end if

    call flip_buffer(comm,port)
  end subroutine comm_send_real64

  subroutine comm_receive_real64( comm, port )
    type(sll_comm_real64), pointer :: comm
    sll_int32, intent(in)          :: port
    sll_int32, dimension(MPI_STATUS_SIZE) :: stat
    sll_int64 :: tag
    sll_int32 :: ierr

    ! error checking...
    if(.not. port_is_busy(comm, port) ) then
       print *, 'comm_receive_real64() error: imbalanced send and receive calls'
       stop
    end if

    call MPI_Wait(mpi_request(comm,port), stat, ierr)
    if(ierr .ne. MPI_SUCCESS) then
       print *, 'comm_receive_real64(), MPI_Wait error'
       stop
    end if

    call MPI_Get_Count(stat, MPI_DOUBLE_PRECISION, count, ierr)
    if(ierr .ne. MPI_SUCCESS) then
       print *, 'comm_receive_real64(), MPI_Get_Count() error'
       stop
    end if

    ! check on the other buffer
    call flip_buffer(comm,port)

    if(.not. port_is_busy(comm, port) ) then
       print *, 'comm_receive_real64() error: imbalanced send and receive calls'
       stop
    end if

    call MPI_Wait(mpi_request(comm,port),MPI_STATUS_IGNORE, ierr)
    if(ierr .ne. MPI_SUCCESS) then
       print *, 'comm_receive_real64(), MPI_Wait error in second call.'
       stop
    end if

    tag = receive_tag( comm%ports(port)%bit, port, comm%ports(port)%other_port)

    call MPI_Irecv( &
         comm%ports(port)%buffer(bit)%data, &
         size, &
         MPI_DOUBLE_PRECISION, &
         comm%ports(port)%other_rank, &
         tag, &
         comm%collective%comm, &
         mpi_request(comm,port), &
         ierr )
    call sll_test_mpi_error(ierr, 'MPI_Irecv error')

    call flip_buffer(comm,port)
  end subroutine comm_receive_real64

  subroutine delete_comm_real64( comm )
    type(sll_comm_real64), pointer :: comm
    ! pending finish...
  end subroutine delete_comm_real64

  function port_num_is_valid( num )
    logical               :: port_num_is_valid
    sll_int32, intent(in) :: num
    if( (num .le. 0) .or. (num .ge. MAX_NUM_PORTS ) ) then 
       port_num_is_valid = .false.
    else
       port_num_is_valid = .true.
    end if
  end function port_num_is_valid

  ! Design choice: the function port_num_is_valid() is tailored to be used
  ! within an SLL_ASSERT() call. As such, it would not be available if the
  ! DEBUG flag is turned off. This call may be critical enough that it could
  ! warrant a subroutine that is called all the time, not only in DEBUG mode.
  ! The following is a possibility:
  !
  !  subroutine port_num_is_valid( num )
  !    sll_int32, intent(in) :: num
  !    if( num .le. 0 ) then
  !       write (*, '(a)') 'comm: Negative number of ports requested'
  !       stop 'comm: bad argument, negative number of ports.'
  !    else if ( num .ge. MAX_NUM_PORTS )
  !       write (*, '(a)') 'comm: Too many ports requested'
  !       stop 'comm: bad argument, too many ports.'
  !    end if
  !  end subroutine port_num_is_valid
  !
  ! Thus, make a choice and delete the one not chosen.




  ! sll_new_comm() creates a ... new comm (what else?) from a collective.
  ! Arguments:
  ! - col:      the parent collective
  ! - n_ports:  tne number of ports that the new collective has
  ! - buf_size: the memory buffer size for each port
  ! The returned comm is uninitialized. The only procedure one can apply
  ! to it is sll_comm_configure() to declare a topology. Afterwards, 
  ! sll_initialize_comm() will initialize the comm for this topology.

  function sll_new_comm( col, n_ports, buf_size )
    type(sll_collective_t), pointer :: col
    type(sll_comm_t), pointer       :: sll_new_comm
    sll_int32, intent(in)           :: n_ports
    sll_int32, intent(in)           :: buf_size
    sll_int32                       :: ierr
    SLL_ASSERT( .true. .eq. port_num_is_valid(n_ports) )
    SLL_ALLOCATE( sll_new_comm, ierr )
    sll_new_comm%parent_collective => col

  end function sll_new_comm


  ! sll_delete_comm() deletes the given comm. This subroutine must be 
  ! called by all the processes in the comm.

  subroutine sll_delete_comm( com )

  end subroutine sll_delete_comm




  ! sll_get_comm_parent() returns the collective which is parent to the comm.

  function sll_get_comm_parent( com )

  end function sll_get_comm_parent



  ! sll_get_comm_num_ports() returns the number of ports for the comm.

  function sll_get_comm_num_ports( com )

  end function sll_get_com_num_ports



  ! sll_get_comm_buffer_size() returns the size of the buffers for this comm.

  function sll_get_comm_buffer_size( com )

  end function sll_get_comm_buffer_size



  ! sll_get_comm_buffer() gets the message buffer for a port.
  ! Arguments:
  ! - com: the comm we want the buffer of.
  ! - n:   the port in the buffer
  ! Returns a pointer to the beginning of the buffer. This may need to be
  ! interfaced if we require specializations to different types in this 
  ! Fortran implementation.

  function sll_get_comm_buffer( com, n )

  end function sll_get_comm_buffer


  ! sll_configure_comm() sets up an edge in the topology of the comm.
  ! Arguments:
  ! - com:         the comm to configure.
  ! - n:           the local port number
  ! - remote_rank: the rank of hte remote process (as given by 
  !                sll_get_comm_parent().
  ! - remote_port: the remote port number.
  ! This must be called on a new, initialized comm.

  subroutine sll_configure_comm( com, n, remote_rank, remote_port )

  end subroutine sll_configure_comm




  ! sll_initialize_comm() must be called by all processes within the comm
  ! before any buffers are examined or any communication started and
  ! after the topology has been set by sll_configure_comm().

  subroutine sll_initialize_comm( com )

  end subroutine sll_initialize_comm


  ! sll_comm_send() initiates a send on a given port, relinquishing the
  ! message buffer. 
  ! Arguments:
  ! - com:    the comm
  ! - n:      the port
  ! - length: the size of the message in the buffer.

  subroutine sll_comm_send( com, n, length )

  end subroutine sll_comm_send



  ! sll_comm_receive() completes a pending receive on a port, obtaining
  ! a message buffer.
  ! Arguments:
  ! com: the comm
  ! n:   the port
  ! length: the length of the message in the newly acquired buffer.

  subroutine sll_comm_receive( com, n, length )

  end subroutine sll_comm_receive


  ! The following special patterns return a configured and initialized comm.
  !
  ! sll_new_comm_fan() creates a topology in the shape of a fan, centered
  ! at root, and with the ports connected as (root,i) <---->(i,0) (for i
  ! different than root). That is to say, the rank i has port 0 connected 
  ! to root's port i.
  !                          >   port 0, process 1
  !                        /  
  !                       /
  !                      /
  !           port 1   <          
  !     root  port 2   <------>  port 0, process 2
  !           port 3   <
  !                      \
  !                       \
  !                        \
  !                          >   port 0, process 3
  !
  ! Arguments:
  ! col:       collective on which to base this topology
  ! root:      the root of the fan
  ! buf_size:  the size fo the message buffer for all ports.
  ! Returns an initialized sll_comm_t with the fan topology.

  function sll_new_comm_fan( col, root, buf_size )

  end function sll_new_comm_fan



  ! sll_new_comm_ring() creates a topology in which each process of rank
  ! r has its port 0 connected to process r-1 and its port 1 connected to 
  ! process r+1 modulo sll_collective_size.

  function sll_new_comm_ring( col, buf_size )

  end function sll_new_comm_ring


#undef SLL_MAX_NUM_PORTS
#undef BUFFER_PADDING


end module sll_comm_module
