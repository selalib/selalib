! ***************************************************************************
! sll_collective is a module that encapsulates our calls to the MPI library.
! Any module that wants to access distributed multiprocessing capabilities 
! must use sll_collective instead of MPI directly. If one were to desire
! some MPI capability that is not included here, such capability should not
! be plugged in directly but rather this module should be expanded.
!
! By centralizing the interaction with the MPI library we achieve several
! things:
! - obtain an explicitely reduced 'parallel' vocabulary to build our 
!   application,
! - centralize argument checking and some error handling, and
! - adjust the interface to our wishes.
!
! ***************************************************************************

module sll_collective
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"
  implicit none
  ! This is the only place in the prototype that should have to include
  ! the mpi header file.
  include 'mpif.h'
  !implicit none

  ! **********************************************************************
  ! The sll_collective_t is essentially a wrapper around an MPI
  ! communicator. We store some of the most frequently accessed values
  ! to save a little on the overhead of the corresponding function calls.
  !
  ! The original idea was to have this as an opaque type. But how to 
  ! do this in Fortran without forward declarations? Any help would be
  ! appreciated... In any case, the only means allowed to interact with 
  ! this type, outside of this module, are the functions/subroutines 
  ! defined herein.
  !
  !***********************************************************************
  type sll_collective_t
     type(sll_collective_t), pointer :: parent=>null()
     sll_int32                       :: comm   ! communicator
     sll_int32                       :: color
     sll_int32                       :: key
     sll_int32                       :: rank
     sll_int32                       :: size
  end type sll_collective_t

  ! **********************************************************************
  ! A few rare global variables.
  ! The sll_world_collective just wraps around MPI_WORLD_COMM. This gets
  ! initialized in the call to sll_boot_collective().
  !
  ! **********************************************************************

  type(sll_collective_t), pointer    :: sll_world_collective


  interface sll_collective_bcast
     module procedure sll_collective_bcast_real
  end interface
  
  
  interface sll_collective_gather
     module procedure sll_collective_gather_real
  end interface
  
  interface sll_collective_allgather
     module procedure sll_collective_allgather_int
  end interface

  interface sll_collective_allgatherv
     module procedure sll_collective_allgatherv_real
  end interface

  interface sll_collective_gatherv
     module procedure sll_collective_gatherv_real
  end interface
  
  interface sll_collective_scatter
     module procedure sll_collective_scatter_real
  end interface

  interface sll_collective_scatterv
     module procedure sll_collective_scatterv_real
  end interface
  
  interface sll_collective_allreduce
     module procedure sll_collective_allreduce_real, &
                      sll_collective_allreduce_logical
  end interface
  
  interface sll_collective_reduce
     module procedure sll_collective_reduce_real
  end interface

  interface sll_collective_alltoall
     module procedure sll_collective_alltoall_int
  end interface

  interface sll_collective_alltoallV
     module procedure sll_collective_alltoallV_int, &
                      sll_collective_alltoallV_real
  end interface

contains !************************** Operations **************************

  ! First a couple of utilities for this module
  subroutine sll_check_collective_ptr( ptr )
    type(sll_collective_t), pointer :: ptr
    if( .not. associated(ptr)) then
       write (*, '(a)') 'sll_check_collective_ptr: non-associated pointer'
       stop 'SLL_CHECK_COLLECTIVE_PTR'
    end if
  end subroutine sll_check_collective_ptr

  subroutine sll_test_mpi_error( ierr, descriptor )
    sll_int32, intent(in)        :: ierr
    character(len=*), intent(in) :: descriptor

    SLL_ASSERT( ierr .eq. MPI_SUCCESS ) ! redundant, make a choice between
                                        ! asserting and this function.
    if( ierr .ne. MPI_SUCCESS ) then
       write (*, '(a, a)') 'MPI error code failure: ', descriptor
    end if
  end subroutine sll_test_mpi_error

  ! Following what is exposed to the users of the module.

  ! sll_boot_collective allocates and initializes the global variable 
  ! sll_world_collective and boots the MPI environment.

  subroutine sll_boot_collective( )
    sll_int32 :: ierr
    call MPI_Init(ierr)
    SLL_ALLOCATE( sll_world_collective, ierr )
    sll_world_collective%comm     = MPI_COMM_WORLD
    sll_world_collective%color    = 0
    sll_world_collective%key      = 0
    call MPI_COMM_RANK( MPI_COMM_WORLD, sll_world_collective%rank, ierr )
    call sll_test_mpi_error( ierr, 'sll_boot_collective(): MPI_COMM_RANK()' )
    call MPI_COMM_SIZE( MPI_COMM_WORLD, sll_world_collective%size, ierr )
    call sll_test_mpi_error( ierr, 'sll_boot_collective(): MPI_COMM_SIZE()')
  end subroutine sll_boot_collective

  subroutine sll_halt_collective( )
    sll_int32 :: ierr
    call MPI_BARRIER( MPI_COMM_WORLD, ierr )
    call sll_test_mpi_error( ierr, &
         'sll_halt_collective(): MPI_BARRIER()' )
    call MPI_Finalize(ierr)
  end subroutine sll_halt_collective

  ! sll_new_collective follows a somewhat similar calling convention as 
  ! MPI_COMM_SPLIT() and is a wrapper around it. Error checking
  ! is done internally. The departure from the standard syntax is to permit 
  ! a call like:
  !
  ! type(sll_collective_t), pointer :: new_col
  !
  ! followed by
  !
  ! new_col = sll_new_collective( parent, color, key )
  !
  ! Note that the library's basic types are always pointers. The alternative
  ! syntax would be:
  !
  ! type(sll_collective_t), pointer :: new_col
  !
  ! call sll_new_collective( parent, color, key, new_col )
  !
  ! The collective stores the information:
  ! - a pointer to the parent collective from which it was split
  ! - the values of color and key used for the splitting
  ! - the other fields, rank and size are pre-populated for convenience.

  function sll_new_collective( parent, color, key )
    type(sll_collective_t), pointer  :: parent
    sll_int32, intent(in)            :: color
    sll_int32, intent(in)            :: key
    sll_int32                        :: ierr
    type(sll_collective_t), pointer  :: sll_new_collective

    call sll_check_collective_ptr( parent )
    SLL_ALLOCATE( sll_new_collective, ierr )
    sll_new_collective%parent => parent
    sll_new_collective%color  = color
    sll_new_collective%key    = key
    call MPI_COMM_SPLIT( parent%comm, color, key, sll_new_collective%comm, &
                         ierr)
    call sll_test_mpi_error( ierr, 'sll_new_collective(): MPI_COMM_SPLIT()')
    ! fill out the rest of the fields.
    call MPI_COMM_RANK(sll_new_collective%comm, sll_new_collective%rank, ierr)
    call sll_test_mpi_error(ierr, 'sll_new_collective(): MPI_COMM_RANK()')
    call MPI_COMM_SIZE(sll_new_collective%comm, sll_new_collective%size, ierr)
    call sll_test_mpi_error(ierr, 'sll_new_collective(): MPI_COMM_RANK()')
  end function sll_new_collective

  subroutine sll_delete_collective( col )
    type(sll_collective_t), pointer :: col
    sll_int32                       :: ierr
    call sll_check_collective_ptr( col )
    SLL_DEALLOCATE( col, ierr )
  end subroutine sll_delete_collective

  function sll_get_collective_comm( col )
    type(sll_collective_t), pointer :: col
    sll_int32                       :: sll_get_collective_comm
    call sll_check_collective_ptr( col )
    sll_get_collective_comm = col%comm
  end function sll_get_collective_comm

  function sll_get_collective_rank( col )
    type(sll_collective_t), pointer :: col
    sll_int32                       :: sll_get_collective_rank
    call sll_check_collective_ptr( col )
    sll_get_collective_rank = col%rank
  end function sll_get_collective_rank

  function sll_get_collective_color( col )
    type(sll_collective_t), pointer :: col
    sll_int32                       :: sll_get_collective_color
    call sll_check_collective_ptr( col )
    sll_get_collective_color = col%color
  end function sll_get_collective_color

  function sll_get_collective_size( col )
    type(sll_collective_t), pointer :: col
    sll_int32                       :: sll_get_collective_size
    call sll_check_collective_ptr( col )
    sll_get_collective_size = col%size
  end function sll_get_collective_size

  function sll_get_collective_parent( col )
    type(sll_collective_t), pointer :: col
    type(sll_collective_t), pointer :: sll_get_collective_parent
    call sll_check_collective_ptr( col )
    sll_get_collective_parent => col%parent
  end function sll_get_collective_parent

  subroutine sll_collective_barrier( col )
    type(sll_collective_t), pointer :: col
    sll_int32                         :: ierr
    call sll_check_collective_ptr( col )
    call MPI_BARRIER( col%comm, ierr )
    call sll_test_mpi_error( ierr, 'sll_collective_barrier(): MPI_BARRIER()' )
  end subroutine sll_collective_barrier


  ! Start with something simple, like a buffer of 'real's...
  subroutine sll_collective_bcast_real( col, buffer, size, root )
    type(sll_collective_t), pointer      :: col
    sll_real32, dimension(:), intent(in) :: buffer ! what would change...
    sll_int32, intent(in)                :: size
    sll_int32, intent(in)                :: root
    sll_int32                            :: ierr
    call MPI_BCAST( buffer, size, MPI_REAL, root, col%comm, ierr )
    call sll_test_mpi_error( ierr, &
         'sll_collective_bcast_real(): MPI_BCAST()' )
  end subroutine sll_collective_bcast_real


  ! Consider joining the next 'gather' interfaces into a single one.

  subroutine sll_collective_gather_real( col, send_buf, send_sz, root, &
       rec_buf )
    type(sll_collective_t), pointer      :: col
    sll_real32, dimension(:), intent(in) :: send_buf ! what would change...
    sll_int32                            :: send_sz
    sll_real32, dimension(:), intent(in) :: rec_buf  ! would also change
    sll_int32, intent(in)                :: root
    sll_int32                            :: ierr
    sll_int32                            :: rec_count ! size of receive buf
    ! FIXME: add some argument checking here
    rec_count = send_sz*col%size
    call MPI_GATHER( send_buf, send_sz, MPI_REAL, rec_buf, rec_count, &
         MPI_REAL, root, col%comm, ierr )
    call sll_test_mpi_error( ierr, &
         'sll_collective_gather_real(): MPI_GATHER()' )
  end subroutine sll_collective_gather_real


  subroutine sll_collective_gatherv_real( col, send_buf, send_sz, &
       recvcnts, displs, root, rec_buf )
    type(sll_collective_t), pointer      :: col
    sll_real32, dimension(:), intent(in) :: send_buf ! what would change...
    sll_int32, intent(in)                :: send_sz
    sll_int32, dimension(:), intent(in)  :: recvcnts
    sll_int32, dimension(:), intent(in)  :: displs
    sll_real32, dimension(:), intent(in) :: rec_buf  ! would also change
    sll_int32, intent(in)                :: root
    sll_int32                            :: ierr
    ! FIXME: Argument checking
    call MPI_GATHERV( send_buf, send_sz, MPI_REAL, rec_buf, recvcnts, &
         displs, MPI_REAL, root, col%comm, ierr )
    call sll_test_mpi_error( ierr, &
         'sll_collective_gatherv_real(): MPI_GATHERV()' )
  end subroutine sll_collective_gatherv_real

  subroutine sll_collective_allgather_int( col, send_buf, send_sz, &
       recv_buf, recv_sz )
    type(sll_collective_t), pointer        :: col
    sll_int32, dimension(:), intent(in)    :: send_buf ! what would change...
    sll_int32, intent(in)                  :: send_sz
    sll_int32, dimension(:), intent(inout) :: recv_buf ! would also change
    sll_int32, intent(in)                  :: recv_sz  
    sll_int32                              :: ierr
    ! FIXME: Argument checking
    call sll_check_collective_ptr( col )
    call MPI_BARRIER( col%comm, ierr )
    call sll_test_mpi_error( ierr, &
         'sll_collective_allgather_int(): MPI_BARRIER()' )
    call MPI_ALLGATHER( send_buf(:), send_sz, MPI_INTEGER, &
                        recv_buf(:), recv_sz, MPI_INTEGER, col%comm, ierr )
    call sll_test_mpi_error( ierr, &
         'sll_collective_allgather_int(): MPI_ALLGATHER()' )
    call MPI_BARRIER( col%comm, ierr )
    call sll_test_mpi_error( ierr, &
         'sll_collective_allgather_int(): MPI_BARRIER()' )
  end subroutine

  subroutine sll_collective_allgatherv_real( col, send_buf, send_sz, &
       sizes, displs, rec_buf )
    type(sll_collective_t), pointer      :: col
    sll_real32, dimension(:), intent(in) :: send_buf ! what would change...
    sll_int32, intent(in)                :: send_sz
    sll_int32, dimension(:), intent(in)  :: sizes
    sll_int32, dimension(:), intent(in)  :: displs
    sll_real32, dimension(:), intent(in) :: rec_buf  ! would also change
    sll_int32                            :: ierr
    ! FIXME: argument checking
    call MPI_ALLGATHERV( send_buf, send_sz, MPI_REAL, rec_buf, sizes, &
         displs, MPI_REAL, col%comm, ierr )
    call sll_test_mpi_error( ierr, &
         'sll_collective_allgatherv_real(): MPI_ALLGATHERV()' )
  end subroutine sll_collective_allgatherv_real

  subroutine sll_collective_scatter_real( col, send_buf, size, root, &
       rec_buf )
    type(sll_collective_t), pointer      :: col
    sll_real32, dimension(:), intent(in) :: send_buf ! what would change...
    sll_int32, intent(in)                :: size
    sll_int32, intent(in)                :: root
    sll_real32, dimension(:), intent(in) :: rec_buf  ! would also change
    sll_int32                            :: ierr
    sll_int32                            :: recvcount
    ! FIXME: argument checking
    recvcount = size/col%size
    call MPI_SCATTER( send_buf, size, MPI_REAL, rec_buf, recvcount, &
         MPI_REAL, root, col%comm, ierr )
    call sll_test_mpi_error( ierr, &
         'sll_collective_scatter_real(): MPI_SCATTER()' )
  end subroutine sll_collective_scatter_real

  subroutine sll_collective_scatterv_real( col, send_buf, sizes, displs, &
       rec_szs, root,rec_buf )
    type(sll_collective_t), pointer      :: col
    sll_real32, dimension(:), intent(in) :: send_buf ! what would change...
    sll_int32, intent(in)                :: sizes
    sll_int32, dimension(:), intent(in)  :: displs
    sll_int32, dimension(:), intent(in)  :: rec_szs
    sll_real32, dimension(:), intent(in) :: rec_buf  ! would also change
    sll_int32, intent(in)                :: root
    sll_int32                            :: ierr
    ! FIXME: ARG CHECKING!
    call MPI_SCATTERV( send_buf, sizes, displs, MPI_REAL, rec_buf, &
         rec_szs, MPI_REAL, root, col%comm, ierr )
    call sll_test_mpi_error( ierr, &
         'sll_collective_scatterv_real(): MPI_SCATTERV()' )
  end subroutine sll_collective_scatterv_real

  subroutine sll_collective_allreduce_real( col, send_buf, count, op, &
       rec_buf )
    type(sll_collective_t), pointer       :: col
    sll_real32, dimension(:), intent(in)  :: send_buf ! what would change...
    sll_int32, intent(in)                 :: count
    sll_real32, intent(in)                :: op
    sll_real32, dimension(:), intent(out) :: rec_buf  ! would also change
    sll_int32                             :: ierr
    ! FIXME: ARG CHECKING!
    call MPI_ALLREDUCE( send_buf, rec_buf, count, MPI_REAL, op, &
         col%comm, ierr )
    call sll_test_mpi_error( ierr, &
         'sll_collective_allreduce_real(): MPI_ALLREDUCE()' )
  end subroutine sll_collective_allreduce_real

  subroutine sll_collective_allreduce_logical( col, send_buf, count, op, &
       rec_buf )
    type(sll_collective_t), pointer       :: col
    logical, dimension(:), intent(in)     :: send_buf ! what would change...
    sll_int32, intent(in)                 :: count
    sll_int32, intent(in)                 :: op
    logical, dimension(:), intent(out)    :: rec_buf  ! would also change
    sll_int32                             :: ierr
    ! FIXME: MORE ARG CHECKING!
    call sll_check_collective_ptr( col )
    call MPI_BARRIER( col%comm, ierr ) 
    call MPI_ALLREDUCE( send_buf, rec_buf, count, MPI_LOGICAL, op, &
         col%comm, ierr )
    call sll_test_mpi_error( ierr, &
         'sll_collective_allreduce_logical(): MPI_ALLREDUCE()' )
    call MPI_BARRIER( col%comm, ierr )
  end subroutine sll_collective_allreduce_logical

  subroutine sll_collective_reduce_real( col, send_buf, size, op, root_rank, &
       rec_buf )
    type(sll_collective_t), pointer      :: col
    sll_real32, dimension(:), intent(in) :: send_buf ! what would change...
    sll_real32, intent(in)               :: size
    sll_real32, intent(in)               :: op
    sll_int32, intent(in)                :: root_rank
    sll_real32, dimension(:), intent(in) :: rec_buf  ! would also change
    
    sll_int32                            :: ierr
    ! FIXME: ARG CHECKING!
    call MPI_REDUCE( send_buf, rec_buf, size, MPI_REAL, op, root_rank, &
         col%comm, ierr )
    call sll_test_mpi_error( ierr, &
         'sll_collective_reduce_real(): MPI_REDUCE()' )
  end subroutine sll_collective_reduce_real
  
  subroutine sll_collective_alltoall_int( send_buf, send_count, &
                                          recv_count, recv_buf, col )
    sll_int32, dimension(:), intent(in)  :: send_buf
    sll_int32, intent(in)                :: send_count
    sll_int32, intent(in)                :: recv_count
    sll_int32, dimension(:), intent(out) :: recv_buf
    type(sll_collective_t), pointer      :: col
    sll_int32                            :: ierr
    call sll_check_collective_ptr( col )
    ! FIXME: MORE ARG CHECKING
    call MPI_BARRIER( col%comm, ierr )
    call sll_test_mpi_error( ierr, &
         'sll_collective_alltoall_int(): MPI_BARRIER()' )
    call MPI_ALLTOALL( send_buf, send_count, MPI_INTEGER, &
                       recv_buf, recv_count, MPI_INTEGER, col%comm, ierr )
    call sll_test_mpi_error( ierr, &
         'sll_collective_alltoall_int(): MPI_ALLTOALLV()' )
    call MPI_BARRIER( col%comm, ierr )
    call sll_test_mpi_error( ierr, &
         'sll_collective_alltoall_int(): MPI_BARRIER()' )
  end subroutine sll_collective_alltoall_int

  ! Explore making this effectively typeless... need a macro implementation;
  ! pointer arguments won't work...
  subroutine sll_collective_alltoallV_real( send_buf, send_cnts, &
                                            send_displs, &
                                            recv_buf, recv_cnts, &
                                            recv_displs, col )
    sll_real32, dimension(:), intent(in) :: send_buf
    sll_int32,  dimension(:), intent(in) :: send_cnts
    sll_int32,  dimension(:), intent(in) :: send_displs
    sll_real32, dimension(:), intent(in) :: recv_buf
    sll_int32,  dimension(:), intent(in) :: recv_cnts
    sll_int32,  dimension(:), intent(in) :: recv_displs
    type(sll_collective_t), pointer      :: col
    sll_int32                            :: ierr
    ! FIXME: ARG CHECKING!
    call sll_check_collective_ptr( col )
    call MPI_BARRIER( col%comm, ierr )
    call sll_test_mpi_error( ierr, &
         'sll_collective_alltoallV_real(): MPI_BARRIER()' )
    call MPI_ALLTOALLV( send_buf, send_cnts, send_displs, MPI_REAL, &
                        recv_buf, recv_cnts, recv_displs, MPI_REAL, &
                        col%comm, ierr )
    call sll_test_mpi_error( ierr, &
         'sll_collective_alltoallV_real(): MPI_ALLTOALLV()' )
    call MPI_BARRIER( col%comm, ierr )
    call sll_test_mpi_error( ierr, &
         'sll_collective_alltoallV_real(): MPI_BARRIER()' )
  end subroutine sll_collective_alltoallV_real

  subroutine sll_collective_alltoallV_int( send_buf, send_cnts, &
                                             send_displs, &
                                             recv_buf, recv_cnts, &
                                             recv_displs, col )
    sll_int32, dimension(:), intent(in) :: send_buf
    sll_int32, dimension(:), intent(in) :: send_cnts
    sll_int32, dimension(:), intent(in) :: send_displs
    sll_int32, dimension(:), intent(in) :: recv_buf
    sll_int32, dimension(:), intent(in) :: recv_cnts
    sll_int32, dimension(:), intent(in) :: recv_displs
    type(sll_collective_t), pointer     :: col
    sll_int32                           :: ierr
    ! FIXME: ARG CHECKING!
    call sll_check_collective_ptr( col )
    call MPI_BARRIER( col%comm, ierr )
    call sll_test_mpi_error( ierr, &
         'sll_collective_alltoallV_int(): MPI_BARRIER()' )
    call MPI_ALLTOALLV( send_buf(:), send_cnts(:), send_displs(:), MPI_INTEGER,&
                        recv_buf(:), recv_cnts(:), recv_displs(:), MPI_INTEGER,&
                        col%comm, ierr )
    call sll_test_mpi_error( ierr, &
         'sll_collective_alltoallV_int(): MPI_ALLTOALLV()' )
    call MPI_BARRIER( col%comm, ierr )
   write (*,'(a, i4)') 'results from rank: ', sll_get_collective_rank(col)
 end subroutine sll_collective_alltoallV_int


end module sll_collective
