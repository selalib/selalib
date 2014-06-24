!**************************************************************
!  Copyright INRIA
!  Authors : 
!     CALVI project team
!  
!  This code SeLaLib (for Semi-Lagrangian-Library) 
!  is a parallel library for simulating the plasma turbulence 
!  in a tokamak.
!  
!  This software is governed by the CeCILL-B license 
!  under French law and abiding by the rules of distribution 
!  of free software.  You can  use, modify and redistribute 
!  the software under the terms of the CeCILL-B license as 
!  circulated by CEA, CNRS and INRIA at the following URL
!  "http://www.cecill.info". 
!**************************************************************

!------------------------------------------------------------------------------
! Selalib
!------------------------------------------------------------------------------
! MODULE: sll_collective
!
! DESCRIPTION:
!> @file sll_collective.F90
!> @namespace sll_collective
!> @author Module Author Name and Affiliation
!> @brief Parallelizing facility.
!> @details Selalib applies the principle of modularization throughout all 
!> levels of abstraction of the library and aims at keeping third-party 
!> library modules as what they are: separate library modules. Therefore, 
!> in its current design, even a library like MPI has a single point of entry 
!> to Selalib. The collective communications module is such point of entry. We 
!> focus thus on the functionality offered by MPI, assign wrappers to its 
!> most desirable functionalities and write wrappers around them. These are 
!> the functions that are actually used throughout
!> the program.  This allows to adjust the exposed interfaces, do additional
!> error-checking and would even permit to completely change the means to
!> parallelize a code, by being able to replace MPI in a single file if this
!> were ever needed.
!>
!> \remark MPI provides five types of collective data movement routines that
!!         come in two variants: "simple" in which all communicated items are
!!         the same size and "vector" where each item can be a different size.
!!         The 'v' at the end indicates vector variant.
!>
!> <b> How to use the sll_collective module :</b> \n
!> *****************************************
!>
!> Include the line
!> \code use sll_collective \endcode
!> \warning Never put the line "use mpi"!
!>
!> Any use of the module's functionalities must be preceeded by calling
!> \code call sll_boot_collective() \endcode
!> and to "turn off" the parallel capabilities, one should finish by a call to:
!> \code call sll_halt_collective() \endcode
!> \warning This \a booting of the parallel environment needs to be done 
!> <b> ONLY ONCE </b> in a program.
!>
!> \n
!> \n
!>
!> <b> Comparaison with MPI module :</b> \n
!> *****************************************
!>
!> <table border="1">
!> <tr>
!> <th>sll_collective</th>
!> <th>MPI</th>
!> </tr>
!> <tr>
!> <td>sll_world_collective</td>
!> <td>MPI_COMM_WORLD</td>
!> </tr>
!> <tr>
!> <td>sll_boot_collective()</td>
!> <td>MPI_INIT(integer :: code)</td>
!> </tr>
!> <tr>
!> <td>sll_halt_collective()</td>
!> <td>MPI_FINALIZE(integer :: code)</td>
!> </tr>
!> <tr>
!> <td>sll_get_collective_rank(type(sll_collective_t), pointer :: col )</td>
!> <td>MPI_COMM_RANK(integer :: comm , integer :: rank , integer :: code)</td>
!> </tr>
!> <tr>
!> <td>sll_get_collective_size(type(sll_collective_t), pointer :: col )</td>
!> <td>MPI_COMM_SIZE(integer :: comm , integer :: size , integer :: code)</td>
!> </tr>
!> <tr>
!> <td>sll_get_collective_color(type(sll_collective_t), pointer :: col )</td>
!> <td></td>
!> </tr>
!> <tr>
!> <td>sll_get_collective_comm(type(sll_collective_t), pointer :: col )</td>
!> <td></td>
!> </tr>
!> <tr>
!> <td>sll_collective_barrier(type(sll_collective_t), pointer :: col )</td>
!> <td>MPI_BARRIER(integer :: comm ,integer :: code)</td>
!> </tr>
!> </table>
!>
!>
! REVISION HISTORY:
! DD Mmm YYYY - Initial Version
! TODO_dd_mmm_yyyy - TODO_describe_appropriate_changes - TODO_name
!------------------------------------------------------------------------------

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

  use mpi
  implicit none
  ! This is the only place in the prototype that should have to include
  ! the mpi header file.
  !include 'mpif.h'
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

  !> @brief Wrapper around the communicator
  type sll_collective_t
     type(sll_collective_t), pointer :: parent=>null() !< Pointer to parent communicator
     sll_int32                       :: comm   !< Communicator
     !> Control of subset assignment. Processes with the same color 
     !! are in the same new communicator
     sll_int32                       :: color 
     sll_int32                       :: key !< Control of rank assigment
     sll_int32                       :: rank !< Rank of the process
     sll_int32                       :: size !< Communicator size
     sll_int32                       :: thread_level_required
     sll_int32                       :: thread_level_provided
  end type sll_collective_t


  ! **********************************************************************
  ! A few rare global variables.
  ! The sll_world_collective just wraps around MPI_WORLD_COMM. This gets
  ! initialized in the call to sll_boot_collective().
  !
  ! **********************************************************************

  !> The Communicator (The same role as MPI_COMM_WORLD)
  type(sll_collective_t), pointer    :: sll_world_collective 


  !> @brief Broadcasts a message from the process with rank
  !>        "root" to all other processes of the communicator.
  interface sll_collective_bcast
     !> @brief Broadcasts a message of real type from the process with 
     !>        rank "root" to all other processes of the communicator.
     module procedure sll_collective_bcast_real
  end interface
  
  !> @brief Gathers together values from a group of processes.
  interface sll_collective_gather
     !> @brief Gathers together values of real type from a group of processes.
     module procedure sll_collective_gather_real32, &
          sll_collective_gather_real64
  end interface
  
  !> @brief Gathers data from all tasks and distribute the combined 
  !!        data to all tasks.
  interface sll_collective_allgather
     module procedure sll_collective_allgather_int, &
          sll_collective_allgather_real64
  end interface

  !> @brief Gathers data from all tasks and deliver the combined
  !!        data to all tasks 
  interface sll_collective_allgatherv
     module procedure sll_collective_allgatherv_real32, &
          sll_collective_allgatherv_real64
  end interface

  !> @brief Gathers into specified locations from all processes in a group.
  interface sll_collective_gatherv
     module procedure sll_collective_gatherv_real, &
          sll_collective_gatherv_real64
  end interface
  
  !> @brief Sends data from one process to all other processes
  !!        in a communicator.
  interface sll_collective_scatter
     module procedure sll_collective_scatter_real
  end interface

  !> @brief Scatters a buffer in parts to all processes in a communicator.
  interface sll_collective_scatterv
     module procedure sll_collective_scatterv_real
  end interface
  
  !> @brief Combines values from all processes and distributes
  !!        the result back to all processes. 
  interface sll_collective_allreduce
     module procedure sll_collective_allreduce_real32, &
                      sll_collective_allreduce_real64, &
                      sll_collective_allreduce_logical
  end interface
  
  !> @brief Reduces values on all processes to a single value.
  interface sll_collective_reduce
     module procedure sll_collective_reduce_real32, &
                      sll_collective_reduce_real64, &
                      sll_collective_reduce_int, &
                      sll_collective_reduce_logical
  end interface

  !> @brief Sends data from all to all processes.
  interface sll_collective_alltoall
     module procedure sll_collective_alltoall_int, &
          sll_collective_alltoall_double, &
          sll_collective_alltoall_complex_double
  end interface

  !> @brief Sends data from all to all processes; each process may send a
  !!        different amount of data and provide displacements for the
  !!        input and output data.
  interface sll_collective_alltoallV
     module procedure sll_collective_alltoallV_int, &
                      sll_collective_alltoallV_real, &
                      sll_collective_alltoallV_double, &
                      sll_collective_alltoallV_complex_double
  end interface




contains !************************** Operations **************************




  ! First a couple of utilities for this module
  !> @brief Checks if the pointer \a ptr is associated to an object.
  !> @param[in] ptr pointer to \a sll_collective_t type
  subroutine sll_check_collective_ptr( ptr )
    type(sll_collective_t), pointer :: ptr
    if( .not. associated(ptr)) then
       write (*, '(a)') 'sll_check_collective_ptr: non-associated pointer'
       stop 'SLL_CHECK_COLLECTIVE_PTR'
    end if
  end subroutine sll_check_collective_ptr

  !> @brief Checks the good execution of collective instruction.
  !> @details Many functions of collective library returns an error code
  !>          for verify the good execution.
  !>          With this function you can check the code and print the
  !>          message \a descriptor if you encountered an error.
  !> @param[in] ierr error's code 
  !> @param[in] descriptor message to print in error case
  subroutine sll_test_mpi_error( ierr, descriptor )
    sll_int32, intent(in)        :: ierr
    character(len=*), intent(in) :: descriptor

    if( ierr .ne. MPI_SUCCESS ) then
       write (*, '(a, a)') 'MPI error code failure: ', descriptor
    end if
  end subroutine sll_test_mpi_error

  ! Following what is exposed to the users of the module.

  ! The following function is a little tricky. The only test of equality we
  ! use is the identifier of the underlying communicator. Can the collective
  ! differ in any other way?
  function collectives_are_same( col1, col2 )
    logical :: collectives_are_same
    type(sll_collective_t), pointer :: col1
    type(sll_collective_t), pointer :: col2
    SLL_ASSERT( associated(col1) )
    SLL_ASSERT( associated(col2) )
    if( col1%comm == col2%comm ) then
       collectives_are_same = .true.
    else
       collectives_are_same = .false.
    end if
  end function collectives_are_same

  ! sll_boot_collective allocates and initializes the global variable 
  ! sll_world_collective and boots the MPI environment.

  !> @brief Starts the paralell environment 
  subroutine sll_boot_collective( )
    sll_int32 :: ierr

    SLL_ALLOCATE( sll_world_collective, ierr )


#ifndef MPI_THREAD_MULTIPLE

    call MPI_Init(ierr)
#else
    
    sll_world_collective%thread_level_required = MPI_THREAD_MULTIPLE
    
    call MPI_Init_Thread(sll_world_collective%thread_level_required, &
                         sll_world_collective%thread_level_provided, &
                         ierr)
#endif
    sll_world_collective%comm     = MPI_COMM_WORLD
    sll_world_collective%color    = 0
    sll_world_collective%key      = 0
    call MPI_COMM_RANK( MPI_COMM_WORLD, sll_world_collective%rank, ierr )
    call sll_test_mpi_error( ierr, 'sll_boot_collective(): MPI_COMM_RANK()' )
    call MPI_COMM_SIZE( MPI_COMM_WORLD, sll_world_collective%size, ierr )
    call sll_test_mpi_error( ierr, 'sll_boot_collective(): MPI_COMM_SIZE()')
  end subroutine sll_boot_collective

  !> @brief Ends the paralell environment 
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

  !> @brief Creates new (wrapper around) communicators based on colors and keys
  !> @param parent a pointer to the parent collective from which it was split
  !> @param[in] color the values of color used for the splitting
  !> @param[in] key the values of key used for the splitting
  !> @return a pointer to the new collective object
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

  !> @brief Marks the communicator object for deallocation
  !> @param col a pointer to the collective object
  subroutine sll_delete_collective( col )
    type(sll_collective_t), pointer :: col
    sll_int32                       :: ierr
    ! Why don't use MPI_COMM_FREE(col%comm,ierr)
    call sll_check_collective_ptr( col )
    SLL_DEALLOCATE( col, ierr )
  end subroutine sll_delete_collective

  !> @brief Gets the id (integer) of the communicator
  !> @param col Pointer to collective object
  function sll_get_collective_comm( col )
    type(sll_collective_t), pointer :: col
    sll_int32                       :: sll_get_collective_comm
    call sll_check_collective_ptr( col )
    sll_get_collective_comm = col%comm
  end function sll_get_collective_comm

 !> @brief Determines the rank of the calling process in the communicator
 !> @param col Collective object
 !> @return rank of the calling process in \a col
 function sll_get_collective_rank( col )
    type(sll_collective_t), pointer :: col
    sll_int32                       :: sll_get_collective_rank
    call sll_check_collective_ptr( col )
    sll_get_collective_rank = col%rank
  end function sll_get_collective_rank

 !> @brief Determines the color of the calling process in the communicator
 !> @param col Collective object
 !> @return color of the calling process in \a col
  function sll_get_collective_color( col )
    type(sll_collective_t), pointer :: col
    sll_int32                       :: sll_get_collective_color
    call sll_check_collective_ptr( col )
    sll_get_collective_color = col%color
  end function sll_get_collective_color

 !> @brief Determines the size of the group associated with a communicator 
 !> @param col Communicator
 !> @return size of the group associated with \a col
  function sll_get_collective_size( col )
    type(sll_collective_t), pointer :: col
    sll_int32                       :: sll_get_collective_size
    call sll_check_collective_ptr( col )
    sll_get_collective_size = col%size
  end function sll_get_collective_size
 
 !> @brief Gets collective parent 
 !> @param col Wrapper around the communicator
 !> @return pointer to collective parent
  function sll_get_collective_parent( col )
    type(sll_collective_t), pointer :: col
    type(sll_collective_t), pointer :: sll_get_collective_parent
    call sll_check_collective_ptr( col )
    sll_get_collective_parent => col%parent
  end function sll_get_collective_parent

  !> @brief Blocks until all processes in the communicator
  !>        have reached this routine.
  !> @param col Communicator
  subroutine sll_collective_barrier( col )
    type(sll_collective_t), pointer :: col
    sll_int32                         :: ierr
    call sll_check_collective_ptr( col )
    call MPI_BARRIER( col%comm, ierr )
    call sll_test_mpi_error( ierr, &
                               'sll_collective_barrier(): MPI_BARRIER()' )
  end subroutine sll_collective_barrier


  ! Start with something simple, like a buffer of 'real's...
  !> @brief Broadcasts a message from the process with rank "root"
  !>        to all other processes of the communicator
  !> @param col Wrapper around the communicator
  !> @param[in] buffer starting address of buffer
  !> @param[in] size number of entries in buffer
  !> @param[in] root rank of broadcast root
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
  !> @brief Gathers together real values from a group of processes
  !> @param[in] col Wrapper around the communicator
  !> @param[in] send_buf starting address of send buffer
  !> @param[in] send_sz number of elements in send buffer
  !> @param[in] root rank of broadcast root
  !> @param[out] rec_buf address of receive buffer
  subroutine sll_collective_gather_real32( col, send_buf, send_sz, root, &
       rec_buf )
    type(sll_collective_t), pointer      :: col
    sll_real32, dimension(:), intent(in) :: send_buf ! what would change...
    sll_int32                            :: send_sz
    sll_real32, dimension(:), intent(in) :: rec_buf  ! would also change
    sll_int32, intent(in)                :: root
    sll_int32                            :: ierr
    !sll_int32                            :: rec_count ! size of receive buf
    ! FIXME: add some argument checking here
    !rec_count = send_sz*col%size
    !Note that the 5th argument at the root indicates the number of items
    !it receives from each task. It is not the total number of items received.
    call MPI_GATHER( send_buf, send_sz, MPI_REAL, rec_buf, send_sz, &
         MPI_REAL, root, col%comm, ierr )
    call sll_test_mpi_error( ierr, &
         'sll_collective_gather_real(): MPI_GATHER()' )
  end subroutine sll_collective_gather_real32

  subroutine sll_collective_gather_real64( col, send_buf, send_sz, root, &
       rec_buf )
    type(sll_collective_t), pointer      :: col
    sll_real64, dimension(:), intent(in) :: send_buf ! what would change...
    sll_int32                            :: send_sz
    sll_real64, dimension(:), intent(in) :: rec_buf  ! would also change
    sll_int32, intent(in)                :: root
    sll_int32                            :: ierr
    !sll_int32                            :: rec_count ! size of receive buf
    ! FIXME: add some argument checking here
    !rec_count = send_sz*col%size
    !Note that the 5th argument at the root indicates the number of items
    !it receives from each task. It is not the total number of items received.
    call MPI_GATHER( send_buf, send_sz, MPI_DOUBLE_PRECISION, rec_buf, send_sz, &
         MPI_DOUBLE_PRECISION, root, col%comm, ierr )
    call sll_test_mpi_error( ierr, &
         'sll_collective_gather_real(): MPI_GATHER()' )
  end subroutine sll_collective_gather_real64

  !> @brief Gathers real values into specified locations from all processes in a group
  !> @param[in] col Wrapper around the communicator
  !> @param[in] send_buf starting address of send buffer
  !> @param[in] send_count number of elements in send buffer
  !> @param[in] recvcnts integer array (of length group size) containing 
  !!            the number of elements that are received from each process 
  !> @param[in] displs integer array. Entry i specifies the displacement 
  !!            relative to rec_buf at which to place the incoming data from process i
  !> @param[in] root rank of receiving process
  !> @param[out] rec_buf address of receive buffer
  subroutine sll_collective_gatherv_real( col, send_buf, send_count, &
       recvcnts, displs, root, rec_buf )
    type(sll_collective_t), pointer      :: col
    sll_real32, dimension(:), intent(in) :: send_buf ! what would change...
    sll_int32, intent(in)                :: send_count
    sll_int32, dimension(:), intent(in)  :: recvcnts
    sll_int32, dimension(:), intent(in)  :: displs
    sll_real32, dimension(:), intent(out) :: rec_buf  ! would also change
    sll_int32, intent(in)                :: root
    sll_int32                            :: ierr
    ! FIXME: Argument checking
    call sll_check_collective_ptr( col )
    ! displs, rec_buf and recvcnts significant only for root
    if (col%rank .eq. root) then 
      SLL_ASSERT( SIZE(recvcnts) .eq. col%size )
      SLL_ASSERT( SIZE(displs) .eq. col%size )
      SLL_ASSERT( SIZE(rec_buf) .eq. SUM(recvcnts) )
    endif
    call MPI_GATHERV( &
      send_buf, &
      send_count, &
      MPI_REAL, &
      rec_buf, &
      recvcnts, &
      displs, &
      MPI_REAL, &
      root, &
      col%comm, &
      ierr )
    call sll_test_mpi_error( ierr, &
         'sll_collective_gatherv_real(): MPI_GATHERV()' )
  end subroutine sll_collective_gatherv_real



  !> @brief Gathers real64 values into specified locations from all processes in a group
  !> @param[in] col Wrapper around the communicator
  !> @param[in] send_buf starting address of send buffer
  !> @param[in] send_count number of elements in send buffer
  !> @param[in] recvcnts integer array (of length group size) containing 
  !!            the number of elements that are received from each process 
  !> @param[in] displs integer array. Entry i specifies the displacement 
  !!            relative to rec_buf at which to place the incoming data from process i
  !> @param[in] root rank of receiving process
  !> @param[out] rec_buf address of receive buffer
  subroutine sll_collective_gatherv_real64( col, send_buf, send_count, &
       recvcnts, displs, root, rec_buf )
    type(sll_collective_t), pointer      :: col
    sll_real64, dimension(:), intent(in) :: send_buf ! what would change...
    sll_int32, intent(in)                :: send_count
    sll_int32, dimension(:), intent(in)  :: recvcnts
    sll_int32, dimension(:), intent(in)  :: displs
    sll_real64, dimension(:), intent(out) :: rec_buf  ! would also change
    sll_int32, intent(in)                :: root
    sll_int32                            :: ierr
    ! FIXME: Argument checking
    call sll_check_collective_ptr( col )
    ! displs, rec_buf and recvcnts significant only for root
    if (col%rank .eq. root) then 
      SLL_ASSERT( SIZE(recvcnts) .eq. col%size )
      SLL_ASSERT( SIZE(displs) .eq. col%size )
      SLL_ASSERT( SIZE(rec_buf) .eq. SUM(recvcnts) )
    endif
    call MPI_GATHERV( send_buf, send_count, MPI_REAL8,rec_buf,recvcnts,&
         displs, MPI_REAL8, root, col%comm, ierr )
    call sll_test_mpi_error( ierr, &
         'sll_collective_gatherv_real64(): MPI_GATHERV()' )
  end subroutine sll_collective_gatherv_real64








  !> @brief Gathers integer data from all tasks and 
  !!        distribute the combined data to all tasks
  !> @param[in] col Wrapper around the communicator
  !> @param[in] send_buf starting address of send buffer
  !> @param[in] send_sz number of elements in send buffer
  !> @param[out] rec_buf address of receive buffer
  !> @param[in] recv_sz number of elements received from any process
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

  subroutine sll_collective_allgather_real64( &
    col, &
    send_buf, &
    send_sz, &
    recv_buf, &
    recv_sz )

    type(sll_collective_t), pointer         :: col
    sll_real64, dimension(:), intent(in)    :: send_buf 
    sll_int32, intent(in)                   :: send_sz
    sll_real64, dimension(:), intent(out)   :: recv_buf ! would change
    sll_int32, dimension(:), intent(in)     :: recv_sz  
    sll_int32                               :: ierr
    ! FIXME: Argument checking
    call sll_check_collective_ptr( col )
    call MPI_BARRIER( col%comm, ierr )
    call sll_test_mpi_error( ierr, &
         'sll_collective_allgather_int(): MPI_BARRIER()' )
    call MPI_ALLGATHER( send_buf(:), send_sz, MPI_DOUBLE_PRECISION, &
                        recv_buf(:), recv_sz, MPI_DOUBLE_PRECISION, &
                        col%comm, ierr )
    call sll_test_mpi_error( ierr, &
         'sll_collective_allgather_int(): MPI_ALLGATHER()' )
    call MPI_BARRIER( col%comm, ierr )
    call sll_test_mpi_error( ierr, &
         'sll_collective_allgather_int(): MPI_BARRIER()' )
  end subroutine sll_collective_allgather_real64


  !> @brief Gathers real data from all tasks and 
  !!        deliver the combined data to all tasks
  !> @param[in] col Wrapper around the communicator
  !> @param[in] send_buf starting address of send buffer
  !> @param[in] send_cnt number of elements in send buffer
  !> @param[in] displs integer array. Entry i specifies the displacement
  !> @param[in] rec_cnt integer array containing the number of elements 
  !!            that are to be received from each process
  !> @param[out] rec_buf address of receive buffer
  subroutine sll_collective_allgatherv_real32( col, send_buf, send_cnt, &
       rec_cnt, displs, rec_buf )
    type(sll_collective_t), pointer      :: col
    sll_real32, dimension(:), intent(in) :: send_buf ! what would change...
    sll_int32, intent(in)                :: send_cnt
    sll_int32, dimension(:), intent(in)  :: rec_cnt
    sll_int32, dimension(:), intent(in)  :: displs
    sll_real32, dimension(:), intent(out) :: rec_buf  ! would also change
    sll_int32                            :: ierr
    ! FIXME: argument checking
    SLL_ASSERT(col%size .eq. SIZE(displs))
    SLL_ASSERT(col%size .eq. SIZE(rec_cnt))
    call MPI_ALLGATHERV( send_buf, send_cnt, MPI_REAL, rec_buf,rec_cnt,&
         displs, MPI_REAL, col%comm, ierr )
    call sll_test_mpi_error( ierr, &
         'sll_collective_allgatherv_real32(): MPI_ALLGATHERV()' )
  end subroutine sll_collective_allgatherv_real32


  subroutine sll_collective_allgatherv_real64( &
    col, &
    send_buf, &
    send_cnt, &
    rec_cnt, &
    displs, &
    rec_buf )

    type(sll_collective_t), pointer       :: col
    sll_real64, dimension(:), intent(in)  :: send_buf ! what would change...
    sll_int32, intent(in)                 :: send_cnt
    sll_int32, dimension(:), intent(in)   :: rec_cnt
    sll_int32, dimension(:), intent(in)   :: displs
    sll_real64, dimension(:), intent(out) :: rec_buf  ! would also change
    sll_int32                             :: ierr
    ! FIXME: argument checking
    SLL_ASSERT(col%size .eq. SIZE(displs))
    SLL_ASSERT(col%size .eq. SIZE(rec_cnt))
    call MPI_ALLGATHERV( &
         send_buf, &
         send_cnt, &
         MPI_DOUBLE_PRECISION, &
         rec_buf, &
         rec_cnt, &
         displs, &
         MPI_DOUBLE_PRECISION, &
         col%comm, &
         ierr )
    call sll_test_mpi_error( ierr, &
         'sll_collective_allgatherv_real64(): MPI_ALLGATHERV()' )
  end subroutine sll_collective_allgatherv_real64

  !> @brief Sends data from one process to all other processes
  !>        in a communicator 
  !> @param col Wrapper around the communicator
  !> @param[in] send_buf address of send buffer
  !> @param[in] send_count number of elements sent to each process
  !> @param[in] root rank of sending process
  !> @param[in] rec_buf address of receive buffer
  subroutine sll_collective_scatter_real( col, send_buf, send_count, root, &
       rec_buf )
    type(sll_collective_t), pointer      :: col
    sll_real32, dimension(:), intent(in) :: send_buf ! what would change...
    sll_int32, intent(in)                :: send_count
    sll_int32, intent(in)                :: root
    sll_real32, dimension(:), intent(in) :: rec_buf  ! would also change
    sll_int32                            :: ierr
    !sll_int32                            :: recvcount
    ! FIXME: argument checking
    ! recvcount = size/col%size
    !send_buf and send_count significant only at root
    if(col%rank .eq. root) then
      SLL_ASSERT( SIZE(send_buf) .eq. (send_count*(col%size)) )
    endif
    call MPI_SCATTER( send_buf, send_count, MPI_REAL, rec_buf, send_count, &
         MPI_REAL, root, col%comm, ierr )
    call sll_test_mpi_error( ierr, &
         'sll_collective_scatter_real(): MPI_SCATTER()' )
  end subroutine sll_collective_scatter_real

  !> @brief Scatters a buffer in parts to all processes in a communicator
  !> @param[in] col wrapper around the communicator
  !> @param[in] send_buf starting address of send buffer
  !> @param[in] send_count integer array (of length group size) specifying
  !>            the number of elements to send to each processor
  !> @param[in] displs integer array (of length group size). Entry i 
  !>            specifies the displacement (relative to sendbuf) from 
  !>            which to take the outgoing data to process i
  !> @param[in] recv_count number of elements in receive buffer
  !> @param[in] root rank of sending process
  !> @param[out] rec_buf starting address of receive buffer
  subroutine sll_collective_scatterv_real( col, send_buf, send_count,&
                                        displs,recv_count, root,rec_buf)
    type(sll_collective_t), pointer      :: col
    sll_real32, dimension(:), intent(in) :: send_buf ! what would change...
    sll_int32, dimension(:), intent(in)  :: send_count
    sll_int32, dimension(:), intent(in)  :: displs
    sll_int32, intent(in)                :: recv_count
    sll_real32, dimension(:), intent(out) :: rec_buf  ! would also change
    sll_int32, intent(in)                :: root
    sll_int32                            :: ierr
    ! FIXME: ARG CHECKING!
    SLL_ASSERT( SIZE(send_count) .eq. col%size )
    SLL_ASSERT( SIZE(displs) .eq. col%size )
    call MPI_SCATTERV( send_buf, send_count, displs, MPI_REAL, rec_buf,&
         recv_count, MPI_REAL, root, col%comm, ierr )
    call sll_test_mpi_error( ierr, &
         'sll_collective_scatterv_real(): MPI_SCATTERV()' )
  end subroutine sll_collective_scatterv_real

  !> @brief Combines real values from all processes and 
  !!        distributes the result back to all processes
  !> @param[in] col wrapper around the communicator
  !> @param[in] send_buf starting address of send buffer
  !> @param[in] count number of elements in send buffer
  !> @param[in] op operation
  !> @param[out] rec_buf starting address of receive buffer
  subroutine sll_collective_allreduce_real32( col, send_buf, count, op, &
       rec_buf )
    type(sll_collective_t), pointer       :: col
    sll_real32, dimension(:), intent(in)  :: send_buf ! what would change...
    sll_int32, intent(in)                 :: count
    sll_int32, intent(in)                :: op
    sll_real32, dimension(:), intent(out) :: rec_buf  ! would also change
    sll_int32                             :: ierr
    ! FIXME: ARG CHECKING!
    call sll_check_collective_ptr( col )
    call MPI_BARRIER( col%comm, ierr ) 
    call MPI_ALLREDUCE( &
      send_buf, &
      rec_buf, &
      count, &
      MPI_REAL, &
      op, &
      col%comm, &
      ierr )
    call sll_test_mpi_error( &
      ierr, &
      'sll_collective_allreduce_real32(): MPI_ALLREDUCE()' )
  end subroutine sll_collective_allreduce_real32


  !> @brief Combines real values from all processes and 
  !!        distributes the result back to all processes
  !> @param[in] col wrapper around the communicator
  !> @param[in] send_buf starting address of send buffer
  !> @param[in] count number of elements in send buffer
  !> @param[in] op operation
  !> @param[out] rec_buf starting address of receive buffer
  subroutine sll_collective_allreduce_real64( col, send_buf, count, op, &
       rec_buf )
    type(sll_collective_t), pointer       :: col
    sll_real64, dimension(:), intent(in)  :: send_buf ! what would change...
    sll_int32, intent(in)                 :: count
    sll_int32, intent(in)                :: op
    sll_real64, dimension(:), intent(out) :: rec_buf  ! would also change
    sll_int32                             :: ierr
    ! FIXME: ARG CHECKING!
    call sll_check_collective_ptr( col )
    call MPI_BARRIER( col%comm, ierr ) 
    call MPI_ALLREDUCE( &
      send_buf, &
      rec_buf, &
      count, &
      MPI_DOUBLE_PRECISION, &
      op, &
      col%comm, &
      ierr )
    call sll_test_mpi_error( &
      ierr, &
      'sll_collective_allreduce_real64(): MPI_ALLREDUCE()' )
  end subroutine sll_collective_allreduce_real64



  !> @brief Combines logical values from all processes and 
  !!        distributes the result back to all processes
  !> @param[in] col wrapper around the communicator
  !> @param[in] send_buf starting address of send buffer
  !> @param[in] count number of elements in send buffer
  !> @param[in] op operation
  !> @param[out] rec_buf starting address of receive buffer
  subroutine sll_collective_allreduce_logical( col, send_buf, count,op,&
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

  !> @brief Reduces integer values on all processes to a single value
  !> @param[in] col wrapper around the communicator
  !> @param[in] send_bu address of send buffer
  !> @param[in] size number of elements in send buffer
  !> @param[in] op reduce operation
  !> @param[in] root_rank rank of root process
  !> @param[out] rec_buf address of receive buffer
  subroutine sll_collective_reduce_int( col, send_buf, size, op, root_rank, &
       rec_buf )
    type(sll_collective_t), pointer      :: col
    sll_int32, dimension(:), intent(in) :: send_buf ! what would change...
    sll_int32, intent(in)               :: size
    sll_int32, intent(in)               :: op
    sll_int32, intent(in)                :: root_rank
    sll_int32, dimension(:), intent(in) :: rec_buf  ! would also change
    
    sll_int32                            :: ierr
    ! FIXME: ARG CHECKING!
    call MPI_REDUCE( send_buf, rec_buf, size, MPI_INTEGER, op, root_rank, &
         col%comm, ierr )
    call sll_test_mpi_error( ierr, &
         'sll_collective_reduce_real(): MPI_REDUCE()' )
  end subroutine sll_collective_reduce_int

  !> @brief Reduces real values on all processes to a single value
  !> @param[in] col wrapper around the communicator
  !> @param[in] send_bu address of send buffer
  !> @param[in] size number of elements in send buffer
  !> @param[in] op reduce operation
  !> @param[in] root_rank rank of root process
  !> @param[out] rec_buf address of receive buffer
  subroutine sll_collective_reduce_real32( col, send_buf, size, op, root_rank, &
       rec_buf )
    type(sll_collective_t), pointer      :: col
    sll_real32, dimension(:), intent(in) :: send_buf ! what would change...
    sll_int32, intent(in)               :: size
    sll_int32, intent(in)               :: op
    sll_int32, intent(in)                :: root_rank
    sll_real32, dimension(:), intent(in) :: rec_buf  ! would also change
    
    sll_int32                            :: ierr
    ! FIXME: ARG CHECKING!
    call MPI_REDUCE( send_buf, rec_buf, size, MPI_REAL, op, root_rank, &
         col%comm, ierr )
    call sll_test_mpi_error( ierr, &
         'sll_collective_reduce_real(): MPI_REDUCE()' )
  end subroutine sll_collective_reduce_real32

  subroutine sll_collective_reduce_real64( col, send_buf, size, op, root_rank, &
       rec_buf )
    type(sll_collective_t), pointer       :: col
    sll_real64, dimension(:), intent(in)  :: send_buf 
    sll_int32, intent(in)                 :: size
    sll_int32, intent(in)                 :: op
    sll_int32, intent(in)                 :: root_rank
    sll_real64, dimension(:), intent(out) :: rec_buf 
    sll_int32                             :: ierr

    ! FIXME: ARG CHECKING!
    call MPI_REDUCE( send_buf, rec_buf, size, MPI_DOUBLE_PRECISION, op, &
         root_rank, col%comm, ierr )
    call sll_test_mpi_error( ierr, &
         'sll_collective_reduce_real(): MPI_REDUCE()' )
  end subroutine sll_collective_reduce_real64

  
  !> @brief Reduces logical values on all processes to a single value
  !> @param[in] col wrapper around the communicator
  !> @param[in] send_bu address of send buffer
  !> @param[in] size number of elements in send buffer
  !> @param[in] op reduce operation
  !> @param[in] root_rank rank of root process
  !> @param[out] rec_buf address of receive buffer
  subroutine sll_collective_reduce_logical( col, send_buf, size, op, root_rank, &
       rec_buf )
    type(sll_collective_t), pointer      :: col
    LOGICAL, DIMENSION(:), intent(in) :: send_buf ! what would change...
    sll_int32, intent(in)               :: size
    sll_int32, intent(in)               :: op
    sll_int32, intent(in)                :: root_rank
    LOGICAL, DIMENSION(:), intent(in) :: rec_buf  ! would also change
    
    sll_int32                            :: ierr
    ! FIXME: ARG CHECKING!
    call MPI_REDUCE( send_buf, rec_buf, size, MPI_LOGICAL, op, root_rank, &
         col%comm, ierr )
    call sll_test_mpi_error( ierr, &
         'sll_collective_reduce_real(): MPI_REDUCE()' )
  end subroutine sll_collective_reduce_logical

  !> @brief Sends integer data from all to all processes
  !> @param[in] send_buf starting address of send buffer
  !> @param[in] send_count number of elements to send to each process
  !> @param[in] recv_count number of elements received from any process
  !> @param[out] recv_buf address of receive buffer
  !> @param[in] col wrapper around the communicator
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

  subroutine sll_collective_alltoall_double( send_buf, send_count, &
                                             recv_count, recv_buf, col )
    sll_real64, dimension(:), intent(in)  :: send_buf
    sll_int32, intent(in)                 :: send_count
    sll_int32, intent(in)                 :: recv_count
    sll_real64, dimension(:), intent(out) :: recv_buf
    type(sll_collective_t), pointer       :: col
    sll_int32                             :: ierr
    call sll_check_collective_ptr( col )
    ! FIXME: MORE ARG CHECKING
    call MPI_BARRIER( col%comm, ierr )
    call sll_test_mpi_error( ierr, &
         'sll_collective_alltoall_double(): MPI_BARRIER()' )
    call MPI_ALLTOALL( send_buf, send_count, MPI_DOUBLE_PRECISION, &
                       recv_buf, recv_count, MPI_DOUBLE_PRECISION, &
                       col%comm, ierr )
    call sll_test_mpi_error( ierr, &
         'sll_collective_alltoall_double(): MPI_ALLTOALL()' )
    call MPI_BARRIER( col%comm, ierr )
    call sll_test_mpi_error( ierr, &
         'sll_collective_alltoall_double(): MPI_BARRIER()' )
  end subroutine sll_collective_alltoall_double

  subroutine sll_collective_alltoall_complex_double( send_buf, send_count, &
                                                     recv_count, recv_buf, col )
    sll_comp64, dimension(:), intent(in)  :: send_buf
    sll_int32, intent(in)                 :: send_count
    sll_int32, intent(in)                 :: recv_count
    sll_comp64, dimension(:), intent(out) :: recv_buf
    type(sll_collective_t), pointer       :: col
    sll_int32                             :: ierr
    call sll_check_collective_ptr( col )
    ! FIXME: MORE ARG CHECKING
    call MPI_BARRIER( col%comm, ierr )
    call sll_test_mpi_error( ierr, &
         'sll_collective_alltoall_complex_double(): MPI_BARRIER()' )
    call MPI_ALLTOALL( send_buf, send_count, MPI_DOUBLE_COMPLEX, &
                       recv_buf, recv_count, MPI_DOUBLE_COMPLEX, &
                       col%comm, ierr )
    call sll_test_mpi_error( ierr, &
         'sll_collective_alltoall_complex_double(): MPI_ALLTOALL()' )
    call MPI_BARRIER( col%comm, ierr )
    call sll_test_mpi_error( ierr, &
         'sll_collective_alltoall_complex_double(): MPI_BARRIER()' )
  end subroutine sll_collective_alltoall_complex_double


  ! Explore making this effectively typeless... need a macro implementation;
  ! pointer arguments won't work...
  !> @brief Sends real data from all to all processes; each process may send a
  !!         different amount of data and provide displacements for the
  !!         input and output data.
  !> @param[in] send_buf starting address of send buffer
  !> @param[in] send_cnts integer array equal to the group size
  !!                      specifying the number of elements to 
  !!                      send to each processor
  !> @param[in] send_displs integer array (of length group size). Entry j 
  !!                        specifies the displacement (relative to send_buf)
  !!                        from which to take the outgoing data destined for 
  !!                        process j
  !> @param[out] recv_buf address of receive buffer
  !> @param[in] recv_cnts integer array equal to the group size specifying the 
  !!                      maximum number of elements that can be received from
  !!                      each processor
  !> @param[in] recv_displs integer array (of length group size). Entry i 
  !!                      specifies the displacement (relative to recvbuf 
  !!                      at which to place the incoming data from process i
  !> @param[in] col wrapper around the communicator
  subroutine sll_collective_alltoallV_real( send_buf, send_cnts, &
                                            send_displs, &
                                            recv_buf, recv_cnts, &
                                            recv_displs, col )
    sll_real32, dimension(:), intent(in) :: send_buf
    sll_int32,  dimension(:), intent(in) :: send_cnts
    sll_int32,  dimension(:), intent(in) :: send_displs
    sll_real32, dimension(:), intent(out) :: recv_buf
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

  subroutine sll_collective_alltoallV_double( send_buf, send_cnts, &
                                              send_displs, &
                                              recv_buf, recv_cnts, &
                                              recv_displs, col )
    sll_real64, dimension(:), intent(in)  :: send_buf
    sll_int32,  dimension(:), intent(in)  :: send_cnts
    sll_int32,  dimension(:), intent(in)  :: send_displs
    sll_real64, dimension(:), intent(out) :: recv_buf
    sll_int32,  dimension(:), intent(in)  :: recv_cnts
    sll_int32,  dimension(:), intent(in)  :: recv_displs
    type(sll_collective_t), pointer       :: col
    sll_int32                             :: ierr
    ! FIXME: ARG CHECKING!
    call sll_check_collective_ptr( col )
    call MPI_BARRIER( col%comm, ierr )
    call sll_test_mpi_error( ierr, &
         'sll_collective_alltoallV_double(): MPI_BARRIER()' )
    call MPI_ALLTOALLV( send_buf, send_cnts, send_displs, &
                        MPI_DOUBLE_PRECISION, &
                        recv_buf, recv_cnts, recv_displs, &
                        MPI_DOUBLE_PRECISION, &
                        col%comm, ierr )
    call sll_test_mpi_error( ierr, &
         'sll_collective_alltoallV_double(): MPI_ALLTOALLV()' )
    call MPI_BARRIER( col%comm, ierr )
    call sll_test_mpi_error( ierr, &
         'sll_collective_alltoallV_double(): MPI_BARRIER()' )
  end subroutine sll_collective_alltoallV_double

  subroutine sll_collective_alltoallV_complex_double( send_buf, send_cnts, &
                                                      send_displs, &
                                                      recv_buf, recv_cnts, &
                                                      recv_displs, col )
    sll_comp64, dimension(:), intent(in)  :: send_buf
    sll_int32,  dimension(:), intent(in)  :: send_cnts
    sll_int32,  dimension(:), intent(in)  :: send_displs
    sll_comp64, dimension(:), intent(out) :: recv_buf
    sll_int32,  dimension(:), intent(in)  :: recv_cnts
    sll_int32,  dimension(:), intent(in)  :: recv_displs
    type(sll_collective_t), pointer       :: col
    sll_int32                             :: ierr
    ! FIXME: ARG CHECKING!
    call sll_check_collective_ptr( col )
    call MPI_BARRIER( col%comm, ierr )
    call sll_test_mpi_error( ierr, &
         'sll_collective_alltoallV_complex_double(): MPI_BARRIER()' )
    call MPI_ALLTOALLV( send_buf, send_cnts, send_displs, &
                        MPI_DOUBLE_COMPLEX, &
                        recv_buf, recv_cnts, recv_displs, &
                        MPI_DOUBLE_COMPLEX, &
                        col%comm, ierr )
    call sll_test_mpi_error( ierr, &
         'sll_collective_alltoallV_complex_double(): MPI_ALLTOALLV()' )
    call MPI_BARRIER( col%comm, ierr )
    call sll_test_mpi_error( ierr, &
         'sll_collective_alltoallV_complex_double(): MPI_BARRIER()' )
  end subroutine sll_collective_alltoallV_complex_double



  !> @brief Sends integer data from all to all processes; each process may 
  !!         send a different amount of data and provide displacements for the
  !!         input and output data.
  !> @param[in] send_buf starting address of send buffer
  !> @param[in] send_cnts integer array equal to the group size
  !!                      specifying the number of elements to 
  !!                      send to each processor
  !> @param[in] send_displs integer array (of length group size). Entry j 
  !!                        specifies the displacement (relative to send_buf)
  !!                        from which to take the outgoing data destined for 
  !!                        process j
  !> @param[out] recv_buf address of receive buffer
  !> @param[in] recv_cnts integer array equal to the group size specifying the 
  !!                      maximum number of elements that can be received from
  !!                      each processor
  !> @param[in] recv_displs integer array (of length group size). Entry i 
  !!                      specifies the displacement (relative to recvbuf 
  !!                      at which to place the incoming data from process i
  !> @param[in] col wrapper around the communicator
  subroutine sll_collective_alltoallV_int( send_buf, send_cnts, &
                                           send_displs, &
                                           recv_buf, recv_cnts, &
                                           recv_displs, col )
    sll_int32, dimension(:), intent(in) :: send_buf
    sll_int32, dimension(:), intent(in) :: send_cnts
    sll_int32, dimension(:), intent(in) :: send_displs
    sll_int32, dimension(:), intent(out) :: recv_buf
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
 end subroutine sll_collective_alltoallV_int

!
!FOR SOME REASON THIS ARGUMENT OF THIS FUNCTION IS DECLARED ALLOCATBLE AND 
!POSING PROBLEM WITH F95 STANDARD. WE SHOULD COME BACK TO THIS
!
!warning sll_collective_alltoallV_int_simple is not fixed

  ! This toutine is a simpler version of the sll_collective_alltoallV_int subroutine
  subroutine sll_collective_alltoallV_int_simple( send_buf, send_cnts, &
       recv_buf,col )
    sll_int32, dimension(:), intent(in) :: send_buf
    sll_int32, dimension(:), intent(in) :: send_cnts
    !sll_int32, ,allocatable dimension(:) :: recv_buf
    sll_int32, pointer, dimension(:) :: recv_buf
    sll_int32, allocatable, dimension(:) :: send_displs
    sll_int32, allocatable, dimension(:) :: recv_cnts
    sll_int32, allocatable, dimension(:) :: recv_displs
    type(sll_collective_t), pointer     :: col
    sll_int32                          :: ierr,size_comm,i,sendcnts_size
    sll_int32                          :: dum
    sendcnts_size = size(send_cnts)
    
    ! FIXME: ARG CHECKING!
    call sll_check_collective_ptr( col )
    call MPI_BARRIER( col%comm, ierr )
    call sll_test_mpi_error( ierr, &
         'sll_collective_alltoallV_int(): MPI_BARRIER()' )
    
    size_comm = sll_get_collective_size(col)
    SLL_ASSERT( sendcnts_size .eq. size_comm )
    
    ! Define RECV_CNTS
    SLL_ALLOCATE(recv_cnts(size_comm),ierr)
    call sll_collective_alltoall_int( send_cnts ,1 ,1, &
         recv_cnts, col)
    
    ! Define RECV_BUF
    dum = SUM(recv_cnts)
    SLL_ALLOCATE(recv_buf(dum),ierr)
    
    ! Define SEND_DISPLS
    SLL_ALLOCATE(send_displs(size_comm),ierr)
    send_displs(1)=0
    do i=2,size_comm
       send_displs(i)=send_displs(i-1)+send_cnts(i-1)
    enddo
    
    ! Define RECV_DISPLS
    SLL_ALLOCATE(recv_displs(size_comm),ierr)
    recv_displs(1)=0
    do i=2,size_comm
       recv_displs(i)=recv_displs(i-1)+recv_cnts(i-1)
    enddo
    
    call MPI_ALLTOALLV( send_buf(:), send_cnts(:), send_displs(:), MPI_INTEGER,&
                       recv_buf(:), recv_cnts(:), recv_displs(:), MPI_INTEGER,&
                       col%comm, ierr )
    call sll_test_mpi_error( ierr, &
         'sll_collective_alltoallV_int(): MPI_ALLTOALLV()' )
    call MPI_BARRIER( col%comm, ierr )
    SLL_DEALLOCATE_ARRAY(recv_displs,ierr)
    SLL_DEALLOCATE_ARRAY(send_displs,ierr)
    SLL_DEALLOCATE_ARRAY(recv_cnts,ierr)
  end subroutine sll_collective_alltoallV_int_simple

  ! Explore if the Irecv calls can be made into collective calls in this module
  !  subroutine sll_collective_Irecv( )

end module sll_collective
