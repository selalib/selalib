!------------------------------------------------------------------------------
! SELALIB
!------------------------------------------------------------------------------
! MODULE:   sll_m_xdmf_light_parallel
!
! DESCRIPTION:
!> @ingroup xdmf_io_parallel
!> @authors Yaman Güçlü - <yaman.guclu@gmail.com>
!> @brief   Construct the XML component of an XDMF database (parallel).
!> @todo    Add detailed description
!------------------------------------------------------------------------------
module sll_m_xdmf_light_parallel

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_working_precision.h"

  use sll_m_collective, only: &
    sll_o_collective_bcast, &
    sll_o_collective_gather, &
    sll_t_collective_t, &
    sll_f_get_collective_rank, &
    sll_f_get_collective_size

  use sll_m_xdmf_light_serial, only: &
    sll_t_xdmf_file

  use sll_mpi, only: &
    mpi_any_source, &
    mpi_any_tag, &
    mpi_character, &
    mpi_integer, &
    mpi_recv, &
    mpi_send, &
    mpi_source, &
    mpi_status_size

  implicit none

  public :: &
    sll_t_xdmf_parallel_file

  private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!==============================================================================

  !> Maximum length of variable-length strings to be passed through MPI
  integer, parameter :: maxlen = 256

  !----------------------------------------------------------------------------
  !> XDMF parallel file
  type :: sll_t_xdmf_parallel_file

    type(sll_t_collective_t), pointer    :: comm => null()
    sll_int32                          :: rank = -1
    type(sll_t_xdmf_file), allocatable :: xdmf_file

  contains
    procedure :: init      => t_xdmf_parallel__init
    procedure :: write     => t_xdmf_parallel__write
    procedure :: delete    => t_xdmf_parallel__delete
    procedure :: add_grid  => t_xdmf_parallel__add_grid
    procedure :: add_field => t_xdmf_parallel__add_field

  end type sll_t_xdmf_parallel_file

!==============================================================================
contains
!==============================================================================

  subroutine t_xdmf_parallel__init( self, time, comm )
    class(sll_t_xdmf_parallel_file), intent(  out) :: self
    sll_real64                     , intent(in   ) :: time
    type(sll_t_collective_t)         , pointer       :: comm

    ! Fill in fields
    self%comm => comm
    self%rank =  sll_f_get_collective_rank( comm )
    
    ! MASTER: Allocate and initialize sequential XDMF file
    if (self%rank == 0) then
      allocate( self%xdmf_file )
      call self%xdmf_file%init( time )
    end if

  end subroutine t_xdmf_parallel__init

  !----------------------------------------------------------------------------
  subroutine t_xdmf_parallel__write( self, fname )
    class(sll_t_xdmf_parallel_file), intent(in) :: self
    character(len=*)               , intent(in) :: fname

    ! MASTER: write XDMF file 
    if (self%rank == 0) then
      call self%xdmf_file%write( fname )
    end if

  end subroutine t_xdmf_parallel__write

  !----------------------------------------------------------------------------
  subroutine t_xdmf_parallel__delete( self )
    class(sll_t_xdmf_parallel_file), intent(inout) :: self

    ! MASTER: delete and deallocate sequential XDMF file
    if (self%rank == 0) then
      if (allocated( self%xdmf_file )) then
        call self%xdmf_file%delete()
        deallocate( self%xdmf_file )
      end if
    end if

    ! Reset internal fields
    self%comm => null()
    self%rank = -1

  end subroutine t_xdmf_parallel__delete

  !----------------------------------------------------------------------------
  subroutine t_xdmf_parallel__add_grid( self, grid_name, x1_path, x2_path, dims, gid )
    class(sll_t_xdmf_parallel_file), intent(inout) :: self
    character(len=*)               , intent(in   ) :: grid_name
    character(len=*)               , intent(in   ) :: x1_path
    character(len=*)               , intent(in   ) :: x2_path
    sll_int32                      , intent(in   ) :: dims(2)
    sll_int32                      , intent(  out) :: gid

    sll_int32 :: buf_gid(1)

    ! MASTER: add grid to sequential XDMF file
    ! TODO: other processes might want to add a grid
    if (self%rank == 0) then
      call self%xdmf_file%add_grid( grid_name, x1_path, x2_path, dims, gid )
      buf_gid(1) = gid
    end if

    ! Broadcast grid ID to all other processes
    call sll_o_collective_bcast( self%comm, buf_gid, 0 )
    gid = buf_gid(1)

  end subroutine t_xdmf_parallel__add_grid

  !----------------------------------------------------------------------------
  subroutine t_xdmf_parallel__add_field( self, &
      grid_id, field_name, field_path, &
      to_file )

    class(sll_t_xdmf_parallel_file), intent(inout) :: self
    sll_int32                      , intent(in   ) :: grid_id
    character(len=*)               , intent(in   ) :: field_name
    character(len=*)               , intent(in   ) :: field_path
    logical                        , intent(in   ) :: to_file

    sll_int32             :: i         ! index for cycles
    sll_int32             :: nprocs    ! number of processors
    sll_int32             :: comm      ! MPI communicator
    sll_int32             :: ierr      ! MPI error
    logical               :: buf(1)    ! send buffer for gather
    logical, allocatable  :: recbuf(:) ! receive buffer for gather
    character(len=11)     :: rank_str  ! rank converted to string
    sll_int32             :: stat(mpi_status_size)
    sll_int32             :: sender_rank

    ! Send/receive buffers (receive for master, send for all others)
    sll_int32             :: buf_gid  ! grid_id
    character(len=maxlen) :: buf_fn   ! field_name
    character(len=maxlen) :: buf_fp   ! field_path

    ! Some info about communicator
    nprocs = sll_f_get_collective_size( self%comm )
    comm   = self%comm%comm

    ! Fill in send buffer (all processes)
    buf(1) = to_file

    ! MASTER ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    if (self%rank == 0) then
      ! Allocate receive buffer for 'to_file' logicals, then gather all of them
      allocate( recbuf(0:nprocs-1) )
      call sll_o_collective_gather( self%comm, buf, 0, recbuf )
      !print *, "PROC #0: gather"

      ! Write array info on XDMF file
      if (to_file) then
        call self%xdmf_file%add_field( grid_id, field_name, field_path )
        !print *, "PROC #0: write array info to xmf"
      end if
      !print *, "PROC #0: recbuf = ", recbuf
      if ( nprocs>1 )then
        ! Receive data from other processes, and write array info to XDMF file
        do i = 1, count( recbuf(1:) )

          ! Receive grid ID from any process, then field name/path from same proc

          !    MPI_RECV( BUF, COUNT, DATATYPE, SOURCE, TAG, COMM, STATUS, IERROR )
          call mpi_recv( buf_gid, 1, mpi_integer, mpi_any_source, &
            mpi_any_tag, comm, stat, ierr )

          sender_rank = stat( mpi_source )

          call mpi_recv( buf_fn, maxlen, mpi_character, sender_rank, &
            mpi_any_tag, comm, stat, ierr )

          call mpi_recv( buf_fp, maxlen, mpi_character, sender_rank, &
            mpi_any_tag, comm, stat, ierr )

          write( rank_str, '(i8)' ) sender_rank
          !print *, "PROC #0: data received from processor "// &
          !  trim( adjustl( rank_str ) )

          ! Add field to sequential XDMF file
          call self%xdmf_file%add_field( buf_gid,trim( buf_fn ),trim( buf_fp ) )
          !print *, "PROC #0: field added to XDMF file"

        end do
      endif
    ! NOT MASTER ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    else
      allocate( recbuf(0) ) ! allocate empty buffer (ifort complains otherwise)
      call sll_o_collective_gather( self%comm, buf, 0, recbuf )
      write( rank_str, '(i8)' ) self%rank
      !print *, "PROC #"//trim( adjustl( rank_str ) )//": gather"

      if (to_file) then
        ! Fill in send buffers
        buf_gid = grid_id
        buf_fn  = adjustl( field_name )
        buf_fp  = adjustl( field_path )

        ! Send data to master
        !    MPI_SEND( BUF, COUNT, DATATYPE, DEST, TAG, COMM, IERROR )
        call mpi_send( buf_gid, 1     , mpi_integer  , 0, 9, comm, ierr )
        call mpi_send( buf_fn , maxlen, mpi_character, 0, 9, comm, ierr )
        call mpi_send( buf_fp , maxlen, mpi_character, 0, 9, comm, ierr )
        !print *, "PROC #"//trim( adjustl( rank_str ) )//": send"
      end if

    end if

  end subroutine t_xdmf_parallel__add_field

end module sll_m_xdmf_light_parallel
