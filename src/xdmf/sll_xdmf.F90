module sll_m_xdmf

#include "sll_working_precision.h"
#include "sll_assert.h"

  use sll_collective, only:  &
    sll_collective_t,        &
    sll_get_collective_size, &
    sll_get_collective_rank, &
    sll_collective_bcast,    &
    sll_collective_gather

  use sll_xml_io, only:    &
    sll_xml_file_create,   &
    sll_xml_file_close,    &
    sll_xml_dataitem

!==============================================================================
! XDMF FILE
!==============================================================================
  type :: sll_t_xdmf_file
    character(len=64)                :: name
    sll_int32                        :: id
    type(sll_collective_t), pointer  :: comm => null()
  contains
    procedure         :: create           => xdmf__create
    procedure         :: new_grid         => xdmf__new_grid
    procedure         :: close_grid       => xdmf__close_grid
    procedure         :: write_array      => xdmf__write_array   ! not parallel
    procedure         :: collect_datasets => xdmf__collect_datasets
    procedure         :: close            => xdmf__close
  end type sll_t_xdmf_file

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
contains
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

!==============================================================================
! UTILITIES
!==============================================================================
  subroutine split( path, head, tail )
    character(len=*), intent(in   ) :: path
    character(len=*), intent(  out) :: head
    character(len=*), intent(  out) :: tail

    sll_int32 :: nc
    sll_int32 :: i

    ! Number of non-blank characters in path string
    nc = len_trim( path, i32 )

    ! If last character is '/', tail is empty
    if (path(nc:nc) == '/') then
      head = path(1:nc)
      tail = ''
    end if

    ! Search backwards (from right to left) for '/' character, and split path
    do i = nc-1,1,-1
      if (path(i:i) == '/') then
        head = path(1:i) 
        tail = path(i+1:nc)
        return
      end if
    end do

    ! If no '/' character was found, head is empty
    head = ''
    tail = path(1:nc)

  end subroutine split

!==============================================================================
! XDMF FILE
!==============================================================================
  subroutine xdmf__create( self, xml_file_name, comm )
    class(sll_t_xdmf_file), intent(  out) :: self
    character(len=*)      , intent(in   ) :: xml_file_name
    type(sll_collective_t), pointer       :: comm

    sll_int32 :: xml_file_id, error
    sll_int32 :: bcast_buf(1)

    ! Master creates XML file and writes header
    if (sll_get_collective_rank( comm ) == 0) then
      call sll_xml_file_create( xml_file_name, xml_file_id, error )
      SLL_ASSERT( error == 0 )
      bcast_buf(1) = xml_file_id
    end if

    ! Master broadcasts file ID to all other processes
    call sll_collective_bcast( comm, bcast_buf, 0 )

    ! Non-master procs read file ID from buffer
    if (sll_get_collective_rank( comm ) == 0) then
      xml_file_id = bcast_buf(1)
    end if

    ! Fill in fields
    self%name =  xml_file_name
    self%id   =  xml_file_id
    self%comm => comm

  end subroutine xdmf__create

!------------------------------------------------------------------------------
  subroutine xdmf__new_grid( self, grid_name, x1_path, x2_path, dims, time )
    class(sll_t_xdmf_file), intent(in) :: self
    character(len=*)      , intent(in) :: grid_name
    character(len=*)      , intent(in) :: x1_path
    character(len=*)      , intent(in) :: x2_path
    sll_int32             , intent(in) :: dims(2)
    sll_real64, optional  , intent(in) :: time

    !----------------------------------------------------
    ! ONLY MASTER WRITES TO FILE
    if (sll_get_collective_rank( self%comm ) /= 0) return
    !----------------------------------------------------

    ! Open 'Grid' block
    write( self%id,"(a)" ) &
        "<Grid Name='"//trim( adjustl( grid_name ) )//"' GridType='Uniform'>"

    ! If time was provided, add 'Time' block
    if (present( time )) then
      write( self%id, "(a13,g15.3,a3)" ) "<Time Value='",time,"'/>"
    end if

    ! Add 'Topology' block
    write( self%id,"(a,2i5,a)" ) &
        "<Topology TopologyType='2DSMesh' NumberOfElements='", &
        dims(2), dims(1), "'/>"

    ! Add 'Geometry' block with two dataitems
    write( self%id,"(a)" ) "<Geometry GeometryType='X_Y'>"
    call sll_xml_dataitem( self%id, trim( adjustl( x1_path ) ), &
      dims(1), dims(2), 'HDF' )
    call sll_xml_dataitem( self%id, trim( adjustl( x2_path ) ), &
      dims(1), dims(2), 'HDF' )
    write( self%id,"(a)" ) "</Geometry>"

  end subroutine xdmf__new_grid

!------------------------------------------------------------------------------
  subroutine xdmf__close_grid( self )
    class(sll_t_xdmf_file), intent(in) :: self

    !----------------------------------------------------
    ! ONLY MASTER WRITES TO FILE
    if (sll_get_collective_rank( self%comm ) /= 0) return
    !----------------------------------------------------

    ! Close grid element
    write( self%id, "(a)" ) "</Grid>"

  end subroutine xdmf__close_grid

!------------------------------------------------------------------------------
  subroutine xdmf__write_array( self, path, dims, filetype )
    class(sll_t_xdmf_file), intent(in) :: self
    character(len=*)      , intent(in) :: path
    sll_int32             , intent(in) :: dims(2)
    character(len=*)      , intent(in) :: filetype

    character(len=256) :: group
    character(len=64)  :: dataset

    !----------------------------------------------------
    ! ONLY MASTER WRITES TO FILE
    if (sll_get_collective_rank( self%comm ) /= 0) return
    !----------------------------------------------------

    SLL_ASSERT( filetype == 'HDF' .or. filetype == 'Binary' )
    call split( path, group, dataset )

    write( self%id,"(a)" ) &
    "<Attribute Name='"//trim(dataset)// &
    "' AttributeType='Scalar' Center='Node'>"
    call sll_xml_dataitem( self%id, path, dims(1), dims(2), filetype )
    write( self%id,"(a)" ) "</Attribute>"

  end subroutine xdmf__write_array

!------------------------------------------------------------------------------
  subroutine xdmf__close( self )
    class(sll_t_xdmf_file), intent(in) :: self

    sll_int32 :: error

    !----------------------------------------------------
    ! ONLY MASTER WRITES TO FILE
    if (sll_get_collective_rank( self%comm ) /= 0) return
    !----------------------------------------------------

    call sll_xml_file_close( self%id, error )
    SLL_ASSERT( error == 0 )

  end subroutine xdmf__close

!------------------------------------------------------------------------------
  subroutine xdmf__collect_datasets( self, &
    dataset, dims, filetype, &
    to_file )

    use mpi

    class(sll_t_xdmf_file), intent(in) :: self
    character(len=*)      , intent(in) :: dataset
    sll_int32             , intent(in) :: dims(2)
    character(len=*)      , intent(in) :: filetype
    logical               , intent(in) :: to_file

    logical              :: buf(1)    ! buffer for send (all processes)
    logical, allocatable :: recbuf(:) ! buffer for receive (only master)
    sll_int32            :: i, comm, np, rank, ierr
    !sll_int32            :: mpi_status(mpi_status_size) !PN change name because
    ! there is a conflict when we use MPICH as MPI library
    sll_int32            :: ompi_status(mpi_status_size)
    character(len=8)     :: rank_str

    ! Some info about communicator
    comm = self%comm%comm
    np   = sll_get_collective_size( self%comm )
    rank = sll_get_collective_rank( self%comm )

    ! Fill in send buffer (all processes)
    buf(1) = to_file

    ! MASTER
    if (rank==0) then
      allocate( recbuf(0:np-1) )
!      recbuf(0) = to_file
!      !    MPI_GATHER( SENDBUF, SENDCOUNT, SENDTYPE, RECVBUF, RECVCOUNT, RECVTYPE, ROOT, COMM, IERROR )
!      call mpi_gather( mpi_in_place, 1, mpi_logical, recbuf, 1, mpi_logical, 0, comm, ierr )
      call sll_collective_gather( self%comm, buf, 0, recbuf )
      print *, "PROC #0: gather"
      ! Write array info on XDMF file
      if (to_file) then
        call self%write_array( dataset, dims, filetype )
        print *, "PROC #0: write array info to xmf"
      end if

      print *, "PROC #0: recbuf = ", recbuf
      do i=1,count( recbuf(1:) )
        !    MPI_RECV( BUF, COUNT, DATATYPE, SOURCE, TAG, COMM, STATUS, IERROR )
        call mpi_recv( dataset, len( dataset ), mpi_character, mpi_any_source, mpi_any_tag, comm, ompi_status, ierr )
        print *, "PROC #0: receive data from processor ", ompi_status( mpi_source )
        call self%write_array( dataset, dims, filetype )
      end do

    ! NOT MASTER
    else
!      allocate( recbuf(1) )
!      !    MPI_GATHER( SENDBUF, SENDCOUNT, SENDTYPE, RECVBUF, RECVCOUNT, RECVTYPE, ROOT, COMM, IERROR )
!      call mpi_gather( buf, 1, mpi_logical, recbuf, 0, 0, 0, comm, ierr )
      call sll_collective_gather( self%comm, buf, 0, recbuf )
      write( rank_str, '(i8)' ) rank
      print *, "PROC #"//trim( adjustl( rank_str ) )//": gather"
      if (buf(1)) then
        !    MPI_SEND( BUF, COUNT, DATATYPE, DEST, TAG, COMM, IERROR )
        call mpi_send( dataset, len( dataset ), mpi_character, 0, 9, comm, ierr )
        print *, "PROC #"//trim( adjustl( rank_str ) )//": send"
      end if

    end if

  end subroutine xdmf__collect_datasets

!==============================================================================

end module sll_m_xdmf
