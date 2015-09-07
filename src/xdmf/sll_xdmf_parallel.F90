module sll_m_xdmf_parallel

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

  use sll_m_xml, only: &
    sll_t_xml_document,    &
    sll_t_xml_element

  implicit none
  private

!==============================================================================
! User-defined types
!==============================================================================

  !----------------------------------------------------------------------------
  !> Pointer to grid
  type :: t_grid_ptr
    type(sll_t_xml_element), pointer :: p
    sll_int32, allocatable           :: dims(:)
  end type t_grid_ptr

  !----------------------------------------------------------------------------
  !> XDMF file
  type, public :: sll_t_xdmf_file
    character(len=64)                :: name
    sll_int32                        :: id
    ! NEW:
    sll_real64                            :: time   =  0.0_f64
    type(sll_collective_t)  , pointer     :: comm   => null()
    type(sll_t_xml_element) , pointer     :: domain => null()
    type(sll_t_xml_document), allocatable :: xml_doc
    type(t_grid_ptr)        , allocatable :: grids(:)
  contains
    procedure         :: create           => xdmf__create
    procedure         :: new_grid         => xdmf__new_grid
    procedure         :: close_grid       => xdmf__close_grid
    procedure         :: write_array      => xdmf__write_array   ! not parallel
    procedure         :: collect_datasets => xdmf__collect_datasets
    procedure         :: close            => xdmf__close
    ! NEW:
    procedure :: init      => t_xdmf__init
    procedure :: write     => t_xdmf__write
    procedure :: delete    => t_xdmf__delete
    procedure :: add_grid  => t_xdmf__add_grid
    procedure :: add_field => t_xdmf__add_field

  end type sll_t_xdmf_file

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
contains
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

!==============================================================================
! NEW
!==============================================================================
  subroutine t_xdmf__init( self, time, comm )
    class(sll_t_xdmf_file), intent(  out) :: self
    sll_real64            , intent(in   ) :: time
    type(sll_collective_t), pointer       :: comm

    type(sll_t_xml_element), pointer :: root

    ! Fill in fields
    self%name =  "Not defined"  ! old
    self%id   =  999            ! old
    self%time =  time
    self%comm => comm

    !----------------------------------------------------
    ! ONLY MASTER WRITES TO FILE
    if (sll_get_collective_rank( self%comm ) /= 0) return
    !----------------------------------------------------

    ! Allocate XML document
    allocate( self%xml_doc )

    ! Add header
    call self%xml_doc%add_header_line( "<?xml version='1.0' ?>" )
    call self%xml_doc%add_header_line( "<!DOCTYPE Xdmf SYSTEM 'Xdmf.dtd' []>" )

    ! Create root element (only one may exist!)
    root => self%xml_doc%new_element( 'Xdmf' )
    call root%add_attribute( 'Version', '2.0' )

    ! Add Domain element to root (only one for our purposes)
    self%domain => root%new_element( 'Domain' )

  end subroutine t_xdmf__init

  !----------------------------------------------------------------------------
  subroutine t_xdmf__write( self, fname )
    class(sll_t_xdmf_file), intent(in) :: self
    character(len=*)      , intent(in) :: fname

    !----------------------------------------------------
    ! ONLY MASTER WRITES TO FILE
    if (sll_get_collective_rank( self%comm ) /= 0) return
    !----------------------------------------------------

    call self%xml_doc%write( fname )

  end subroutine t_xdmf__write

  !----------------------------------------------------------------------------
  subroutine t_xdmf__delete( self )
    class(sll_t_xdmf_file), intent(inout) :: self

    self%time   =  0.0_f64
    self%comm   => null()
    self%domain => null()

    if (allocated( self%xml_doc )) then
      call self%xml_doc%delete()
      deallocate( self%xml_doc )
    end if

    if (allocated( self%grids )) then
      deallocate( self%grids )
    end if

  end subroutine t_xdmf__delete

  !----------------------------------------------------------------------------
  subroutine t_xdmf__add_grid( self, grid_name, x1_path, x2_path, dims, gid )
    class(sll_t_xdmf_file), intent(inout) :: self
    character(len=*)      , intent(in   ) :: grid_name
    character(len=*)      , intent(in   ) :: x1_path
    character(len=*)      , intent(in   ) :: x2_path
    sll_int32             , intent(in   ) :: dims(2)
    sll_int32             , intent(  out) :: gid

    sll_int32                        :: ng
    character(len=64)                :: time_str
    character(len=:), allocatable    :: dims_str
    type(t_grid_ptr), allocatable    :: tmp(:)
    type(sll_t_xml_element), pointer :: grid, geometry
    type(sll_t_xml_element), pointer :: time, topology, dataitem

    !----------------------------------------------------
    ! ONLY MASTER WRITES TO FILE
    if (sll_get_collective_rank( self%comm ) /= 0) return
    !----------------------------------------------------

    ! Prepare strings with data
    write( time_str, '(f3.1)' ) self%time ! TODO: investigate format options
    call ints_to_string( dims, dims_str )

    ! Add new grid to domain
    grid => self%domain%new_element( 'Grid' )
    call grid%add_attribute( 'Name', trim( grid_name ) )
    call grid%add_attribute( 'GridType', 'Uniform' ) ! only option for now

    ! Add time to grid
    time => grid%new_element( 'Time' )
    call time%add_attribute( 'Value', trim( adjustl( time_str ) ) )

    ! Add topology to grid
    topology => grid%new_element( 'Topology' )
    call topology%add_attribute( 'TopologyType', '2DSMesh' ) ! only option now
    call topology%add_attribute( 'NumberOfElements', dims_str )
    
    ! Add geometry to grid
    geometry => grid%new_element( 'Geometry' )
    call geometry%add_attribute( 'GeometryType', 'X_Y' ) ! only option for now

    ! Add X axis to geometry
    dataitem => geometry%new_element( 'DataItem' )
    call dataitem%add_attribute( 'Dimensions', dims_str )
    call dataitem%add_attribute( 'NumberType', 'Float'  )
    call dataitem%add_attribute( 'Precision' , '8'      )
    call dataitem%add_attribute( 'Format'    , 'HDF'    ) ! only option for now
    call dataitem%add_chardata( trim( x1_path ) )        ! only option for now

    ! Add Y axis to geometry
    dataitem => geometry%new_element( 'DataItem' )
    call dataitem%add_attribute( 'Dimensions', dims_str )
    call dataitem%add_attribute( 'NumberType', 'Float'  )
    call dataitem%add_attribute( 'Precision' , '8'      )
    call dataitem%add_attribute( 'Format'    , 'HDF'    ) ! only option for now
    call dataitem%add_chardata( trim( x2_path ) )        ! only option for now

    ! Determine size of 'self%grids' array
    if (allocated( self%grids )) then
      ng = size( self%grids )
    else
      ng = 0
    end if

    ! Extend 'self%grids' array and store pointer to new grid
    allocate( tmp(ng+1) )
    tmp(1:ng) = self%grids(1:ng)
    tmp(ng+1)%p    => grid
    tmp(ng+1)%dims =  dims
    call move_alloc( from=tmp, to=self%grids )

    ! Return the grid id, i.e. the last index in the 'self%grids' array
    gid = ng+1

  end subroutine t_xdmf__add_grid

  !----------------------------------------------------------------------------
  subroutine t_xdmf__add_field( self, grid_id, field_name, field_path )
    class(sll_t_xdmf_file), intent(inout) :: self
    sll_int32             , intent(in   ) :: grid_id
    character(len=*)      , intent(in   ) :: field_name
    character(len=*)      , intent(in   ) :: field_path

    character(len=:), allocatable    :: dims_str
    type(sll_t_xml_element), pointer :: grid, field, dataitem

    !----------------------------------------------------
    ! ONLY MASTER WRITES TO FILE
    if (sll_get_collective_rank( self%comm ) /= 0) return
    !----------------------------------------------------

    ! Prepare strings with data
    call ints_to_string( self%grids( grid_id )%dims, dims_str )

    ! Create new field (scalar, nodal)
    field => self%grids( grid_id )%p%new_element( 'Attribute' )
    call field%add_attribute( 'Name'         , trim( field_name ) )
    call field%add_attribute( 'AttributeType', 'Scalar' ) ! only option for now
    call field%add_attribute( 'Center'       , 'Node'   ) ! only option for now

    ! Add new dataitem
    dataitem => field%new_element( 'DataItem' )
    call dataitem%add_attribute( 'Dimensions', dims_str )
    call dataitem%add_attribute( 'NumberType', 'Float'  )
    call dataitem%add_attribute( 'Precision' , '8'      )
    call dataitem%add_attribute( 'Format'    , 'HDF'    ) ! only option for now
    call dataitem%add_chardata( trim( adjustl( field_path ) ) )

  end subroutine t_xdmf__add_field

!==============================================================================
! UTILITIES
!==============================================================================

  subroutine ints_to_string( ints, str )
    sll_int32,                     intent(in   ) :: ints(:)
    character(len=:), allocatable, intent(  out) :: str

    sll_int32                      :: i, ni, nc, lc, str_len
    character(len=11), allocatable :: tmp(:)
    sll_int32        , allocatable :: ints_len(:)

    ! Allocate an homogeneous array of character
    ni = size( ints )
    allocate( tmp(ni) )
    allocate( ints_len(ni) )

    ! Write integers to an omogeneous array of character,
    ! as left-justified strings, and store length of each trimmed string
    do i = 1, ni
      write( tmp(i), '(i11)' ) ints(i)
      tmp(i)      = adjustl ( tmp(i) )
      ints_len(i) = len_trim( tmp(i), i32 )
    end do

    ! Allocate single string with minimum length
    str_len = sum( ints_len ) + ni - 1
    allocate( character(len=str_len) :: str )

    ! Write trimmed strings to single string, separated by blank space
    lc = 0
    do i = 1, ni
      nc = ints_len(i)
      str(lc+1:lc+nc) = trim( tmp(i) )
      lc = lc+nc
      if (i /= ni) then
        str(lc+1:lc+1) = ' '
        lc = lc+1
      end if
    end do

  end subroutine ints_to_string

  !----------------------------------------------------------------------------
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
    if (sll_get_collective_rank( comm ) /= 0) then
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

end module sll_m_xdmf_parallel
