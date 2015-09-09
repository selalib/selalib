module sll_m_xdmf

#include "sll_working_precision.h"

  use sll_m_xml, only: &
    sll_t_xml_document,    &
    sll_t_xml_element

  implicit none
  private

  !----------------------------------------------------------------------------
  !> Pointer to grid
  type :: t_xdmf_grid
    type(sll_t_xml_element), pointer :: xml_grid
    sll_int32, allocatable           :: dims(:)
  end type t_xdmf_grid

  !----------------------------------------------------------------------------
  !> XDMF file, sequential
  type, public :: sll_t_xdmf_file

    sll_real64                            :: time =  0.0_f64
    type(sll_t_xml_document), allocatable :: xml_doc
    type(sll_t_xml_element) , pointer     :: xml_domain => null()
    type(t_xdmf_grid)       , allocatable :: grids(:)

  contains

    procedure :: init      => t_xdmf__init
    procedure :: write     => t_xdmf__write
    procedure :: delete    => t_xdmf__delete
    procedure :: add_grid  => t_xdmf__add_grid
    procedure :: add_field => t_xdmf__add_field

  end type sll_t_xdmf_file

!==============================================================================
contains
!==============================================================================

  !----------------------------------------------------------------------------
  !> Initialize XDMF file (set time, allocate storage, store pointer to domain)
  subroutine t_xdmf__init( self, time )
    class(sll_t_xdmf_file), intent(  out) :: self
    sll_real64            , intent(in   ) :: time

    type(sll_t_xml_element), pointer :: root

    ! Fill in fields
    self%time =  time

    ! Allocate XML document
    allocate( self%xml_doc )

    ! Add header
    call self%xml_doc%add_header_line( "<?xml version='1.0' ?>" )
    call self%xml_doc%add_header_line( "<!DOCTYPE Xdmf SYSTEM 'Xdmf.dtd' []>" )

    ! Create root element (only one may exist!)
    root => self%xml_doc%new_element( 'Xdmf' )
    call root%add_attribute( 'Version', '2.0' )

    ! Add Domain element to root (only one for our purposes)
    self%xml_domain => root%new_element( 'Domain' )

  end subroutine t_xdmf__init

  !----------------------------------------------------------------------------
  !> Write XML file
  subroutine t_xdmf__write( self, fname )
    class(sll_t_xdmf_file), intent(in) :: self
    character(len=*)      , intent(in) :: fname

    call self%xml_doc%write( fname )

  end subroutine t_xdmf__write

  !----------------------------------------------------------------------------
  !> Delete XDMF file (deallocate storage, nullify pointers)
  subroutine t_xdmf__delete( self )
    class(sll_t_xdmf_file), intent(inout) :: self

    self%time       =  0.0_f64
    self%xml_domain => null()

    if (allocated( self%xml_doc )) then
      call self%xml_doc%delete()
      deallocate( self%xml_doc )
    end if

    if (allocated( self%grids )) then
      deallocate( self%grids )
    end if

  end subroutine t_xdmf__delete

  !----------------------------------------------------------------------------
  !> Add grid to file (new grid ID is returned)
  subroutine t_xdmf__add_grid( self, grid_name, x1_path, x2_path, dims, gid )
    class(sll_t_xdmf_file), intent(inout) :: self
    character(len=*)      , intent(in   ) :: grid_name
    character(len=*)      , intent(in   ) :: x1_path
    character(len=*)      , intent(in   ) :: x2_path
    sll_int32             , intent(in   ) :: dims(2)
    sll_int32             , intent(  out) :: gid

    sll_int32                        :: ng
    character(len=64)                :: time_str
    character(len= :), allocatable   :: dims_str
    type(t_xdmf_grid), allocatable   :: tmp(:)
    type(sll_t_xml_element), pointer :: grid, geometry
    type(sll_t_xml_element), pointer :: time, topology, dataitem

    ! Prepare strings with data
    write( time_str, '(f3.1)' ) self%time ! TODO: investigate format options
    call ints_to_string( dims, dims_str )

    ! Add new grid to domain
    grid => self%xml_domain%new_element( 'Grid' )
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
    tmp(ng+1)%xml_grid => grid
    tmp(ng+1)%dims     =  dims
    call move_alloc( from=tmp, to=self%grids )

    ! Return the grid id, i.e. the last index in the 'self%grids' array
    gid = ng+1

  end subroutine t_xdmf__add_grid

  !----------------------------------------------------------------------------
  !> Add field to grid (selected by its grid ID)
  subroutine t_xdmf__add_field( self, grid_id, field_name, field_path )
    class(sll_t_xdmf_file), intent(inout) :: self
    sll_int32             , intent(in   ) :: grid_id
    character(len=*)      , intent(in   ) :: field_name
    character(len=*)      , intent(in   ) :: field_path

    character(len=:), allocatable    :: dims_str
    type(sll_t_xml_element), pointer :: grid, field, dataitem

    ! Prepare strings with data
    call ints_to_string( self%grids( grid_id )%dims, dims_str )

    ! Create new field (scalar, nodal)
    field => self%grids( grid_id )%xml_grid%new_element( 'Attribute' )
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

  !----------------------------------------------------------------------------
  !> Write an array of integers to a single string:
  !>   . Numbers are separated by a blank space;  
  !>   . String is allocated here with minimum length.
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
  !> Split path into head (directory) and tail (file)
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

end module sll_m_xdmf
