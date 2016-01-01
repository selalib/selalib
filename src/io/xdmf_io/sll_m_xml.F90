!------------------------------------------------------------------------------
! SELALIB
!------------------------------------------------------------------------------
! MODULE:   sll_m_xml
!
! DESCRIPTION:
!> @ingroup xdmf_io
!> @authors Yaman Güçlü    - <yaman.guclu@gmail.com>
!> @authors Marco Restelli - <marco.restelli@gmail.com>
!> @brief   Facilities for constructing an XML tree and printing it to file.
!> @todo    Add detailed description
!------------------------------------------------------------------------------
module sll_m_xml

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  implicit none

  public :: &
    sll_t_xml_document, &
    sll_t_xml_element

  private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 
!==============================================================================

  character, parameter :: nl = achar(10)  ! newline character
  integer  , parameter :: maxlen = 256    ! max string length
 
  !============================================================================
  ! USER-DEFINED TYPES
  !============================================================================
 
  !----------------------------------------------------------------------------
  !> Base class for all the XML entities which can appear in content
  type, abstract :: c_xml_item
  contains
   procedure(i_xml_item__write ), deferred :: write
   procedure(i_xml_item__delete), deferred :: delete
  end type c_xml_item
 
  !----------------------------------------------------------------------------
  !> XML type: linked list of XML entities, used in element content
  type, extends(c_xml_item) :: t_xml_content
    class(c_xml_item)   , pointer     :: item
    type (t_xml_content), pointer     :: next => null()
  contains
    procedure :: write       => t_xml_content__write
    procedure :: delete      => t_xml_content__delete
    procedure :: new_content => t_xml_content__new_content
  end type
  
  !----------------------------------------------------------------------------
  !> XML attribute type
  type :: t_xml_attribute
!    character(len=:), allocatable :: name
    character(len=maxlen) :: name
!   Notice: it could be useful introducing a type for the attvalue,
!   representing the formal definition
!   AttValue ::=  '"' ([^<&"] | Reference)* '"'
!            |  "'" ([^<&'] | Reference)* "'"
!    character(len=:), allocatable :: attvalue
    character(len=maxlen) :: attvalue
  contains
    procedure :: to_string => attribute_name_val_string
  end type t_xml_attribute

  !----------------------------------------------------------------------------
  !> XML element type
  type, extends(c_xml_item) :: sll_t_xml_element
!    character(len=:)     ,  allocatable :: name
    character(len=maxlen)              :: name
    type(t_xml_attribute), allocatable :: attributes(:)
    type(t_xml_content)  , allocatable :: content
  contains
    procedure :: write                => t_xml_element__write
    procedure :: delete               => t_xml_element__delete
    procedure :: new_element          => t_xml_element__new_element
    procedure :: add_attribute        => t_xml_element__add_attribute
    procedure :: add_chardata_string  => t_xml_element__add_chardata_string
    procedure :: add_chardata_printer => t_xml_element__add_chardata_printer
    generic   :: add_chardata   => add_chardata_string, add_chardata_printer
  end type sll_t_xml_element

  !----------------------------------------------------------------------------
  !> XML document type
  type :: sll_t_xml_document
    character(len=maxlen)  , allocatable :: header_lines(:)
    type(sll_t_xml_element), allocatable :: root
  contains
    procedure :: add_header_line => t_xml_document__add_header_line
    procedure :: new_element     => t_xml_document__new_element
    procedure :: write           => t_xml_document__write
    procedure :: delete          => t_xml_document__delete
  end type sll_t_xml_document

  !----------------------------------------------------------------------------
  !> XML abstract class: generic printer for writing chardata to file
  type, abstract :: c_text_data_printer
  contains
    procedure(i_print_text    ), deferred, pass(self) :: print_text
    procedure(i_delete_printer), deferred, pass(self) :: delete
  end type c_text_data_printer
 
  !----------------------------------------------------------------------------
  !> XML type: default printer for writing chardata
  type, extends(c_text_data_printer) :: t_default_text_data_printer
!    character(len=:), allocatable :: text
    character(len=maxlen) :: text
  contains
    procedure, pass(self) :: print_text => default_print_text
    procedure, pass(self) :: delete     => default_delete
  end type t_default_text_data_printer

  !----------------------------------------------------------------------------
  !> XML type: chardata
  type, extends(c_xml_item) :: t_xml_chardata
    class(c_text_data_printer), allocatable :: chardata
  contains
    procedure :: write  => t_xml_chardata__write
    procedure :: delete => t_xml_chardata__delete
  end type t_xml_chardata

  !============================================================================
  ! ABSTRACT INTERFACES
  !============================================================================

  !----------------------------------------------------------------------------
  !> Write XML content to file
  abstract interface
   subroutine i_xml_item__write( self, indent, fid )
    import :: c_xml_item
    implicit none
    class(c_xml_item), intent(in) :: self
    integer          , intent(in) :: indent
    integer          , intent(in) :: fid
   end subroutine i_xml_item__write
  end interface
 
  !----------------------------------------------------------------------------
  !> Delete XML item (deallocate everything)
  abstract interface
   subroutine i_xml_item__delete( self )
    import :: c_xml_item
    implicit none
    class(c_xml_item), intent(inout) :: self
   end subroutine i_xml_item__delete
  end interface
 
  !----------------------------------------------------------------------------
  !> Write chardata to file
  abstract interface
   subroutine i_print_text( self, indent, fid )
    import :: c_text_data_printer
    implicit none
    class(c_text_data_printer), intent(in) :: self
    integer                   , intent(in) :: indent
    integer                   , intent(in) :: fid
   end subroutine i_print_text
  end interface
 
  !----------------------------------------------------------------------------
  !> Delete chardata printer
  abstract interface
   subroutine i_delete_printer( self )
    import :: c_text_data_printer
    implicit none
    class(c_text_data_printer), intent(inout) :: self
   end subroutine i_delete_printer
  end interface

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
contains
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

!------------------------------------------------------------------------------
  recursive subroutine t_xml_content__write( self, indent, fid )
    class(t_xml_content), intent(in) :: self
    integer             , intent(in) :: indent
    integer             , intent(in) :: fid

    if (.not. associated ( self%item )) then
      ! TODO: Give error, item should always exist
    end if

    ! 1) Write local item
    call self%item%write( indent, fid )

    ! 2) Write following container: recursive step
    if (associated( self%next )) call self%next%write( indent, fid )

  end subroutine t_xml_content__write

!------------------------------------------------------------------------------
  recursive subroutine t_xml_content__delete( self )
    class(t_xml_content), intent(inout) :: self

    if (.not. associated ( self%item )) then
      ! TODO: Give error, item should always exist
    end if
      
    ! 1) Delete (and deallocate) following container: recursive step
    if (associated( self%next )) then
      call self%next%delete()
      deallocate( self%next )
    end if

    ! 2) Delete (and deallocate) local item
    call self%item%delete()
    deallocate( self%item )

  end subroutine t_xml_content__delete

!------------------------------------------------------------------------------
  recursive function t_xml_content__new_content( self ) result(new_cont)
    class(t_xml_content), intent(inout) :: self

    type(t_xml_content), pointer :: new_cont

    if (.not. associated ( self%item )) then
      ! TODO: Give error, item should always exist
    end if

    if (associated ( self%next )) then
      new_cont => self%next%new_content()
    else
      allocate( t_xml_content :: self%next )
      new_cont => self%next
    end if

  end function t_xml_content__new_content

!------------------------------------------------------------------------------
  subroutine t_xml_document__add_header_line( self, line )
    class(sll_t_xml_document), intent(inout) :: self
    character(len=*)         , intent(in   ) :: line

    integer :: nl
    character(len=maxlen), allocatable :: tmp(:)

    if (allocated( self%header_lines )) then
      nl = size( self%header_lines )
      allocate( tmp(nl+1) )
      tmp(1:nl) = self%header_lines(1:nl)
    else
      nl = 0
      allocate( tmp(1) )
    end if

    tmp(nl+1) = trim( line )

    call move_alloc( from=tmp, to=self%header_lines )

  end subroutine t_xml_document__add_header_line

!------------------------------------------------------------------------------
  function t_xml_document__new_element( self, name ) result( new_root )
    class(sll_t_xml_document), target, intent(inout) :: self
    character(len=*)                 , intent(in   ) :: name

    type(sll_t_xml_element), pointer :: new_root

    if (allocated( self%root )) then
      ! TODO: give error!!! Only one root element may exist
    else
      allocate( self%root )
      self%root%name = trim( name )
    end if

    new_root => self%root

  end function t_xml_document__new_element

!------------------------------------------------------------------------------
  subroutine t_xml_document__write( self, fname )
    class(sll_t_xml_document), intent(in) :: self
    character(len=*)         , intent(in) :: fname

    integer :: fid, i

    ! Create a new file using default properties
    open( file=fname, status='replace', form='formatted', newunit=fid )

    ! Print header lines
    if (allocated( self%header_lines )) then
      do i = lbound( self%header_lines, 1 ), ubound( self%header_lines, 1 )
        write(fid,'(a)') trim( self%header_lines(i) )
      end do
    end if

    ! Print root element, recursively
    if (allocated( self%root )) then
      call self%root%write( 0, fid )
    end if

    ! Close file 
    close( fid )
  
  end subroutine t_xml_document__write

!------------------------------------------------------------------------------
  subroutine t_xml_document__delete( self )
    class(sll_t_xml_document), intent(inout) :: self

    ! Remove header lines
    if (allocated( self%header_lines )) deallocate( self%header_lines )

    ! Remove root element (and therefore the whole XML tree, recursively)
    if (allocated( self%root )) then
      call self%root%delete()
      deallocate( self%root )
    end if

  end subroutine t_xml_document__delete

!------------------------------------------------------------------------------
  recursive subroutine t_xml_element__write( self, indent, fid )
    class(sll_t_xml_element), intent(in) :: self
    integer                 , intent(in) :: indent
    integer                 , intent(in) :: fid

    logical :: empty_element
    integer :: na, i
    character(len=2) :: closing

    ! Does element have content?
    empty_element = .not. allocated( self%content )

    ! How many attributes in element?
    if (allocated( self%attributes )) then
      na = size( self%attributes )
    else
      na = 0
    endif

    ! Select closing tag for first line
    if (empty_element) then
      closing = '/>'
    else
      closing = '>'
    end if

    ! Write first line (only line if element is empty)
    write(fid,'(*(a))') &
      repeat( ' ', indent )//'<'//trim( self%name ), & ! element name
      (' '//self%attributes(i)%to_string(), i=1,na), & ! attributes
      trim( closing )                                  ! closing tag

    ! Write contents and last line (do nothing if element is empty)
    if (.not. empty_element) then
      call self%content%write( indent+2, fid )
      write(fid,'(a)') repeat( ' ',indent )//'</'//trim( self%name )//'>'
    endif

  end subroutine t_xml_element__write

!------------------------------------------------------------------------------
  recursive subroutine t_xml_element__delete( self )
    class(sll_t_xml_element), intent(inout) :: self

    ! Deallocate attributes
    if (allocated( self%attributes )) then
      deallocate( self%attributes )
    end if

    ! Delete and deallocate content
    if (allocated( self%content )) then
      call self%content%delete()
      deallocate( self%content )
    end if

  end subroutine t_xml_element__delete

!------------------------------------------------------------------------------
  subroutine t_xml_element__add_attribute( self, name, attvalue )
    class(sll_t_xml_element), intent(inout) :: self
    character(len=*)        , intent(in   ) :: name
    character(len=*)        , intent(in   ) :: attvalue

    integer :: na
    type( t_xml_attribute ), allocatable :: tmp(:)

    if (allocated( self%attributes )) then
      na = size( self%attributes )
      allocate( tmp(na+1) )
      tmp(1:na) = self%attributes(1:na)
    else
      na = 0
      allocate( tmp(1) )
    end if

    tmp(na+1)%name     = name
    tmp(na+1)%attvalue = attvalue 

    call move_alloc( from=tmp, to=self%attributes )

  end subroutine t_xml_element__add_attribute

!------------------------------------------------------------------------------
  subroutine t_xml_element__add_chardata_string( self, string )
    class(sll_t_xml_element), target, intent(inout) :: self
    character(len=*),                 intent(in   ) :: string

    type(t_xml_content), pointer :: new_cont  ! local variable

    ! Allocate new content (container)
    if (.not. allocated( self%content )) then
      allocate( self%content )
      new_cont => self%content
    else
      new_cont => self%content%new_content()
    end if

    ! Create new chardata inside the container
    allocate( t_xml_chardata :: new_cont%item )
    select type( new_item => new_cont%item ); type is( t_xml_chardata )

    ! Create new default data printer inside the chardata
    allocate( t_default_text_data_printer :: new_item%chardata )
    select type( new_cdata => new_item%chardata )
           type is( t_default_text_data_printer )

    ! Add text to printer
    new_cdata%text = adjustl( string )

    end select
    end select

  end subroutine t_xml_element__add_chardata_string

!------------------------------------------------------------------------------
  subroutine t_xml_element__add_chardata_printer( self, printer )
    class(sll_t_xml_element), target, intent(inout) :: self
    class(c_text_data_printer),       intent(in   ) :: printer

    type(t_xml_content), pointer :: new_cont  ! local variable

    ! Allocate new content (container)
    if (.not. allocated( self%content )) then
      allocate( self%content )
      new_cont => self%content
    else
      new_cont => self%content%new_content()
    end if

    ! Create new chardata inside container
    allocate( t_xml_chardata :: new_cont%item )
    select type( new_item => new_cont%item ); type is( t_xml_chardata )

    ! Copy user-defined printer inside chardata
    allocate( new_item%chardata, source=printer )

    end select

  end subroutine t_xml_element__add_chardata_printer

!------------------------------------------------------------------------------
  function t_xml_element__new_element( self, name ) result( new_elem )
    class(sll_t_xml_element), target, intent(inout) :: self
    character(len=*),                 intent(in   ) :: name

    type(sll_t_xml_element), pointer :: new_elem  ! output argument
    type(    t_xml_content), pointer :: new_cont  ! local variable

    ! Allocate new content (container)
    if (.not. allocated( self%content )) then
      allocate( self%content )
      new_cont => self%content
    else
      new_cont => self%content%new_content()
    end if

    ! Create new element inside the container
    allocate( sll_t_xml_element :: new_cont%item )
    select type( new_item => new_cont%item ); type is( sll_t_xml_element )
      new_elem => new_item
      new_elem%name = name
    end select

  end function t_xml_element__new_element

!------------------------------------------------------------------------------
  subroutine t_xml_chardata__write( self, indent, fid )
    class(t_xml_chardata), intent(in) :: self
    integer              , intent(in) :: indent
    integer              , intent(in) :: fid

    call self%chardata%print_text( indent, fid )

  end subroutine t_xml_chardata__write

!------------------------------------------------------------------------------
  subroutine t_xml_chardata__delete( self )
    class(t_xml_chardata), intent(inout) :: self

    if (allocated( self%chardata )) then
      call self%chardata%delete()
      deallocate( self%chardata )
    end if

  end subroutine

!------------------------------------------------------------------------------
  subroutine default_print_text( self, indent, fid )
    class(t_default_text_data_printer), intent(in) :: self
    integer                           , intent(in) :: indent
    integer                           , intent(in) :: fid
     
    integer :: i

    write(fid,'(*(a))') (' ', i=1,indent), trim( self%text )

  end subroutine default_print_text

!------------------------------------------------------------------------------
  subroutine default_delete( self )
    class(t_default_text_data_printer), intent(inout) :: self

    return
    print*, storage_size(self)
    ! Do nothing for now
  end subroutine

!------------------------------------------------------------------------------
  pure function attribute_name_val_string( self ) result( s )
    class(t_xml_attribute), intent(in) :: self
    !character( len = len(self%name)+1+len(self%attvalue) ) :: s
    character(len=:), allocatable :: s

    integer :: char_len

    char_len = len_trim( self%name ) + len_trim( self%attvalue ) + 3
    allocate( character(len=char_len) :: s )

    s = trim( self%name )//"='"//trim( self%attvalue )//"'"

  end function attribute_name_val_string

!------------------------------------------------------------------------------

end module sll_m_xml

