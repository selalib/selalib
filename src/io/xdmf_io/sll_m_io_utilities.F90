!------------------------------------------------------------------------------
! SELALIB
!------------------------------------------------------------------------------
! MODULE:   sll_m_io_utilities
!
! DESCRIPTION:
!> @ingroup xdmf_io
!> @authors Yaman Güçlü - <yaman.guclu@gmail.com>
!> @brief   Collection of functions/subroutines operating on files and strings
!> @todo    Add detailed description
!------------------------------------------------------------------------------
module sll_m_io_utilities

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_working_precision.h"

  implicit none

  public :: &
    sll_f_check_empty_file, &
    sll_f_check_equal_files, &
    sll_s_ints_to_string, &
    sll_s_read_file, &
    sll_s_remove_file, &
    sll_s_split_path

  private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!==============================================================================
contains
!==============================================================================

  !----------------------------------------------------------------------------
  !> Remove file (dangerous function!)
  subroutine sll_s_remove_file( filename )
    character(len=*), intent(in) :: filename

    integer :: iunit

    open( newunit=iunit, file=filename, status='OLD' )
    close( iunit, status='DELETE' )

  end subroutine sll_s_remove_file

  !----------------------------------------------------------------------------
  !> Read whole file content into allocatable string
  subroutine sll_s_read_file( filename, str )
    character(len=*),              intent(in   ) :: filename
    character(len=:), allocatable, intent(  out) :: str

    integer :: iunit, istat, filesize
    character(len=1) :: c

    ! Open required file
    open( newunit=iunit, file=filename, status='OLD', &
      form='UNFORMATTED', access='STREAM' )

    ! How many characters in file
    inquire( file=filename, size=filesize )
    if (filesize < 0) then
      write(*,*) 'ERROR: negative file size'
      stop
    end if

    ! Read whole file into one string
    allocate( character(len=filesize) :: str )
    read( iunit, pos=1 ) str

    ! Make sure it was all read by trying to read more
    read( iunit, pos=filesize+1, iostat=istat ) c
    if (.not. IS_IOSTAT_END(istat)) then
      write(*,*) 'Error: file was not completely read'
      stop
    end if

    ! Close file
    close( iunit, iostat=istat )

  end subroutine sll_s_read_file

  !----------------------------------------------------------------------------
  ! Return True if two files are identical, False otherwise
  function sll_f_check_equal_files( filename1, filename2 ) result( equal )
    character(len=*), intent(in   ) :: filename1
    character(len=*), intent(in   ) :: filename2
    logical :: equal

    character(len=:), allocatable :: str1
    character(len=:), allocatable :: str2

    call sll_s_read_file( filename1, str1 )
    call sll_s_read_file( filename2, str2 )

    equal = (str1 == str2)
    
  end function sll_f_check_equal_files

  !----------------------------------------------------------------------------
  !> Return True if file is empty, False otherwise
  function sll_f_check_empty_file( filename ) result( empty )
    character(len=*), intent(in) :: filename
    logical :: empty

    integer :: filesize

    inquire( file=filename, size=filesize )
    empty = (filesize == 0)

  end function sll_f_check_empty_file

  !----------------------------------------------------------------------------
  !> Write an array of integers to a single string:
  !>   . Numbers are separated by a blank space;  
  !>   . String is allocated here with minimum length.
  subroutine sll_s_ints_to_string( ints, str )
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

  end subroutine sll_s_ints_to_string

  !----------------------------------------------------------------------------
  !> Split path into head (directory) and tail (file)
  subroutine sll_s_split_path( path, head, tail )
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

  end subroutine sll_s_split_path

!==============================================================================
end module sll_m_io_utilities
