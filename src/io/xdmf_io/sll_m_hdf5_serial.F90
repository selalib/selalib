!------------------------------------------------------------------------------
! SELALIB
!------------------------------------------------------------------------------
! MODULE:   sll_m_hdf5_serial
!
! DESCRIPTION:
!> @ingroup xdmf_io
!> @authors Yaman Güçlü - <yaman.guclu@gmail.com>
!> @brief   Simple object-oriented wrapper to Pierre's "sll_m_hdf5_io_serial".
!> @todo    Add detailed description
!------------------------------------------------------------------------------
module sll_m_hdf5_serial

#ifndef NOHDF5

#include "sll_errors.h"
#include "sll_working_precision.h"

  use sll_m_hdf5_io_serial, only: &
    sll_hdf5_file_create, &
    sll_hdf5_file_open, &
    sll_hdf5_file_close, &
    sll_hdf5_write_array

  implicit none
  public :: sll_t_hdf5_serial
  private

  !> Maximum length of non-allocatable strings in module
  integer, parameter :: maxlen = 256

!==============================================================================

  type :: sll_t_hdf5_serial

    character(len=maxlen) :: filename = ""  !< file name
    integer               :: file_id  = -1  !< unit number
    logical               :: is_open  = .false.

  contains
    procedure :: init                => t_hdf5_serial__init
    procedure :: create              => t_hdf5_serial__create
    procedure :: open                => t_hdf5_serial__open
    procedure :: close               => t_hdf5_serial__close
    procedure :: delete              => t_hdf5_serial__delete

    procedure :: write_dble_array_1d => t_hdf5_serial__write_dble_array_1d
    procedure :: write_dble_array_2d => t_hdf5_serial__write_dble_array_2d
    procedure :: write_dble_array_3d => t_hdf5_serial__write_dble_array_3d
    procedure :: write_int_array_1d  => t_hdf5_serial__write_int_array_1d
    procedure :: write_int_array_2d  => t_hdf5_serial__write_int_array_2d
    procedure :: write_int_array_3d  => t_hdf5_serial__write_int_array_3d
    generic   :: write_array         => write_dble_array_1d, &
                                        write_dble_array_2d, &
                                        write_dble_array_3d, &
                                        write_int_array_1d, &
                                        write_int_array_2d, &
                                        write_int_array_3d
  end type sll_t_hdf5_serial

!==============================================================================
contains
!==============================================================================

  subroutine t_hdf5_serial__init( self, filename )
    class(sll_t_hdf5_serial), intent(  out) :: self
    character(len=*)        , intent(in   ) :: filename

    self%filename = filename

  end subroutine t_hdf5_serial__init

  !----------------------------------------------------------------------------
  subroutine t_hdf5_serial__create( self )
    class(sll_t_hdf5_serial), intent(inout) :: self

    character(len=*), parameter   :: this_sub_name = "t_hdf5_serial__create"
    !character(len=:), allocatable :: err_msg
    !integer                       :: file_id
    integer                       :: error

    call sll_hdf5_file_create( self%filename, self%file_id, error )

    if (error /= 0) then
      SLL_ERROR( this_sub_name, "Cannot create HDF5 file "//trim(self%filename))
    end if

    self%is_open = .true.

  end subroutine t_hdf5_serial__create

  !----------------------------------------------------------------------------
  subroutine t_hdf5_serial__open( self )
    class(sll_t_hdf5_serial), intent(inout) :: self

    character(len=*), parameter   :: this_sub_name = "t_hdf5_serial__open"
    !character(len=:), allocatable :: err_msg
    !integer                       :: file_id
    integer                       :: error

    if (self%is_open) then
      SLL_ERROR( this_sub_name, "File is already open")
    end if

    call sll_hdf5_file_open( self%filename, self%file_id, error )

    if (error /= 0) then
      SLL_ERROR( this_sub_name, "Cannot open HDF5 file " // trim( self%filename ))
    end if

    self%is_open = .true.

  end subroutine t_hdf5_serial__open

  !----------------------------------------------------------------------------
  subroutine t_hdf5_serial__close( self )
    class(sll_t_hdf5_serial), intent(inout) :: self

    character(len=*), parameter   :: this_sub_name = "t_hdf5_serial__close"
    !character(len=:), allocatable :: err_msg
    integer                       :: error

    if (.not. self%is_open) then
      SLL_ERROR( this_sub_name, "File is already closed")
    end if

    call sll_hdf5_file_close( self%file_id, error )

    if (error /= 0) then
      SLL_ERROR( this_sub_name, "Cannot close HDF5 file " // trim( self%filename ))
    end if

    self%is_open = .false.

  end subroutine t_hdf5_serial__close

  !----------------------------------------------------------------------------
  subroutine t_hdf5_serial__delete( self )
    class(sll_t_hdf5_serial), intent(inout) :: self

    if (self%is_open) call self%close()

    self%filename = ""
    self%file_id  = -1
    self%is_open  = .false.

  end subroutine t_hdf5_serial__delete

  !----------------------------------------------------------------------------
  subroutine t_hdf5_serial__write_dble_array_1d( self, array, dsetname )
    class(sll_t_hdf5_serial), intent(inout) :: self
    sll_real64              , intent(in   ) :: array(:)
    character(len=*)        , intent(in   ) :: dsetname

    character(len=*), parameter   :: this_sub_name = "t_hdf5_serial__write_dble_array_1d"
    !character(len=:), allocatable :: err_msg
    integer                       :: error

    if (.not. self%is_open) then
      SLL_ERROR( this_sub_name, "File is closed")
    end if

    call sll_hdf5_write_array( self%file_id, array, dsetname, error )

    if (error /= 0) then
      SLL_ERROR( this_sub_name, "Cannot write to file " // trim( self%filename ))
    end if

  end subroutine t_hdf5_serial__write_dble_array_1d

  !----------------------------------------------------------------------------
  subroutine t_hdf5_serial__write_dble_array_2d( self, array, dsetname )
    class(sll_t_hdf5_serial), intent(inout) :: self
    sll_real64              , intent(in   ) :: array(:,:)
    character(len=*)        , intent(in   ) :: dsetname

    character(len=*), parameter   :: this_sub_name = "t_hdf5_serial__write_dble_array_2d"
    !character(len=:), allocatable :: err_msg
    integer                       :: error

    if (.not. self%is_open) then
      SLL_ERROR( this_sub_name, "File is closed")
    end if

    call sll_hdf5_write_array( self%file_id, array, dsetname, error )

    if (error /= 0) then
      SLL_ERROR( this_sub_name, "Cannot write to file " // trim( self%filename ))
    end if

  end subroutine t_hdf5_serial__write_dble_array_2d

  !----------------------------------------------------------------------------
  subroutine t_hdf5_serial__write_dble_array_3d( self, array, dsetname )
    class(sll_t_hdf5_serial), intent(inout) :: self
    sll_real64              , intent(in   ) :: array(:,:,:)
    character(len=*)        , intent(in   ) :: dsetname

    character(len=*), parameter   :: this_sub_name = "t_hdf5_serial__write_dble_array_3d"
    !character(len=:), allocatable :: err_msg
    integer                       :: error

    if (.not. self%is_open) then
      SLL_ERROR( this_sub_name, "File is closed")
    end if

    call sll_hdf5_write_array( self%file_id, array, dsetname, error )

    if (error /= 0) then
      SLL_ERROR( this_sub_name, "Cannot write to file " // trim( self%filename ))
    end if

  end subroutine t_hdf5_serial__write_dble_array_3d

  !----------------------------------------------------------------------------
  subroutine t_hdf5_serial__write_int_array_1d( self, array, dsetname )
    class(sll_t_hdf5_serial), intent(inout) :: self
    sll_int32               , intent(in   ) :: array(:)
    character(len=*)        , intent(in   ) :: dsetname

    character(len=*), parameter   :: this_sub_name = "t_hdf5_serial__write_int_array_1d"
    !character(len=:), allocatable :: err_msg
    integer                       :: error

    if (.not. self%is_open) then
      SLL_ERROR( this_sub_name, "File is closed")
    end if

    call sll_hdf5_write_array( self%file_id, array, dsetname, error )

    if (error /= 0) then
      SLL_ERROR( this_sub_name, "Cannot write to file " // trim( self%filename ))
    end if

  end subroutine t_hdf5_serial__write_int_array_1d

  !----------------------------------------------------------------------------
  subroutine t_hdf5_serial__write_int_array_2d( self, array, dsetname )
    class(sll_t_hdf5_serial), intent(inout) :: self
    sll_int32               , intent(in   ) :: array(:,:)
    character(len=*)        , intent(in   ) :: dsetname

    character(len=*), parameter   :: this_sub_name = "t_hdf5_serial__write_int_array_2d"
    !character(len=:), allocatable :: err_msg
    integer                       :: error

    if (.not. self%is_open) then
      SLL_ERROR( this_sub_name, "File is closed")
    end if

    call sll_hdf5_write_array( self%file_id, array, dsetname, error )

    if (error /= 0) then
      SLL_ERROR( this_sub_name, "Cannot write to file " // trim( self%filename ))
    end if

  end subroutine t_hdf5_serial__write_int_array_2d

  !----------------------------------------------------------------------------
  subroutine t_hdf5_serial__write_int_array_3d( self, array, dsetname )
    class(sll_t_hdf5_serial), intent(inout) :: self
    sll_int32               , intent(in   ) :: array(:,:,:)
    character(len=*)        , intent(in   ) :: dsetname

    character(len=*), parameter   :: this_sub_name = "t_hdf5_serial__write_int_array_3d"
    !character(len=:), allocatable :: err_msg
    integer                       :: error

    if (.not. self%is_open) then
      SLL_ERROR( this_sub_name, "File is closed")
    end if

    call sll_hdf5_write_array( self%file_id, array, dsetname, error )

    if (error /= 0) then
      SLL_ERROR( this_sub_name, "Cannot write to file " // trim( self%filename ))
    end if


  end subroutine t_hdf5_serial__write_int_array_3d

#endif /* NOHDF5 */

end module sll_m_hdf5_serial
