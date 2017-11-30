!> @ingroup file_io_serial
!> @brief   Test serial HDF5 routines
!> @author  Pierre Navaro, INRIA
!> @author  Yaman Güçlü, IPP Garching
!> @details Test serial HDF5 write/read for 2D array of double precision floats,
!>          as well as write/read of scalar attributes (real64 and integer)
!>          attached to dataset of aforementioned array.

program test_hdf5_io_serial
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_assert.h"
#include "sll_working_precision.h"

  use sll_m_hdf5_io_serial, only:   &
    sll_t_hdf5_ser_handle         , &
    sll_s_hdf5_ser_file_create    , &
    sll_s_hdf5_ser_file_close     , &
    sll_s_hdf5_ser_file_open      , &
    sll_o_hdf5_ser_write_array    , &
    sll_o_hdf5_ser_read_array     , &
    sll_o_hdf5_ser_write_attribute, &
    sll_o_hdf5_ser_read_attribute

  implicit none
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  character(len=*), parameter :: fname     = "test_hdf5_io_serial.h5" ! File name
  character(len=*), parameter :: dsetname  = "real_array"      !     Dataset name
  character(len=*), parameter :: attr1name = "scalar_1"        ! Attribute-1 name
  character(len=*), parameter :: attr2name = "scalar_2"        ! Attribute-2 name
  real(f64)       , parameter :: attr1_in  = -1.67_f64 ! Value of attribute-1 to write
  integer         , parameter :: attr2_in  = 19        ! Value of attribute-2 to write

  type(sll_t_hdf5_ser_handle) :: fid
  integer                     :: ferror
  real(f64), allocatable      :: a(:,:)     !  Local data to write
  real(f64), allocatable      :: b(:,:)     ! Global data to read
  real(f64)                   :: attr1_out  ! Value of attribute-1 to read
  integer                     :: attr2_out  ! Value of attribute-2 to read
  integer                     :: i, j
  logical                     :: equal(0:2)

  !-----------------------------------------------------------------------------
  ! Check array WRITE/READ
  !-----------------------------------------------------------------------------

  write(*,'(/a)') "------------------------------------------------------------"
  write(*, '(a)') "TESTING WRITE/READ OF MULTI-DIMENSIONAL ARRAYS:"
  write(*, '(a)') "------------------------------------------------------------"

  ! Create rectangular matrix (2D array) with 7x11 elements and a(i,j) = i-j
  allocate( a(7,11) )
  do j = 1,11
    do i = 1,7
      a(i,j) = real( i-j, f64 )
    end do
  end do

  ! Create HDF5 file
  call sll_s_hdf5_ser_file_create( fname, fid, ferror )
  SLL_ASSERT( ferror == 0 )

  ! Write array to file
  call sll_o_hdf5_ser_write_array( fid, a, dsetname, ferror )
  SLL_ASSERT( ferror == 0 )

  ! Close file
  call sll_s_hdf5_ser_file_close( fid, ferror )
  SLL_ASSERT( ferror == 0 )

  ! Allocate memory for reading data from file
  allocate( b(7,11) )
  b(:,:) = 0.0_f64

  ! Read array from file
  call sll_s_hdf5_ser_file_open( fname, fid, ferror )
  call sll_o_hdf5_ser_read_array( fid, b, dsetname, ferror )
  call sll_s_hdf5_ser_file_close( fid, ferror )

  ! Check correctness of file
  equal(0) = all( a==b )

  if (equal(0)) then
    write(*,'(/a)') "OK: a and b arrays are identical"
  else
    write(*,'(/a)') "ERROR: a and b arrays differ"
  end if

  !-----------------------------------------------------------------------------
  ! Check attribute WRITE/READ
  !-----------------------------------------------------------------------------

  write(*,'(/a)') "------------------------------------------------------------"
  write(*, '(a)') "TESTING WRITE/READ OF SCALAR ATTRIBUTES:"
  write(*, '(a)') "------------------------------------------------------------"
  write(*,'(5a12/)') "type", "name", "value[in]", "value[out]", "status"

  ! Attach attributes to array
  call sll_s_hdf5_ser_file_open( fname, fid, ferror )
  call sll_o_hdf5_ser_write_attribute( fid, dsetname, attr1name, attr1_in, ferror )
  call sll_o_hdf5_ser_write_attribute( fid, dsetname, attr2name, attr2_in, ferror )
  call sll_s_hdf5_ser_file_close( fid, ferror )

  ! Read attributes
  call sll_s_hdf5_ser_file_open( fname, fid, ferror )
  call sll_o_hdf5_ser_read_attribute( fid, dsetname, attr1name, attr1_out, ferror )
  call sll_o_hdf5_ser_read_attribute( fid, dsetname, attr2name, attr2_out, ferror )
  call sll_s_hdf5_ser_file_close( fid, ferror )

  ! Check correctness of attributes
  equal(1) = (attr1_in == attr1_out)
  equal(2) = (attr2_in == attr2_out)

  write(*,'(a12,a12,g12.3,g12.3,a12)') "real(f64)", attr1name, attr1_in, attr1_out, status( equal(1) )
  write(*,'(a12,a12,i12,i12,a12/)')    "integer"  , attr2name, attr2_in, attr2_out, status( equal(2) )

  if (all( equal(1:2) )) then
    write(*,'(a)') "OK: all attributes are identical"
  else
    write(*,'(a)') "ERROR: attributes differ"
  end if

  !-----------------------------------------------------------------------------
  ! Output for CTest
  !-----------------------------------------------------------------------------

  if (all( equal )) then
    write(*,'(/a/)') ">>>>>>>>>>>>>>>>>  UNIT TEST HAS PASSED!  <<<<<<<<<<<<<<<<<<"
  else
    write(*,'(/a/)') ">>>>>>>>>>>>>>>>>  UNIT TEST HAS FAILED!  <<<<<<<<<<<<<<<<<<"
  end if

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
contains
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  function status( equal ) result( msg )
    logical, intent(in) :: equal
    character(len=4) :: msg
    msg = merge( " OK ", "FAIL", equal )
  end function


end program test_hdf5_io_serial
