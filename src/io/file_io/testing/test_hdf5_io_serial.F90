!
!> @ingroup file_io_serial
!> @brief   Test serial HDF5 write/read for 2D array of double precision floats
!
program test_hdf5_io_serial
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_assert.h"
#include "sll_working_precision.h"

  use sll_m_hdf5_io_serial, only: &
    sll_t_hdf5_ser_handle,      &
    sll_s_hdf5_ser_file_create, &
    sll_s_hdf5_ser_file_close,  &
    sll_s_hdf5_ser_file_open,   &
    sll_o_hdf5_ser_write_array, &
    sll_o_hdf5_ser_read_array

  implicit none
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  type(sll_t_hdf5_ser_handle) :: fid
  integer                     :: ferror
  real(f64), allocatable      :: a(:,:) !  Local data to write
  real(f64), allocatable      :: b(:,:) ! Global data to read
  character(len=*), parameter :: fname = "test_hdf5_io_serial.h5" ! File name
  character(len=*), parameter :: dsetname = "real_array"          ! Dataset name
  integer                     :: i, j
  logical                     :: equal

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
  equal = .true.
  do j = 1,11
    do i = 1,7
      if (a(i,j) /= b(i,j)) then
        equal = .false.
      end if
    end do
  end do

  ! Output for CTest
  if (equal) then
    print *, "PASSED"
  else
    print *, "ERROR: a and b arrays not identical"
  end if

end program test_hdf5_io_serial
