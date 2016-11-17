!
!> @ingroup file_io_parallel
!> @brief   Test parallel write/read for 2D array of double precision floats
!> @author  Yaman Güçlü, IPP Garching
!
program test_hdf5_io_parallel
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_errors.h"
#include "sll_assert.h"
#include "sll_working_precision.h"

  use sll_m_hdf5_io_parallel, only: &
    sll_o_hdf5_file_create, &
    sll_o_hdf5_file_close,  &
    sll_o_hdf5_file_open,  &
    sll_o_hdf5_write_array, &
    sll_o_hdf5_read_array

  use sll_mpi, only: &
    mpi_comm_world, &
    mpi_info_null, &
    mpi_init, &
    mpi_comm_size, &
    mpi_comm_rank, &
    mpi_finalize

  use hdf5, only: &
    hsize_t, &
    hssize_t

  implicit none
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  ! HDF5 variables (to be converted into 'hsize_t' and 'hssize_t')
  integer :: dset_dims(2) ! Dataset dimensions in the file (global)
!  integer :: chnk_dims(2) ! Chunk dimensions (global)
  integer :: block (2)    ! Dimensions of local array
  integer :: offset(2)    ! Offset of local array within global dataset

  ! MPI variables
  integer :: perror
  integer :: pcomm, pinfo
  integer :: psize, prank

  ! SELALIB variables
  integer(i32) :: fid
  integer(i32) :: ferror
  real(f64), allocatable :: a(:,:) !  Local data to write
  real(f64), allocatable :: b(:,:) ! Global data to read
  character(len=*), parameter :: fname = "test_hdf5_io_parallel.h5" ! File name
  character(len=*), parameter :: dsetname = "real_array"         ! Dataset name
  integer :: i, j
  
  pcomm = mpi_comm_world
  pinfo = mpi_info_null

  call mpi_init( perror )
  call mpi_comm_size( pcomm, psize, perror )
  call mpi_comm_rank( pcomm, prank, perror )

  ! Quit if mpi size is not 4
  if (psize /= 4) then
    SLL_ERROR( "program 'test_hdf5_io_parallel'", "Test requires 4 processes" )
  endif

  !------------------------------------------------
  ! Partitioning of 8x10 domain across 4 processes:
  !     ___________________
  !    |       |           |
  !  3 |  (0)  |    (1)    |
  !    |_______|___________|
  !    |       |           |
  !    |       |           |
  !  5 |  (2)  |    (3)    |
  !    |       |           |
  !    |_______|___________|
  !        4         6
  !
  ! NOTE: the HDF5 library interprets and stores
  !       the data as row major, therefore the
  !       resulting array will appear transposed 
  !------------------------------------------------

  ! Global array dimensions
  dset_dims = [8,10]

  ! Create data space for local array
  select case( prank )
  case(0)
    block  = [3,4]
    offset = [0,0]
  case(1)
    block  = [3,6]
    offset = [0,4]
  case(2)
    block  = [5,4]
    offset = [3,0]
  case(3)
    block  = [5,6]
    offset = [3,4]
  end select
  allocate( a(block(1),block(2)) )

  ! Set all values in local array = mpi rank
  a(:,:) = real( prank, f64 )

  ! Create parallel HDF5 file
  call sll_o_hdf5_file_create( fname, pcomm, fid, ferror )
  SLL_ASSERT( ferror == 0 )

  ! Parallel write array to file
  call sll_o_hdf5_write_array( &
    file_id     = fid, &
    global_size = int( dset_dims,  hsize_t ), &
    offset      = int(    offset, hssize_t ), &
    array       = a, &
    dsetname    = dsetname, &
    error       = ferror )
  SLL_ASSERT( ferror == 0 )

  ! Close file
  call sll_o_hdf5_file_close( fid, ferror )
  SLL_ASSERT( ferror == 0 )

  ! Allocate memory for reading data from file
  allocate( b(block(1),block(2)) )
  b(:,:) = -1.0_f64

  ! Read file in parallel
!    call sll_o_hdf5_file_open( fname, pcomm, fid, ferror )
  call sll_o_hdf5_file_open( fid, fname, pcomm, ferror ) !TODO: fix input order
  call sll_o_hdf5_read_array( &
    file_id     = fid, &
    global_size = int( dset_dims,  hsize_t ), &
    offset      = int( offset, hssize_t ), &
    array       = b, &
    dsetname    = dsetname, &
    error       = ferror )
  call sll_o_hdf5_file_close( fid, ferror )

  ! Check correctness of file
  do i = 1, block(1)
    do j = 1, block(2)
      if (a(i,j) /= b(i,j)) then
        print *, "ERROR"
      end if
    end do
  end do

  ! Turn off parallel environment
  call mpi_finalize( perror )

end program test_hdf5_io_parallel
