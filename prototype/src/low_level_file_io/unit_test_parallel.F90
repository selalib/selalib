program test_io_parallel

use mpi
use sll_collective
use sll_xml_io

#include "sll_remap.h"
#include "sll_memory.h"
#include "sll_working_precision.h"
#include "misc_utils.h"

implicit none

! Take a 2D array of dimensions ni*nj where ni, nj are the dimensions of
! the full array.
integer , parameter                       :: nx = 512
integer , parameter                       :: ny = 256
! Local sizes
integer                                   :: loc_sz_i_init
integer                                   :: loc_sz_j_init

! the process mesh
integer                                   :: npi
integer                                   :: npj
sll_int32                                 :: gi, gj
integer                                   :: ierr
integer                                   :: myrank
sll_int64                                 :: colsz        ! collective size
type(layout_2D), pointer                  :: layout

sll_real64                                :: rand_real
integer                                   :: i, j
sll_int32, dimension(2)                   :: global_indices
logical                                   :: test_passed

integer                  :: mx, my

!Parameters
real(8)                  :: tcpu1
real(8)                  :: tcpu2

real(8), dimension(:,:), allocatable :: xdata, ydata, zdata
character(len=8), parameter :: xfile = "xdata.h5" ! File name
character(len=8), parameter :: yfile = "ydata.h5" ! File name
character(len=8), parameter :: zfile = "zdata.h5" ! File name
character(len=8), parameter :: xdset = "xdataset" ! Dataset name
character(len=8), parameter :: ydset = "ydataset" ! Dataset name
character(len=8), parameter :: zdset = "zdataset" ! Dataset name

sll_int32 :: file_id, error

test_passed = .true.

! Boot parallel environment
call sll_boot_collective()

colsz  = sll_get_collective_size(sll_world_collective)
myrank = sll_get_collective_rank(sll_world_collective)

if( myrank .eq. 0) then
   print *, ' '
   print *, '--------------- HDF5 parallel test ---------------------'
   print *, ' '
   print"('Running a test on ',i4,' processes')", colsz
   call flush()
end if

if (.not. is_power_of_two(colsz)) then     
   print *, 'This test needs to run in a number of processes which is ',&
        'a power of 2.'
   call sll_halt_collective()
   stop
end if

layout => new_layout_2D( sll_world_collective )        
call factorize_in_random_2powers_2d(colsz, npi, npj)

if( myrank .eq. 0 ) then
   print *, 'source configuration: ', npi, npj
end if

call initialize_layout_with_distributed_2D_array( &
     nx, ny, npi, npj, layout )
     
call compute_local_sizes_2d( layout, loc_sz_i_init, loc_sz_j_init)        

! initialize the local data    
do j=1,loc_sz_j_init 
   do i=1,loc_sz_i_init
   end do
end do
     
call sll_collective_barrier(sll_world_collective)

tcpu1 = MPI_WTIME()

call compute_local_sizes_2d( layout, mx, my)        

SLL_ALLOCATE(xdata(mx,my),ierr)
SLL_ALLOCATE(ydata(mx,my),ierr)
SLL_ALLOCATE(zdata(mx,my),ierr)

do j = 1, my
   do i = 1, mx
      global_indices =  local_to_global_2D( layout, (/i, j/) )
      gi = global_indices(1)
      gj = global_indices(2)
      xdata(i,j) = float(gi-1)
      ydata(i,j) = float(gj-1)
      zdata(i,j) = myrank * xdata(i,j) * ydata(i,j)
   end do
end do

call write_hdf5_dataset(xfile,xdset,xdata)
call write_hdf5_dataset(yfile,ydset,ydata)
call write_hdf5_dataset(zfile,zdset,zdata)

call sll_xml_file_create("parallel.xmf",file_id,error)
call sll_xml_grid_geometry(file_id, xfile, xdset, nx, yfile, ydset, ny )
call sll_xml_field(file_id,'Z', "zdata.h5:/zdataset",nx,ny,'HDF','Node')
call sll_xml_file_close(file_id,error)

tcpu2 = MPI_WTIME()
if (myrank == 0) &
   write(*,"(//10x,' Temps CPU = ', G15.3, ' sec' )") (tcpu2-tcpu1)*colsz

call delete_layout_2D( layout )
call sll_halt_collective()

if( myrank .eq. 0) print *, 'PASSED'
  
contains

  subroutine factorize_in_random_2powers_2d(n, n1, n2)
    sll_int64, intent(in) :: n
    integer, intent(out)  ::n1, n2
    integer   :: expo, expo1, expo2
    if (.not.is_power_of_two(n)) then   
       print*, 'The number of processors must be a power of 2'
       stop
    endif 
    expo = int(log(real(n))/log(2.))  
    call random_number(rand_real)
    expo1 = int(rand_real*expo)
    expo2 = expo - expo1
    n1 = 2**expo1
    n2 = 2**expo2
  end subroutine factorize_in_random_2powers_2d

  subroutine compute_local_sizes_2d( layout, loc_sz_i, loc_sz_j )
    type(layout_2D), pointer :: layout
    sll_int32, intent(out) :: loc_sz_i
    sll_int32, intent(out) :: loc_sz_j
    sll_int32 :: i_min
    sll_int32 :: i_max
    sll_int32 :: j_min
    sll_int32 :: j_max
    sll_int32 :: my_rank
    if( .not. associated(layout) ) then
       print *, 'not-associated layout passed to new_distributed_mesh_2D'
       print *, 'Exiting...'
       STOP
    end if
    my_rank = sll_get_collective_rank(get_layout_2D_collective(layout))
    i_min = get_layout_2D_i_min( layout, my_rank )
    i_max = get_layout_2D_i_max( layout, my_rank )
    j_min = get_layout_2D_j_min( layout, my_rank )
    j_max = get_layout_2D_j_max( layout, my_rank )
    loc_sz_i = i_max - i_min + 1
    loc_sz_j = j_max - j_min + 1
  end subroutine compute_local_sizes_2d

  subroutine write_hdf5_dataset(filename, dsetname, fdata)

    use hdf5 ! This module contains all necessary modules 
        
    implicit none

    character(len=*), intent(in) :: filename 
    character(len=*), intent(in) :: dsetname 
    sll_real64, intent(in), dimension(:,:) :: fdata

    integer(HID_T) :: file_id       ! File identifier 
    integer(HID_T) :: dset_id       ! Dataset identifier 
    integer(HID_T) :: filespace     ! Dataspace identifier in file 
    integer(HID_T) :: plist_id      ! Property list identifier 
    integer(HID_T) :: memspace      ! Dataspace identifier in memory

    integer(HSIZE_T), dimension(2) :: dimsfi = (/nx,ny/)
    integer(HSIZE_T), dimension(2) :: datadims = (/nx,ny/)
    integer(HSIZE_T), dimension(2) :: chunk_dims ! Chunks dimensions
    
    integer(HSIZE_T),  dimension(2) :: count  
    integer(HSSIZE_T), dimension(2) :: offset 
    integer(HSIZE_T),  dimension(2) :: stride
    integer(HSIZE_T),  dimension(2) :: block

    integer :: rank = 2 ! Dataset rank 

    integer :: error  ! Error flags

    ! MPI definitions and calls.
    integer :: mpierror       ! MPI error flag
    integer :: comm, info
    integer :: mpi_size, mpi_rank

    comm = MPI_COMM_WORLD
    info = MPI_INFO_NULL

    call MPI_COMM_SIZE(comm, mpi_size, mpierror)
    call MPI_COMM_RANK(comm, mpi_rank, mpierror) 

    chunk_dims = (/mx,my/) ! Chunks dimensions


    ! Initialize HDF5 library and Fortran interfaces.
    call h5open_f(error) 

    ! Setup file access property list with parallel I/O access.
    call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, error)
    call h5pset_fapl_mpio_f(plist_id, comm, info, error)

    ! Create the file collectively.
    call h5fcreate_f(filename, H5F_ACC_TRUNC_F, file_id, error, access_prp = plist_id)
    call h5pclose_f(plist_id, error)
    
    ! Create the data space for the  dataset. 
    call h5screate_simple_f(rank, datadims, filespace, error)
    call h5screate_simple_f(rank, chunk_dims, memspace, error)

    ! Create chunked dataset.
    !
    call h5pcreate_f(H5P_DATASET_CREATE_F, plist_id, error)
    call h5pset_chunk_f(plist_id, rank, chunk_dims, error)
    call h5dcreate_f(file_id, dsetname, H5T_NATIVE_DOUBLE, filespace, &
                     dset_id, error, plist_id)
    call h5sclose_f(filespace, error)

    ! Each process defines dataset in memory and writes it to the hyperslab
    ! in the file. 
    stride(:) = 1
    count(:)  = 1
    block(1:2)  = (/mx,my/)

    offset(1) =  get_layout_2D_i_min( layout, mpi_rank ) - 1
    offset(2) =  get_layout_2D_j_min( layout, mpi_rank ) - 1
    ! 
    ! Select hyperslab in the file.
    !
    call h5dget_space_f(dset_id, filespace, error)
    call h5sselect_hyperslab_f (filespace, H5S_SELECT_SET_F, offset, count, error, &
                            stride, block)
    !
    ! Create property list for collective dataset write
    !
    call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error) 
    call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)
    !
    ! Write the dataset collectively. 
    !
    call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, fdata, dimsfi, error,   &
                    file_space_id = filespace, mem_space_id = memspace, &
                    xfer_prp = plist_id)
    !
    ! Write the dataset independently. 
    !
    !    call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, data, dimsfi,error, &
    !                     file_space_id = filespace, mem_space_id = memspace)
    !
    !
    ! Close dataspaces.
    !
    call h5sclose_f(filespace, error)
    call h5sclose_f(memspace, error)
    !
    ! Close the dataset.
    !
    call h5dclose_f(dset_id, error)
    !
    ! Close the property list.
    !
    call h5pclose_f(plist_id, error)
    !
    ! Close the file.
    !
    call h5fclose_f(file_id, error)
    !
    ! Close FORTRAN interfaces and HDF5 library.
    !
    call h5close_f(error)
    
  end subroutine write_hdf5_dataset

end program test_io_parallel
