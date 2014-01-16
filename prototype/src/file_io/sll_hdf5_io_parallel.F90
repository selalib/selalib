!**************************************************************
!  Copyright INRIA
!  Authors : 
!     Pierre Navaro 
!  
!  This code SeLaLib (for Semi-Lagrangian-Library) 
!  is a parallel library for simulating the plasma turbulence 
!  in a tokamak.
!  
!  This software is governed by the CeCILL-B license 
!  under French law and abiding by the rules of distribution 
!  of free software.  You can  use, modify and redistribute 
!  the software under the terms of the CeCILL-B license as 
!  circulated by CEA, CNRS and INRIA at the following URL
!  "http://www.cecill.info". 
!**************************************************************

!> @brief
!> Implements the functions to write hdf5 file to store heavy data
!> @details
!> With HDF5 you can store several datasets in a single file.
!> - HDF5 file (http://www.hdfgroup.org/HDF5/)
module sll_hdf5_io_parallel
#include "sll_working_precision.h"
  
#ifndef NOHDF5

  use hdf5
  use sll_collective

  implicit none

  !> Write array in hdf5 file
  interface sll_hdf5_write_array
     module procedure sll_hdf5_write_array_1d
     module procedure sll_hdf5_write_array_2d
     module procedure sll_hdf5_write_array_3d
  end interface

contains
  
  !> Create HDF5 file
  !>    - Initialize fortran interface
  !>    - Create a new file using default properties
  subroutine sll_hdf5_file_create(filename, file_id, error)

    character(len=*), intent(in)  :: filename   !< file name
    integer,          intent(out) :: error      !< error code
    integer(hid_t),   intent(out) :: file_id    !< file unit number
    integer(hid_t)                :: plist_id   
    integer                       :: comm, info
    integer                       :: mpi_size
    integer                       :: mpi_rank

    comm     = sll_world_collective%comm
    info     = MPI_INFO_NULL
    mpi_size = sll_get_collective_size(sll_world_collective)
    mpi_rank = sll_get_collective_rank(sll_world_collective)

    call h5open_f(error) 
    call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, error)
    call h5pset_fapl_mpio_f(plist_id, comm, info, error)
    call h5fcreate_f(filename, H5F_ACC_TRUNC_F, file_id, error, access_prp = plist_id)
    call h5pclose_f(plist_id, error)

  end subroutine sll_hdf5_file_create

  !> Open HDF5 file
  !>    - Initialize fortran interface
  !>    - Open a HDF5 file
  subroutine sll_hdf5_file_open(file_id,filename,error)

    character(len=*) , intent(in)  :: filename    !< file name
    integer(hid_t)                 :: file_id     !< file unit number
    integer,           intent(out) :: error       !< error code
    integer(hid_t)                 :: plist_id    
    integer                        :: comm, info
    integer                        :: mpi_size
    integer                        :: mpi_rank

    comm     = sll_world_collective%comm
    info     = MPI_INFO_NULL
    mpi_size = sll_get_collective_size(sll_world_collective)
    mpi_rank = sll_get_collective_rank(sll_world_collective)
    
    call h5open_f(error) 
    call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, error)
    call h5pset_fapl_mpio_f(plist_id, comm, info, error)
    call h5fopen_f(filename, H5F_ACC_RDONLY_F, file_id, error, access_prp = plist_id)
    call h5pclose_f(plist_id, error)

  end subroutine sll_hdf5_file_open

  !> Close HDF5 file 
  subroutine sll_hdf5_file_close(file_id,error)
    integer(hid_t), intent(in)     :: file_id   !< file unit number
    integer, intent(out)           :: error     !< error code

    !
    ! Close property list and the file.
    !
    call h5fclose_f(file_id, error)
    !
    ! Close FORTRAN interface
    !
    call h5close_f(error)
  end subroutine sll_hdf5_file_close

  !> Write a 1D array of float in double precision in a HDF5 file
  !> - Create a dataspace with 1 dimensions
  !> - Write the dataset
  !> - Close dataset and dataspace
  subroutine sll_hdf5_write_array_1d(file_id,global_size,offset,array,dsetname,error)
    integer, parameter                     :: dspace_dims=1
    character(len=*), intent(in)           :: dsetname
    integer(hsize_t), intent(in)           :: global_size(dspace_dims)
    sll_int32, intent(out)                 :: error
    sll_real64, intent(in)                 :: array(:)
    integer(hid_t)                         :: file_id
    integer(hid_t)                         :: plist_id
    sll_int32                              :: rank, i
    integer(hsize_t)                       :: array_dims(dspace_dims)
    integer(hid_t)                         :: dset_id
    integer(hid_t)                         :: memspace
    integer(hid_t)                         :: filespace
    integer(HSIZE_T)                       :: dimsfi(dspace_dims)
    integer(HSIZE_T)                       :: count(dspace_dims)
    integer(HSSIZE_T)                      :: offset(dspace_dims)
    integer(HSIZE_T)                       :: stride(dspace_dims)
    integer(HSIZE_T)                       :: block(dspace_dims)
    rank = dspace_dims
    do i = 1, rank
      array_dims(i) = size(array,i)
    end do
    dimsfi = global_size
    call h5screate_simple_f(rank, global_size, filespace, error)
    call h5screate_simple_f(rank, array_dims, memspace, error)
    call h5pcreate_f(H5P_DATASET_CREATE_F, plist_id, error)
    call h5pset_chunk_f(plist_id, rank, array_dims, error)
    call h5dcreate_f(file_id, dsetname, H5T_NATIVE_DOUBLE, filespace,  &
                     dset_id, error, plist_id)
    call h5pclose_f(plist_id, error)
    call h5sclose_f(filespace, error)
    stride = 1
    count  = 1
    block  = array_dims
    call h5dget_space_f(dset_id, filespace, error)
    call h5sselect_hyperslab_f (filespace, H5S_SELECT_SET_F,           &
                                offset, count, error, stride, block)
    call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
    call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)
    call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, array, dimsfi, error,  &
                    file_space_id = filespace, mem_space_id = memspace,&
                    xfer_prp = plist_id)
    call h5pclose_f(plist_id, error)
    call h5sclose_f(filespace, error)
    call h5sclose_f(memspace, error)
    call h5dclose_f(dset_id, error)
  end subroutine 

  !> Write a 2D array of float in double precision in a HDF5 file
  !> - Create a dataspace with 2 dimensions
  !> - Write the dataset
  !> - Close dataset and dataspace
  subroutine sll_hdf5_write_array_2d(file_id,global_size,offset,array,dsetname,error)
    integer, parameter                     :: dspace_dims=2
    character(len=*), intent(in)           :: dsetname
    integer(hsize_t), intent(in)           :: global_size(dspace_dims)
    sll_int32, intent(out)                 :: error
    sll_real64, intent(in)                 :: array(:,:)
    integer(hid_t)                         :: file_id
    integer(hid_t)                         :: plist_id
    sll_int32                              :: rank, i
    integer(hsize_t)                       :: array_dims(dspace_dims)
    integer(hid_t)                         :: dset_id
    integer(hid_t)                         :: memspace
    integer(hid_t)                         :: filespace
    integer(HSIZE_T)                       :: dimsfi(dspace_dims)
    integer(HSIZE_T)                       :: count(dspace_dims)
    integer(HSSIZE_T)                      :: offset(dspace_dims)
    integer(HSIZE_T)                       :: stride(dspace_dims)
    integer(HSIZE_T)                       :: block(dspace_dims)
    rank = dspace_dims
    do i = 1, rank
      array_dims(i) = size(array,i)
    end do
    dimsfi = global_size
    call h5screate_simple_f(rank, global_size, filespace, error)
    call h5screate_simple_f(rank, array_dims, memspace, error)
    call h5pcreate_f(H5P_DATASET_CREATE_F, plist_id, error)
    call h5pset_chunk_f(plist_id, rank, array_dims, error)
    call h5dcreate_f(file_id, dsetname, H5T_NATIVE_DOUBLE, filespace,  &
                     dset_id, error, plist_id)
    call h5pclose_f(plist_id, error)
    call h5sclose_f(filespace, error)
    stride = 1
    count  = 1
    block  = array_dims
    call h5dget_space_f(dset_id, filespace, error)
    call h5sselect_hyperslab_f (filespace, H5S_SELECT_SET_F,           &
                                offset, count, error, stride, block)
    call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
    call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)
    call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, array, dimsfi, error,  &
                    file_space_id = filespace, mem_space_id = memspace,&
                    xfer_prp = plist_id)
    call h5pclose_f(plist_id, error)
    call h5sclose_f(filespace, error)
    call h5sclose_f(memspace, error)
    call h5dclose_f(dset_id, error)
  end subroutine 

  !> Write a 3D array of float in double precision in a HDF5 file
  !> - Create a dataspace with 3 dimensions
  !> - Write the dataset
  !> - Close dataset and dataspace
  subroutine sll_hdf5_write_array_3d(file_id,global_size,offset,array,dsetname,error)
    integer, parameter                     :: dspace_dims=3
    character(len=*), intent(in)           :: dsetname
    integer(hsize_t), intent(in)           :: global_size(dspace_dims)
    sll_int32, intent(out)                 :: error
    sll_real64, intent(in)                 :: array(:,:,:)
    integer(hid_t)                         :: file_id
    integer(hid_t)                         :: plist_id
    sll_int32                              :: rank, i
    integer(hsize_t)                       :: array_dims(dspace_dims)
    integer(hid_t)                         :: dset_id
    integer(hid_t)                         :: memspace
    integer(hid_t)                         :: filespace
    integer(HSIZE_T)                       :: dimsfi(dspace_dims)
    integer(HSIZE_T)                       :: count(dspace_dims)
    integer(HSSIZE_T)                      :: offset(dspace_dims)
    integer(HSIZE_T)                       :: stride(dspace_dims)
    integer(HSIZE_T)                       :: block(dspace_dims)
    rank = dspace_dims
    do i = 1, rank
      array_dims(i) = size(array,i)
    end do
    dimsfi = global_size
    call h5screate_simple_f(rank, global_size, filespace, error)
    call h5screate_simple_f(rank, array_dims, memspace, error)
    call h5pcreate_f(H5P_DATASET_CREATE_F, plist_id, error)
    call h5pset_chunk_f(plist_id, rank, array_dims, error)
    call h5dcreate_f(file_id, dsetname, H5T_NATIVE_DOUBLE, filespace,  &
                     dset_id, error, plist_id)
    call h5pclose_f(plist_id, error)
    call h5sclose_f(filespace, error)
    stride = 1
    count  = 1
    block  = array_dims
    call h5dget_space_f(dset_id, filespace, error)
    call h5sselect_hyperslab_f (filespace, H5S_SELECT_SET_F,           &
                                offset, count, error, stride, block)
    call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
    call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)
    call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, array, dimsfi, error,  &
                    file_space_id = filespace, mem_space_id = memspace,&
                    xfer_prp = plist_id)
    call h5pclose_f(plist_id, error)
    call h5sclose_f(filespace, error)
    call h5sclose_f(memspace, error)
    call h5dclose_f(dset_id, error)
  end subroutine 
!Gysela functions that can be useful for future
!  ! HDF5 saving for an integer
!  subroutine HDF5_integer_saving(file_id,int,dsetname)
!  ! HDF5 saving for a real double
!  subroutine HDF5_real_saving(file_id,rd,dsetname)
!  ! HDF5 saving for a 1D array of integer
!  subroutine HDF5_array1D_saving_int(file_id,array1D,dim1,dsetname)
!  ! gzip HDF5 saving for a 1D array of real*4
!  subroutine HDF5_array1D_saving_r4(file_id,array1D,dim1,dsetname)
!  ! gzip HDF5 saving for a 1D array of real*8
!  subroutine HDF5_array1D_saving(file_id,array1D,dim1,dsetname)
!  ! gzip HDF5 saving for a 2D array double
!  subroutine HDF5_array2D_saving(file_id,array2D,dim1,dim2,dsetname)
!  ! gzip HDF5 saving for a 3D array double
!  subroutine HDF5_array3D_saving(file_id,array3D,dim1,dim2,dim3,dsetname)
!  ! gzip HDF5 saving for a 4D array
!  subroutine HDF5_array4D_saving(file_id,array4d,dim1,dim2,dim3,dim4,dsetname)
!  ! gzip HDF5 saving for a 5D array
!  subroutine HDF5_array5D_saving(file_id,array5d,dim1,dim2,dim3,dim4,dim5,dsetname)
!  
!  ! HDF5 reading for an integer 
!  subroutine HDF5_integer_reading(file_id,itg,dsetname)
!  ! HDF5 reading for a real double
!  subroutine HDF5_real_reading(file_id,rd,dsetname)
!  ! HDF5 reading for an array 1D double
!  subroutine HDF5_array1D_reading(file_id,array1D,dsetname)
!  ! HDF5 reading for an array 2DS double
!  subroutine HDF5_array2D_reading(file_id,array2D,dsetname)
!  ! HDF5 reading for an array 3D double
!  subroutine HDF5_array3D_reading(file_id,array3D,dsetname)
!  ! HDF5 reading for an array 4D
!  subroutine HDF5_array4D_reading(file_id,array4D,dsetname,error)
!  ! HDF5 reading for an array 5D
!  subroutine HDF5_array5D_reading(file_id,array5D,dsetname)

#endif

end module sll_hdf5_io_parallel

