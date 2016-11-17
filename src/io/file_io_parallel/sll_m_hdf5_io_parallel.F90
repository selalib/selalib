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

!> @ingroup file_io_parallel
!> @brief
!> Parallel version of sll_hdf5_io
!> @details
!> With HDF5 you can store several datasets in a single file.
!> - HDF5 file (http://www.hdfgroup.org/HDF5/)
module sll_m_hdf5_io_parallel
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#ifndef NOHDF5

#include "sll_assert.h"
#include "sll_working_precision.h"

  use hdf5, only: &
    h5dclose_f, &
    h5dcreate_f, &
    h5dget_space_f, &
    h5dopen_f, &
    h5dread_f, &
    h5dwrite_f, &
    h5f_acc_rdonly_f, &
    h5f_acc_trunc_f, &
    h5fcreate_f, &
    h5fd_mpio_collective_f, &
    h5fopen_f, &
    h5open_f, &
    h5p_dataset_create_f, &
    h5p_dataset_xfer_f, &
    h5p_file_access_f, &
    h5pclose_f, &
    h5pcreate_f, &
    h5pset_chunk_f, &
    h5pset_dxpl_mpio_f, &
    h5pset_fapl_mpio_f, &
    h5s_select_and_f, &
    h5s_select_set_f, &
    h5sclose_f, &
    h5screate_simple_f, &
    h5sselect_hyperslab_f, &
    h5t_native_double, &
    hid_t, &
    hsize_t, &
    hssize_t

  use sll_m_hdf5_io_serial, only: &
    sll_o_hdf5_file_close

  use sll_mpi, only: &
    mpi_info_null

  implicit none

  public :: &
    sll_o_hdf5_file_close, &
    sll_o_hdf5_file_create, &
    sll_o_hdf5_write_array, &
    sll_o_hdf5_file_open, &
    sll_o_hdf5_read_array

  private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  !> Create new HDF5 file
  interface sll_o_hdf5_file_create
    module procedure sll_hdf5_par_file_create
  end interface

  !> Open existing HDF5 file
  interface sll_o_hdf5_file_open
    module procedure sll_hdf5_par_file_open
  end interface

  !> Write array in HDF5 file
  interface sll_o_hdf5_write_array
     module procedure sll_hdf5_par_write_dble_array_1d
     module procedure sll_hdf5_par_write_dble_array_2d
     module procedure sll_hdf5_par_write_dble_array_3d
     module procedure sll_hdf5_par_write_dble_array_4d
     module procedure sll_hdf5_par_write_dble_array_5d
     module procedure sll_hdf5_par_write_dble_array_6d
  end interface

  !> Read array form HDF5 file
  interface sll_o_hdf5_read_array
     module procedure sll_hdf5_par_read_array_6d
  end interface sll_o_hdf5_read_array

contains

  !-----------------------------------------------------------------------------
  !> Create HDF5 file
  !>    - Initialize fortran interface
  !>    - Create a new file using default properties
  !-----------------------------------------------------------------------------
  subroutine sll_hdf5_par_file_create( filename, comm, file_id, error )

    character(len=*), intent(in   ) :: filename   !< file name
    integer         , intent(in   ) :: comm       !< MPI comm
    integer(hid_t)  , intent(  out) :: file_id    !< file unit number
    integer         , intent(  out) :: error      !< error code

    integer(hid_t) :: plist_id
    integer        :: info
    info = MPI_INFO_NULL

    call h5open_f(error) 
    SLL_ASSERT(error==0)
    call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, error)
    SLL_ASSERT(error==0)
    call h5pset_fapl_mpio_f(plist_id, comm, info, error)
    SLL_ASSERT(error==0)
    call h5fcreate_f(filename, H5F_ACC_TRUNC_F, file_id, error, access_prp = plist_id)
    SLL_ASSERT(error==0)
    call h5pclose_f(plist_id, error)
    SLL_ASSERT(error==0)

  end subroutine sll_hdf5_par_file_create

  !-----------------------------------------------------------------------------
  !> Open HDF5 file
  !>    - Initialize fortran interface
  !>    - Open a HDF5 file
  !-----------------------------------------------------------------------------
  subroutine sll_hdf5_par_file_open( file_id, filename, comm, error )

    integer(hid_t)  , intent(  out) :: file_id     !< file unit number
    character(len=*), intent(in   ) :: filename    !< file name
    integer         , intent(in   ) :: comm        !< error code
    integer         , intent(  out) :: error       !< error code

    integer(hid_t) :: plist_id
    integer        :: info
    info = MPI_INFO_NULL
    
    call h5open_f(error) 
    SLL_ASSERT(error==0)
    call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, error)
    SLL_ASSERT(error==0)
    call h5pset_fapl_mpio_f(plist_id, comm, info, error)
    SLL_ASSERT(error==0)
    call h5fopen_f(filename, H5F_ACC_RDONLY_F, file_id, error, access_prp = plist_id)
    SLL_ASSERT(error==0)
    call h5pclose_f(plist_id, error)
    SLL_ASSERT(error==0)

  end subroutine sll_hdf5_par_file_open

!  !> Close HDF5 file
!  subroutine sll_o_hdf5_file_close(file_id,error)
!    integer(hid_t), intent(in)     :: file_id   !< file unit number
!    integer, intent(out)           :: error     !< error code
!
!    !
!    ! Close property list and the file.
!    !
!    call h5fclose_f(file_id, error)
!    SLL_ASSERT(error==0)
!    !
!    ! Close FORTRAN interface
!    !
!    call h5close_f(error)
!    SLL_ASSERT(error==0)
!  end subroutine sll_o_hdf5_file_close

  !-----------------------------------------------------------------------------
  !> Write a 1D array of float in double precision in a HDF5 file
  !> - Create a dataspace with 1 dimensions
  !> - Write the dataset
  !> - Close dataset and dataspace
  !-----------------------------------------------------------------------------
  subroutine sll_hdf5_par_write_dble_array_1d( &
      file_id, global_size, offset, array, dsetname, error )
    integer, parameter               :: dspace_dims = 1
    integer(   hid_t), intent(in   ) :: file_id
    integer( hsize_t), intent(in   ) :: global_size(dspace_dims)
    integer(hssize_t), intent(in   ) :: offset     (dspace_dims)
    sll_real64       , intent(in   ) :: array(:)
    character(len=*) , intent(in   ) :: dsetname
    sll_int32        , intent(  out) :: error

#define  DATATYPE  H5T_NATIVE_DOUBLE
#include "sll_k_hdf5_write_array.f90"

  end subroutine sll_hdf5_par_write_dble_array_1d

  !-----------------------------------------------------------------------------
  !> Write a 2D array of float in double precision in a HDF5 file
  !> - Create a dataspace with 2 dimensions
  !> - Write the dataset
  !> - Close dataset and dataspace
  !-----------------------------------------------------------------------------
  subroutine sll_hdf5_par_write_dble_array_2d( &
      file_id, global_size, offset, array, dsetname, error )
    integer, parameter               :: dspace_dims = 2
    integer(   hid_t), intent(in   ) :: file_id
    integer( hsize_t), intent(in   ) :: global_size(dspace_dims)
    integer(hssize_t), intent(in   ) :: offset     (dspace_dims)
    sll_real64       , intent(in   ) :: array(:,:)
    character(len=*) , intent(in   ) :: dsetname
    sll_int32        , intent(  out) :: error

#define  DATATYPE  H5T_NATIVE_DOUBLE
#include "sll_k_hdf5_write_array.f90"

  end subroutine sll_hdf5_par_write_dble_array_2d

  !-----------------------------------------------------------------------------
  !> Write a 3D array of float in double precision in a HDF5 file
  !> - Create a dataspace with 3 dimensions
  !> - Write the dataset
  !> - Close dataset and dataspace
  !-----------------------------------------------------------------------------
  subroutine sll_hdf5_par_write_dble_array_3d( &
      file_id, global_size, offset, array, dsetname, error )
    integer, parameter               :: dspace_dims = 3
    integer(   hid_t), intent(in   ) :: file_id
    integer( hsize_t), intent(in   ) :: global_size(dspace_dims)
    integer(hssize_t), intent(in   ) :: offset     (dspace_dims)
    sll_real64       , intent(in   ) :: array(:,:,:)
    character(len=*) , intent(in   ) :: dsetname
    sll_int32        , intent(  out) :: error

#define  DATATYPE  H5T_NATIVE_DOUBLE
#include "sll_k_hdf5_write_array.f90"

  end subroutine sll_hdf5_par_write_dble_array_3d

  !-----------------------------------------------------------------------------
  !> Write a 4D array of float in double precision in a HDF5 file
  !> - Create a dataspace with 4 dimensions
  !> - Write the dataset
  !> - Close dataset and dataspace
  !-----------------------------------------------------------------------------
  subroutine sll_hdf5_par_write_dble_array_4d( &
      file_id, global_size, offset, array, dsetname, error )
    integer, parameter               :: dspace_dims = 4
    integer(   hid_t), intent(in   ) :: file_id
    integer( hsize_t), intent(in   ) :: global_size(dspace_dims)
    integer(hssize_t), intent(in   ) :: offset     (dspace_dims)
    sll_real64       , intent(in   ) :: array(:,:,:,:)
    character(len=*) , intent(in   ) :: dsetname
    sll_int32        , intent(  out) :: error

#define  DATATYPE  H5T_NATIVE_DOUBLE
#include "sll_k_hdf5_write_array.f90"

  end subroutine sll_hdf5_par_write_dble_array_4d

  !-----------------------------------------------------------------------------
  !> Write a 5D array of float in double precision in a HDF5 file
  !> - Create a dataspace with 5 dimensions
  !> - Write the dataset
  !> - Close dataset and dataspace
  !-----------------------------------------------------------------------------
  subroutine sll_hdf5_par_write_dble_array_5d( &
      file_id, global_size, offset, array, dsetname, error )
    integer, parameter               :: dspace_dims = 5
    integer(   hid_t), intent(in   ) :: file_id
    integer( hsize_t), intent(in   ) :: global_size(dspace_dims)
    integer(hssize_t), intent(in   ) :: offset     (dspace_dims)
    sll_real64       , intent(in   ) :: array(:,:,:,:,:)
    character(len=*) , intent(in   ) :: dsetname
    sll_int32        , intent(  out) :: error

#define  DATATYPE  H5T_NATIVE_DOUBLE
#include "sll_k_hdf5_write_array.f90"

  end subroutine sll_hdf5_par_write_dble_array_5d

  !-----------------------------------------------------------------------------
  !> Write a 6D array of float in double precision in a HDF5 file
  !> - Create a dataspace with 6 dimensions
  !> - Write the dataset
  !> - Close dataset and dataspace
  !-----------------------------------------------------------------------------
  subroutine sll_hdf5_par_write_dble_array_6d( &
      file_id, global_size, offset, array, dsetname, error )
    integer, parameter               :: dspace_dims = 6
    integer(   hid_t), intent(in   ) :: file_id
    integer( hsize_t), intent(in   ) :: global_size(dspace_dims)
    integer(hssize_t), intent(in   ) :: offset     (dspace_dims)
    sll_real64       , intent(in   ) :: array(:,:,:,:,:,:)
    character(len=*) , intent(in   ) :: dsetname
    sll_int32        , intent(  out) :: error

#define  DATATYPE  H5T_NATIVE_DOUBLE
#include "sll_k_hdf5_write_array.f90"

  end subroutine sll_hdf5_par_write_dble_array_6d

  !-----------------------------------------------------------------------------
  !> Read from a 6D array of float in double precision from a HDF5 file
  !> into a 6D array
  !-----------------------------------------------------------------------------
  subroutine sll_hdf5_par_read_array_6d( &
      file_id, global_size, offset, array, dsetname, error )
    integer, parameter           :: dspace_dims=6
    character(len=*), intent(in) :: dsetname
    integer(hsize_t), intent(in) :: global_size(dspace_dims)
    sll_real64, intent(inout)    :: array(:,:,:,:,:,:)
    integer(hid_t)               :: file_id
    integer(hid_t)               :: plist_id
    sll_int32                    :: rank, i
    integer(hsize_t)             :: array_dims(dspace_dims)
    integer(hid_t)               :: dset_id
    integer(hid_t)               :: memspace
    integer(hid_t)               :: filespace
    integer(HSIZE_T)             :: dimsfi(dspace_dims)
    integer(HSIZE_T)             :: count(dspace_dims)
    integer(HSSIZE_T)            :: offset(dspace_dims)
    integer(HSIZE_T)             :: stride(dspace_dims)
    integer(HSIZE_T)             :: block(dspace_dims)
    sll_int32                    :: error

    rank = dspace_dims
    do i = 1, rank
      array_dims(i) = int(size(array,i),HSIZE_T)
    end do
    dimsfi = global_size
    call h5screate_simple_f(rank, global_size, filespace, error)
    SLL_ASSERT(error==0)
    call h5screate_simple_f(rank, array_dims, memspace, error)
    SLL_ASSERT(error==0)
    call h5pcreate_f(H5P_DATASET_CREATE_F, plist_id, error)
    SLL_ASSERT(error==0)
    call h5pset_chunk_f(plist_id, rank, array_dims, error)
    SLL_ASSERT(error==0)
    !call h5dcreate_f(file_id, dsetname, H5T_NATIVE_DOUBLE, filespace,  &
    !                 dset_id, error, plist_id)
    call h5dopen_f(file_id,dsetname,dset_id,error)
    SLL_ASSERT(error==0)
    call h5pclose_f(plist_id, error)
    SLL_ASSERT(error==0)
    call h5sclose_f(filespace, error)
    SLL_ASSERT(error==0)
    stride = 1
    count  = 1
    block  = array_dims
    call h5dget_space_f(dset_id, filespace, error)
    SLL_ASSERT(error==0)
    call h5sselect_hyperslab_f (filespace, H5S_SELECT_SET_F,           &
                                offset, count, error, stride, block)
    SLL_ASSERT(error==0)
    call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
    SLL_ASSERT(error==0)
    call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)
    SLL_ASSERT(error==0)
    call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, array, dimsfi, error,  &
                    file_space_id = filespace, mem_space_id = memspace,&
                    xfer_prp = plist_id)
    SLL_ASSERT(error==0)
    call h5pclose_f(plist_id, error)
    SLL_ASSERT(error==0)
    call h5sclose_f(filespace, error)
    SLL_ASSERT(error==0)
    call h5sclose_f(memspace, error)
    SLL_ASSERT(error==0)
    call h5dclose_f(dset_id, error)
    SLL_ASSERT(error==0)

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

end module sll_m_hdf5_io_parallel

