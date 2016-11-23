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
!> @brief   Parallel version of sll_hdf5_io
!> @author  Pierre Navaro, INRIA
!> @author  Yaman Güçlü, IPP Garching
!> @details
!> With HDF5 you can store several datasets in a single file.
!> - HDF5 file (http://www.hdfgroup.org/HDF5/)
module sll_m_hdf5_io_parallel
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#ifndef NOHDF5

#include "sll_assert.h"
#include "sll_working_precision.h"

  use hdf5, only: &
    h5close_f, &
    h5dclose_f, &
    h5dcreate_f, &
    h5dget_space_f, &
    h5dopen_f, &
    h5dread_f, &
    h5dwrite_f, &
    h5f_acc_rdonly_f, &
    h5f_acc_trunc_f, &
    h5fclose_f, &
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

  use sll_mpi, only: &
    mpi_info_null

  implicit none

  public :: &
    sll_t_hdf5_par_handle,      &
    sll_s_hdf5_par_file_create, &
    sll_s_hdf5_par_file_open,   &
    sll_s_hdf5_par_file_close,  &
    sll_o_hdf5_par_write_array, &
    sll_o_hdf5_par_read_array

  private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  !> Opaque object around HDF5 file id
  type :: sll_t_hdf5_par_handle
    integer(hid_t), private :: file_id
  end type

  !-----------------------------------------------------------------------------
  !> @brief
  !> Collectively write distributed nD array into HDF5 file
  !>
  !> @detail
  !> Write distributed nD Fortran array of real(f64) or integer(i32) into HDF5 file
  !> Each process writes its own data block into a global HDF5 dataset
  !>
  !> @param[in]  handle       parallel file handle
  !> @param[in]  global_size  global shape of distributed nD array
  !> @param[in]  offset       offset of local data block within global array
  !> @param[in]  array        local data block (nD array) written by one process
  !> @param[in]  dsetname     HDF5 dataset name
  !> @param[out] error        HDF5 error code
  !> @param[in]  chunk_dims   shape of HDF5 chunks (CHUNKED storage layout)
  !-----------------------------------------------------------------------------
  interface sll_o_hdf5_par_write_array
    module procedure sll_hdf5_par_write_dble_array_1d
    module procedure sll_hdf5_par_write_dble_array_2d
    module procedure sll_hdf5_par_write_dble_array_3d
    module procedure sll_hdf5_par_write_dble_array_4d
    module procedure sll_hdf5_par_write_dble_array_5d
    module procedure sll_hdf5_par_write_dble_array_6d
  end interface

  !-----------------------------------------------------------------------------
  !> @brief
  !> Collectively read distributed nD array from HDF5 file
  !>
  !> @detail
  !> Read distributed nD Fortran array of real(f64) or integer(i32) from HDF5 file
  !> Each process reads its own data block from a global HDF5 dataset
  !>
  !> @param[in]  handle       parallel file handle
  !> @param[in]  global_size  global shape of distributed nD array
  !> @param[in]  offset       offset of local data block within global array
  !> @param[out] array        local data block (nD array) read by one process
  !> @param[in]  dsetname     HDF5 dataset name
  !> @param[out] error        HDF5 error code
  !-----------------------------------------------------------------------------
  interface sll_o_hdf5_par_read_array
    module procedure sll_hdf5_par_read_dble_array_1d
    module procedure sll_hdf5_par_read_dble_array_2d
    module procedure sll_hdf5_par_read_dble_array_3d
    module procedure sll_hdf5_par_read_dble_array_4d
    module procedure sll_hdf5_par_read_dble_array_5d
    module procedure sll_hdf5_par_read_dble_array_6d
  end interface sll_o_hdf5_par_read_array

contains

  !-----------------------------------------------------------------------------
  !> Create new HDF5 file
  !>    - Initialize fortran interface
  !>    - Create a new file using default properties
  !>
  !> @param[in]  filename  file name
  !> @param[in]  comm      MPI communicator
  !> @param[out] handle    parallel file handle
  !> @param[out] error     HDF5 error code
  !-----------------------------------------------------------------------------
  subroutine sll_s_hdf5_par_file_create( filename, comm, handle, error )
    character(len=*)           , intent(in   ) :: filename
    integer                    , intent(in   ) :: comm
    type(sll_t_hdf5_par_handle), intent(  out) :: handle
    integer                    , intent(  out) :: error

    integer(hid_t) :: plist_id
    integer        :: info
    info = mpi_info_null

    call h5open_f(error) 
    SLL_ASSERT(error==0)
    call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, error)
    SLL_ASSERT(error==0)
    call h5pset_fapl_mpio_f(plist_id, comm, info, error)
    SLL_ASSERT(error==0)
    call h5fcreate_f(filename, H5F_ACC_TRUNC_F, handle%file_id, error, access_prp = plist_id)
    SLL_ASSERT(error==0)
    call h5pclose_f(plist_id, error)
    SLL_ASSERT(error==0)

  end subroutine sll_s_hdf5_par_file_create

  !-----------------------------------------------------------------------------
  !> Open existing HDF5 file
  !>    - Initialize fortran interface
  !>    - Open a HDF5 file
  !>
  !> @param[in]  filename  file name
  !> @param[in]  comm      MPI communicator
  !> @param[out] handle    parallel file handle
  !> @param[out] error     HDF5 error code
  !-----------------------------------------------------------------------------
  subroutine sll_s_hdf5_par_file_open( filename, comm, handle, error )
    character(len=*)           , intent(in   ) :: filename
    integer                    , intent(in   ) :: comm
    type(sll_t_hdf5_par_handle), intent(  out) :: handle
    integer                    , intent(  out) :: error

    integer(hid_t) :: plist_id
    integer        :: info
    info = mpi_info_null

    call h5open_f(error) 
    SLL_ASSERT(error==0)
    call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, error)
    SLL_ASSERT(error==0)
    call h5pset_fapl_mpio_f(plist_id, comm, info, error)
    SLL_ASSERT(error==0)
    call h5fopen_f(filename, H5F_ACC_RDONLY_F, handle%file_id, error, access_prp = plist_id)
    SLL_ASSERT(error==0)
    call h5pclose_f(plist_id, error)
    SLL_ASSERT(error==0)

  end subroutine sll_s_hdf5_par_file_open

  !-----------------------------------------------------------------------------
  !> Close existing HDF5 file
  !>
  !> @param[in]  handle  parallel file handle
  !> @param[out] error   HDF5 error code
  !-----------------------------------------------------------------------------
  subroutine sll_s_hdf5_par_file_close( handle, error )
    type(sll_t_hdf5_par_handle), intent(in   ) :: handle
    integer                    , intent(  out) :: error

    call h5fclose_f( handle%file_id, error ) ! Close property list and file
    SLL_ASSERT(error==0)
    call h5close_f( error )                  ! Close FORTRAN interface
    SLL_ASSERT(error==0)
  end subroutine sll_s_hdf5_par_file_close

  !-----------------------------------------------------------------------------
  !> Write 1D array of double precision floats into HDF5 file
  !-----------------------------------------------------------------------------
  subroutine sll_hdf5_par_write_dble_array_1d( &
      handle, global_size, offset, array, dsetname, error, chunk_dims )
    integer, parameter                         :: dspace_dims = 1
    type(sll_t_hdf5_par_handle), intent(in   ) :: handle
    integer(i64)               , intent(in   ) :: global_size(dspace_dims)
    integer(i64)               , intent(in   ) :: offset     (dspace_dims)
    real(f64)                  , intent(in   ) :: array(:)
    character(len=*)           , intent(in   ) :: dsetname
    integer                    , intent(  out) :: error
    integer(i64),     optional , intent(in   ) :: chunk_dims(dspace_dims)

#define  DATATYPE  H5T_NATIVE_DOUBLE
#include "sll_k_hdf5_par_write_array.F90"
#undef   DATATYPE

  end subroutine sll_hdf5_par_write_dble_array_1d

  !-----------------------------------------------------------------------------
  !> Write 2D array of double precision floats into HDF5 file
  !-----------------------------------------------------------------------------
  subroutine sll_hdf5_par_write_dble_array_2d( &
      handle, global_size, offset, array, dsetname, error, chunk_dims )
    integer, parameter                         :: dspace_dims = 2
    type(sll_t_hdf5_par_handle), intent(in   ) :: handle
    integer(i64)               , intent(in   ) :: global_size(dspace_dims)
    integer(i64)               , intent(in   ) :: offset     (dspace_dims)
    real(f64)                  , intent(in   ) :: array(:,:)
    character(len=*)           , intent(in   ) :: dsetname
    integer                    , intent(  out) :: error
    integer(i64),     optional , intent(in   ) :: chunk_dims(dspace_dims)

#define  DATATYPE  H5T_NATIVE_DOUBLE
#include "sll_k_hdf5_par_write_array.F90"
#undef   DATATYPE

  end subroutine sll_hdf5_par_write_dble_array_2d

  !-----------------------------------------------------------------------------
  !> Write 3D array of double precision floats into HDF5 file
  !-----------------------------------------------------------------------------
  subroutine sll_hdf5_par_write_dble_array_3d( &
      handle, global_size, offset, array, dsetname, error, chunk_dims )
    integer, parameter                         :: dspace_dims = 3
    type(sll_t_hdf5_par_handle), intent(in   ) :: handle
    integer(i64)               , intent(in   ) :: global_size(dspace_dims)
    integer(i64)               , intent(in   ) :: offset     (dspace_dims)
    real(f64)                  , intent(in   ) :: array(:,:,:)
    character(len=*)           , intent(in   ) :: dsetname
    integer                    , intent(  out) :: error
    integer(i64),     optional , intent(in   ) :: chunk_dims(dspace_dims)

#define  DATATYPE  H5T_NATIVE_DOUBLE
#include "sll_k_hdf5_par_write_array.F90"
#undef   DATATYPE

  end subroutine sll_hdf5_par_write_dble_array_3d

  !-----------------------------------------------------------------------------
  !> Write 4D array of double precision floats into HDF5 file
  !-----------------------------------------------------------------------------
  subroutine sll_hdf5_par_write_dble_array_4d( &
      handle, global_size, offset, array, dsetname, error, chunk_dims )
    integer, parameter                         :: dspace_dims = 4
    type(sll_t_hdf5_par_handle), intent(in   ) :: handle
    integer(i64)               , intent(in   ) :: global_size(dspace_dims)
    integer(i64)               , intent(in   ) :: offset     (dspace_dims)
    real(f64)                  , intent(in   ) :: array(:,:,:,:)
    character(len=*)           , intent(in   ) :: dsetname
    integer                    , intent(  out) :: error
    integer(i64),     optional , intent(in   ) :: chunk_dims(dspace_dims)

#define  DATATYPE  H5T_NATIVE_DOUBLE
#include "sll_k_hdf5_par_write_array.F90"
#undef   DATATYPE

  end subroutine sll_hdf5_par_write_dble_array_4d

  !-----------------------------------------------------------------------------
  !> Write 5D array of double precision floats into HDF5 file
  !-----------------------------------------------------------------------------
  subroutine sll_hdf5_par_write_dble_array_5d( &
      handle, global_size, offset, array, dsetname, error, chunk_dims )
    integer, parameter                         :: dspace_dims = 5
    type(sll_t_hdf5_par_handle), intent(in   ) :: handle
    integer(i64)               , intent(in   ) :: global_size(dspace_dims)
    integer(i64)               , intent(in   ) :: offset     (dspace_dims)
    real(f64)                  , intent(in   ) :: array(:,:,:,:,:)
    character(len=*)           , intent(in   ) :: dsetname
    integer                    , intent(  out) :: error
    integer(i64),     optional , intent(in   ) :: chunk_dims(dspace_dims)

#define  DATATYPE  H5T_NATIVE_DOUBLE
#include "sll_k_hdf5_par_write_array.F90"
#undef   DATATYPE

  end subroutine sll_hdf5_par_write_dble_array_5d

  !-----------------------------------------------------------------------------
  !> Write 6D array of double precision floats into HDF5 file
  !-----------------------------------------------------------------------------
  subroutine sll_hdf5_par_write_dble_array_6d( &
      handle, global_size, offset, array, dsetname, error, chunk_dims )
    integer, parameter                         :: dspace_dims = 6
    type(sll_t_hdf5_par_handle), intent(in   ) :: handle
    integer(i64)               , intent(in   ) :: global_size(dspace_dims)
    integer(i64)               , intent(in   ) :: offset     (dspace_dims)
    real(f64)                  , intent(in   ) :: array(:,:,:,:,:,:)
    character(len=*)           , intent(in   ) :: dsetname
    integer                    , intent(  out) :: error
    integer(i64),     optional , intent(in   ) :: chunk_dims(dspace_dims)

#define  DATATYPE  H5T_NATIVE_DOUBLE
#include "sll_k_hdf5_par_write_array.F90"
#undef   DATATYPE

  end subroutine sll_hdf5_par_write_dble_array_6d

  !-----------------------------------------------------------------------------
  !> Read 1D array of double precision floats from HDF5 file
  !-----------------------------------------------------------------------------
  subroutine sll_hdf5_par_read_dble_array_1d( &
      handle, global_size, offset, array, dsetname, error )
    integer, parameter                         :: dspace_dims = 1
    type(sll_t_hdf5_par_handle), intent(in   ) :: handle
    integer(i64)               , intent(in   ) :: global_size(dspace_dims)
    integer(i64)               , intent(in   ) :: offset(dspace_dims)
    real(f64)                  , intent(  out) :: array(:)
    character(len=*)           , intent(in   ) :: dsetname
    integer                    , intent(  out) :: error

#define  DATATYPE  H5T_NATIVE_DOUBLE
#include "sll_k_hdf5_par_read_array.F90"
#undef   DATATYPE

  end subroutine sll_hdf5_par_read_dble_array_1d

  !-----------------------------------------------------------------------------
  !> Read 2D array of double precision floats from HDF5 file
  !-----------------------------------------------------------------------------
  subroutine sll_hdf5_par_read_dble_array_2d( &
      handle, global_size, offset, array, dsetname, error )
    integer, parameter                         :: dspace_dims = 2
    type(sll_t_hdf5_par_handle), intent(in   ) :: handle
    integer(i64)               , intent(in   ) :: global_size(dspace_dims)
    integer(i64)               , intent(in   ) :: offset(dspace_dims)
    real(f64)                  , intent(  out) :: array(:,:)
    character(len=*)           , intent(in   ) :: dsetname
    integer                    , intent(  out) :: error

#define  DATATYPE  H5T_NATIVE_DOUBLE
#include "sll_k_hdf5_par_read_array.F90"
#undef   DATATYPE

  end subroutine sll_hdf5_par_read_dble_array_2d

  !-----------------------------------------------------------------------------
  !> Read 3D array of double precision floats from HDF5 file
  !-----------------------------------------------------------------------------
  subroutine sll_hdf5_par_read_dble_array_3d( &
      handle, global_size, offset, array, dsetname, error )
    integer, parameter                         :: dspace_dims = 3
    type(sll_t_hdf5_par_handle), intent(in   ) :: handle
    integer(i64)               , intent(in   ) :: global_size(dspace_dims)
    integer(i64)               , intent(in   ) :: offset(dspace_dims)
    real(f64)                  , intent(  out) :: array(:,:,:)
    character(len=*)           , intent(in   ) :: dsetname
    integer                    , intent(  out) :: error

#define  DATATYPE  H5T_NATIVE_DOUBLE
#include "sll_k_hdf5_par_read_array.F90"
#undef   DATATYPE

  end subroutine sll_hdf5_par_read_dble_array_3d

  !-----------------------------------------------------------------------------
  !> Read 4D array of double precision floats from HDF5 file
  !-----------------------------------------------------------------------------
  subroutine sll_hdf5_par_read_dble_array_4d( &
      handle, global_size, offset, array, dsetname, error )
    integer, parameter                         :: dspace_dims = 4
    type(sll_t_hdf5_par_handle), intent(in   ) :: handle
    integer(i64)               , intent(in   ) :: global_size(dspace_dims)
    integer(i64)               , intent(in   ) :: offset(dspace_dims)
    real(f64)                  , intent(  out) :: array(:,:,:,:)
    character(len=*)           , intent(in   ) :: dsetname
    integer                    , intent(  out) :: error

#define  DATATYPE  H5T_NATIVE_DOUBLE
#include "sll_k_hdf5_par_read_array.F90"
#undef   DATATYPE

  end subroutine sll_hdf5_par_read_dble_array_4d

  !-----------------------------------------------------------------------------
  !> Read 5D array of double precision floats from HDF5 file
  !-----------------------------------------------------------------------------
  subroutine sll_hdf5_par_read_dble_array_5d( &
      handle, global_size, offset, array, dsetname, error )
    integer, parameter                         :: dspace_dims = 5
    type(sll_t_hdf5_par_handle), intent(in   ) :: handle
    integer(i64)               , intent(in   ) :: global_size(dspace_dims)
    integer(i64)               , intent(in   ) :: offset(dspace_dims)
    real(f64)                  , intent(  out) :: array(:,:,:,:,:)
    character(len=*)           , intent(in   ) :: dsetname
    integer                    , intent(  out) :: error

#define  DATATYPE  H5T_NATIVE_DOUBLE
#include "sll_k_hdf5_par_read_array.F90"
#undef   DATATYPE

  end subroutine sll_hdf5_par_read_dble_array_5d

  !-----------------------------------------------------------------------------
  !> Read 6D array of double precision floats from HDF5 file
  !-----------------------------------------------------------------------------
  subroutine sll_hdf5_par_read_dble_array_6d( &
      handle, global_size, offset, array, dsetname, error )
    integer, parameter                         :: dspace_dims = 6
    type(sll_t_hdf5_par_handle), intent(in   ) :: handle
    integer(i64)               , intent(in   ) :: global_size(dspace_dims)
    integer(i64)               , intent(in   ) :: offset(dspace_dims)
    real(f64)                  , intent(  out) :: array(:,:,:,:,:,:)
    character(len=*)           , intent(in   ) :: dsetname
    integer                    , intent(  out) :: error

#define  DATATYPE  H5T_NATIVE_DOUBLE
#include "sll_k_hdf5_par_read_array.F90"
#undef   DATATYPE

  end subroutine sll_hdf5_par_read_dble_array_6d


#endif

end module sll_m_hdf5_io_parallel
