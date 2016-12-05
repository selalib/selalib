!**************************************************************
!  Copyright INRIA
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

!> @ingroup file_io
!> @brief
!> Implements the functions to write hdf5 file to store heavy data.
!> @author  Pierre Navaro, INRIA
!> @author  Yaman Güçlü, IPP Garching
!> @details
!> With HDF5 you can store several datasets in a single file.
!> @note 
!> link with sll_file_io
module sll_m_hdf5_io_serial
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#ifndef NOHDF5

#include "sll_assert.h"
#include "sll_working_precision.h"

  use hdf5, only: &
    h5close_f, &
    h5dclose_f, &
    h5dcreate_f, &
    h5dopen_f, &
    h5dread_f, &
    h5dwrite_f, &
    h5f_acc_rdonly_f, &
    h5f_acc_rdwr_f, &
    h5f_acc_trunc_f, &
    h5fclose_f, &
    h5fcreate_f, &
    h5fopen_f, &
    h5open_f, &
    h5sclose_f, &
    h5screate_f, &
    h5screate_simple_f, &
    h5t_native_double, &
    h5t_native_integer, &
    h5t_native_character, &
    hid_t, &
    hsize_t

  implicit none

  public :: &
    sll_t_hdf5_ser_handle,      &
    sll_s_hdf5_ser_file_create, &
    sll_s_hdf5_ser_file_open,   &
    sll_s_hdf5_ser_file_close,  &
    sll_o_hdf5_ser_write_array, &
    sll_o_hdf5_ser_read_array,  &
    sll_s_hdf5_ser_write_file

  private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  !> Opaque object around HDF5 file id
  type :: sll_t_hdf5_ser_handle
    integer(hid_t), private :: file_id
  end type

  !-----------------------------------------------------------------------------
  !> @brief
  !> Write nD array of double precision floats or integers into HDF5 file
  !>
  !> @param[in]  handle    file handle
  !> @param[in]  array     multi-dimensional array
  !> @param[in]  dsetname  HDF5 dataset name
  !> @param[out] error     HDF5 error code
  !-----------------------------------------------------------------------------
  interface sll_o_hdf5_ser_write_array
     module procedure sll_hdf5_ser_write_dble_array_1d
     module procedure sll_hdf5_ser_write_dble_array_2d
     module procedure sll_hdf5_ser_write_dble_array_3d
     module procedure sll_hdf5_ser_write_int_array_1d
     module procedure sll_hdf5_ser_write_int_array_2d
     module procedure sll_hdf5_ser_write_int_array_3d
  end interface

  !-----------------------------------------------------------------------------
  !> @brief
  !> Read nD array of double precision floats or integers from HDF5 file
  !>
  !> @param[in]  handle    file handle
  !> @param[out] array     multi-dimensional array
  !> @param[in]  dsetname  HDF5 dataset name
  !> @param[out] error     HDF5 error code
  !-----------------------------------------------------------------------------
  interface sll_o_hdf5_ser_read_array
    module procedure sll_hdf5_ser_read_dble_array_1d
    module procedure sll_hdf5_ser_read_dble_array_2d
    module procedure sll_hdf5_ser_read_dble_array_3d
  end interface

contains
  
  !-----------------------------------------------------------------------------
  !> @brief Create new HDF5 file
  !> @details To use this subroutine HDF5_ENABLE should be set to ON 
  !> in CMake configuration
  !>
  !> @param[in]  filename  file name
  !> @param[out] handle    file handle
  !> @param[out] error     HDF5 error code
  !-----------------------------------------------------------------------------
  subroutine sll_s_hdf5_ser_file_create( filename, handle, error )
    character(len=*)           , intent(in   ) :: filename  !< file name
    type(sll_t_hdf5_ser_handle), intent(  out) :: handle    !< file handle
    integer                    , intent(  out) :: error     !< error code
    
    call h5open_f(error)
    SLL_ASSERT(error==0)
    call h5fcreate_f( filename, H5F_ACC_TRUNC_F, handle%file_id, error )
    SLL_ASSERT(error==0)
  end subroutine sll_s_hdf5_ser_file_create

  !-----------------------------------------------------------------------------
  !> Open existing HDF5 file
  !>
  !> @param[in]  filename  file name
  !> @param[out] handle    file handle
  !> @param[out] error     HDF5 error code
  !-----------------------------------------------------------------------------
  subroutine sll_s_hdf5_ser_file_open( filename, handle, error )
    character(len=*)           , intent(in   ) :: filename  !< file name
    type(sll_t_hdf5_ser_handle), intent(  out) :: handle    !< file handle
    integer                    , intent(  out) :: error     !< error code
    
    call h5open_f(error)
    SLL_ASSERT(error==0)
    call h5fopen_f( filename, H5F_ACC_RDWR_F, handle%file_id, error )
    SLL_ASSERT(error==0)
  end subroutine sll_s_hdf5_ser_file_open

  !-----------------------------------------------------------------------------
  !> Close existing HDF5 file
  !>
  !> @param[in]  handle  file handle
  !> @param[out] error   HDF5 error code
  !-----------------------------------------------------------------------------
  subroutine sll_s_hdf5_ser_file_close( handle, error )
    type(sll_t_hdf5_ser_handle), intent(in   ) :: handle
    integer                    , intent(  out) :: error

    call h5fclose_f( handle%file_id, error ) ! Close property list and file
    SLL_ASSERT(error==0)
    call h5close_f( error )                  ! Close FORTRAN interface
    SLL_ASSERT(error==0)
  end subroutine sll_s_hdf5_ser_file_close

  !-----------------------------------------------------------------------------
  !> Write 1D array of int32 into HDF5 file
  !-----------------------------------------------------------------------------
  subroutine sll_hdf5_ser_write_int_array_1d( handle, array, dsetname, error )
    integer, parameter                         :: rank = 1
    type(sll_t_hdf5_ser_handle), intent(in   ) :: handle
    integer(i32)               , intent(in   ) :: array(:)
    character(len=*)           , intent(in   ) :: dsetname
    integer                    , intent(  out) :: error

#define  DATATYPE  H5T_NATIVE_INTEGER
#include "sll_k_hdf5_ser_write_array.F90"
#undef   DATATYPE

  end subroutine sll_hdf5_ser_write_int_array_1d

  !-----------------------------------------------------------------------------
  !> Write 2D array of int32 into HDF5 file
  !-----------------------------------------------------------------------------
  subroutine sll_hdf5_ser_write_int_array_2d( handle, array, dsetname, error )
    integer, parameter                         :: rank = 2
    type(sll_t_hdf5_ser_handle), intent(in   ) :: handle
    integer(i32)               , intent(in   ) :: array(:,:)
    character(len=*)           , intent(in   ) :: dsetname
    integer                    , intent(  out) :: error

#define  DATATYPE  H5T_NATIVE_INTEGER
#include "sll_k_hdf5_ser_write_array.F90"
#undef   DATATYPE

  end subroutine sll_hdf5_ser_write_int_array_2d

  !-----------------------------------------------------------------------------
  !> Write 3D array of int32 into HDF5 file
  !-----------------------------------------------------------------------------
  subroutine sll_hdf5_ser_write_int_array_3d( handle, array, dsetname, error )
    integer, parameter                         :: rank = 3
    type(sll_t_hdf5_ser_handle), intent(in   ) :: handle
    integer(i32)               , intent(in   ) :: array(:,:,:)
    character(len=*)           , intent(in   ) :: dsetname
    integer                    , intent(  out) :: error

#define  DATATYPE  H5T_NATIVE_INTEGER
#include "sll_k_hdf5_ser_write_array.F90"
#undef   DATATYPE

  end subroutine sll_hdf5_ser_write_int_array_3d

  !-----------------------------------------------------------------------------
  !> Write 1D array of float64 into HDF5 file
  !-----------------------------------------------------------------------------
  subroutine sll_hdf5_ser_write_dble_array_1d( handle, array, dsetname, error )
    integer, parameter                         :: rank = 1
    type(sll_t_hdf5_ser_handle), intent(in   ) :: handle
    real(f64)                  , intent(in   ) :: array(:)
    character(len=*)           , intent(in   ) :: dsetname
    integer                    , intent(  out) :: error

#define  DATATYPE  H5T_NATIVE_DOUBLE
#include "sll_k_hdf5_ser_write_array.F90"
#undef   DATATYPE

  end subroutine sll_hdf5_ser_write_dble_array_1d

  !-----------------------------------------------------------------------------
  !> Write 2D array of float64 into HDF5 file
  !-----------------------------------------------------------------------------
  subroutine sll_hdf5_ser_write_dble_array_2d( handle, array, dsetname, error )
    integer, parameter                         :: rank = 2
    type(sll_t_hdf5_ser_handle), intent(in   ) :: handle
    real(f64)                  , intent(in   ) :: array(:,:)
    character(len=*)           , intent(in   ) :: dsetname
    integer                    , intent(  out) :: error

#define  DATATYPE  H5T_NATIVE_DOUBLE
#include "sll_k_hdf5_ser_write_array.F90"
#undef   DATATYPE

  end subroutine sll_hdf5_ser_write_dble_array_2d

  !-----------------------------------------------------------------------------
  !> Write 3D array of float64 into HDF5 file
  !-----------------------------------------------------------------------------
  subroutine sll_hdf5_ser_write_dble_array_3d( handle, array, dsetname, error )
    integer, parameter                         :: rank = 3
    type(sll_t_hdf5_ser_handle), intent(in   ) :: handle
    real(f64)                  , intent(in   ) :: array(:,:,:)
    character(len=*)           , intent(in   ) :: dsetname
    integer                    , intent(  out) :: error

#define  DATATYPE  H5T_NATIVE_DOUBLE
#include "sll_k_hdf5_ser_write_array.F90"
#undef   DATATYPE

  end subroutine sll_hdf5_ser_write_dble_array_3d

  !-----------------------------------------------------------------------------
  !> Read 1D array of float64 from HDF5 file
  !-----------------------------------------------------------------------------
  subroutine sll_hdf5_ser_read_dble_array_1d( handle, array, dsetname, error )
    integer, parameter                         :: rank = 1
    type(sll_t_hdf5_ser_handle), intent(in   ) :: handle
    real(f64)                  , intent(  out) :: array(:)
    character(len=*)           , intent(in   ) :: dsetname
    integer                    , intent(  out) :: error

#define  DATATYPE  H5T_NATIVE_DOUBLE
#include "sll_k_hdf5_ser_read_array.F90"
#undef   DATATYPE

  end subroutine sll_hdf5_ser_read_dble_array_1d

  !-----------------------------------------------------------------------------
  !> Read 2D array of float64 from HDF5 file
  !-----------------------------------------------------------------------------
  subroutine sll_hdf5_ser_read_dble_array_2d( handle, array, dsetname, error )
    integer, parameter                         :: rank = 2
    type(sll_t_hdf5_ser_handle), intent(in   ) :: handle
    real(f64)                  , intent(  out) :: array(:,:)
    character(len=*)           , intent(in   ) :: dsetname
    integer                    , intent(  out) :: error

#define  DATATYPE  H5T_NATIVE_DOUBLE
#include "sll_k_hdf5_ser_read_array.F90"
#undef   DATATYPE

  end subroutine sll_hdf5_ser_read_dble_array_2d

  !-----------------------------------------------------------------------------
  !> Read 3D array of float64 from HDF5 file
  !-----------------------------------------------------------------------------
  subroutine sll_hdf5_ser_read_dble_array_3d( handle, array, dsetname, error )
    integer, parameter                         :: rank = 3
    type(sll_t_hdf5_ser_handle), intent(in   ) :: handle
    real(f64)                  , intent(  out) :: array(:,:,:)
    character(len=*)           , intent(in   ) :: dsetname
    integer                    , intent(  out) :: error

#define  DATATYPE  H5T_NATIVE_DOUBLE
#include "sll_k_hdf5_ser_read_array.F90"
#undef   DATATYPE

  end subroutine sll_hdf5_ser_read_dble_array_3d

  !-----------------------------------------------------------------------------
  !> Write Fortran string to HDF5 file as 1D array of characters
  !>
  !> @param[in]  handle    file handle
  !> @param[in]  string    data string
  !> @param[in]  dsetname  HDF5 dataset name
  !> @param[out] error     HDF5 error code
  !-----------------------------------------------------------------------------
  subroutine sll_hdf5_ser_write_char_array( handle, string, dsetname, error )
    type(sll_t_hdf5_ser_handle), intent(in   ) :: handle
    character(len=*)           , intent(in   ) :: string
    character(len=*)           , intent(in   ) :: dsetname
    integer                    , intent(  out) :: error

    integer(hsize_t) :: array_dim(1)
    integer(hid_t)   :: dataset_id
    integer(hid_t)   :: dataspace_id

    array_dim = int( len(string), hsize_t )
 
    ! Create dataspace
    call h5screate_simple_f( 1, array_dim, dataspace_id, error )
    SLL_ASSERT(error==0)

    ! Create dataset
    call h5dcreate_f( handle%file_id, dsetname, H5T_NATIVE_CHARACTER, &
                      dataspace_id, dataset_id, error )
    SLL_ASSERT(error==0)

    ! Write dataset
    call h5dwrite_f( dataset_id, H5T_NATIVE_CHARACTER, string, array_dim, error )
    SLL_ASSERT(error==0)

    ! Close dataspace
    call h5sclose_f( dataspace_id, error )
    SLL_ASSERT(error==0)

    ! Close dataset
    call h5dclose_f( dataset_id, error )
    SLL_ASSERT(error==0)
  
  end subroutine sll_hdf5_ser_write_char_array

  !-----------------------------------------------------------------------------
  !> Read complete file to string and store it into HDF5 data structure
  !>
  !> @param[in]  handle    file handle
  !> @param[in]  filename  name of file to be copied
  !> @param[in]  dsetname  HDF5 dataset name
  !> @param[out] error     HDF5 error code
  !-----------------------------------------------------------------------------
  subroutine sll_s_hdf5_ser_write_file( handle, filename, dsetname, error )
    type(sll_t_hdf5_ser_handle), intent(in   ) :: handle
    character(len=*)           , intent(in   ) :: filename
    character(len=*)           , intent(in   ) :: dsetname
    integer                    , intent(  out) :: error

    integer                       :: fileid, filesize
    character(len=:), allocatable :: content
       
    open( newunit=fileid, file=filename,&
            form='UNFORMATTED', access='STREAM', iostat=error )

    !get size of file
    inquire( file=filename, size=filesize )
      if (filesize>0) then

        ! Allocate memory
        allocate( character(len=filesize) :: content )

        ! Read the whole file into one string
        read( fileid, pos=1, iostat=error ) content 

        if (error==0) then
          !try to read more in case not everything was read
          read( fileid, pos=filesize+1, iostat=error ) content
          if (.not. IS_IOSTAT_END(error)) then
            write(*,*) "ERROR: file ", filename, "was not completely read."
          end if
        else
          write(*,*) 'ERROR reading file: ', filename
        end if

      else
        write(*,*) 'Error getting size of file: ', filename
      end if

    close( fileid )
    call sll_hdf5_ser_write_char_array( handle, content, dsetname, error )

  end subroutine sll_s_hdf5_ser_write_file  
  

#endif

end module sll_m_hdf5_io_serial
