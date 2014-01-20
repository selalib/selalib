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
!> Implements the functions to write hdf5 file to store heavy data.
!> @details
!> With HDF5 you can store several datasets in a single file.
!> @note 
!> link with sll_file_io
module sll_hdf5_io
#include "sll_working_precision.h"
  
#ifndef NOHDF5

use hdf5, only: hid_t, hsize_t, H5Screate_simple_f, h5t_native_double, &
                h5fopen_f, H5F_ACC_RDONLY_F, H5F_ACC_TRUNC_F, &
                h5dcreate_f, h5dwrite_f,h5dclose_f, h5open_f,&
                h5screate_f, h5sclose_f, h5fclose_f, h5fcreate_f
  
  implicit none
  
  !> Write a nD array of float in double precision in a HDF5 file 
  !>\param[in]  file_id file unit number
  !>\param[in]  array array
  !>\param[in]  dsetname dataset name
  !>\param[out] error dataset error code
  interface sll_hdf5_write_array
     module procedure sll_hdf5_write_array_1d
     module procedure sll_hdf5_write_array_2d
     module procedure sll_hdf5_write_array_3d
  end interface

contains
  
  !> Create HDF5 file
  subroutine sll_hdf5_file_create(filename,file_id,error)
    character(len=*) , intent(in)  :: filename  !< file name
    integer(hid_t)   , intent(out) :: file_id   !< unit number
    integer,           intent(out) :: error     !< error code
    
    call H5open_f(error)
    call H5Fcreate_f(filename,H5F_ACC_TRUNC_F,file_id,error)
  end subroutine sll_hdf5_file_create

  !> Open HDF5 file
  subroutine sll_hdf5_file_open(filename,file_id,error)
    character(len=*) , intent(in)  :: filename  !< file name
    integer(hid_t)   , intent(out) :: file_id   !< unit number
    integer,           intent(out) :: error     !< error code
    
    call H5open_f(error)
    call H5Fopen_f(filename,H5F_ACC_RDONLY_F,file_id,error)
  end subroutine sll_hdf5_file_open

  !> Close HDF5 file 
  subroutine sll_hdf5_file_close(file_id,error)
    integer(hid_t), intent(in) :: file_id  !< file unit number
    integer, intent(out) :: error          !< error code

    call H5Fclose_f(file_id,error)
  end subroutine sll_hdf5_file_close

  !> Write a 1d array in hdf5 file
  subroutine sll_hdf5_write_array_1d(file_id,array,dsetname,error)
    integer(hid_t)  , intent(in) :: file_id      !< file unit number
    character(len=*), intent(in) :: dsetname     !< hdf5 dataset name
    sll_int32, intent(out)       :: error        !< error code
    sll_int32                    :: rank, i
    sll_real64, intent(in)       :: array(:)     !< data array
    integer(hsize_t)             :: array_dims(1)
    integer(hid_t)               :: dataset_id
    integer(hid_t)               :: dataspace_id
    rank = 1      
    do i = 1, rank
      array_dims(i) = size(array,i)
    end do
    call H5Screate_simple_f(rank,array_dims,dataspace_id,error)
    call H5Dcreate_f(file_id,                                            &
                     dsetname,                                           &
                     H5T_NATIVE_DOUBLE,                                  &
                     dataspace_id,                                       &
                     dataset_id,                                         &
                     error)
    call H5Dwrite_f(dataset_id,H5T_NATIVE_DOUBLE,array,array_dims,error)
    call H5Sclose_f(dataspace_id,error)
    call H5Dclose_f(dataset_id,error)
  end subroutine sll_hdf5_write_array_1d

  !> Write a 2d array in hdf5 file
  subroutine sll_hdf5_write_array_2d(file_id,array,dsetname,error)
    integer(hid_t)  , intent(in) :: file_id      !< file unit number
    character(len=*), intent(in) :: dsetname     !< hdf5 dataset name
    sll_int32, intent(out)       :: error        !< error code
    sll_int32                    :: rank, i
    sll_real64, intent(in)       :: array(:,:)   !< data array
    integer(hsize_t)             :: array_dims(2)
    integer(hid_t)               :: dataset_id
    integer(hid_t)               :: dataspace_id
    rank = 2
    do i = 1, rank
      array_dims(i) = size(array,i)
    end do
    call H5Screate_simple_f(rank,array_dims,dataspace_id,error)
    call H5Dcreate_f(file_id,                                            &
                     dsetname,                                           &
                     H5T_NATIVE_DOUBLE,                                  &
                     dataspace_id,                                       &
                     dataset_id,                                         &
                     error)
    call H5Dwrite_f(dataset_id,H5T_NATIVE_DOUBLE,array,array_dims,error)
    call H5Sclose_f(dataspace_id,error)
    call H5Dclose_f(dataset_id,error)
  end subroutine sll_hdf5_write_array_2d

  !> Write a 3d array in hdf5 file
  subroutine sll_hdf5_write_array_3d(file_id,array,dsetname,error)
    integer(hid_t)  , intent(in) :: file_id      !< file unit number
    character(len=*), intent(in) :: dsetname     !< hdf5 dataset name
    sll_int32, intent(out)       :: error        !< error code
    sll_int32                    :: rank, i
    sll_real64, intent(in)       :: array(:,:,:) !< data array
    integer(hsize_t)             :: array_dims(3)
    integer(hid_t)               :: dataset_id
    integer(hid_t)               :: dataspace_id
    rank = 3
    do i = 1, rank
      array_dims(i) = size(array,i)
    end do
    call H5Screate_simple_f(rank,array_dims,dataspace_id,error)
    call H5Dcreate_f(file_id,                                            &
                     dsetname,                                           &
                     H5T_NATIVE_DOUBLE,                                  &
                     dataspace_id,                                       &
                     dataset_id,                                         &
                     error)
    call H5Dwrite_f(dataset_id,H5T_NATIVE_DOUBLE,array,array_dims,error)
    call H5Sclose_f(dataspace_id,error)
    call H5Dclose_f(dataset_id,error)
  end subroutine sll_hdf5_write_array_3d

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

end module sll_hdf5_io

