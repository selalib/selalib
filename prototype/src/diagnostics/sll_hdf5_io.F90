module sll_hdf5_io

#ifndef NOHDF5
use hdf5
#include "sll_working_precision.h"

implicit none

contains
  
! Create HDF5 file 
subroutine sll_hdf5_file_create(filename,file_id,error)
character(len=*) , intent(in)  :: filename  
integer(hid_t)   , intent(out) :: file_id   
integer,           intent(out) :: error

! Initialize fortran interface
call H5open_f(error)

! Create a new file using default properties
call H5Fcreate_f(filename,H5F_ACC_TRUNC_F,file_id,error)
end subroutine sll_hdf5_file_create

! Open HDF5 file 
subroutine sll_hdf5_file_open(filename,file_id,error)
character(len=*) , intent(in)  :: filename  
integer(hid_t)   , intent(out) :: file_id   
integer,           intent(out) :: error

! Initialize fortran interface
call H5open_f(error)

! Open the HDF5 file
call H5Fopen_f(filename,H5F_ACC_RDONLY_F,file_id,error)
end subroutine sll_hdf5_file_open

! Close HDF5 file 
subroutine sll_hdf5_file_close(file_id,error)
integer(hid_t), intent(in) :: file_id
integer :: error

call H5Fclose_f(file_id,error)
end subroutine sll_hdf5_file_close

#define NEW_HDF5_FUNCTION(func_name, dataspace_dimension, array_name_and_dims)      \
subroutine func_name(file_id,array,dsetname,error);                                 \
integer(hid_t)  , intent(in)           :: file_id;                                  \
character(len=*), intent(in)           :: dsetname;                                 \
sll_int32, intent(out)                 :: error;                                    \
sll_int32                              :: rank, i;                                  \
sll_real64, intent(in)                 :: array_name_and_dims;                      \
integer(hsize_t)                       :: array_dims(dataspace_dimension);          \
integer(hid_t)                         :: dataset_id;                               \
integer(hid_t)                         :: dataspace_id;                             \
rank = dataspace_dimension;                                                         \
do i = 1, rank;                                                                     \
   array_dims(i) = size(array,i);                                                   \
end do;                                                                             \
call H5Screate_simple_f(rank,array_dims,dataspace_id,error);                        \
call H5Dcreate_f(file_id,dsetname,H5T_NATIVE_DOUBLE,dataspace_id,dataset_id,error); \
call H5Dwrite_f(dataset_id,H5T_NATIVE_DOUBLE,array,array_dims,error);               \
call H5Sclose_f(dataspace_id,error);                                                \
call H5Dclose_f(dataset_id,error);                                                  \
end subroutine func_name

NEW_HDF5_FUNCTION(sll_hdf5_write_array_1d, 1, array(:))
NEW_HDF5_FUNCTION(sll_hdf5_write_array_2d, 2, array(:,:))
NEW_HDF5_FUNCTION(sll_hdf5_write_array_3d, 3, array(:,:,:))

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

