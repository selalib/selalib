program diagnostics_tester
#include "sll_working_precision.h"
  use sll_diagnostics
  implicit none

  sll_int32 :: nx, nv, i, j
  sll_real64, allocatable, dimension(:,:) :: x
  sll_real64, allocatable, dimension(:,:) :: v
  sll_real64, allocatable, dimension(:,:) :: f
  sll_real64 :: angle, xt, vt, pi, R
  
  pi = 4.*atan(1.)

  nx = 128
  nv = 64

  ! Create the coordinate data.
  allocate(x(nx,nv))
  allocate(v(nx,nv))

  do j = 1, nv
     vt = real(j-1)/(nv-1)
     angle = vt * 2. * pi
     do i = 1, nx
        xt = real(i-1) / float(nx-1)
        R = (1.-xt)*2. + xt*5.
        x(i,j) = R * cos(angle)
        v(i,j) = R * sin(angle)
     end do
  end do 

  ! Create the scalar data.
  allocate(f(nx-1,nv-1))
  do j = 1, nv-1
     do i = 1, nx-1
        f(i,j) = i*sin(float(j-1))
        f(i,j) = i*sin(float(j-1))
     end do
  end do
 
  call write_vec1d( f, nx, nv, "vec1d", "mesh")
  call write_vec2d( x, v, nx, nv, "vec2d", "mesh")
  call write_mesh(x, v, nx, nv, "mesh")

  call write_hdf5_data()

contains

  subroutine write_hdf5_data()

     use hdf5 
     implicit none
     integer :: k, error, i, j
     real :: xt, yt, pi, R, angle
     integer(hid_t)   ::  file_id, dataset_id, dataspace_id
     integer(hsize_t) ::  data_dims(2)
     character(len=2) :: coordNames(2)
     integer, parameter :: nx = 30, ny = 20
     real ,dimension(:), allocatable :: x, y
     real ,dimension(:), allocatable :: cells_values, nodes_values

     pi = 4. * atan(1.)
     !call write_xdmf_xml(nx, ny)

     ! Create the coordinate data.
     allocate(x((nx+1)*(ny+1)))
     allocate(y((nx+1)*(ny+1)))
     k = 1
     do j = 1, ny+1
        yt = float(j-1)/ny
        angle = yt * 2. * pi
        do i = 1, nx+1
           xt = float(i-1) / float(nx)
           R = (1.-xt)*2. + xt*5.
           x(k) = R * cos(angle)
           y(k) = R * sin(angle)
           k=k+1
        end do
     end do 

     ! Create the scalar data.
     allocate(cells_values(nx*ny))
     k = 0
     do j = 1, ny-1
        do i = 1, nx-1
           k = k+1
           cells_values(k) = float(j-1)
        end do
     end do
 
     allocate(nodes_values((nx+1)*(ny+1)))
     k = 0
     do j = 1, ny+1
        do i = 1, nx+1
           k = k+1
           nodes_values(k) = float(i-1)
        end do
     end do
 
     !Initialize FORTRAN interface.
     call h5open_f (error)
     !Create a new file using default properties.
     call H5Fcreate_f("xdmf2d.h5", H5F_ACC_TRUNC_F, file_id, error)

     !Write the data file.
     coordNames(1) = "/X"
     coordNames(2) = "/Y"
 
     !Write separate coordinate arrays for the x and y coordinates.
     data_dims(1) = nx+1
     data_dims(2) = ny+1
     call H5Screate_simple_f(2, data_dims, dataspace_id, error)
     call H5Dcreate_f(file_id, coordnames(1), H5T_NATIVE_REAL, &
                      dataspace_id, dataset_id, error)
     call H5Dwrite_f(dataset_id, H5T_NATIVE_REAL, x, data_dims, error)
     call H5Dclose_f(dataset_id,error);
     call H5Sclose_f(dataspace_id,error);

     data_dims(1) = nx+1
     data_dims(2) = ny+1
     call H5Screate_simple_f(2, data_dims, dataspace_id, error)
     call H5Dcreate_f(file_id, coordnames(2), H5T_NATIVE_REAL, &
                      dataspace_id, dataset_id, error)
     call H5Dwrite_f(dataset_id, H5T_NATIVE_REAL, y, data_dims, error)
     call H5Dclose_f(dataset_id,error);
     call H5Sclose_f(dataspace_id,error);
 
     !Write the scalar data.
     data_dims(1) = nx;
     data_dims(2) = ny;
     call H5Screate_simple_f(2, data_dims, dataspace_id, error);
     call H5Dcreate_f(file_id, "/Pressure", H5T_NATIVE_REAL, &
                      dataspace_id, dataset_id, error);
     call H5Dwrite_f(dataset_id, H5T_NATIVE_REAL, cells_values, data_dims, error);
     call H5Dclose_f(dataset_id,error);
     call H5Sclose_f(dataspace_id,error);
 
     data_dims(1) = nx+1
     data_dims(2) = ny+1
     call H5Screate_simple_f(2, data_dims, dataspace_id, error);
     call H5Dcreate_f(file_id, "/VelocityX", H5T_NATIVE_REAL, &
                      dataspace_id, dataset_id, error);
     call H5Dwrite_f(dataset_id, H5T_NATIVE_REAL, nodes_values, data_dims, error);
     call H5Dclose_f(dataset_id,error);
     call H5Sclose_f(dataspace_id,error);
 
     !Free the data.
     deallocate(x);
     deallocate(y);
     deallocate(cells_values);
     deallocate(nodes_values);
 
     !Terminate access to the file.
     call H5Fclose_f(file_id, error)

     !Close FORTRAN interface.
     call H5close_f(error) 

  end subroutine write_hdf5_data 

end program diagnostics_tester
