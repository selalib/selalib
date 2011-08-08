module sll_diagnostics
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"
use sll_working_precision
use sll_memory
use sll_assertion

implicit none

contains  ! ****************************************************************


  subroutine write_fxv( f, x, v)
  use hdf5
  sll_real64, dimension(:,:) :: x, v, f
  sll_int32 :: nx, nv
  integer(hid_t)   :: file_id, dataset_id, dataspace_id
  integer(hsize_t) :: data_dims(2)
  character(len=2) :: coordNames(2)
  character(len=10) :: fieldNames(4)
  integer :: error

  nx = size(x,1)
  nv = size(x,2)

  SLL_ASSERT(size(v,1) == nx)
  SLL_ASSERT(size(v,2) == nv)
  SLL_ASSERT(size(f,1) == nx)
  SLL_ASSERT(size(f,2) == nv)

  !Initialize FORTRAN interface.
  call h5open_f (error)
  !Create a new file using default properties.
  call H5Fcreate_f("fxv.h5", H5F_ACC_TRUNC_F, file_id, error)

  !Write the data file.
  coordNames(1) = "/x"
  coordNames(2) = "/v"
  fieldNames(1) = "/fxv_cells"
  fieldNames(2) = "/fxv_nodes"
 
  !Write separate coordinate arrays for the x and v coordinates.

  data_dims(1) = nx
  data_dims(2) = nv

  call H5Screate_simple_f(2, data_dims, dataspace_id, error)
  call H5Dcreate_f(file_id, coordnames(1), H5T_NATIVE_REAL, &
                   dataspace_id, dataset_id, error)
  call H5Dwrite_f(dataset_id, H5T_NATIVE_REAL, sngl(x), data_dims, error)
  call H5Dclose_f(dataset_id,error);
  call H5Sclose_f(dataspace_id,error);

  call H5Screate_simple_f(2, data_dims, dataspace_id, error)
  call H5Dcreate_f(file_id, coordnames(2), H5T_NATIVE_REAL, &
                   dataspace_id, dataset_id, error)
  call H5Dwrite_f(dataset_id, H5T_NATIVE_REAL, sngl(v), data_dims, error)
  call H5Dclose_f(dataset_id,error);
  call H5Sclose_f(dataspace_id,error);
 
  !Write the scalar data.
  data_dims(1) = nx-1
  data_dims(2) = nv-1
  call H5Screate_simple_f(2, data_dims, dataspace_id, error);
  call H5Dcreate_f(file_id, fieldNames(1), H5T_NATIVE_REAL, &
                   dataspace_id, dataset_id, error);
  call H5Dwrite_f(dataset_id, H5T_NATIVE_REAL, sngl(f(1:nx-1,1:nv-1)), data_dims, error);
  call H5Dclose_f(dataset_id,error);
  call H5Sclose_f(dataspace_id,error);

  data_dims(1) = nx
  data_dims(2) = nv
  call H5Screate_simple_f(2, data_dims, dataspace_id, error);
  call H5Dcreate_f(file_id, fieldNames(2), H5T_NATIVE_REAL, &
                   dataspace_id, dataset_id, error);
  call H5Dwrite_f(dataset_id, H5T_NATIVE_REAL, sngl(f), data_dims, error);
  call H5Dclose_f(dataset_id,error);
  call H5Sclose_f(dataspace_id,error);
 
  !Terminate access to the file.
  call H5Fclose_f(file_id, error)

  !write the xmf file readable by VisIt
  open(10, file="fxv.xmf")
  write(10,"(a)")"<?xml version='1.0' ?>"
  write(10,"(a)")"<!DOCTYPE Xdmf SYSTEM 'Xdmf.dtd' []>"
  write(10,"(a)")"<Xdmf Version='2.0'>"
  write(10,"(a)")"<Domain>"
  write(10,"(a)")"<Grid Name='mesh' GridType='Uniform'>"
  write(10,"(a,2i4,a)")"<Topology TopologyType='2DSMesh' NumberOfElements='",nv,nx,"'/>"
  write(10,"(a)")"<Geometry GeometryType='X_Y'>"
  write(10,"(a,2i4,a)")"<DataItem Dimensions='",nv,nx,"' NumberType='Float' Precision='4' Format='HDF'>"
  write(10,"(a)")"fxv.h5:"//coordNames(1)
  write(10,"(a)")"</DataItem>"
  write(10,"(a,2i4,a)")"<DataItem Dimensions='",nv,nx,"' NumberType='Float' Precision='4' Format='HDF'>"
  write(10,"(a)")"fxv.h5:"//coordNames(2)
  write(10,"(a)")"</DataItem>"
  write(10,"(a)")"</Geometry>"
  write(10,"(a)")"<Attribute Name='CellsValues' AttributeType='Scalar' Center='Cell'>"
  write(10,"(a,2i4,a)")"<DataItem Dimensions='",nv-1,nx-1,"' NumberType='Float' Precision='4' Format='HDF'>"
  write(10,"(a)")"fxv.h5:"//fieldNames(1)
  write(10,"(a)")"</DataItem>"
  write(10,"(a)")"</Attribute>"
  write(10,"(a)")"<Attribute Name='NodesValues' AttributeType='Scalar' Center='Node'>"
  write(10,"(a,2i4,a)")"<DataItem Dimensions='",nv,nx,"' NumberType='Float' Precision='4' Format='HDF'>"
  write(10,"(a)")"fxv.h5:"//fieldNames(2)
  write(10,"(a)")"</DataItem>"
  write(10,"(a)")"</Attribute>"
  write(10,"(a)")"</Grid>"
  write(10,"(a)")"</Domain>"
  write(10,"(a)")"</Xdmf>"
  close(10)

  end subroutine

end module sll_diagnostics
