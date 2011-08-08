module sll_diagnostics
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"
use sll_working_precision
use sll_memory
use sll_assertion

implicit none

contains  ! ****************************************************************


  subroutine write_mesh( x, v, nx, nv, mname)
  use hdf5
  character(len=*), intent(in) :: mname ! mesh file prefix
  sll_real64, dimension(:,:), intent(in) :: v
  sll_real64, dimension(:,:), intent(in) :: x
  sll_int32, intent(in) :: nx
  sll_int32, intent(in) :: nv

  integer(hid_t)   :: file_id, dataset_id, dataspace_id
  integer(hsize_t) :: data_dims(2)
  character(len=2) :: coordNames(2)
  integer :: error

  SLL_ASSERT(size(x,1) == nx)
  SLL_ASSERT(size(v,1) == nx)
  SLL_ASSERT(size(x,2) == nv)
  SLL_ASSERT(size(v,2) == nv)

  !Initialize FORTRAN interface.
  call h5open_f (error)
  !Create a new file using default properties.
  call H5Fcreate_f(trim(mname)//".h5", H5F_ACC_TRUNC_F, file_id, error)

  !Write the data file.
  coordNames(1) = "/x"
  coordNames(2) = "/v"
 
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
 
  !Terminate access to the file.
  call H5Fclose_f(file_id, error)

  !write the xmf file readable by VisIt
  open(10, file=trim(mname)//".xmf")
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
  write(10,"(a)")"</Grid>"
  write(10,"(a)")"</Domain>"
  write(10,"(a)")"</Xdmf>"
  close(10)

  end subroutine write_mesh

  subroutine write_vec1d( f, nx, nv, fname, mname)
  use hdf5
  character(len=*), intent(in) :: fname !prefix for field file
  character(len=*), intent(in) :: mname !prefix for mesh file
  sll_real64, dimension(:,:), intent(in) :: f
  sll_int32, intent(in) :: nx
  sll_int32, intent(in) :: nv

  integer(hid_t)   :: file_id, dataset_id, dataspace_id
  integer(hsize_t) :: data_dims(2)
  character(len=2) :: coordNames(2)
  integer :: error
  logical :: flag

  inquire(file=trim(mname)//".h5", exist=flag) 

  if (.not. flag) then
     call errout(6,"W","sll_diagnostics:write_vec1d", "Mesh file does not exist'" )
  end if

  !Initialize FORTRAN interface.
  call h5open_f (error)
  !Create a new file using default properties.
  call H5Fcreate_f(trim(fname)//".h5", H5F_ACC_TRUNC_F, file_id, error)

  !Write the data file.
  coordNames(1) = "/x"
  coordNames(2) = "/v"
 
  !Write the scalar data.
  data_dims(1) = nx-1
  data_dims(2) = nv-1
  call H5Screate_simple_f(2, data_dims, dataspace_id, error);
  call H5Dcreate_f(file_id, trim(fname)//"_cells", H5T_NATIVE_REAL, &
                   dataspace_id, dataset_id, error);
  call H5Dwrite_f(dataset_id, H5T_NATIVE_REAL, sngl(f(1:nx-1,1:nv-1)), data_dims, error);
  call H5Dclose_f(dataset_id,error);
  call H5Sclose_f(dataspace_id,error);

  data_dims(1) = nx
  data_dims(2) = nv
  call H5Screate_simple_f(2, data_dims, dataspace_id, error);
  call H5Dcreate_f(file_id, trim(fname)//"_nodes", H5T_NATIVE_REAL, &
                   dataspace_id, dataset_id, error);
  call H5Dwrite_f(dataset_id, H5T_NATIVE_REAL, sngl(f), data_dims, error);
  call H5Dclose_f(dataset_id,error);
  call H5Sclose_f(dataspace_id,error);
 
  !Terminate access to the file.
  call H5Fclose_f(file_id, error)

  !write the xmf file readable by VisIt
  open(10, file=trim(fname)//".xmf")
  write(10,"(a)")"<?xml version='1.0' ?>"
  write(10,"(a)")"<!DOCTYPE Xdmf SYSTEM 'Xdmf.dtd' []>"
  write(10,"(a)")"<Xdmf Version='2.0'>"
  write(10,"(a)")"<Domain>"
  write(10,"(a)")"<Grid Name='mesh' GridType='Uniform'>"
  write(10,"(a,2i4,a)")"<Topology TopologyType='2DSMesh' NumberOfElements='",nv,nx,"'/>"
  write(10,"(a)")"<Geometry GeometryType='X_Y'>"
  write(10,"(a,2i4,a)")"<DataItem Dimensions='",nv,nx,"' NumberType='Float' Precision='4' Format='HDF'>"
  write(10,"(a)")trim(mname)//".h5:"//coordNames(1)
  write(10,"(a)")"</DataItem>"
  write(10,"(a,2i4,a)")"<DataItem Dimensions='",nv,nx,"' NumberType='Float' Precision='4' Format='HDF'>"
  write(10,"(a)")trim(mname)//".h5:"//coordNames(2)
  write(10,"(a)")"</DataItem>"
  write(10,"(a)")"</Geometry>"
  write(10,"(a)")"<Attribute Name='CellsValues' AttributeType='Scalar' Center='Cell'>"
  write(10,"(a,2i4,a)")"<DataItem Dimensions='",nv-1,nx-1,"' NumberType='Float' Precision='4' Format='HDF'>"
  write(10,"(a)")trim(fname)//".h5:"//trim(fname)//"_cells"
  write(10,"(a)")"</DataItem>"
  write(10,"(a)")"</Attribute>"
  write(10,"(a)")"<Attribute Name='NodesValues' AttributeType='Scalar' Center='Node'>"
  write(10,"(a,2i4,a)")"<DataItem Dimensions='",nv,nx,"' NumberType='Float' Precision='4' Format='HDF'>"
  write(10,"(a)")trim(fname)//".h5:"//trim(fname)//"_nodes"
  write(10,"(a)")"</DataItem>"
  write(10,"(a)")"</Attribute>"
  write(10,"(a)")"</Grid>"
  write(10,"(a)")"</Domain>"
  write(10,"(a)")"</Xdmf>"
  close(10)

  end subroutine write_vec1d

  subroutine write_vec2d( fx, fv, nx, nv, fname, mname)
  use hdf5
  character(len=*), intent(in) :: fname !prefix for field file
  character(len=*), intent(in) :: mname !prefix for mesh file
  sll_real64, dimension(:,:), intent(in) :: fx
  sll_real64, dimension(:,:), intent(in) :: fv
  sll_int32, intent(in) :: nx
  sll_int32, intent(in) :: nv

  integer(hid_t)   :: file_id, dataset_id, dataspace_id
  integer(hsize_t) :: data_dims(2)
  character(len=2) :: coordNames(2)
  integer :: error
  logical :: flag

  inquire(file=trim(mname)//".h5", exist=flag) 

  if (.not. flag) then
     call errout(6,"W","sll_diagnostics:write_vec2d", "Mesh file does not exist'" )
  end if

  !Initialize FORTRAN interface.
  call h5open_f (error)
  !Create a new file using default properties.
  call H5Fcreate_f(trim(fname)//".h5", H5F_ACC_TRUNC_F, file_id, error)

  !Write the data file.
  coordNames(1) = "/x"
  coordNames(2) = "/v"
 
  !Write the scalar data.
  data_dims(1) = nx-1
  data_dims(2) = nv-1
  call H5Screate_simple_f(2, data_dims, dataspace_id, error);
  call H5Dcreate_f(file_id, trim(fname)//"x_cells", H5T_NATIVE_REAL, &
                   dataspace_id, dataset_id, error);
  call H5Dwrite_f(dataset_id, H5T_NATIVE_REAL, sngl(fx(1:nx-1,1:nv-1)), data_dims, error);
  call H5Dclose_f(dataset_id,error);
  call H5Sclose_f(dataspace_id,error);

  data_dims(1) = nx
  data_dims(2) = nv
  call H5Screate_simple_f(2, data_dims, dataspace_id, error);
  call H5Dcreate_f(file_id, trim(fname)//"x_nodes", H5T_NATIVE_REAL, &
                   dataspace_id, dataset_id, error);
  call H5Dwrite_f(dataset_id, H5T_NATIVE_REAL, sngl(fx), data_dims, error);
  call H5Dclose_f(dataset_id,error);
  call H5Sclose_f(dataspace_id,error);

  data_dims(1) = nx-1
  data_dims(2) = nv-1
  call H5Screate_simple_f(2, data_dims, dataspace_id, error);
  call H5Dcreate_f(file_id, trim(fname)//"v_cells", H5T_NATIVE_REAL, &
                   dataspace_id, dataset_id, error);
  call H5Dwrite_f(dataset_id, H5T_NATIVE_REAL, sngl(fv(1:nx-1,1:nv-1)), data_dims, error);
  call H5Dclose_f(dataset_id,error);
  call H5Sclose_f(dataspace_id,error);

  data_dims(1) = nx
  data_dims(2) = nv
  call H5Screate_simple_f(2, data_dims, dataspace_id, error);
  call H5Dcreate_f(file_id, trim(fname)//"v_nodes", H5T_NATIVE_REAL, &
                   dataspace_id, dataset_id, error);
  call H5Dwrite_f(dataset_id, H5T_NATIVE_REAL, sngl(fv), data_dims, error);
  call H5Dclose_f(dataset_id,error);
  call H5Sclose_f(dataspace_id,error);
 
  !Terminate access to the file.
 
  !Terminate access to the file.
  call H5Fclose_f(file_id, error)

  !write the xmf file readable by VisIt
  open(10, file=trim(fname)//".xmf")
  write(10,"(a)")"<?xml version='1.0' ?>"
  write(10,"(a)")"<!DOCTYPE Xdmf SYSTEM 'Xdmf.dtd' []>"
  write(10,"(a)")"<Xdmf Version='2.0'>"
  write(10,"(a)")"<Domain>"
  write(10,"(a)")"<Grid Name='mesh' GridType='Uniform'>"
  write(10,"(a,2i4,a)")"<Topology TopologyType='2DSMesh' NumberOfElements='",nv,nx,"'/>"
  write(10,"(a)")"<Geometry GeometryType='X_Y'>"
  write(10,"(a,2i4,a)")"<DataItem Dimensions='",nv,nx,"' NumberType='Float' Precision='4' Format='HDF'>"
  write(10,"(a)")trim(mname)//".h5:"//coordNames(1)
  write(10,"(a)")"</DataItem>"
  write(10,"(a,2i4,a)")"<DataItem Dimensions='",nv,nx,"' NumberType='Float' Precision='4' Format='HDF'>"
  write(10,"(a)")trim(mname)//".h5:"//coordNames(2)
  write(10,"(a)")"</DataItem>"
  write(10,"(a)")"</Geometry>"
  write(10,"(a)")"<Attribute Name='x_CellsValues' AttributeType='Scalar' Center='Cell'>"
  write(10,"(a,2i4,a)")"<DataItem Dimensions='",nv-1,nx-1,"' NumberType='Float' Precision='4' Format='HDF'>"
  write(10,"(a)")trim(fname)//".h5:"//trim(fname)//"x_cells"
  write(10,"(a)")"</DataItem>"
  write(10,"(a)")"</Attribute>"
  write(10,"(a)")"<Attribute Name='x_NodesValues' AttributeType='Scalar' Center='Node'>"
  write(10,"(a,2i4,a)")"<DataItem Dimensions='",nv,nx,"' NumberType='Float' Precision='4' Format='HDF'>"
  write(10,"(a)")trim(fname)//".h5:"//trim(fname)//"x_nodes"
  write(10,"(a)")"</DataItem>"
  write(10,"(a)")"</Attribute>"
  write(10,"(a)")"<Attribute Name='v_CellsValues' AttributeType='Scalar' Center='Cell'>"
  write(10,"(a,2i4,a)")"<DataItem Dimensions='",nv-1,nx-1,"' NumberType='Float' Precision='4' Format='HDF'>"
  write(10,"(a)")trim(fname)//".h5:"//trim(fname)//"v_cells"
  write(10,"(a)")"</DataItem>"
  write(10,"(a)")"</Attribute>"
  write(10,"(a)")"<Attribute Name='v_NodesValues' AttributeType='Scalar' Center='Node'>"
  write(10,"(a,2i4,a)")"<DataItem Dimensions='",nv,nx,"' NumberType='Float' Precision='4' Format='HDF'>"
  write(10,"(a)")trim(fname)//".h5:"//trim(fname)//"v_nodes"
  write(10,"(a)")"</DataItem>"
  write(10,"(a)")"</Attribute>"
  write(10,"(a)")"</Grid>"
  write(10,"(a)")"</Domain>"
  write(10,"(a)")"</Xdmf>"
  close(10)

  end subroutine write_vec2d


  !Title: Subroutine Errout
  
  !Subroutine: errout
  !     Outputs an error message
  !     PRTFIL - unit number for print-out
  !     SEVRTY - 'W' - Warning 'F' - Fatal
  !     WHERE  - in which program or subroutine
  !     ErrMsg - error message
  subroutine errout( prtfil, sevrty, lwhere, ErrMsg )

    implicit none
    integer, intent(in) ::  prtfil
    character(len=1) :: sevrty 
    character(len=*) :: lwhere , ErrMsg

    write( prtfil, * )

    select case ( sevrty )  !     *** Severity ***
    case ( 'W' )
    write(prtfil,"(/10x,a)") '*** WARNING ***'
    case ( 'F' )
    write(prtfil,"(/10x,a)") '*** FATAL ERROR ***'
    case default
    write(prtfil,"(/10x,a)") '*** FATAL ERROR ***'
    write(prtfil,"(/10x,a)") 'Error handler (ERROUT) called with unknown severity level: ', SEVRTY
    end select

    write( prtfil,"(/10x,a)") 'Generated by program or subroutine: ', trim(lwhere)

    write( prtfil,"(/10x,a)") trim(ErrMsg)
    write( prtfil,"(/10x,a)")

    ! return or stop depending on severity

    if ( sevrty == 'W' ) then
       return
    else
       stop 'Fatal Error: See print file for details'
    end if

end subroutine errout


end module sll_diagnostics
