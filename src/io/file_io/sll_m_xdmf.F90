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
!> Implements the functions to write xdmf file plotable by VisIt
!> @details
!> In <b> XDMF </b> (eXtensible Data Model and Format) the description of the 
!> data is separate from the values themselves. Light data is stored using XML, 
!> Heavy data is stored using HDF5 or Binary files. \n
!>
module sll_m_xdmf
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_assert.h"
#include "sll_errors.h"
#include "sll_memory.h"
#include "sll_working_precision.h"

  use sll_m_ascii_io, only: &
    sll_o_ascii_write_array

  use sll_m_utilities, only: &
    sll_s_int2string

  use sll_m_xml_io, only: &
    sll_o_xml_field, &
    sll_s_xml_file_close, &
    sll_s_xml_file_create, &
    sll_o_xml_grid_geometry

#ifdef NOHDF5
  use sll_m_binary_io, only: &
    sll_s_binary_file_create, &
    sll_s_binary_file_close, &
    sll_o_binary_write_array

#else
  use hdf5, only: hid_t
  use sll_m_hdf5_io_serial, only: &
    sll_s_hdf5_ser_file_create, &
    sll_s_hdf5_ser_file_close, &
    sll_o_hdf5_ser_write_array

#endif
  implicit none

  public :: &
    sll_s_plot_f, &
    sll_s_plot_f_cartesian, &
    sll_s_xdmf_close, &
    sll_s_xdmf_corect2d_nodes, &
    sll_s_xdmf_corect3d_nodes, &
    sll_s_xdmf_curv2d_nodes, &
    sll_s_xdmf_curv3d_nodes, &
    sll_o_xdmf_open, &
    sll_s_xdmf_rect2d_nodes, &
    sll_s_xdmf_rect3d_nodes, &
    sll_o_xdmf_write_array

  private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
!> Create a xmf file 
interface sll_o_xdmf_open
   module procedure sll_xdmf_open_2d
   module procedure sll_xdmf_open_3d
end interface

!> Write the field in xdmf format
interface sll_o_xdmf_write_array
   module procedure sll_xdmf_array_2d
   module procedure sll_xdmf_array_3d
end interface
  
contains  
  
!> Add the the good value of time in VisIt plot.
subroutine sll_xdmf_set_time(file_id, time)
sll_real64, intent(in) :: time       !< input time
sll_int32, intent(in)  :: file_id    !< file unit number
logical                :: i_opened

inquire (file_id, opened=i_opened)

if (i_opened) then
   write(file_id,"(a13,g15.3,a3)") "<Time Value='",time,"'/>"
else
   SLL_ERROR( "sll_xdmf_set_time", "This xdmf file is not not open." )
end if

end subroutine sll_xdmf_set_time


!>Open a XDMF format file for a 2d plot
subroutine sll_xdmf_open_2d(file_name,mesh_name,nnodes_x1,nnodes_x2,file_id,error)

character(len=*), intent(in) :: file_name  !< xmf file name
character(len=*), intent(in) :: mesh_name  !< mesh file name
sll_int32, intent(in)        :: nnodes_x1  !< x node number
sll_int32, intent(in)        :: nnodes_x2  !< y node number
sll_int32, intent(out)       :: file_id    !< file unit number
sll_int32, intent(out)       :: error      !< error code
    
call sll_s_xml_file_create(trim(file_name),file_id,error)
call sll_o_xml_grid_geometry(file_id, trim(mesh_name), &
                           nnodes_x1, nnodes_x2)

end subroutine sll_xdmf_open_2d

!>Open a XDMF format file for a 3d plot
subroutine sll_xdmf_open_3d( file_name, &
                             mesh_name, &
                             nnodes_x1, &
                             nnodes_x2, &
                             nnodes_x3, &
                             file_id,   &
                             error)
    
character(len=*), intent(in) :: file_name  !< xmf file name
character(len=*), intent(in) :: mesh_name  !< mesh name
sll_int32, intent(out)       :: file_id    !< file unit number
sll_int32, intent(out)       :: error      !< error code
sll_int32, intent(in)        :: nnodes_x1  !< x nodes number
sll_int32, intent(in)        :: nnodes_x2  !< y nodes number
sll_int32, intent(in)        :: nnodes_x3  !< z nodes number

call sll_s_xml_file_create(trim(file_name),file_id,error)
call sll_o_xml_grid_geometry(file_id,trim(mesh_name), &
                           nnodes_x1,nnodes_x2,nnodes_x3)

end subroutine sll_xdmf_open_3d
  
!>Write 2d array in binary or hdf5 file and the matching line in XDMF file
subroutine sll_xdmf_array_2d(mesh_name,array,array_name,error,xmffile_id,center)

character(len=*), intent(in)    :: mesh_name   !< mesh name
sll_real64, intent(in)          :: array(:,:)  !< data array
character(len=*), intent(in)    :: array_name  !< array name (hdf5 dataset)
sll_int32, intent(out)          :: error       !< error code
sll_int32                       :: npoints_x1  !< x nodes number
sll_int32                       :: npoints_x2  !< y nodes number
sll_int32, intent(in), optional :: xmffile_id  !< xmf file unit number
character(len=4), optional      :: center      !< "Node" or "Cell"

#ifndef NOHDF5
integer(hid_t) :: hfile_id
#else
sll_int32                       :: file_id     !< hdf5 file unit number
#endif


npoints_x1 = size(array,1)
npoints_x2 = size(array,2)

#ifdef NOHDF5
call sll_s_binary_file_create(trim(mesh_name)//"-"//trim(array_name)//".bin", &
                            file_id,error)
call sll_o_binary_write_array(file_id,array,error)
call sll_s_binary_file_close(file_id,error)
#else
call sll_s_hdf5_ser_file_create( trim(mesh_name)//"-"//trim(array_name)//".h5", &
                          hfile_id, error )
call sll_o_hdf5_ser_write_array( hfile_id, array, "/"//trim(array_name), error )
call sll_s_hdf5_ser_file_close( hfile_id, error )
#endif

if ( present(xmffile_id) .and. present(center)) then
#ifdef NOHDF5
   call sll_o_xml_field(xmffile_id,trim(array_name), &
                      trim(mesh_name)//"-"//trim(array_name)//".bin", &
                      npoints_x1,npoints_x2,'Binary',center)
#else
   call sll_o_xml_field( &
        xmffile_id, &
        trim(array_name), &
        trim(mesh_name)//"-"//trim(array_name)//".h5:/"//trim(array_name), &
        npoints_x1,&
        npoints_x2,&
        'HDF',&
        center)
#endif
end if
end subroutine sll_xdmf_array_2d

!>Write 3d array in binary or hdf5 file and the matching line in XDMF file
subroutine sll_xdmf_array_3d(mesh_name,array,array_name,error,xmffile_id,center)

character(len=*), intent(in)    :: mesh_name    !< mesh name
sll_real64, intent(in)          :: array(:,:,:) !< data array
character(len=*), intent(in)    :: array_name   !< hdf5 dataset name
sll_int32, intent(out)          :: error        !< error code
sll_int32, intent(in), optional :: xmffile_id   !< xmf file unit number
character(len=4), optional      :: center       !< "Node" or "Cell"
sll_int32                       :: npoints_x1   !< x nodes number
sll_int32                       :: npoints_x2   !< y nodes number
sll_int32                       :: npoints_x3   !< z nodes number
    
#ifndef NOHDF5
integer(hid_t) :: hfile_id
#else
sll_int32                       :: file_id      !< hdf5 file unit number
#endif

npoints_x1 = size(array,1)
npoints_x2 = size(array,2)
npoints_x3 = size(array,3)

#ifdef NOHDF5
call sll_s_binary_file_create(trim(mesh_name)//"-"//trim(array_name)//".bin", &
                            file_id,error)
call sll_o_binary_write_array(file_id,array,error)
call sll_s_binary_file_close(file_id,error)
#else
call sll_s_hdf5_ser_file_create( trim(mesh_name)//"-"//trim(array_name)//".h5", &
                          hfile_id, error )
call sll_o_hdf5_ser_write_array( hfile_id, array, "/"//trim(array_name), error )
call sll_s_hdf5_ser_file_close( hfile_id, error )
#endif

if ( present(xmffile_id) .and. present(center)) then

#ifdef NOHDF5
   call sll_o_xml_field(xmffile_id,trim(array_name), &
                      trim(mesh_name)//"-"//trim(array_name)//".bin", &
                      npoints_x1,npoints_x2,npoints_x3,'Binary',center)
#else
   call sll_o_xml_field( &
         xmffile_id, &
         trim(array_name), &
         trim(mesh_name)//"-"//trim(array_name)//".h5:/"//trim(array_name), &
         npoints_x1,&
         npoints_x2,&
         npoints_x3,&
         'HDF', &
         center)
#endif

end if

end subroutine sll_xdmf_array_3d


!>Subroutine to write a 2D array in xdmf format
!>The field is describe on a cartesian mesh
!>Axis are perpendicular and spacing is constant
subroutine sll_s_xdmf_corect2d_nodes( file_name,  &
                                      array,      &
                                      array_name, &
                                      eta1_min,   &
                                      delta_eta1, &
                                      eta2_min,   &
                                      delta_eta2, &
                                      file_format,&
                                      iplot,      & 
                                      time        ) 

sll_real64, intent(in)       :: array(:,:)  !< data array
character(len=*), intent(in) :: file_name   !< xmf file name
character(len=*), intent(in) :: array_name  !< field name
sll_int32                    :: error       !< error code
sll_real64                   :: eta1_min    !< x min
sll_real64                   :: eta2_min    !< y min
sll_real64                   :: delta_eta1  !< dx
sll_real64                   :: delta_eta2  !< dy
sll_int32                    :: file_id     !< xmf file unit number
sll_int32                    :: nx1         !< x nodes number
sll_int32                    :: nx2         !< y nodes number
character(len=4), optional   :: file_format !< "HDF5" or "Binary"
sll_int32       , optional   :: iplot       !< plot index
sll_real64      , optional   :: time        !< time value

character(len=4)             :: cplot     
#ifndef NOHDF5
integer(hid_t)               :: hfile_id    !< h5 file unit number
#endif
    
nx1 = size(array,1)
nx2 = size(array,2)

if( present(iplot)) then
  call sll_s_int2string(iplot, cplot)
  call sll_s_xml_file_create(file_name//cplot//".xmf",file_id,error)
else
  call sll_s_xml_file_create(file_name//".xmf",file_id,error)
end if

write(file_id,"(a)")"<Grid Name='mesh' GridType='Uniform'>"
if (present(time)) then
  write(file_id,"(a,f12.6,a)") "<Time Value='",time,"'/>"
end if
write(file_id,"(a,2i5,a)")"<Topology TopologyType='2DCoRectMesh' NumberOfElements='", &
                          nx2,nx1,"'/>"
write(file_id,"(a)")"<Geometry GeometryType='ORIGIN_DXDY'>"
write(file_id,"(a)")"<DataItem Dimensions='2' NumberType='Float' Format='XML'>"
write(file_id,"(2f12.5)") eta1_min, eta2_min
write(file_id,"(a)")"</DataItem>"
write(file_id,"(a)")"<DataItem Dimensions='2' NumberType='Float' Format='XML'>"
write(file_id,"(2f12.5)") delta_eta1, delta_eta2
write(file_id,"(a)")"</DataItem>"
write(file_id,"(a)")"</Geometry>"
write(file_id,"(a)")"<Attribute Name='"//array_name//"' AttributeType='Scalar' Center='Node'>"

if(present(file_format)) then
  if( file_format == "HDF5") then
    write(file_id,"(a,2i5,a)")"<DataItem Dimensions='",nx2,nx1, &
                              "' NumberType='Float' Precision='8' Format='HDF'>"
    if (present(iplot)) then
      write(file_id,"(a)")array_name//cplot//".h5:/node_values"
    else
      write(file_id,"(a)")array_name//".h5:/node_values"
    end if
#ifndef NOHDF5
    if (present(iplot)) then
      call sll_s_hdf5_ser_file_create( array_name//cplot//".h5", hfile_id, error )
    else
      call sll_s_hdf5_ser_file_create( array_name//".h5", hfile_id, error )
    end if
    call sll_o_hdf5_ser_write_array( hfile_id, array, "/node_values", error )
    call sll_s_hdf5_ser_file_close( hfile_id, error )
#endif
  end if
else
  write(file_id,"(a,2i5,a)")"<DataItem Dimensions='",nx2,nx1, &
                             "' NumberType='Float' Precision='8' Format='XML'>"
  call sll_o_ascii_write_array(file_id,array,error)
end if
write(file_id,"(a)")"</DataItem>"
write(file_id,"(a)")"</Attribute>"
call sll_s_xml_file_close(file_id,error)

end subroutine sll_s_xdmf_corect2d_nodes

!>Subroutine to write a 3D array in xdmf format
!>The field is describe on a cartesian mesh
!>Axis are perpendicular and spacing is constant
subroutine sll_s_xdmf_corect3d_nodes( file_name,  &
                                      array,      &
                                      array_name, &
                                      eta1_min,   &
                                      delta_eta1, &
                                      eta2_min,   &
                                      delta_eta2, &
                                      eta3_min,   &
                                      delta_eta3, &
                                      file_format,&
                                      iplot       ) 

sll_real64, intent(in)       :: array(:,:,:)!< data array
character(len=*), intent(in) :: file_name   !< xmf file name
character(len=*), intent(in) :: array_name  !< field name
sll_int32                    :: error       !< error code
sll_real64                   :: eta1_min    !< x min
sll_real64                   :: eta2_min    !< y min
sll_real64                   :: eta3_min    !< z min
sll_real64                   :: delta_eta1  !< dx
sll_real64                   :: delta_eta2  !< dy
sll_real64                   :: delta_eta3  !< dz
sll_int32                    :: file_id     !< xmf file unit number
sll_int32                    :: nx1         !< x nodes number
sll_int32                    :: nx2         !< y nodes number
sll_int32                    :: nx3         !< z nodes number
character(len=4), optional   :: file_format !< "HDF5" or "Binary"
sll_int32       , optional   :: iplot       !< plot index

character(len=4)             :: cplot
#ifndef NOHDF5
integer(hid_t)               :: hfile_id    !< h5 file unit number
#endif
    
nx1 = size(array,1)
nx2 = size(array,2)
nx3 = size(array,3)

if (present(iplot)) then
  call sll_s_int2string(iplot, cplot)
  call sll_s_xml_file_create(file_name//cplot//".xmf",file_id,error)
else
  call sll_s_xml_file_create(file_name//".xmf",file_id,error)
end if

write(file_id,"(a)")"<Grid Name='mesh' GridType='Uniform'>"
write(file_id,"(a,3i5,a)")"<Topology TopologyType='3DCoRectMesh' NumberOfElements='", &
                          nx3,nx2,nx1,"'/>"
write(file_id,"(a)")"<Geometry GeometryType='ORIGIN_DXDYDZ'>"
write(file_id,"(a)")"<DataItem Dimensions='3' NumberType='Float' Format='XML'>"
write(file_id,"(3f12.5)") eta1_min, eta2_min, eta3_min
write(file_id,"(a)")"</DataItem>"
write(file_id,"(a)")"<DataItem Dimensions='3' NumberType='Float' Format='XML'>"
write(file_id,"(3f12.5)") delta_eta1, delta_eta2, delta_eta3
write(file_id,"(a)")"</DataItem>"
write(file_id,"(a)")"</Geometry>"
write(file_id,"(a)")"<Attribute Name='"//array_name//"' AttributeType='Scalar' Center='Node'>"
if(present(file_format) .and. file_format == "HDF5") then
  write(file_id,"(a,3i5,a)")"<DataItem Dimensions='",nx3,nx2,nx1, &
                            "' NumberType='Float' Precision='8' Format='HDF'>"
  if (present(iplot)) then
    write(file_id,"(a)")array_name//cplot//".h5:/node_values"
  else
    write(file_id,"(a)")array_name//".h5:/node_values"
  end if
#ifndef NOHDF5
  if (present(iplot)) then
    call sll_s_hdf5_ser_file_create( array_name//cplot//".h5", hfile_id, error )
  else
    call sll_s_hdf5_ser_file_create( array_name//".h5", hfile_id, error )
  end if
  call sll_o_hdf5_ser_write_array( hfile_id, array, "/node_values", error )
  call sll_s_hdf5_ser_file_close( hfile_id, error )
#endif
else
  write(file_id,"(a,3i5,a)")"<DataItem Dimensions='",nx3,nx2,nx1, &
                            "' NumberType='Float' Precision='8' Format='XML'>"
  call sll_o_ascii_write_array(file_id,array,error)
end if
write(file_id,"(a)")"</DataItem>"
write(file_id,"(a)")"</Attribute>"
call sll_s_xml_file_close(file_id,error)

end subroutine sll_s_xdmf_corect3d_nodes

!>Subroutine to write a 2D array in xdmf format.
!>The field is describe on a cartesian mesh.
!>Axis are perpendicular and spacing is define by eta1 and eta2 arrays.
subroutine sll_s_xdmf_rect2d_nodes( file_name,   &
                                    array,       &
                                    array_name,  &
                                    eta1,        &
                                    eta2,        &
                                    file_format, &
                                    iplot,       & 
                                    time         ) 

sll_real64, intent(in)       :: array(:,:) !< data array
sll_real64, intent(in)       :: eta1(:)    !< x data
sll_real64, intent(in)       :: eta2(:)    !< y data
character(len=*), intent(in) :: file_name  !< xmf file name
character(len=*), intent(in) :: array_name !< array name
sll_int32                    :: error      !< error code
sll_int32                    :: file_id    !< xmf file unit number
sll_int32                    :: nx1        !< x nodes number
sll_int32                    :: nx2        !< y nodes number
character(len=4), optional   :: file_format!< file format "HDF5" or "Binary"
sll_int32,  optional         :: iplot      !< plot index
sll_real64, optional         :: time       !< time
sll_int32                    :: i, j
character(len=4)             :: cplot
#ifndef NOHDF5
integer(hid_t)               :: hfile_id    !< h5 file unit number
#endif
    
nx1 = size(array,1)
nx2 = size(array,2)

SLL_ASSERT(nx1 == size(eta1))
SLL_ASSERT(nx2 == size(eta2))

if (present(iplot)) then
  call sll_s_int2string(iplot, cplot)
  call sll_s_xml_file_create(file_name//cplot//".xmf",file_id,error)
else
  call sll_s_xml_file_create(file_name//".xmf",file_id,error)
end if
write(file_id,"(a)")"<Grid Name='mesh' GridType='Uniform'>"
if (present(time)) then
  write(file_id,"(a,f12.6,a)") "<Time Value='",time,"'/>"
end if
write(file_id,"(a,2i5,a)")"<Topology TopologyType='2DRectMesh' NumberOfElements='", &
                          nx2,nx1,"'/>"
write(file_id,"(a)")"<Geometry GeometryType='VXVY'>"
write(file_id,"(a,i5,a)")"<DataItem Dimensions='",nx1, &
                         "' NumberType='Float' Format='XML'>"
write(file_id,*) (eta1(i),i=1,nx1)
write(file_id,"(a)")"</DataItem>"
write(file_id,"(a,i5,a)")"<DataItem Dimensions='",nx2, &
                             "' NumberType='Float' Format='XML'>"
write(file_id,*) (eta2(j),j=1,nx2)
write(file_id,"(a)")"</DataItem>"
write(file_id,"(a)")"</Geometry>"
write(file_id,"(a)")"<Attribute Name='"//array_name//"' AttributeType='Scalar' Center='Node'>"
if(present(file_format) .and. file_format == "HDF5") then
   write(file_id,"(a,2i5,a)")"<DataItem Dimensions='",nx2,nx1, &
                             "' NumberType='Float' Precision='8' Format='HDF'>"
  if (present(iplot)) then
     write(file_id,"(a)")array_name//cplot//".h5:/node_values"
  else
     write(file_id,"(a)")array_name//".h5:/node_values"
  end if
#ifndef NOHDF5
  if (present(iplot)) then
    call sll_s_hdf5_ser_file_create( array_name//cplot//".h5", hfile_id, error )
  else
    call sll_s_hdf5_ser_file_create( array_name//".h5", hfile_id, error )
  end if
  call sll_o_hdf5_ser_write_array( hfile_id, array, "/node_values", error )
  call sll_s_hdf5_ser_file_close( hfile_id, error )
#endif
else
  write(file_id,"(a,2i5,a)")"<DataItem Dimensions='",nx2,nx1, &
                            "' NumberType='Float' Precision='4' Format='XML'>"
  call sll_o_ascii_write_array(file_id,array,error)
end if
write(file_id,"(a)")"</DataItem>"
write(file_id,"(a)")"</Attribute>"
call sll_s_xml_file_close(file_id,error)

end subroutine sll_s_xdmf_rect2d_nodes

!>Subroutine to write a 3D array in xdmf format.
!>The field is describe on a cartesian mesh.
!>Axis are perpendicular and spacing is define by eta1, eta2 and eta3 arrays.
subroutine sll_s_xdmf_rect3d_nodes( file_name,   &
                                    array,       &
                                    array_name,  &
                                    eta1,        &
                                    eta2,        &
                                    eta3,        &
                                    file_format, &
                                    iplot        ) 

sll_real64, intent(in)       :: array(:,:,:) !< data array
sll_real64, intent(in)       :: eta1(:)      !< x data
sll_real64, intent(in)       :: eta2(:)      !< y data
sll_real64, intent(in)       :: eta3(:)      !< z data
character(len=*), intent(in) :: file_name    !< xmf file name
character(len=*), intent(in) :: array_name   !< array name
sll_int32                    :: error        !< error code
sll_int32                    :: file_id      !< xmf file unit number
sll_int32                    :: nx1          !< x nodes number
sll_int32                    :: nx2          !< y nodes number
sll_int32                    :: nx3          !< z nodes number
character(len=4), optional   :: file_format  !< file format "HDF5" or "Binary"
sll_int32       , optional   :: iplot        !< plot index
sll_int32                    :: i, j, k
#ifndef NOHDF5
integer(hid_t)               :: hfile_id    !< h5 file unit number
#endif
character(len=4)             :: cplot
    
nx1 = size(array,1)
nx2 = size(array,2)
nx3 = size(array,3)

SLL_ASSERT(nx1 == size(eta1))
SLL_ASSERT(nx2 == size(eta2))
SLL_ASSERT(nx3 == size(eta3))

if (present(iplot)) then
  call sll_s_int2string(iplot, cplot)
  call sll_s_xml_file_create(file_name//cplot//".xmf",file_id,error)
else
  call sll_s_xml_file_create(file_name//".xmf",file_id,error)
end if

write(file_id,"(a)")"<Grid Name='mesh' GridType='Uniform'>"
write(file_id,"(a,3i5,a)")"<Topology TopologyType='3DRectMesh' NumberOfElements='", &
                          nx3,nx2,nx1,"'/>"
write(file_id,"(a)")"<Geometry GeometryType='VXVYVZ'>"
write(file_id,"(a,i5,a)")"<DataItem Dimensions='",nx1, &
                         "' NumberType='Float' Format='XML'>"
write(file_id,*) (eta1(i),i=1,nx1)
write(file_id,"(a)")"</DataItem>"
write(file_id,"(a,i5,a)")"<DataItem Dimensions='",nx2, &
                             "' NumberType='Float' Format='XML'>"
write(file_id,*) (eta2(j),j=1,nx2)
write(file_id,"(a)")"</DataItem>"
write(file_id,"(a,i5,a)")"<DataItem Dimensions='",nx3, &
                             "' NumberType='Float' Format='XML'>"
write(file_id,*) (eta3(k),k=1,nx3)
write(file_id,"(a)")"</DataItem>"
write(file_id,"(a)")"</Geometry>"
write(file_id,"(a)")"<Attribute Name='"//array_name//"' AttributeType='Scalar' Center='Node'>"
if(present(file_format) .and. file_format == "HDF5") then
  write(file_id,"(a,3i5,a)")"<DataItem Dimensions='",nx3,nx2,nx1, &
                                 "' NumberType='Float' Precision='8' Format='HDF'>"
  if (present(iplot)) then
    write(file_id,"(a)")array_name//cplot//".h5:/node_values"
  else
    write(file_id,"(a)")array_name//".h5:/node_values"
  end if
#ifndef NOHDF5
  if (present(iplot)) then
    call sll_s_hdf5_ser_file_create( array_name//cplot//".h5", hfile_id, error )
  else
    call sll_s_hdf5_ser_file_create( array_name//".h5", hfile_id, error )
  end if
  call sll_o_hdf5_ser_write_array( hfile_id, array, "/node_values", error )
  call sll_s_hdf5_ser_file_close( hfile_id, error )
#endif
else
  write(file_id,"(a,3i5,a)")"<DataItem Dimensions='",nx3,nx2,nx1, &
                            "' NumberType='Float' Precision='4' Format='XML'>"
  call sll_o_ascii_write_array(file_id,array,error)
end if
write(file_id,"(a)")"</DataItem>"
write(file_id,"(a)")"</Attribute>"
call sll_s_xml_file_close(file_id,error)

end subroutine sll_s_xdmf_rect3d_nodes


!>Subroutine to write a 2D array in xdmf format.
!>The field is describe on a cartesian mesh.
!>Nodes coordinates are defined by x and y (2d arrays).
subroutine sll_s_xdmf_curv2d_nodes( file_name,   &
                                    array,       &
                                    array_name,  &
                                    eta1,        &
                                    eta2,        &
                                    file_format, &
                                    iplot        ) 

sll_real64, intent(in)       :: array(:,:)  !< data array
sll_real64, intent(in)       :: eta1(:,:)   !< x data
sll_real64, intent(in)       :: eta2(:,:)   !< y data
character(len=*), intent(in) :: file_name   !< xmf file name
character(len=*), intent(in) :: array_name  !< array name
sll_int32                    :: error       !< error code
sll_int32                    :: file_id     !< xmf file unit number
sll_int32                    :: nx1         !< x nodes number
sll_int32                    :: nx2         !< y nodes number
character(len=4), optional   :: file_format !< file format "HDF5" or "Binary"
sll_int32,        optional   :: iplot       !< plot index
character(len=4)             :: cplot 
#ifndef NOHDF5
integer(hid_t)               :: hfile_id    !< h5 file unit number
#endif

nx1 = size(array,1)
nx2 = size(array,2)

SLL_ASSERT(nx1 == size(eta1,1))
SLL_ASSERT(nx2 == size(eta2,2))

if (present(iplot)) then
  call sll_s_int2string(iplot, cplot)
  call sll_s_xml_file_create(file_name//cplot//".xmf",file_id,error)
else
  call sll_s_xml_file_create(file_name//".xmf",file_id,error)
end if

write(file_id,"(a)")"<Grid Name='mesh' GridType='Uniform'>"
write(file_id,"(a,2i5,a)")"<Topology TopologyType='2DSMesh' NumberOfElements='", &
                          nx2,nx1,"'/>"
write(file_id,"(a)")"<Geometry GeometryType='X_Y'>"

if(present(file_format) .and. file_format == "HDF5") then

  write(file_id,"(a,2i5,a)")"<DataItem Dimensions='",nx2,nx1, &
                            "' NumberType='Float' Precision='8' Format='HDF'>"
  write(file_id,"(a)")array_name//".h5:/x1_values"
  write(file_id,"(a)")"</DataItem>"
  write(file_id,"(a,2i5,a)")"<DataItem Dimensions='",nx2,nx1, &
                            "' NumberType='Float' Precision='8' Format='HDF'>"
  if (present(iplot)) then
    write(file_id,"(a)")array_name//cplot//".h5:/x2_values"
  else
    write(file_id,"(a)")array_name//".h5:/x2_values"
  end if
  write(file_id,"(a)")"</DataItem>"

#ifndef NOHDF5
  call sll_s_hdf5_ser_file_create( file_name//".h5", hfile_id, error )
  call sll_o_hdf5_ser_write_array( hfile_id, eta1, "/x1_values", error )
  call sll_o_hdf5_ser_write_array( hfile_id, eta2, "/x2_values", error )
#endif

else

  write(file_id,"(a,2i5,a)")"<DataItem Dimensions='",nx2,nx1, &
                            "' NumberType='Float' Precision='4' Format='XML'>"
  call sll_o_ascii_write_array(file_id,eta1,error)
  write(file_id,"(a)")"</DataItem>"
  write(file_id,"(a,2i5,a)")"<DataItem Dimensions='",nx2,nx1, &
                            "' NumberType='Float' Precision='4' Format='XML'>"
  call sll_o_ascii_write_array(file_id,eta2,error)
  write(file_id,"(a)")"</DataItem>"

end if

write(file_id,"(a)")"</Geometry>"
write(file_id,"(a)") &
"<Attribute Name='"//array_name//"' AttributeType='Scalar' Center='Node'>"

#ifndef NOHDF5
if(present(file_format) .and. file_format == "HDF5") then

  write(file_id,"(a,2i5,a)")"<DataItem Dimensions='",nx2,nx1, &
                            "' NumberType='Float' Precision='8' Format='HDF'>"
  if (present(iplot)) then
    write(file_id,"(a)")array_name//cplot//".h5:/node_values"
  else
    write(file_id,"(a)")array_name//".h5:/node_values"
  end if

  call sll_o_hdf5_ser_write_array( hfile_id, array, "/node_values", error )
  call sll_s_hdf5_ser_file_close( hfile_id, error )

else
#endif

  write(file_id,"(a,2i5,a)")"<DataItem Dimensions='",nx2,nx1, &
                            "' NumberType='Float' Precision='4' Format='XML'>"
  call sll_o_ascii_write_array(file_id,array,error)

#ifndef NOHDF5
end if
#endif

write(file_id,"(a)")"</DataItem>"
write(file_id,"(a)")"</Attribute>"
call sll_s_xml_file_close(file_id,error)

end subroutine sll_s_xdmf_curv2d_nodes

!>@brief
!>Subroutine to write a 3D array in xdmf format.
!>The field is describe on a cartesian mesh.
!>Nodes coordinates are defined by x,y,z (3d arrays).
!>@details
!> Example:
!>@code
!>call sll_o_xdmf_open(file_name,mesh_name,nnodes_x1,nnodes_x2,nnodes_x3,file_id,error)
!>call sll_o_xdmf_write_array(mesh_name,x1,'x1',error)
!>call sll_o_xdmf_write_array(mesh_name,x2,'x2',error)
!>call sll_o_xdmf_write_array(mesh_name,x3,'x3',error)
!>call sll_o_xdmf_write_array("field3d",df,"NodeVal",error,file_id,"Node")
!>call sll_o_xdmf_write_array("field3d",df(1:ncells_x1,1:ncells_x2,1:ncells_x3), &
!>                          "CellVal",error,file_id,"Cell")
!>call sll_s_xdmf_close(file_id,error)
!>@endcode

subroutine sll_s_xdmf_curv3d_nodes( file_name,   &
                                    array,       &
                                    array_name,  &
                                    eta1,        &
                                    eta2,        &
                                    eta3,        &
                                    file_format, & 
                                    iplot        )

sll_real64, intent(in)       :: array(:,:,:)  !< data array
sll_real64, intent(in)       :: eta1(:,:,:)   !< x data
sll_real64, intent(in)       :: eta2(:,:,:)   !< y data
sll_real64, intent(in)       :: eta3(:,:,:)   !< z data
character(len=*), intent(in) :: file_name     !< xmf file name
character(len=*), intent(in) :: array_name    !< array name
sll_int32                    :: error         !< error code
sll_int32                    :: file_id       !< xmf file unit number
sll_int32                    :: nx1           !< x nodes number
sll_int32                    :: nx2           !< y nodes number
sll_int32                    :: nx3           !< z nodes number
character(len=4), optional   :: file_format   !< file format "HDF5" or "Binary"
sll_int32       , optional   :: iplot         !< plot index

character(len=4)             :: cplot         
#ifndef NOHDF5
integer(hid_t)               :: hfile_id    !< h5 file unit number
#endif

nx1 = size(array,1)
nx2 = size(array,2)
nx3 = size(array,3)

SLL_ASSERT(nx1 == size(eta1,1))
SLL_ASSERT(nx2 == size(eta2,2))
SLL_ASSERT(nx3 == size(eta3,3))

if (present(iplot)) then
  call sll_s_int2string(iplot, cplot)
  call sll_s_xml_file_create(file_name//cplot//".xmf",file_id,error)
else
  call sll_s_xml_file_create(file_name//".xmf",file_id,error)
end if

write(file_id,"(a)")"<Grid Name='mesh' GridType='Uniform'>"
write(file_id,"(a,3i5,a)")"<Topology TopologyType='3DSMesh' NumberOfElements='", &
                          nx3,nx2,nx1,"'/>"
write(file_id,"(a)")"<Geometry GeometryType='X_Y_Z'>"

if(present(file_format) .and. file_format == "HDF5") then

  write(file_id,"(a,3i5,a)")"<DataItem Dimensions='",nx3,nx2,nx1, &
                            "' NumberType='Float' Precision='8' Format='HDF'>"
  write(file_id,"(a)")array_name//".h5:/x1_values"
  write(file_id,"(a)")"</DataItem>"
  write(file_id,"(a,3i5,a)")"<DataItem Dimensions='",nx3,nx2,nx1, &
                            "' NumberType='Float' Precision='8' Format='HDF'>"
  write(file_id,"(a)")array_name//".h5:/x2_values"
  write(file_id,"(a)")"</DataItem>"
  write(file_id,"(a,3i5,a)")"<DataItem Dimensions='",nx3,nx2,nx1, &
                            "' NumberType='Float' Precision='8' Format='HDF'>"
  write(file_id,"(a)")array_name//".h5:/x3_values"
  write(file_id,"(a)")"</DataItem>"

#ifndef NOHDF5
  call sll_s_hdf5_ser_file_create( array_name//".h5", hfile_id, error )
  call sll_o_hdf5_ser_write_array( hfile_id, eta1, "/x1_values", error )
  call sll_o_hdf5_ser_write_array( hfile_id, eta2, "/x2_values", error )
  call sll_o_hdf5_ser_write_array( hfile_id, eta3, "/x3_values", error )
#endif

else

  write(file_id,"(a,3i5,a)")"<DataItem Dimensions='",nx3,nx2,nx1, &
                            "' NumberType='Float' Precision='4' Format='XML'>"
  call sll_o_ascii_write_array(file_id,eta1,error)
  write(file_id,"(a)")"</DataItem>"
  write(file_id,"(a,3i5,a)")"<DataItem Dimensions='",nx3,nx2,nx1, &
                            "' NumberType='Float' Precision='4' Format='XML'>"
  call sll_o_ascii_write_array(file_id,eta2,error)
  write(file_id,"(a)")"</DataItem>"
  write(file_id,"(a,3i5,a)")"<DataItem Dimensions='",nx3,nx2,nx1, &
                             "' NumberType='Float' Precision='4' Format='XML'>"
  call sll_o_ascii_write_array(file_id,eta3,error)
  write(file_id,"(a)")"</DataItem>"

end if

write(file_id,"(a)")"</Geometry>"
if (present(iplot)) then
  write(file_id,"(a)") &
  "<Attribute Name='"//array_name//cplot//"' AttributeType='Scalar' Center='Node'>"
else
  write(file_id,"(a)") &
  "<Attribute Name='"//array_name//"' AttributeType='Scalar' Center='Node'>"
end if

#ifndef NOHDF5
if(present(file_format) .and. file_format == "HDF5") then

   write(file_id,"(a,3i5,a)")"<DataItem Dimensions='",nx3,nx2,nx1, &
                             "' NumberType='Float' Precision='8' Format='HDF'>"
  if (present(iplot)) then
    write(file_id,"(a)")array_name//cplot//".h5:/node_values"
  else
    write(file_id,"(a)")array_name//".h5:/node_values"
  end if

   call sll_o_hdf5_ser_write_array( hfile_id, array," /node_values", error )
   call sll_s_hdf5_ser_file_close( hfile_id, error )

else
#endif

   write(file_id,"(a,3i5,a)")"<DataItem Dimensions='",nx3,nx2,nx1, &
                             "' NumberType='Float' Precision='4' Format='XML'>"
   call sll_o_ascii_write_array(file_id,array,error)

#ifndef NOHDF5
end if
#endif

write(file_id,"(a)")"</DataItem>"
write(file_id,"(a)")"</Attribute>"
call sll_s_xml_file_close(file_id,error)

end subroutine sll_s_xdmf_curv3d_nodes

!> Close the XML file and finish to write last lines.
subroutine sll_s_xdmf_close(file_id,error)
sll_int32, intent(in)  :: file_id !< xml file unit number
sll_int32, intent(out) :: error   !< error code
       
write(file_id,"(a)")"</Grid>"
write(file_id,"(a)")"</Domain>"
write(file_id,"(a)")"</Xdmf>"
close(file_id)
error = 0
end subroutine sll_s_xdmf_close

!> Plot 2d distribution function for VisIt
!> @param[in]  iplot      plot counter.
!> @param[in]  f          function values .
!> @param[in]  vec_x1     node positions on x1.
!> @param[in]  nnodes_x1  nodes number on x1.
!> @param[in]  vec_x2     node positions on x2.
!> @param[in]  nnodes_x2  nodes number on x2.
!> @param[in]  array_name file name.
!> @param[in]  time       Add time value on plot title.
!> @details
!> This routine will create a file named array_name_[iplot].xmf
subroutine sll_s_plot_f_cartesian( iplot,      &
                                   f,          &
                                   vec_x1,     &
                                   nnodes_x1,  &
                                   vec_x2,     &
                                   nnodes_x2,  &
                                   array_name, &
                                   time)    

  sll_int32,                  intent(in) :: iplot
  sll_real64, dimension(:,:), intent(in) :: f
  sll_real64, dimension(:),   intent(in) :: vec_x1
  sll_int32,                  intent(in) :: nnodes_x1
  sll_real64, dimension(:),   intent(in) :: vec_x2    
  sll_int32,                  intent(in) :: nnodes_x2
  character(len=*),           intent(in) :: array_name 
  sll_real64                             :: time

  sll_int32                               :: file_id
  sll_int32                               :: error
  sll_real64, dimension(:,:), allocatable :: x1
  sll_real64, dimension(:,:), allocatable :: x2
  sll_int32                               :: i
  sll_int32                               :: j
  character(len=4)                        :: cplot

#ifndef NOHDF5
integer(hid_t) :: hfile_id
#endif
  
  SLL_ASSERT(iplot > 0)
  if (iplot == 1) then

    SLL_ALLOCATE(x1(nnodes_x1,nnodes_x2), error)
    SLL_ALLOCATE(x2(nnodes_x1,nnodes_x2), error)
    do j = 1,nnodes_x2
      do i = 1,nnodes_x1
        x1(i,j) = vec_x1(i) !x1_min+real(i-1,f32)*dx1
        x2(i,j) = vec_x2(j) !x2_min+real(j-1,f32)*dx2
      end do
    end do

#ifndef NOHDF5
    call sll_s_hdf5_ser_file_create( "cartesian_mesh-x1.h5", hfile_id, error )
    call sll_o_hdf5_ser_write_array( hfile_id, x1, "/x1", error )
    call sll_s_hdf5_ser_file_close( hfile_id, error )
    call sll_s_hdf5_ser_file_create( "cartesian_mesh-x2.h5", hfile_id, error )
    call sll_o_hdf5_ser_write_array( hfile_id, x2, "/x2", error )
    call sll_s_hdf5_ser_file_close( hfile_id, error )
#endif

    deallocate(x1)
    deallocate(x2)

  end if

  call sll_s_int2string(iplot,cplot)
  call sll_o_xdmf_open(trim(array_name)//cplot//".xmf","cartesian_mesh", &
                     nnodes_x1,nnodes_x2,file_id,error)

  write(file_id,"(a,f8.3,a)") "<Time Value='",time,"'/>"

  call sll_o_xdmf_write_array(trim(array_name)//cplot,f,"values", &
                            error,file_id,"Node")
  call sll_s_xdmf_close(file_id,error)

end subroutine sll_s_plot_f_cartesian


!> Plot 2d distribution function for VisIt
!> @param[in]  iplot      plot counter.
!> @param[in]  f          function values .
!> @param[in]  nnodes_x1  nodes number on x1.
!> @param[in]  nnodes_x2  nodes number on x2.
!> @param[in]  array_name file name.
!> @param[in]  mesh_name  file name for mesh
!> @param[in]  time       Add time value on plot title.
!> @param[in]  x1         2d-array for x1 (create mesh if provided)
!> @param[in]  x2         2d-array for x2 (create mesh if provided)
!> @details
!> This routine will create a file named array_name_[iplot].xmf
!> TODO suggestion: merge sll_s_plot_f and sll_s_plot_f_cartesian
 
subroutine sll_s_plot_f( &
  iplot, &
  f, &  
  nnodes_x1, &
  nnodes_x2,  &
  array_name, &
  mesh_name, &
  time, &
  x1, &
  x2)    

  sll_int32,  intent(in) :: iplot
  sll_real64, intent(in) :: f(:,:)
  sll_int32,  intent(in) :: nnodes_x1
  sll_int32,  intent(in) :: nnodes_x2
  character(len=*), intent(in) :: array_name 
  character(len=*), intent(in) :: mesh_name 
  sll_real64, intent(in) :: time
  sll_real64, intent(in), optional :: x1(:,:)
  sll_real64, intent(in), optional :: x2(:,:)

  sll_int32 :: file_id
  sll_int32 :: ierr
  character(len=4) :: cplot

#ifndef NOHDF5
  integer(hid_t) :: hfile_id
  
  if (present(x1).and.present(x2)) then


    call sll_s_hdf5_ser_file_create( trim(mesh_name)//"-x1.h5", hfile_id, ierr )
    call sll_o_hdf5_ser_write_array( hfile_id, x1, "/x1", ierr )
    call sll_s_hdf5_ser_file_close( hfile_id, ierr )
    call sll_s_hdf5_ser_file_create( trim(mesh_name)//"-x2.h5", hfile_id, ierr )
    call sll_o_hdf5_ser_write_array( hfile_id, x2, "/x2", ierr )
    call sll_s_hdf5_ser_file_close( hfile_id, ierr )

  end if

#endif

  call sll_s_int2string(iplot,cplot)
  call sll_o_xdmf_open( &
    trim(array_name)//cplot//".xmf", &
    trim(mesh_name), &
    nnodes_x1, &
    nnodes_x2, &
    file_id, &
    ierr)

  write(file_id,"(a,f8.3,a)") "<Time Value='",time,"'/>"

  call sll_o_xdmf_write_array( &
    trim(array_name)//cplot, &
    f, &
    "values", &
    ierr, &
    file_id, &
    "Node")
    
  call sll_s_xdmf_close(file_id,ierr)

end subroutine sll_s_plot_f



end module sll_m_xdmf
