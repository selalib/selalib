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

!> @author Pierre Navaro
!> @brief
!> Implements the functions to write xml file to store light data
!> @details
!> With XDMF file you can describe data to plot them with VisIt
module sll_xml_io
#include "sll_working_precision.h"
#include "sll_assert.h"
  
  !> write a data item in the xml file
  interface sll_xml_dataitem
     module procedure sll_xml_dataitem_2d
     module procedure sll_xml_dataitem_3d
  end interface
  
  !> write a data attribute in the xml file
  interface sll_xml_field
     module procedure sll_xml_field_2d
     module procedure sll_xml_field_3d
  end interface
  
  !> write grid description in the xml file
  interface sll_xml_grid_geometry
     module procedure sll_xml_grid_geometry_2d_high_level
     module procedure sll_xml_grid_geometry_2d_low_level
     module procedure sll_xml_grid_geometry_3d_high_level
     module procedure sll_xml_grid_geometry_3d_low_level
  end interface
  
contains
  
  !> Create the XML file and begin to write first lines.
  !> You get the file unit number.
  subroutine sll_xml_file_create(filename,file_id,error)
    character(len=*) , intent(in)  :: filename   !< file name
    sll_int32        , intent(out) :: file_id    !< file unit number
    sll_int32        , intent(out) :: error      !< error code
    logical                        :: lopen
    
    error=0
    do 100 file_id=20,99
       inquire(unit=file_id,opened=lopen)
       if(lopen) then
          cycle
       else
          open(file_id,status='SCRATCH',err=100)
          close(file_id,status='DELETE',err=100)
          goto 200
       end if
       
100    continue
       error=1
200    continue
       error=0
       
       ! Create a new file using default properties
       inquire(file=filename,opened=lopen)
       SLL_ASSERT(.not. lopen)
       
       open(file_id,FILE=filename,FORM='FORMATTED',IOSTAT=error)
       rewind(file_id)
       
       write(file_id,"(a)")"<?xml version='1.0' ?>"
       write(file_id,"(a)")"<!DOCTYPE Xdmf SYSTEM 'Xdmf.dtd' []>"
       write(file_id,"(a)")"<Xdmf Version='2.0'>"
       write(file_id,"(a)")"<Domain>"
       
     end subroutine sll_xml_file_create
     
     !> Close the XML file and finish to write last lines.
     !> You give the file unit number.
     !> \param[in]  file_id - the unit number or your xml file
     !> \param[out] error   - error parameter
     subroutine sll_xml_file_close(file_id,error)
       sll_int32, intent(in)  :: file_id !< xml file unit number
       sll_int32, intent(out) :: error   !< error code
       
       write(file_id,"(a)")"</Grid>"
       write(file_id,"(a)")"</Domain>"
       write(file_id,"(a)")"</Xdmf>"
       close(file_id)
       error = 0
     end subroutine sll_xml_file_close
     
     !> Write the description of a scalar field on a 2D mesh.
     !> \param[in] file_id is the unit number or your xml file
     !> \param[in] filename is the file name where the heavy data are 
     !> (bin or h5)
     !> \param[in] nnodes_x1 - nodes number along direction 1
     !> \param[in] nnodes_x2 - nodes number along direction 2
     !> \param[in] filetype  - heavy data format 'HDF' or 'Binary'
     !>
     !> The file named filename must exist.
     !>
     subroutine sll_xml_dataitem_2d( file_id,   &
                                     filename,  &
                                     nnodes_x1, &
                                     nnodes_x2, &
                                     filetype )
       
       sll_int32, intent(in)        :: file_id   !< file unit number
       character(len=*), intent(in) :: filename  !< xmf file name
       character(len=*), intent(in) :: filetype  !< data file format
       sll_int32, intent(in)        :: nnodes_x1 !< x nodes number
       sll_int32, intent(in)        :: nnodes_x2 !< y nodes number
       
       SLL_ASSERT(filetype == 'HDF' .or. filetype == 'Binary')
       write(file_id,"(a,2i5,a)")"<DataItem Dimensions='",nnodes_x2,nnodes_x1, &
            "' NumberType='Float' Precision='8' Format='"//trim(filetype)//"'>"
       write(file_id,"(a)")trim(filename)
       write(file_id,"(a)")"</DataItem>"
     end subroutine sll_xml_dataitem_2d
     
     !> Write the description of a scalar field on a 3D mesh.
     subroutine sll_xml_dataitem_3d( file_id,   &
                                     filename,  &
                                     nnodes_x1, &
                                     nnodes_x2, &
                                     nnodes_x3, &
                                     filetype)
       
       sll_int32, intent(in)        :: file_id   !< file unit number
       character(len=*), intent(in) :: filename  !< xmf file name
       character(len=*), intent(in) :: filetype  !< "HDF" or "Binary"
       sll_int32, intent(in)        :: nnodes_x1 !< x nodes number
       sll_int32, intent(in)        :: nnodes_x2 !< y nodes number
       sll_int32, intent(in)        :: nnodes_x3 !< z nodes number
       
       SLL_ASSERT(filetype == 'HDF' .or. filetype == 'Binary')
       write(file_id,"(a,3i5,a)")"<DataItem Dimensions='",nnodes_x3, &
            nnodes_x2,nnodes_x1, &
            "' NumberType='Float' Precision='8' Format='"//trim(filetype)//"'>"
       write(file_id,"(a)")trim(filename)
       write(file_id,"(a)")"</DataItem>"
     end subroutine sll_xml_dataitem_3d
     
     !> Write the description of a scalar field on a 2D mesh.
     !> \param[in] fieldname the dataset name where the heavy data are 
     !> (hdf5 case)
     !> \param[in] filename  the file name where the heavy data are 
     !> (bin or h5)
     !> \param[in] npoints_1 nodes or cells number along direction 1
     !> \param[in] npoints_2 nodes or cells number along direction 2
     !> \param[in] center    values are centered on nodes or cells 
     !>
     !> The file named filename-fieldname.bin must exist in case of binary 
     !> output.
     !> The file named filename.h5 with dataset fieldname must exist in case 
     !> of hdf5 output.
     subroutine sll_xml_field_2d( file_id,   &
                                  fieldname, &
                                  filename,  &
                                  npoints_1, &
                                  npoints_2, &
                                  filetype,  &
                                  center)

       sll_int32,        intent(in) :: file_id   !< the unit number or your xml file
       character(len=*), intent(in) :: filename 
       character(len=*), intent(in) :: fieldname
       character(len=*), intent(in) :: center
       character(len=*), intent(in) :: filetype   !< "HDF" or "Binary"
       sll_int32,        intent(in) :: npoints_1
       sll_int32,        intent(in) :: npoints_2
       
       write(file_id,"(a)") &
       "<Attribute Name='"//trim(fieldname)//"' AttributeType='Scalar' Center='"//center//"'>"
       call sll_xml_dataitem_2d(file_id,filename,npoints_1,npoints_2,filetype)
       write(file_id,"(a)")"</Attribute>"
     end subroutine sll_xml_field_2d
     
     !> Write the description of a scalar field on a 3D mesh.
     !> \param[in] file_id   the unit number or your xml file
     !> \param[in] fieldname the dataset name where the heavy data are 
     !> (hdf5 case)
     !> \param[in] filename  the file name where the heavy data are 
     !> (bin or h5)
     !> \param[in] npoints_1 nodes or cells number along direction 1
     !> \param[in] npoints_2 nodes or cells number along direction 2
     !> \param[in] npoints_3 nodes or cells number along direction 3
     !> \param[in] center    values are centered on nodes or cells 
     !>
     !> The file named filename-fieldname.bin must exist in case of binary 
     !> output.
     !> The file named filename.h5 with dataset fieldname must exist in case 
     !> of hdf5 output.
     subroutine sll_xml_field_3d( file_id,   &
                                  fieldname, &
                                  filename,  &
                                  npoints_1, &
                                  npoints_2, &
                                  npoints_3, &
                                  filetype,  &
                                  center)

       sll_int32,        intent(in) :: file_id
       character(len=*), intent(in) :: filename
       character(len=*), intent(in) :: fieldname
       character(len=*), intent(in) :: center
       character(len=*), intent(in) :: filetype   !< "HDF" or "Binary"
       sll_int32,        intent(in) :: npoints_1
       sll_int32,        intent(in) :: npoints_2
       sll_int32,        intent(in) :: npoints_3
       
       write(file_id,"(a)") &
       "<Attribute Name='"//trim(fieldname)//"' AttributeType='Scalar' Center='"//center//"'>"
       call sll_xml_dataitem_3d( file_id, &
            filename, &
            npoints_1, &
            npoints_2, &
            npoints_3, &
            filetype)
       write(file_id,"(a)")"</Attribute>"
     end subroutine sll_xml_field_3d
     
     !> Write the description of a 2D strutured grid
     !> mesh with its nodes coordinates contains in filename-x1 and filename-x2.
     !> \param[in] file_id is the unit number or your xml file
     !> \param[in] filename is the file name where the coordinates data are (bin or h5)
     !> \param[in] nnodes_x1 - nodes number along direction 1
     !> \param[in] nnodes_x2 - nodes number along direction 2
     !>
     !> The file named filename-x1.bin and filename-x2.bin must exist in case of binary output.
     !> The file named filename.h5 with dataset x1 and x2 must exist in case of hdf5 output.
     !>
     subroutine sll_xml_grid_geometry_2d_high_level( &
       file_id, &
       filename, &
       nnodes_x1, &
       nnodes_x2 )

       sll_int32, intent(in) :: file_id
       character(len=*), intent(in) :: filename
       sll_int32, intent(in) :: nnodes_x1
       sll_int32, intent(in) :: nnodes_x2
       
#ifdef NOHDF5
       call sll_xml_grid_geometry_2d_low_level( file_id, &
            trim(filename)//"-x1.bin", nnodes_x1, &
            trim(filename)//"-x2.bin", nnodes_x2 )
#else
       call sll_xml_grid_geometry_2d_low_level( file_id, &
            trim(filename)//"-x1.h5", nnodes_x1, &
            trim(filename)//"-x2.h5", nnodes_x2, "x1", "x2" )
#endif
       
       
     end subroutine sll_xml_grid_geometry_2d_high_level

     !> Write the description of a 2D strutured grid
     !> mesh with its nodes coordinates contains in filename-x1 and filename-x2.
     !> \param[in] file_id is the unit number or your xml file
     !> \param[in] x1filename is the file name where the coordinates x1 are (bin or h5)
     !> \param[in] x2filename is the file name where the coordinates x2 are (bin or h5)
     !> \param[in] x1dsetname is the dataset name of coordinates x1 are (bin or h5)
     !> \param[in] x2dsetname is the dataset name of coordinates x2 are (bin or h5)
     !> \param[in] nnodes_x1 - nodes number along direction 1
     !> \param[in] nnodes_x2 - nodes number along direction 2
     !>
     !> The file named x*filename-x*dsetname.bin must exists
     !> The file named x*filename-x*dsetname.h5 with dataset x*dsetname must exists. \n
     !> Low level version where you have to set dataset names in hdf5 files
     subroutine sll_xml_grid_geometry_2d_low_level( file_id,    &
                                                    x1filename, &
                                                    nnodes_x1,  &
                                                    x2filename, &
                                                    nnodes_x2,  &
                                                    x1dsetname, &
                                                    x2dsetname) 

       sll_int32, intent(in)        :: file_id
       character(len=*), intent(in) :: x1filename
       character(len=*), intent(in) :: x2filename
       sll_int32, intent(in)        :: nnodes_x1
       sll_int32, intent(in)        :: nnodes_x2
       character(len=*), optional   :: x1dsetname
       character(len=*), optional   :: x2dsetname
       
       write(file_id,"(a)")"<Grid Name='mesh' GridType='Uniform'>"
       write(file_id, &
            "(a,2i5,a)")"<Topology TopologyType='2DSMesh' NumberOfElements='", &
            nnodes_x2,nnodes_x1,"'/>"
       write(file_id,"(a)")"<Geometry GeometryType='X_Y'>"
       
#ifdef NOHDF5
       call sll_xml_dataitem_2d(file_id, &
                                trim(x1filename), &
                                nnodes_x1,nnodes_x2,'Binary')
       call sll_xml_dataitem_2d(file_id, &
                                trim(x2filename), &
                                nnodes_x1,nnodes_x2,'Binary')
#else
       call sll_xml_dataitem_2d(file_id, &
                                trim(x1filename)//":/"//x1dsetname, &
                                nnodes_x1,nnodes_x2,'HDF')
       call sll_xml_dataitem_2d(file_id, &
                                trim(x2filename)//":/"//x2dsetname, &
                                nnodes_x1,nnodes_x2,'HDF')
#endif
       
       write(file_id,"(a)")"</Geometry>"
     end subroutine sll_xml_grid_geometry_2d_low_level
     

     !> Write the description of a 3D structured curvilinear grid
     !> mesh with its nodes coordinates contains in filename-x1 and filename-x2.
     !> High level version where dataset names in hdf5 files are set automatically
     subroutine sll_xml_grid_geometry_3d_high_level(file_id, filename,  &
       nnodes_x1, nnodes_x2, nnodes_x3)
       
       sll_int32, intent(in)        :: file_id    !< xmf file unit number
       character(len=*), intent(in) :: filename   !< xmf file name
       sll_int32, intent(in)        :: nnodes_x1  !< x nodes number
       sll_int32, intent(in)        :: nnodes_x2  !< y nodes number
       sll_int32, intent(in)        :: nnodes_x3  !< z nodes number
       
       write(file_id,"(a)")"<Grid Name='mesh' GridType='Uniform'>"
       write(file_id,"(a,3i5,a)")"<Topology TopologyType='3DSMesh' NumberOfElements='", &
            nnodes_x3,nnodes_x2,nnodes_x1,"'/>"
       write(file_id,"(a)")"<Geometry GeometryType='X_Y_Z'>"
       
#ifdef NOHDF5
       
       call sll_xml_dataitem_3d(file_id,trim(filename)//"-x1.bin", &
                                nnodes_x1,nnodes_x2,nnodes_x3,'Binary')
       call sll_xml_dataitem_3d(file_id,trim(filename)//"-x2.bin", &
                                nnodes_x1,nnodes_x2,nnodes_x3,'Binary')
       call sll_xml_dataitem_3d(file_id,trim(filename)//"-x3.bin", &
                                nnodes_x1,nnodes_x2,nnodes_x3,'Binary')
#else
       
       call sll_xml_dataitem_3d(file_id,trim(filename)//"-x1.h5:/x1", &
                                nnodes_x1,nnodes_x2,nnodes_x3,'HDF')
       call sll_xml_dataitem_3d(file_id,trim(filename)//"-x2.h5:/x2", &
                                nnodes_x1,nnodes_x2,nnodes_x3,'HDF')
       call sll_xml_dataitem_3d(file_id,trim(filename)//"-x3.h5:/x3", &
                                nnodes_x1,nnodes_x2,nnodes_x3,'HDF')
       
#endif
       
       write(file_id,"(a)")"</Geometry>"
       
     end subroutine sll_xml_grid_geometry_3d_high_level

     !> Write the description of a 3D structured curvilinear grid
     !> mesh with its nodes coordinates contains in filename-x1 and filename-x2.
     !> Low level version where dataset names in hdf5 files must be set.
     subroutine sll_xml_grid_geometry_3d_low_level(file_id,       &
                             x1filename, nnodes_x1,   &
                             x2filename, nnodes_x2,   &
                             x3filename, nnodes_x3,   &
                             x1dsetname, x2dsetname, x3dsetname  )
       
       sll_int32, intent(in)                  :: file_id    !< xmf unif file number
       character(len=*), intent(in)           :: x1filename !< x data file name
       character(len=*), intent(in)           :: x2filename !< y data file name
       character(len=*), intent(in)           :: x3filename !< x datz file name
       character(len=*), intent(in), optional :: x1dsetname !< x dataset name
       character(len=*), intent(in), optional :: x2dsetname !< y dataset name
       character(len=*), intent(in), optional :: x3dsetname !< z dataset name
       sll_int32, intent(in)                  :: nnodes_x1  !< x nodes number
       sll_int32, intent(in)                  :: nnodes_x2  !< y nodes number
       sll_int32, intent(in)                  :: nnodes_x3  !< z nodes number
       
       write(file_id,"(a)")"<Grid Name='mesh' GridType='Uniform'>"
       write(file_id,"(a,3i5,a)") &
          "<Topology TopologyType='3DSMesh' NumberOfElements='", &
          nnodes_x3,nnodes_x2,nnodes_x1,"'/>"
       write(file_id,"(a)")"<Geometry GeometryType='X_Y_Z'>"
       
#ifdef NOHDF5
       
       call sll_xml_dataitem_3d(file_id,x1filename, &
                                nnodes_x1,nnodes_x2,nnodes_x3,'Binary')
       call sll_xml_dataitem_3d(file_id,x2filename, &
                                nnodes_x1,nnodes_x2,nnodes_x3,'Binary')
       call sll_xml_dataitem_3d(file_id,x3filename, &
                                nnodes_x1,nnodes_x2,nnodes_x3,'Binary')
#else
       
       call sll_xml_dataitem_3d(file_id,x1filename//":/"//x1dsetname, &
                                nnodes_x1,nnodes_x2,nnodes_x3,'HDF')
       call sll_xml_dataitem_3d(file_id,x2filename//":/"//x2dsetname, &
                                nnodes_x1,nnodes_x2,nnodes_x3,'HDF')
       call sll_xml_dataitem_3d(file_id,x3filename//":/"//x3dsetname, &
                                nnodes_x1,nnodes_x2,nnodes_x3,'HDF')
       
#endif
       
       write(file_id,"(a)")"</Geometry>"
       
     end subroutine sll_xml_grid_geometry_3d_low_level
     
     
end module sll_xml_io
   
