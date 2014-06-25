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
!> Implements the functions to write xdmf file plotable by VisIt
!> @details
!> In <b> XDMF </b> (eXtensible Data Model and Format) the description of the 
!> data is separate from the values themselves. Light data is stored using XML, 
!> Heavy data is stored using Parallel HDF5. These files are readable by 
!> Paraview.
module sll_xdmf_serial_blocks
#include "sll_working_precision.h"
#include "sll_assert.h"
  
  use sll_collective
  use sll_utilities
#ifndef NOHDF5
  use hdf5
  use sll_hdf5_io_serial
#endif
  use sll_ascii_io
  use sll_xml_io
  
  implicit none
  
  interface sll_xdmf_open
     module procedure sll_xdmf_open_serial_blocks
  end interface

  !> Write and array in an xmf file
  interface sll_xdmf_write_array
     module procedure sll_xdmf_array_2d_serial_blocks
  end interface
  
  !> Close the xdmf file
  interface sll_xdmf_close
     module procedure sll_xdmf_close_serial_blocks
  end interface

contains  
  
  !>Open a XDMF format file for a 2d plot
  subroutine sll_xdmf_open_serial_blocks( file_name,  &
                                          file_id,    &
                                          error)

    character(len=*), intent(in) :: file_name !< xmf file name 
    sll_int32                    :: file_id   !< file unit number
    sll_int32                    :: error     !< error code
    sll_int32                    :: prank     !< processor number id
    
    prank = sll_get_collective_rank(sll_world_collective)
    if (prank == 0) then
       call sll_xml_file_create(trim(file_name),file_id,error)
       write(file_id,'(a)')  &
   "<Grid Name=""AllDomain"" GridType=""Collection"" CollectionType=""Spatial"">"
    end if

  end subroutine sll_xdmf_open_serial_blocks

  !>Write 2d array in parallel hdf5 file and the matching line in XDMF file
  subroutine sll_xdmf_array_2d_serial_blocks(mesh_name, &
                               array,array_name,error,&
                               xmf,center)

    character(len=*), intent(in)     :: mesh_name      !< file with mesh coordinates
    sll_real64, intent(in)           :: array(:,:)     !< data array
    character(len=*), intent(in)     :: array_name     !< name of the field
    integer(HID_T)                   :: file_id        !< data file unit number
!    sll_int32                        :: npoints_x1     !< nodes number x
!    sll_int32                        :: npoints_x2     !< nodes number y
    sll_int32, intent(in), optional  :: xmf            !< xml file unit number
    character(len=4), optional       :: center         !< "Node" or "Cell"
    sll_int32, intent(out)           :: error          !< error code
    sll_int32                        :: prank
    sll_int32                        :: psize
    character(len=4)                 :: crank
!    sll_int32                        :: iproc
!    character(len=4)                 :: cproc
!    character(len=4)                 :: cplot

    prank = sll_get_collective_rank(sll_world_collective)
    psize = sll_get_collective_size(sll_world_collective)
    call int2string(prank,crank)

#ifndef NOHDF5
    call sll_hdf5_file_create( &
       trim(mesh_name)//"-"//trim(array_name)//"-"//crank//".h5", &
       file_id,error)
    call sll_hdf5_write_array(file_id,array,"/"//trim(array_name),error)
    call sll_hdf5_file_close(file_id, error)
#endif

!   if (prank == 0) then
!       do iproc = 0, psize-1
!          call int2string(iproc,cproc)
!          write(xmf,'(a)')"<Grid Name=""SubDomain"" GridType=""Uniform"">"
!          write(xmf,'(a,2i6,a)') &
!          "<Topology TopologyType=""2DSMesh"" NumberOfElements='",npoints_x2,npoints_x1,"'/>"
!          write(xmf,'(a)')"<Geometry GeometryType=""X_Y"">"
!          write(xmf,'(a,2i6,a)')"<DataItem Dimensions='",npoints_x2,npoints_x1, &
!          "' NumberType=""Float"" Precision=""4"" Format=""HDF"">"
!          write(xmf,'(a)')mesh_name//cproc//".h5:"//coordnames(1)
!          write(xmf,'(a)')"</DataItem>"
!          write(xmf,'(a,2i6,a)')"<DataItem Dimensions='",npoints_x2,npoints_x1, &
!          "' NumberType=""Float"" Precision=""4"" Format=""HDF"">"
!          write(xmf,'(a)')mesh_name//cproc//".h5:"//coordnames(2)
!          write(xmf,'(a)')"</DataItem>"
!          write(xmf,'(a)')"</Geometry>"
!          write(xmf,'(a)') &
!          "<Attribute Name=""Z"" AttributeType=""Scalar"" Center=""Node"">"
!          write(xmf,'(a,2i6,a)')"<DataItem Dimensions='",npoints_x2,npoints_x1, &
!          "' NumberType=""Float"" Precision=""8"" Format=""HDF"">"
!          write(xmf,'(a)')array_name//cproc//"-"//cplot//".h5:/Z"
!          write(xmf,'(a)')"</DataItem>"
!          write(xmf,'(a)')"</Attribute>"
!          write(xmf,'(a)')"</Grid>"
!       end do
!    end if
   

  end subroutine sll_xdmf_array_2d_serial_blocks

!> Close the XML file and finish to write last lines.
  subroutine sll_xdmf_close_serial_blocks(file_id,error)
  sll_int32, intent(in) :: file_id !< file unit number
  sll_int32, intent(out) :: error  !< error code
  sll_int32 :: prank
  prank = sll_get_collective_rank(sll_world_collective)
  if (prank==0) then
     call sll_xml_file_close(file_id,error)
  end if
  end subroutine sll_xdmf_close_serial_blocks


end module sll_xdmf_serial_blocks
