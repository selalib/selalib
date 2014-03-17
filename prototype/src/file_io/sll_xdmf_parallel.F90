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
module sll_xdmf_parallel
#include "sll_working_precision.h"
#include "sll_assert.h"
  
  use sll_collective
  use sll_hdf5_io_parallel
  use sll_ascii_io
  use sll_xml_io
  use hdf5
  
  implicit none
  
  !>Create the xdmf file
  interface sll_xdmf_open
     module procedure sll_xdmf_open_2d_parallel
     module procedure sll_xdmf_open_3d_parallel
  end interface

  !> Write and array in an xmf file
  interface sll_xdmf_write_array
     module procedure sll_xdmf_array_2d_parallel
     module procedure sll_xdmf_array_3d_parallel
  end interface
  
  !> Close the xdmf file
  interface sll_xdmf_close
     module procedure sll_xdmf_close_parallel
  end interface
contains  
  
  !>Open a XDMF format file for a 2d plot
  subroutine sll_xdmf_open_2d_parallel(rank,file_name,mesh_name,nnodes_x1,nnodes_x2,file_id,error)

    sll_int32, intent(in)        :: rank      !< processor number id
    character(len=*), intent(in) :: file_name !< xmf file name 
    character(len=*), intent(in) :: mesh_name !< file name that contains mesh coordinates
    sll_int32                    :: file_id   !< file unit number
    sll_int32                    :: error     !< error code
    sll_int32                    :: nnodes_x1 !< nodes number x
    sll_int32                    :: nnodes_x2 !< nodes number y
    
    if (rank == 0) then
       call sll_xml_file_create(trim(file_name),file_id,error)
       call sll_xml_grid_geometry(file_id, trim(mesh_name), nnodes_x1, nnodes_x2)
    end if

  end subroutine sll_xdmf_open_2d_parallel

  !>Open a XDMF format file for a 3d plot
  subroutine sll_xdmf_open_3d_parallel( &
    rank,      &
    file_name, &
    mesh_name, &
    nnodes_x1, &
    nnodes_x2, &
    nnodes_x3, &
    file_id,   &
    error)
    
    sll_int32, intent(in)        :: rank       !< processor number id
    character(len=*), intent(in) :: file_name  !< xml file name
    character(len=*), intent(in) :: mesh_name  !< file name that contains data coordinates
    sll_int32                    :: nnodes_x1  !< nodes number x
    sll_int32                    :: nnodes_x2  !< nodes number y
    sll_int32                    :: nnodes_x3  !< nodes number z
    sll_int32, intent(out)       :: file_id    !< file unit number
    sll_int32, intent(out)       :: error      !< error code
    
    if (rank == 0) then
       call sll_xml_file_create(trim(file_name),file_id,error)
       call sll_xml_grid_geometry(file_id, trim(mesh_name),  &
                                  nnodes_x1, nnodes_x2, nnodes_x3)
    end if

  end subroutine sll_xdmf_open_3d_parallel
  
  !>Write 2d array in parallel hdf5 file and the matching line in XDMF file
  subroutine sll_xdmf_array_2d_parallel(mesh_name,global_dims, offset,&
                               array,array_name,error,&
                               xmffile_id,center)

    character(len=*), intent(in)     :: mesh_name      !< file with mesh coordinates
    sll_real64, intent(in)           :: array(:,:)     !< data array
    character(len=*), intent(in)     :: array_name     !< name of the field
    integer(HSSIZE_T)                :: offset(2)      !< offset
    integer(HSIZE_T)                 :: global_dims(2) !< global dimensions
    integer(HID_T)                   :: file_id        !< data file unit number
    sll_int32                        :: npoints_x1     !< nodes number x
    sll_int32                        :: npoints_x2     !< nodes number y
    sll_int32, intent(in), optional  :: xmffile_id     !< xml file unit number
    character(len=4), optional       :: center         !< "Node" or "Cell"
    sll_int32, intent(out)           :: error          !< error code
    sll_int32                        :: myrank

    call sll_hdf5_file_create(trim(mesh_name)//"-"//trim(array_name)//".h5", &
                              file_id,error)
    call sll_hdf5_write_array(file_id,global_dims,offset, &
                              array,"/"//trim(array_name),error)
    call sll_hdf5_file_close(file_id, error)

    myrank = sll_get_collective_rank(sll_world_collective)

    if ( present(xmffile_id) .and. present(center) .and. myrank==0) then
       npoints_x1 = int(global_dims(1),4)
       npoints_x2 = int(global_dims(2),4)

#ifndef NOHDF5
       call sll_xml_field( &
            xmffile_id, &
            trim(array_name), &
            trim(mesh_name)//"-"//trim(array_name)//".h5:/"//trim(array_name), &
            npoints_x1,&
            npoints_x2,&
            'HDF', &
            center)
#else
       call sll_xml_field(xmffile_id,trim(array_name), &
                          trim(mesh_name)//"-"//trim(array_name)//".bin", &
                          npoints_x1,npoints_x2,'Binary',center)
#endif
    end if

  end subroutine sll_xdmf_array_2d_parallel

  !>Write 3d array in binary or hdf5 file and the matching line in XDMF file
  subroutine sll_xdmf_array_3d_parallel(mesh_name,global_dims,offset, &
                               array,array_name,error,xmffile_id,center)

    character(len=*), intent(in)    :: mesh_name      !< file with mesh coordinates
    sll_real64, intent(in)          :: array(:,:,:)   !< data array
    character(len=*), intent(in)    :: array_name     !< name of the field
    integer(HSSIZE_T)               :: offset(3)      !< offset
    integer(HSIZE_T)                :: global_dims(3) !< global dimensions
    integer(HID_T)                  :: file_id        !< data file unit number
    sll_int32                       :: npoints_x1     !< nodes number x
    sll_int32                       :: npoints_x2     !< nodes number y
    sll_int32                       :: npoints_x3     !< nodes number z
    sll_int32, intent(in), optional :: xmffile_id     !< xml file unit number
    character(len=4), optional      :: center         !< "Node" or "Cell"
    sll_int32, intent(out)          :: error          !< error code
    sll_int32                       :: myrank
    
    call sll_hdf5_file_create(trim(mesh_name)//"-"//trim(array_name)//".h5", &
                             file_id,error)
    call sll_hdf5_write_array(file_id,global_dims,offset,array, &
                              "/"//trim(array_name),error)
    call sll_hdf5_file_close(file_id, error)

    myrank = sll_get_collective_rank(sll_world_collective)
    if ( present(xmffile_id) .and. present(center) .and. myrank==0) then
       npoints_x1 = int(global_dims(1),4)
       npoints_x2 = int(global_dims(2),4)
       npoints_x3 = int(global_dims(3),4)

#ifndef NOHDF5
       call sll_xml_field( &
            xmffile_id, &
            trim(array_name), &
            trim(mesh_name)//"-"//trim(array_name)//".h5:/"//trim(array_name), &
            npoints_x1,&
            npoints_x2,&
            npoints_x3,&
            'HDF', &
            center)
#else
       call sll_xml_field(xmffile_id,trim(array_name), &
                          trim(mesh_name)//"-"//trim(array_name)//".bin", &
                          npoints_x1,npoints_x2,npoints_x3,'Binary',center)
#endif
    end if

  end subroutine sll_xdmf_array_3d_parallel

!> Close the XML file and finish to write last lines.
  subroutine sll_xdmf_close_parallel(file_id,error)
  sll_int32, intent(in) :: file_id !< file unit number
  sll_int32, intent(out) :: error  !< error code
  sll_int32 :: myrank
  myrank = sll_get_collective_rank(sll_world_collective)
  if (myrank==0) then
     call sll_xml_file_close(file_id,error)
  end if
  end subroutine sll_xdmf_close_parallel


end module sll_xdmf_parallel
