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

!> @ingroup file_io_parallel
!> @brief
!> Implements the functions to write xdmf file plotable by VisIt
!> @details
!> In <b> XDMF </b> (eXtensible Data Model and Format) the description of the 
!> data is separate from the values themselves. Light data is stored using XML, 
!> Heavy data is stored using Parallel HDF5. These files are readable by 
!> Paraview.
module sll_m_xdmf_serial_blocks
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_assert.h"
#include "sll_working_precision.h"

  use sll_m_utilities, only: &
    sll_s_int2string

  use sll_m_xml_io, only: &
    sll_s_xml_file_close, &
    sll_s_xml_file_create

  use sll_mpi, only: &
    mpi_comm_rank, &
    mpi_comm_size, &
    mpi_comm_world, &
    mpi_integer, &
    mpi_recv, &
    mpi_send, &
    mpi_status_size

#ifndef NOHDF5
  use hdf5, only: &
    hid_t

  use sll_m_hdf5_io_serial, only: &
    sll_s_hdf5_ser_file_create, &
    sll_s_hdf5_ser_file_close, &
    sll_o_hdf5_ser_write_array

#endif
  implicit none

  public :: &
    sll_s_xdmf_array_2d_serial_blocks, &
    sll_s_xdmf_close_serial_blocks, &
    sll_s_xdmf_open_serial_blocks

  private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
  !> Create a new xmdf file to plot parallel array using hdf5 serial library
  interface sll_xdmf_open
     module procedure sll_s_xdmf_open_serial_blocks
  end interface

  !> Write and array in an xmf file
  interface sll_xdmf_write_array
     module procedure sll_s_xdmf_array_2d_serial_blocks
  end interface
  
  !> Close the xdmf file
  interface sll_xdmf_close
     module procedure sll_s_xdmf_close_serial_blocks
  end interface

contains  
  
  !>Open a XDMF format file for a 2d plot
  subroutine sll_s_xdmf_open_serial_blocks( prefix,  &
                                          x1,      &
                                          x2,      &
                                          xmf,     &
                                          error)

    character(len=*), intent(in) :: prefix  !< xmf file name 
    sll_real64, intent(in)       :: x1(:,:) !< x1 coordinates
    sll_real64, intent(in)       :: x2(:,:) !< x2 coordinates
    sll_int32, intent(out)       :: xmf     !< unit file number
    sll_int32                    :: error   !< error code
    sll_int32                    :: prank   !< processor number id
    character(len=4)             :: crank
#ifndef NOHDF5
    integer(HID_T)               :: h5file_id
#endif
    
    call MPI_COMM_RANK(MPI_COMM_WORLD,prank,error)
    call sll_s_int2string(prank,crank)
    if(prank == 0) then
     call sll_s_xml_file_create(trim(prefix)//".xmf",xmf,error)
    end if
#ifndef NOHDF5
    call sll_s_hdf5_ser_file_create(trim(prefix)//"-mesh-"//crank//".h5",h5file_id,error)
    call sll_o_hdf5_ser_write_array( h5file_id, x1, "/x1", error )
    call sll_o_hdf5_ser_write_array( h5file_id, x2, "/x2", error )
    call sll_s_hdf5_ser_file_close( h5file_id, error )
#endif

  end subroutine sll_s_xdmf_open_serial_blocks

  !>Write 2d array in parallel hdf5 file and the matching line in XDMF file
  subroutine sll_s_xdmf_array_2d_serial_blocks(xmf, prefix, &
                               array,array_name,error,&
                               center)

   sll_int32, intent(in)           :: xmf        !< unit file number
   character(len=*), intent(in)    :: prefix     !< file with mesh coordinates
   sll_real64, intent(in)          :: array(:,:) !< data array
   character(len=*), intent(in)    :: array_name !< name of the field
#ifndef NOHDF5
   integer(HID_T)                  :: file_id    !< data file unit number
#endif
   sll_int32                       :: npts_x1    !< nodes number x
   sll_int32                       :: npts_x2    !< nodes number y
   sll_int32, intent(out)          :: error      !< error code
   character(len=4), optional      :: center     !< "Node" or "Cell"
   sll_int32                       :: prank
   sll_int32                       :: psize
   character(len=4)                :: crank
   sll_int32                       :: iproc
   character(len=4)                :: cproc
   sll_int32, parameter            :: tag=1111
   sll_int32                       :: statut(MPI_STATUS_SIZE)

   call MPI_COMM_RANK(MPI_COMM_WORLD,prank,error)
   call MPI_COMM_SIZE(MPI_COMM_WORLD,psize,error)
   call sll_s_int2string(prank,crank)

#ifndef NOHDF5
   call sll_s_hdf5_ser_file_create( trim(prefix)//"-"//crank//".h5", file_id, error )
   call sll_o_hdf5_ser_write_array( file_id, array, "/"//trim(array_name), error )
   call sll_s_hdf5_ser_file_close( file_id, error )
#endif

   npts_x1 = size(array,1)
   npts_x2 = size(array,2)

   if (prank == 0) then
     write(xmf,'(a)')  &
   "<Grid Name=""AllDomain"" GridType=""Collection"" CollectionType=""Spatial"">"
     do iproc = 0, psize-1
       call sll_s_int2string(iproc,cproc)
       if (iproc > 0) then
          call MPI_RECV(npts_x1,1,MPI_INTEGER,iproc,tag,MPI_COMM_WORLD,statut,error)
          call MPI_RECV(npts_x2,1,MPI_INTEGER,iproc,tag,MPI_COMM_WORLD,statut,error)
       end if

       write(xmf,'(a)')"<Grid Name=""SubDomain"" GridType=""Uniform"">"
       write(xmf,'(a,2i6,a)') &
         "<Topology TopologyType=""2DSMesh"" NumberOfElements='",npts_x2,npts_x1,"'/>"
       write(xmf,'(a)')"<Geometry GeometryType=""X_Y"">"
       write(xmf,'(a,2i6,a)')"<DataItem Dimensions='",npts_x2,npts_x1, &
         "' NumberType=""Float"" Precision=""4"" Format=""HDF"">"
       write(xmf,'(a)')trim(prefix)//"-mesh-"//cproc//".h5:/x1"
       write(xmf,'(a)')"</DataItem>"
       write(xmf,'(a,2i6,a)')"<DataItem Dimensions='",npts_x2,npts_x1, &
         "' NumberType=""Float"" Precision=""4"" Format=""HDF"">"
       write(xmf,'(a)')trim(prefix)//"-mesh-"//cproc//".h5:/x2"
       write(xmf,'(a)')"</DataItem>"
       write(xmf,'(a)')"</Geometry>"
       write(xmf,'(a)') &
         "<Attribute Name='"//trim(array_name)//"' AttributeType=""Scalar"" Center=""Node"">"
       write(xmf,'(a,2i6,a)')"<DataItem Dimensions='",npts_x2,npts_x1, &
         "' NumberType=""Float"" Precision=""8"" Format=""HDF"">"
       write(xmf,'(a)') trim(prefix)//"-"//cproc//".h5:/"//trim(array_name)
       write(xmf,'(a)')"</DataItem>"
       write(xmf,'(a)')"</Attribute>"
       write(xmf,'(a)')"</Grid>"
    end do
  else
    call MPI_SEND(npts_x1, 1, MPI_INTEGER, 0, tag, MPI_COMM_WORLD, error)
    call MPI_SEND(npts_x2, 1, MPI_INTEGER, 0, tag, MPI_COMM_WORLD, error)
  end if

  return
  SLL_ASSERT(present(center))

  end subroutine sll_s_xdmf_array_2d_serial_blocks

!> Close the XML file and finish to write last lines.
  subroutine sll_s_xdmf_close_serial_blocks(file_id,error)
  sll_int32, intent(in) :: file_id !< file unit number
  sll_int32, intent(out) :: error  !< error code
  sll_int32 :: prank
  logical   :: i_opened, i_exist
  character(len=255) :: i_name

  call MPI_COMM_RANK(MPI_COMM_WORLD,prank,error)
  if (prank==0) then
     INQUIRE (file_id, OPENED=I_OPENED, NAME=I_NAME, EXIST=I_EXIST) 
     call sll_s_xml_file_close(file_id,error)
  end if
  end subroutine sll_s_xdmf_close_serial_blocks


end module sll_m_xdmf_serial_blocks
