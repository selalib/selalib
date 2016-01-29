!**************************************************************
!  Copyright INRIA
!  Authors :
!     CALVI project team
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

!> @ingroup particle_methods
!> \brief
!> This module provides some routines for plotting fields and particles
!> during PIC simulations.

module sll_m_pic_viewer

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_assert.h"
#include "sll_memory.h"
#include "sll_working_precision.h"

use sll_m_ascii_io, only: &
    sll_s_ascii_file_create

use sll_m_gnuplot, only: &
    sll_o_gnuplot_2d

use sll_m_coordinate_transformation_2d_base, only: &
    sll_p_io_gnuplot, &
    sll_p_io_xdmf

use sll_m_utilities, only: &
    sll_s_int2string, &
    sll_s_new_file_id

use sll_m_xdmf, only: &
    sll_s_xdmf_corect2d_nodes

use sll_m_cartesian_meshes, only: &
    sll_t_cartesian_mesh_2d

use sll_m_pic_visu, only: &
    sll_o_particles_center_gnuplot

#ifndef NOHDF5
  use hdf5, only: &
    hid_t

  use sll_m_hdf5_io_serial, only: &
    sll_o_hdf5_file_close, &
    sll_o_hdf5_file_create, &
    sll_o_hdf5_write_array

#endif


implicit none

public sll_f_new_pic_viewer_2d, &
       sll_s_pic_viewer_set_format, &
       sll_o_pic_viewer_write

private

type, public :: sll_t_pic_viewer_2d

  type(sll_t_cartesian_mesh_2d), pointer :: mesh
  sll_int32                              :: output_format
  character(len=:), allocatable          :: label

end type

interface sll_o_pic_viewer_write

  module procedure write_2d_field
  module procedure write_2d_particles
  module procedure write_2d_field_and_particles

end interface sll_o_pic_viewer_write



contains

  function sll_f_new_pic_viewer_2d( mesh, &
                                    label &
    ) result(viewer)

    type(sll_t_cartesian_mesh_2d), pointer :: mesh
    type(sll_t_pic_viewer_2d),     pointer :: viewer
    character(len=*)                       :: label

    sll_int32 :: ierr

    SLL_ALLOCATE(viewer, ierr)

    viewer%label = label
    call initialize_pic_viewer_2d( &
         viewer,                   &
         mesh )

  end function sll_f_new_pic_viewer_2d

  subroutine initialize_pic_viewer_2d( &
    viewer,                                &
    mesh ) 

    type(sll_t_pic_viewer_2d),     pointer :: viewer
    type(sll_t_cartesian_mesh_2d), pointer :: mesh

    viewer%mesh => mesh

  end subroutine initialize_pic_viewer_2d

  subroutine sll_s_pic_viewer_set_format( viewer, output_format )

    type(sll_t_pic_viewer_2d) :: viewer
    sll_int32                  :: output_format

    viewer%output_format = output_format

  end subroutine sll_s_pic_viewer_set_format

  subroutine write_2d_field( viewer, field, iplot )

    type(sll_t_pic_viewer_2d)            :: viewer
    sll_real64,                intent(in) :: field(:,:)
    sll_int32,                 intent(in) :: iplot

    sll_real64                            :: xmin
    sll_real64                            :: ymin
    sll_real64                            :: xmax
    sll_real64                            :: ymax
    sll_real64                            :: dx
    sll_real64                            :: dy
    sll_int32                             :: nx
    sll_int32                             :: ny
    sll_int32                             :: error
    sll_int32                             :: file_id

#ifndef NOHDF5
    integer(hid_t) :: hfile_id  
#endif

    xmin = viewer%mesh%eta1_min
    xmax = viewer%mesh%eta1_max
    ymin = viewer%mesh%eta2_min
    ymax = viewer%mesh%eta2_max
    dx   = viewer%mesh%delta_eta1
    dy   = viewer%mesh%delta_eta2

    nx = viewer%mesh%num_cells1
    ny = viewer%mesh%num_cells2

    select case (viewer%output_format)

    case(SLL_P_IO_GNUPLOT)
 
      call sll_o_gnuplot_2d(xmin, xmax, nx,  &
                            ymin, ymax, ny,  &
                            field, viewer%label, &
                            iplot, error)

    case(SLL_P_IO_XDMF)

      open(newunit = file_id, &
           file    = viewer%label//".xmf" )
      write(file_id,"(a)")"<?xml version=""1.0"" ?> <!DOCTYPE Xdmf SYSTEM ""Xdmf.dtd"" []>"
      write(file_id,"(a)")"<Xdmf xmlns:xi=""http://www.w3.org/2003/XInclude"" Version=""2.2"">"
      write(file_id,"(a)")"<Domain>"
      write(file_id,"(a)")"<Grid Name='mesh' GridType='Uniform'>"
      write(file_id,"(a,2i5,a)")"<Topology TopologyType='2DCoRectMesh' NumberOfElements='",ny,nx,"'/>"
      write(file_id,"(a)")"<Geometry GeometryType='ORIGIN_DXDY'>"
      write(file_id,"(a)")"<DataItem Dimensions='2' NumberType='Float' Format='XML'>"
      write(file_id,"(2f10.5)") xmin, ymin
      write(file_id,"(a)")"</DataItem>"
      write(file_id,"(a)")"<DataItem Dimensions='2' NumberType='Float' Format='XML'>"
      write(file_id,"(2f10.5)") dx, dy
      write(file_id,"(a)")"</DataItem>"
      write(file_id,"(a)")"</Geometry>"
      write(file_id,"(a)")"<Attribute Name='field' AttributeType='Scalar' Center='Node'>"
      write(file_id,"(a,2i6,a)")"<DataItem Dimensions='", ny, nx,"' Format='HDF'>"
      write(file_id,"(a)")viewer%label//".h5:/field"
      write(file_id,"(a)")"</DataItem>"
      write(file_id,"(a)")"</Attribute>"
      write(file_id,"(a)")"</Grid>"
      write(file_id,"(a)")"</Domain>"
      write(file_id,"(a)")"</Xdmf>"
      close(file_id)

      call sll_o_hdf5_file_create(viewer%label//".h5",hfile_id,error)
      call sll_o_hdf5_write_array(hfile_id,field,"/field",error)
      call sll_o_hdf5_file_close(hfile_id, error)

    end select

  end subroutine write_2d_field

  subroutine write_2d_particles( viewer, xp, yp, iplot, time )

    type(sll_t_pic_viewer_2d)            :: viewer
    sll_real64,                intent(in) :: xp(:)
    sll_real64,                intent(in) :: yp(:)
    sll_int32,                 intent(in) :: iplot
    sll_real64,                optional   :: time

    sll_real64                            :: xmin
    sll_real64                            :: ymin
    sll_real64                            :: xmax
    sll_real64                            :: ymax
    sll_int32                             :: error
    sll_int32                             :: file_id

#ifndef NOHDF5
    integer(hid_t) :: hfile_id  
#endif

    xmin = viewer%mesh%eta1_min
    xmax = viewer%mesh%eta1_max
    ymin = viewer%mesh%eta2_min
    ymax = viewer%mesh%eta2_max

    select case (viewer%output_format)

    case(SLL_P_IO_GNUPLOT)
 
      call sll_o_particles_center_gnuplot( viewer%label, &
           xp, yp, xmin, xmax, ymin, ymax, iplot, time )

    case(SLL_P_IO_XDMF)

      open(newunit = file_id, &
           file    = viewer%label//".xmf" )
      write(file_id,"(a)")"<?xml version=""1.0"" ?> <!DOCTYPE Xdmf SYSTEM ""Xdmf.dtd"" []>"
      write(file_id,"(a)")"<Xdmf xmlns:xi=""http://www.w3.org/2003/XInclude"" Version=""2.2"">"
      write(file_id,"(a)")"<Domain>"
      write(file_id,"(a)")"<Grid Name=""Particles"" Type=""Uniform"">"
      write(file_id,"(a,i6,a)")"<Topology TopologyType=""Polyvertex"" NumberOfElements=""",size(xp),"""/>"
      write(file_id,"(a)")"<Geometry Type=""X_Y"">"
      write(file_id,"(a,i6,a)")"<DataItem Format=""HDF"" Dimensions=""",size(xp),""">"
      write(file_id,"(a)")viewer%label//".h5:/xp"
      write(file_id,"(a)")"</DataItem>"
      write(file_id,"(a,i6,a)")"<DataItem Format=""HDF"" Dimensions=""",size(yp),""">"
      write(file_id,"(a)")viewer%label//".h5:/yp"
      write(file_id,"(a)")"</DataItem>"
      write(file_id,"(a)")"</Geometry>"
      write(file_id,"(a)")"</Grid>"
      write(file_id,"(a)")"</Domain>"
      write(file_id,"(a)")"</Xdmf>"
      close(file_id)

      call sll_o_hdf5_file_create(viewer%label//".h5",hfile_id,error)
      call sll_o_hdf5_write_array(hfile_id,xp,"/xp",error)
      call sll_o_hdf5_write_array(hfile_id,yp,"/yp",error)
      call sll_o_hdf5_file_close(hfile_id, error)

    end select

  end subroutine write_2d_particles

  subroutine write_2d_field_and_particles( viewer, xp, yp, op, fp, iplot, time )

    type(sll_t_pic_viewer_2d)            :: viewer
    sll_real64,                intent(in) :: xp(:)
    sll_real64,                intent(in) :: yp(:)
    sll_real64,                intent(in) :: op(:)
    sll_real64,                intent(in) :: fp(:,:)
    sll_int32,                 intent(in) :: iplot
    sll_real64,                optional   :: time

    sll_real64                            :: xmin
    sll_real64                            :: ymin
    sll_real64                            :: xmax
    sll_real64                            :: ymax
    sll_real64                            :: x, dx
    sll_real64                            :: y, dy
    sll_int32                             :: i, j, k, nx, ny
    sll_int32                             :: error
    sll_int32                             :: file_id
    character(len=4)                      :: cplot

#ifndef NOHDF5
    integer(hid_t) :: hfile_id  
#endif

    xmin = viewer%mesh%eta1_min
    xmax = viewer%mesh%eta1_max
    ymin = viewer%mesh%eta2_min
    ymax = viewer%mesh%eta2_max
    dx   = viewer%mesh%delta_eta1
    dy   = viewer%mesh%delta_eta2

    nx = viewer%mesh%num_cells1
    ny = viewer%mesh%num_cells2

    call sll_s_int2string( iplot, cplot )

    select case (viewer%output_format)

    case(SLL_P_IO_GNUPLOT)
 
      SLL_ASSERT(size(xp) == size(yp))
      SLL_ASSERT(size(op) == size(yp))

      open(newunit = file_id, &
           file    = viewer%label//"_particles_"//cplot//'.dat' )
      do k=1, size(xp)
         write(file_id,*) sngl(xp(k)),sngl(yp(k)),sngl(op(k))
      end do
      close(file_id)

      open( file    = viewer%label//'_field_'//cplot//'.dat', &
            status  = 'replace',   &
            form    = 'formatted', &
            newunit = file_id,     &
            iostat  = error )

      x = xmin
      do i=1,nx
         y = ymin
         do j=1,ny
            write(file_id,*) sngl(x),sngl(y),sngl(fp(i,j))
            y = y + dy
         end do
         x = x + dx
         write(file_id,*)
      enddo
      close(file_id)

      open( file     = viewer%label//'.gnu', &
            position = 'append',    &
            newunit  = file_id,     &
            iostat   = error )
      if (iplot == 1) then
       rewind(file_id)
      end if

      if ( present(time)) then
        write(file_id,"(A18,G10.3,A1)")"set title 'Time = ",time,"'"
      end if
      write(file_id,"(a)") &
      & "splot '"//viewer%label//"_field_"//cplot//".dat' w l &
      &  , '"//viewer%label//"_particles_"//cplot//".dat' w p "
      close(file_id)

    case(SLL_P_IO_XDMF)

      open(newunit = file_id, &
           file    = viewer%label//".xmf" )
      write(file_id,"(a)")"<?xml version=""1.0"" ?> <!DOCTYPE Xdmf SYSTEM ""Xdmf.dtd"" []>"
      write(file_id,"(a)")"<Xdmf xmlns:xi=""http://www.w3.org/2003/XInclude"" Version=""2.2"">"
      write(file_id,"(a)")"<Domain>"
      write(file_id,"(a)")"<Grid Name=""Points"" Type=""Uniform"">"
      write(file_id,"(a,i6,a)")"<Topology TopologyType=""Polyvertex"" NumberOfElements=""",size(xp),"""/>"
      write(file_id,"(a)")"<Geometry Type=""X_Y"">"
      write(file_id,"(a,i6,a)")"<DataItem Format=""HDF"" Dimensions=""",size(xp),""">"
      write(file_id,"(a)")viewer%label//".h5:/xp"
      write(file_id,"(a)")"</DataItem>"
      write(file_id,"(a,i6,a)")"<DataItem Format=""HDF"" Dimensions=""",size(yp),""">"
      write(file_id,"(a)")viewer%label//".h5:/yp"
      write(file_id,"(a)")"</DataItem>"
      write(file_id,"(a)")"</Geometry>"
      write(file_id,"(a)")"<Attribute Name='omega' AttributeType='Scalar' Center='Node'>"
      write(file_id,"(a,i6,a)")"<DataItem Format=""HDF"" Dimensions=""",size(op),""">"
      write(file_id,"(a)")viewer%label//".h5:/op"
      write(file_id,"(a)")"</DataItem>"
      write(file_id,"(a)")"</Attribute>"
      write(file_id,"(a)")"</Grid>"
      write(file_id,"(a)")"<Grid Name='mesh' GridType='Uniform'>"
      write(file_id,"(a,2i5,a)")"<Topology TopologyType='2DCoRectMesh' NumberOfElements='",ny,nx,"'/>"
      write(file_id,"(a)")"<Geometry GeometryType='ORIGIN_DXDY'>"
      write(file_id,"(a)")"<DataItem Dimensions='2' NumberType='Float' Format='XML'>"
      write(file_id,"(2f10.5)") xmin, ymin
      write(file_id,"(a)")"</DataItem>"
      write(file_id,"(a)")"<DataItem Dimensions='2' NumberType='Float' Format='XML'>"
      write(file_id,"(2f10.5)") dx, dy
      write(file_id,"(a)")"</DataItem>"
      write(file_id,"(a)")"</Geometry>"
      write(file_id,"(a)")"<Attribute Name='omega' AttributeType='Scalar' Center='Node'>"
      write(file_id,"(a,2i6,a)")"<DataItem Dimensions='", ny, nx,"' Format='HDF'>"
      write(file_id,"(a)")viewer%label//".h5:/fp"
      write(file_id,"(a)")"</DataItem>"
      write(file_id,"(a)")"</Attribute>"
      write(file_id,"(a)")"</Grid>"
      write(file_id,"(a)")"</Domain>"
      write(file_id,"(a)")"</Xdmf>"
      close(file_id)

      call sll_o_hdf5_file_create(viewer%label//".h5",hfile_id,error)
      call sll_o_hdf5_write_array(hfile_id,xp,"/xp",error)
      call sll_o_hdf5_write_array(hfile_id,yp,"/yp",error)
      call sll_o_hdf5_write_array(hfile_id,op,"/op",error)
      call sll_o_hdf5_write_array(hfile_id,fp,"/fp",error)
      call sll_o_hdf5_file_close(hfile_id, error)

    end select

  end subroutine write_2d_field_and_particles

end module sll_m_pic_viewer

