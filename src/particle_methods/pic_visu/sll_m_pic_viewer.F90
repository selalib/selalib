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

use hdf5, only: hid_t

use sll_m_hdf5_io_serial, only: &
  sll_s_hdf5_ser_file_create, &
  sll_s_hdf5_ser_file_close, &
  sll_o_hdf5_ser_write_array

#endif

use sll_m_xml_io


implicit none

public sll_f_new_pic_viewer_2d, &
       sll_o_pic_viewer_write

private

type, public :: sll_t_pic_viewer_2d

  type(sll_t_cartesian_mesh_2d), pointer :: mesh
  character(len=72)                      :: label

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
    viewer,                            &
    mesh ) 

    type(sll_t_pic_viewer_2d),     pointer :: viewer
    type(sll_t_cartesian_mesh_2d), pointer :: mesh

    viewer%mesh => mesh

  end subroutine initialize_pic_viewer_2d

  subroutine write_2d_field( viewer, field, iplot )

    type(sll_t_pic_viewer_2d)             :: viewer
    sll_real64,                intent(in) :: field(:,:)
    sll_int32,                 intent(in) :: iplot

    sll_int32                             :: nx
    sll_int32                             :: ny
    sll_int32                             :: error
    sll_int32                             :: file_id
    character(len=4)                      :: cplot

#ifndef NOHDF5
    integer(hid_t) :: hfile_id  
#endif

    call sll_s_int2string(iplot, cplot)
    nx = viewer%mesh%num_cells1
    ny = viewer%mesh%num_cells2

    !Write the light data
    call sll_s_xml_file_create( trim(viewer%label)//cplot//".xmf", file_id, error)
    call sll_o_xml_field( file_id, 'omega', trim(viewer%label)//cplot//".h5:/fp", &
                          nx, ny, "HDF", "Node")
    call sll_s_xml_file_close( file_id, error)

    !Write the heavy data
#ifndef NOHDF5
    call sll_s_hdf5_ser_file_create( trim(viewer%label)//cplot//".h5", hfile_id, error )
    call sll_o_hdf5_ser_write_array( hfile_id, field, "/fp", error )
    call sll_s_hdf5_ser_file_close( hfile_id, error )
#endif

  end subroutine write_2d_field

  subroutine write_2d_particles( viewer, xp, yp, op, iplot, time )

    type(sll_t_pic_viewer_2d)             :: viewer
    sll_real64,                intent(in) :: xp(:)
    sll_real64,                intent(in) :: yp(:)
    sll_real64,                intent(in) :: op(:)
    sll_int32,                 intent(in) :: iplot
    sll_real64,                optional   :: time

    sll_real64                            :: xmin
    sll_real64                            :: ymin
    sll_real64                            :: xmax
    sll_real64                            :: ymax
    sll_int32                             :: error
    sll_int32                             :: file_id
    character(len=4)                      :: cplot

#ifndef NOHDF5
    integer(hid_t) :: hfile_id  
#endif

    call sll_s_int2string(iplot, cplot)

    xmin = viewer%mesh%eta1_min
    xmax = viewer%mesh%eta1_max
    ymin = viewer%mesh%eta2_min
    ymax = viewer%mesh%eta2_max

    call sll_s_xml_file_create( trim(viewer%label)//cplot//".xmf", file_id, error)
    call write_particles_2d_xml( file_id, trim(viewer%label)//cplot, xp, yp )
    call sll_o_xml_field( file_id, 'omega', trim(viewer%label)//cplot//".h5:/op", &
                          size(op), "HDF", "Node")
    call sll_s_xml_file_close( file_id, error)

#ifndef NOHDF5
    call sll_s_hdf5_ser_file_create( trim(viewer%label)//cplot//".h5", hfile_id, error )
    call sll_o_hdf5_ser_write_array( hfile_id, xp, "/xp", error )
    call sll_o_hdf5_ser_write_array( hfile_id, yp, "/yp", error )
    call sll_o_hdf5_ser_write_array( hfile_id, op, "/op", error )
    call sll_s_hdf5_ser_file_close( hfile_id, error )
#endif

  end subroutine write_2d_particles

  subroutine write_2d_field_and_particles( viewer, xp, yp, op, fp, iplot, time )

    type(sll_t_pic_viewer_2d)             :: viewer
    sll_real64,                intent(in) :: xp(:)
    sll_real64,                intent(in) :: yp(:)
    sll_real64,                intent(in) :: op(:)
    sll_real64,                intent(in) :: fp(:,:)
    sll_int32,                 intent(in) :: iplot
    sll_real64,                optional   :: time

    sll_int32                             :: error
    sll_int32                             :: nx
    sll_int32                             :: ny
    sll_int32                             :: file_id
    character(len=4)                      :: cplot

#ifndef NOHDF5
    integer(hid_t) :: hfile_id  
#endif

    nx = viewer%mesh%num_cells1
    ny = viewer%mesh%num_cells2

    call sll_s_int2string( iplot, cplot )

    call sll_s_xml_file_create( trim(viewer%label)//cplot//".xmf", file_id, error)
    call write_particles_2d_xml( file_id, trim(viewer%label)//cplot, xp, yp )
    call sll_o_xml_field( file_id, 'omega', trim(viewer%label)//cplot//".h5:/op", &
                          size(op), "HDF", "Node")
    write(file_id,"(a)")"</Grid>"
    call write_grid_2d_xml( file_id, viewer )
    call sll_o_xml_field( file_id, 'omega', trim(viewer%label)//cplot//".h5:/fp", &
                          nx, ny, "HDF", "Node")
    call sll_s_xml_file_close( file_id, error)

#ifndef NOHDF5
    call sll_s_hdf5_ser_file_create( trim(viewer%label)//cplot//".h5", hfile_id, error )
    call sll_o_hdf5_ser_write_array( hfile_id, xp, "/xp", error )
    call sll_o_hdf5_ser_write_array( hfile_id, yp, "/yp", error )
    call sll_o_hdf5_ser_write_array( hfile_id, op, "/op", error )
    call sll_o_hdf5_ser_write_array( hfile_id, fp, "/fp", error )
    call sll_s_hdf5_ser_file_close( hfile_id, error )
#endif

  end subroutine write_2d_field_and_particles

  subroutine write_grid_2d_xml( file_id, viewer )

    type(sll_t_pic_viewer_2d)            :: viewer
    sll_int32,                intent(in) :: file_id
    sll_real64                           :: xmin
    sll_real64                           :: ymin
    sll_real64                           :: dx
    sll_real64                           :: dy
    sll_int32                            :: nx
    sll_int32                            :: ny

    xmin = viewer%mesh%eta1_min
    ymin = viewer%mesh%eta2_min
    dx   = viewer%mesh%delta_eta1
    dy   = viewer%mesh%delta_eta2

    nx = viewer%mesh%num_cells1
    ny = viewer%mesh%num_cells2

    write(file_id,"(a)")"<Grid Name='Grid' GridType='Uniform'>"
    write(file_id,"(a,2i5,a)")"<Topology TopologyType='2DCoRectMesh' NumberOfElements='",ny,nx,"'/>"
    write(file_id,"(a)")"<Geometry GeometryType='ORIGIN_DXDY'>"
    write(file_id,"(a)")"<DataItem Dimensions='2' NumberType='Float' Format='XML'>"
    write(file_id,"(2f10.5)") xmin, ymin
    write(file_id,"(a)")"</DataItem>"
    write(file_id,"(a)")"<DataItem Dimensions='2' NumberType='Float' Format='XML'>"
    write(file_id,"(2f10.5)") dx, dy
    write(file_id,"(a)")"</DataItem>"
    write(file_id,"(a)")"</Geometry>"

  end subroutine write_grid_2d_xml

  subroutine write_particles_2d_xml( file_id, prefix, xp, yp )

    sll_int32,                intent(in) :: file_id
    character(len=*),         intent(in) :: prefix
    sll_real64,               intent(in) :: xp(:)
    sll_real64,               intent(in) :: yp(:)

    sll_int32                            :: n

    n = size(xp)
    SLL_ASSERT( size(yp) == n)

    write(file_id,"(a)")"<Grid Name=""Particles"" Type=""Uniform"">"
    write(file_id,"(a,i6,a)")"<Topology TopologyType=""Polyvertex"" NumberOfElements=""",n,"""/>"
    write(file_id,"(a)")"<Geometry Type=""X_Y"">"
    write(file_id,"(a,i6,a)")"<DataItem Format=""HDF"" Dimensions=""",n,""">"
    write(file_id,"(a)")prefix//".h5:/xp"
    write(file_id,"(a)")"</DataItem>"
    write(file_id,"(a,i6,a)")"<DataItem Format=""HDF"" Dimensions=""",n,""">"
    write(file_id,"(a)")prefix//".h5:/yp"
    write(file_id,"(a)")"</DataItem>"
    write(file_id,"(a)")"</Geometry>"

  end subroutine write_particles_2d_xml

end module sll_m_pic_viewer

