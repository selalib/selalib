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
    sll_p_io_gnuplot

use sll_m_utilities, only: &
    sll_s_int2string, &
    sll_s_new_file_id

use sll_m_xdmf, only: &
    sll_s_xdmf_corect2d_nodes

use sll_m_cartesian_meshes, only: &
    sll_t_cartesian_mesh_2d

use sll_m_pic_visu, only: &
    sll_o_particles_center_gnuplot

implicit none

type :: sll_t_pic_viewer_2d

  type(sll_t_cartesian_mesh_2d), pointer :: mesh
  sll_int32                              :: output_format
  character(len=:), allocatable          :: label

  contains
  
  procedure :: set_format
  procedure :: write_field
  procedure :: write_particles
  !procedure :: write_field_and_particles

end type

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

  subroutine set_format( viewer, output_format )

    class(sll_t_pic_viewer_2d) :: viewer
    sll_int32                  :: output_format

    viewer%output_format = output_format

  end subroutine set_format

  subroutine write_field( viewer, field, iplot )

    class(sll_t_pic_viewer_2d)            :: viewer
    sll_real64,                intent(in) :: field(:,:)
    sll_int32,                 intent(in) :: iplot

    sll_real64                            :: xmin
    sll_real64                            :: ymin
    sll_real64                            :: xmax
    sll_real64                            :: ymax
    sll_int32                             :: nx
    sll_int32                             :: ny
    sll_int32                             :: error

    xmin = viewer%mesh%eta1_min
    xmax = viewer%mesh%eta1_max
    ymin = viewer%mesh%eta2_min
    ymax = viewer%mesh%eta2_max

    nx = viewer%mesh%num_cells1
    ny = viewer%mesh%num_cells2

    if (viewer%output_format == SLL_P_IO_GNUPLOT) then
 
      call sll_o_gnuplot_2d(xmin, xmax, nx,  &
                            ymin, ymax, ny,  &
                            field, viewer%label, &
                            iplot, error)
    end if

  end subroutine write_field

  subroutine write_particles( viewer, xp, yp, iplot, time )

    class(sll_t_pic_viewer_2d)            :: viewer
    sll_real64,                intent(in) :: xp(:)
    sll_real64,                intent(in) :: yp(:)
    sll_int32,                 intent(in) :: iplot
    sll_real64,                optional   :: time

    sll_real64                            :: xmin
    sll_real64                            :: ymin
    sll_real64                            :: xmax
    sll_real64                            :: ymax
    sll_int32                             :: error

    xmin = viewer%mesh%eta1_min
    xmax = viewer%mesh%eta1_max
    ymin = viewer%mesh%eta2_min
    ymax = viewer%mesh%eta2_max

    if (viewer%output_format == SLL_P_IO_GNUPLOT) then
 
      call sll_o_particles_center_gnuplot( viewer%label, &
           xp, yp, xmin, xmax, ymin, ymax, iplot, time )

    end if

  end subroutine write_particles

end module sll_m_pic_viewer

