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

use sll_m_utilities, only: &
    sll_s_int2string, &
    sll_s_new_file_id

use sll_m_xdmf, only: &
    sll_s_xdmf_corect2d_nodes

use sll_m_cartesian_meshes, only: &
    sll_t_cartesian_mesh_2d

implicit none

type :: sll_t_pic_viewer_2d
  type(sll_t_cartesian_mesh_2d), pointer :: mesh
contains

end type

contains

  function sll_f_new_pic_viewer_2d( mesh ) result(viewer)

    type(sll_t_cartesian_mesh_2d), pointer :: mesh
    type(sll_t_pic_viewer_2d),     pointer :: viewer

    sll_int32 :: ierr

    SLL_ALLOCATE(viewer, ierr)

    call initialize_pic_viewer_2d( &
         viewer,                   &
         mesh )

  end function sll_f_new_pic_viewer_2d

  subroutine initialize_pic_viewer_2d( &
    viewer,                                &
    mesh ) 

    type(sll_t_pic_viewer_2d),     pointer :: viewer
    type(sll_t_cartesian_mesh_2d), pointer :: mesh

  end subroutine initialize_pic_viewer_2d

end module sll_m_pic_viewer
