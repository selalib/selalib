!**************************************************************
!  Copyright INRIA
!  Authors : MCP,ALH
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


module sll_m_pic_visualizer_interface

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_assert.h"
#include "sll_errors.h"
! #include "sll_memory.h"
#include "sll_working_precision.h"

  use sll_m_particle_group_base, only: &
    sll_c_particle_group_base

  use sll_m_pic_lbfr_4d_group, only: &
    sll_t_pic_lbfr_4d_group

  implicit none

  public :: &
    sll_t_pic_visualizer_interface

  private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


  type sll_t_pic_visualizer_interface

      !   no member fields for the moment

    contains

      procedure :: visualize_particle_group

  end type sll_t_pic_visualizer_interface


contains


  subroutine visualize_particle_group( &
          self, &
          particle_group,   &
          array_name,   &
          plot_np_x,  &
          plot_np_y,  &
          plot_np_vx,   &
          plot_np_vy,   &
          iplot )

    class(sll_t_pic_visualizer_interface),     intent( inout ) :: self
    class(sll_c_particle_group_base), pointer, intent( inout ) :: particle_group
    character(len=*),                          intent(in)      :: array_name        !< field name
    sll_int32,                                 intent(in)      :: plot_np_x         !< nb of points in the x  plotting grid
    sll_int32,                                 intent(in)      :: plot_np_y         !< nb of points in the y  plotting grid
    sll_int32,                                 intent(in)      :: plot_np_vx        !< nb of points in the vx plotting grid
    sll_int32,                                 intent(in)      :: plot_np_vy        !< nb of points in the vy plotting grid
    sll_int32,                                 intent(in)      :: iplot             !< plot counter

    select type ( particle_group )

    type is ( sll_t_pic_lbfr_4d_group )

      call particle_group%pic_lbfr_4d_visualize_f_slice_x_vx(array_name, plot_np_x, plot_np_y, plot_np_vx, plot_np_vy, iplot)

    class default
      SLL_ERROR("visualize_particle_group", "procedure not implemented for this type of particle group")

    end select

  end subroutine visualize_particle_group


end module  sll_m_pic_visualizer_interface
