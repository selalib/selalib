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

!> @ingroup particle_methods
!> @author MCP
!> @brief Interface routines for visualizing particle groups.

module sll_m_particle_visualization_interface

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_assert.h"
#include "sll_errors.h"
#include "sll_working_precision.h"

  use sll_m_particle_group_base, only: &
    sll_c_particle_group_base

  use sll_m_pic_lbfr_4d_group, only: &
    sll_t_pic_lbfr_4d_group

  implicit none

  public :: &
    sll_t_plotting_params_2d, &
    sll_s_visualize_particle_group

  private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  !> type used to enforce some conservation properties in the sampling -- (there may be more than just one scalar)
  type sll_t_plotting_params_2d
      character(len=256) :: field_name
      sll_int32          :: plot_np_x         !< nb of points in the x  plotting grid
      sll_int32          :: plot_np_y         !< nb of points in the y  plotting grid
      sll_int32          :: plot_np_vx        !< nb of points in the vx plotting grid
      sll_int32          :: plot_np_vy        !< nb of points in the vy plotting grid
    contains
      procedure :: set_params
  end type sll_t_plotting_params_2d

contains

  subroutine set_params(self, &
          field_name,   &
          plot_np_x,  &
          plot_np_y,  &
          plot_np_vx,   &
          plot_np_vy  &
    )
    class(sll_t_plotting_params_2d),  intent( inout )      :: self
    character(len=*),                 intent( in    )      :: field_name
    sll_int32,                        intent( in    )      :: plot_np_x         !< nb of points in the x  plotting grid
    sll_int32,                        intent( in    )      :: plot_np_y         !< nb of points in the y  plotting grid
    sll_int32,                        intent( in    )      :: plot_np_vx        !< nb of points in the vx plotting grid
    sll_int32,                        intent( in    )      :: plot_np_vy        !< nb of points in the vy plotting grid

    self%field_name = field_name
    self%plot_np_x  = plot_np_x
    self%plot_np_y  = plot_np_y
    self%plot_np_vx = plot_np_vx
    self%plot_np_vy = plot_np_vy

  end subroutine

  !> visualization interface
  subroutine sll_s_visualize_particle_group( &
          particle_group,   &
          plotting_params_2d, &
          iplot )

    class(sll_c_particle_group_base), pointer, intent( inout ) :: particle_group
    class(sll_t_plotting_params_2d),           intent( in    ) :: plotting_params_2d
    sll_int32,                                 intent( in    ) :: iplot             !< plot counter

    select type ( particle_group )

    type is ( sll_t_pic_lbfr_4d_group )
      call particle_group%pic_lbfr_4d_visualize_f_slice_x_vx( &
          trim(plotting_params_2d%field_name),  &
          plotting_params_2d%plot_np_x,   &
          plotting_params_2d%plot_np_y,   &
          plotting_params_2d%plot_np_vx,  &
          plotting_params_2d%plot_np_vy,  &
          iplot )

    class default
      SLL_ERROR("sll_s_visualize_particle_group", "procedure not implemented for this type of particle group")
    end select

  end subroutine sll_s_visualize_particle_group

end module sll_m_particle_visualization_interface