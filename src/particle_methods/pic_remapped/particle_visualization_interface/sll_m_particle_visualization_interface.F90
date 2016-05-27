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
#include "sll_memory.h"
#include "sll_errors.h"
#include "sll_working_precision.h"

  use sll_m_cartesian_meshes, only: &
    sll_f_new_cartesian_mesh_4d, &
    sll_t_cartesian_mesh_4d

  use sll_m_particle_group_base, only: &
    sll_c_particle_group_base

  use sll_m_pic_lbfr_4d_group, only: &
    sll_t_pic_lbfr_4d_group,  &
    SLL_PIC_LBFR_GIVEN_GRID

  use sll_m_gnuplot, only: &
    sll_o_gnuplot_2d

  implicit none

  public :: &
    sll_t_plotting_params_2d, &
    sll_s_visualize_particle_group

  private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  type sll_t_plotting_params_2d
      character(len=256) :: field_name        !< may need to be trimmed when used
      sll_int32          :: plot_count
      sll_int32          :: plot_np_x         !< nb of points in the x  plotting grid  todo: rename plot_np_1
      sll_int32          :: plot_np_y         !< nb of points in the y  plotting grid  todo: discard, really use slices
      sll_int32          :: plot_np_vx        !< nb of points in the vx plotting grid  todo: rename plot_np_2
      sll_int32          :: plot_np_vy        !< nb of points in the vy plotting grid  todo: discard, really use slices
      sll_real64         :: slice_y           ! todo: use this
      sll_real64         :: slice_vy          ! todo: use this
    contains
      procedure :: reset_params
  end type sll_t_plotting_params_2d

contains

  subroutine reset_params(self, &
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

    self%field_name = trim(field_name)
    self%plot_np_x  = plot_np_x
    self%plot_np_y  = plot_np_y
    self%plot_np_vx = plot_np_vx
    self%plot_np_vy = plot_np_vy
    self%plot_count = 0

  end subroutine

  !> visualization interface
  subroutine sll_s_visualize_particle_group( &
          particle_group,   &
          plotting_params_2d, &
          iplot )

    class(sll_c_particle_group_base), pointer, intent( inout ) :: particle_group
    class(sll_t_plotting_params_2d),           intent( inout ) :: plotting_params_2d
    sll_int32,                                 intent( in    ) :: iplot             !< plot counter

    ! sll_real64, dimension(:,:),       pointer :: x_vx_grid_values
    sll_real64, dimension(:,:),                pointer         :: x_vx_grid_values  !< point values of f on the x_vx grid
    type(sll_t_cartesian_mesh_4d),             pointer         :: plotting_grid_4d
    sll_int32     :: reconstruction_set_type
    sll_real64    :: dummy_total_charge
    logical       :: enforce_total_charge
    sll_int32     :: file_id
    sll_int32     :: ierr

    select type ( particle_group )

    type is ( sll_t_pic_lbfr_4d_group )

    SLL_ALLOCATE( x_vx_grid_values(plotting_params_2d%plot_np_x, plotting_params_2d%plot_np_vx), ierr)

    plotting_grid_4d => sll_f_new_cartesian_mesh_4d(  &
        plotting_params_2d%plot_np_x  - 1,  &
        plotting_params_2d%plot_np_y  - 1,  &
        plotting_params_2d%plot_np_vx - 1,  &
        plotting_params_2d%plot_np_vy - 1,  &
        particle_group%remapping_grid_eta_min(1), &
        particle_group%remapping_grid_eta_max(1), &
        particle_group%remapping_grid_eta_min(2), &
        particle_group%remapping_grid_eta_max(2), &
        particle_group%remapping_grid_eta_min(3), &
        particle_group%remapping_grid_eta_max(3), &
        particle_group%remapping_grid_eta_min(4), &
        particle_group%remapping_grid_eta_max(4)  &
    )

    reconstruction_set_type = SLL_PIC_LBFR_GIVEN_GRID
    dummy_total_charge = 0.0_f64
    enforce_total_charge = .false.

    call particle_group%pic_lbfr_4d_reconstruct_f(  &
        reconstruction_set_type,   &
        plotting_grid_4d,          &
        x_vx_grid_values,          &
        dummy_total_charge,        &
        enforce_total_charge       &
    )

    if( plotting_params_2d%plot_count == 0 )then

      ! Open Gnuplot script file: a new ASCII file will be created (replaced if already existing)
      open( file= trim(plotting_params_2d%field_name)//'.gnu', &
        status  = 'replace',   &
        form    = 'formatted', &
        position= 'append',    &
        newunit = file_id,     &
        iostat  = ierr )

      ! Write Gnuplot instructions for plotting f, then close file
      write(file_id,*) "set view 0,0"
      write(file_id,*) "set pm3d"
      write(file_id,*) "set hid"
      close(file_id)
    end if

    call sll_o_gnuplot_2d(  &
        particle_group%remapping_grid_eta_min(1), &
        particle_group%remapping_grid_eta_max(1), &
        plotting_params_2d%plot_np_x,   &                ! (note: this is indeed the nb of plotted points, not 'cells')
        particle_group%remapping_grid_eta_min(3), &
        particle_group%remapping_grid_eta_max(3), &
        plotting_params_2d%plot_np_vx,  &                ! (same comment)
        x_vx_grid_values,   &
        trim(plotting_params_2d%field_name),  &
        iplot,      &
        ierr   &
        )
        ! todo: use this: force_keep_gnu_file=.true. )        ! do not replace existing file, no matter what iplot is

    plotting_params_2d%plot_count = plotting_params_2d%plot_count + 1

    class default
      SLL_ERROR("sll_s_visualize_particle_group", "procedure not implemented for this type of particle group")
    end select

  end subroutine sll_s_visualize_particle_group

end module sll_m_particle_visualization_interface