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
!> @author MCP
!> @brief Interface routines for visualizing particle groups.

module sll_m_particle_visualization_interface

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_assert.h"
#include "sll_errors.h"
#include "sll_memory.h"
#include "sll_working_precision.h"

  use sll_m_cartesian_meshes, only: &
    sll_f_new_cartesian_mesh_4d, &
    sll_t_cartesian_mesh_4d

  use sll_m_particle_group_base, only: &
    sll_c_particle_group_base

  use sll_m_particle_group_2d2v_lbf, only: &
    sll_t_particle_group_2d2v_lbf

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
      sll_int32          :: plot_np_x         !< nb of points in the x  plotting grid
      sll_int32          :: plot_np_vx        !< nb of points in the vx plotting grid
      sll_real64         :: slice_y           ! position of the slice in y dimension
      sll_real64         :: slice_vy          ! position of the slice in vy dimension
    contains
      procedure :: reset_params
  end type sll_t_plotting_params_2d

contains

  subroutine reset_params(self, &
          field_name,   &
          plot_np_x,  &
          plot_np_vx,   &
          slice_y, &
          slice_vy &
    )
    class(sll_t_plotting_params_2d),  intent( inout )      :: self
    character(len=*),                 intent( in    )      :: field_name
    sll_int32,                        intent( in    )      :: plot_np_x
    sll_int32,                        intent( in    )      :: plot_np_vx
    sll_real64,                       intent( in    )      :: slice_y
    sll_real64,                       intent( in    )      :: slice_vy

    self%field_name = trim(field_name)
    self%plot_np_x  = plot_np_x
    self%plot_np_vx = plot_np_vx
    self%slice_y = slice_y
    self%slice_vy = slice_vy
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
    sll_real64, dimension(:,:,:,:),            pointer         :: grid_values_4d    !< point values of f on the 4d grid
    type(sll_t_cartesian_mesh_4d),             pointer         :: plotting_grid_4d
    sll_int32     :: reconstruction_set_type
    sll_real64    :: dummy_total_charge
    logical       :: enforce_total_charge
    logical       :: reconstruct_f_on_last_node(4)
    sll_int32     :: file_id
    sll_int32     :: ierr

    select type ( particle_group )

    type is ( sll_t_particle_group_2d2v_lbf )

    SLL_ALLOCATE( grid_values_4d(plotting_params_2d%plot_np_x, 1, plotting_params_2d%plot_np_vx, 1), ierr)

    SLL_ASSERT( plotting_params_2d%plot_np_x  > 1 )
    SLL_ASSERT( plotting_params_2d%plot_np_vx > 1 )
    SLL_ASSERT( plotting_params_2d%slice_y  >= particle_group%lbf_grid%eta2_min )
    SLL_ASSERT( plotting_params_2d%slice_y  <  particle_group%lbf_grid%eta2_max )
    SLL_ASSERT( plotting_params_2d%slice_vy >= particle_group%lbf_grid%eta4_min )
    SLL_ASSERT( plotting_params_2d%slice_vy <  particle_group%lbf_grid%eta4_max )

    plotting_grid_4d => sll_f_new_cartesian_mesh_4d(  &
        plotting_params_2d%plot_np_x  - 1,  &
        1,  &
        plotting_params_2d%plot_np_vx - 1,  &
        1,  &
        particle_group%lbf_grid%eta1_min, &
        particle_group%lbf_grid%eta1_max, &
        plotting_params_2d%slice_y, &
        particle_group%lbf_grid%eta2_max, &   ! max value along y does not matter because upper bound will not be plotted
        particle_group%lbf_grid%eta3_min, &
        particle_group%lbf_grid%eta3_max, &
        plotting_params_2d%slice_vy, &
        particle_group%lbf_grid%eta4_max  &   ! max value along vy does not matter because upper bound will not be plotted
    )
    reconstruct_f_on_last_node(1) = .true.
    reconstruct_f_on_last_node(2) = .false.
    reconstruct_f_on_last_node(1) = .true.
    reconstruct_f_on_last_node(4) = .false.

    dummy_total_charge = 0.0_f64
    enforce_total_charge = .false.

    call particle_group%reconstruct_f_lbf_on_given_grid(  &
        plotting_grid_4d,          &
        grid_values_4d,            &
        reconstruct_f_on_last_node,&
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
      write(file_id,*) '# run this script with $ gnuplot ' // trim(plotting_params_2d%field_name) // '.gnu' // ' -persist'
      write(file_id,*) 'set view 0,0'
      write(file_id,*) 'set pm3d'
      write(file_id,*) 'set hid'
      close(file_id)
    end if

    call sll_o_gnuplot_2d(  &
        particle_group%lbf_grid%eta1_min, &
        particle_group%lbf_grid%eta1_max, &
        plotting_params_2d%plot_np_x,   &                ! (note: this is indeed the nb of plotted points, not 'cells')
        particle_group%lbf_grid%eta3_min, &
        particle_group%lbf_grid%eta3_max, &
        plotting_params_2d%plot_np_vx,  &                ! (same comment)
        grid_values_4d(:,1,:,1),   &
        trim(plotting_params_2d%field_name),  &
        iplot,      &
        ierr   &
        )
        ! todo: use this: force_keep_gnu_file=.true. )     ! optional argument do avoid replacing existing file even if iplot=1

    plotting_params_2d%plot_count = plotting_params_2d%plot_count + 1

    class default
      SLL_ERROR("sll_s_visualize_particle_group", "procedure not implemented for this type of particle group")
    end select

  end subroutine sll_s_visualize_particle_group

end module sll_m_particle_visualization_interface