!**************************************************************
!  Copyright INRIA
!  Authors : 
!     MCP ALH
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

!> @author MCP ALH

!> @brief Module for groups of particles for particle methods with linearized-backward-flow (lbf) resamplings

! <<<<particle_lbf_group>>>> AAA_ALH_HERE Same method as [[file:../pic_lbfr/sll_m_pic_lbfr_4d_group.F90]] which much less options to allow
! for more optimizations

module sll_m_particle_group_lbf_2d2v

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_assert.h"
#include "sll_errors.h"
#include "sll_memory.h"
#include "sll_working_precision.h"

  use sll_m_particle_lbf_2d2v, only: &
    sll_t_particle_lbf_2d2v

  use sll_m_cartesian_meshes, only: &
    sll_f_new_cartesian_mesh_4d, &
    sll_t_cartesian_mesh_4d


  implicit none

  public :: &
    sll_t_particle_group_lbf_2d2v, &
    sll_s_new_particle_group_lbf_2d2v_ptr

  private

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  !> Group of @ref sll_t_particle_lbf_4d_marker
  type, extends(sll_c_particle_group_base) :: sll_t_particle_group_lbf_2d2v

   !> @name The particles
   !> @{
    sll_int32                                                   :: n_particles_x
    sll_int32                                                   :: n_particles_y
    sll_int32                                                   :: n_particles_vx
    sll_int32                                                   :: n_particles_vy

    type(sll_t_particle_lbf_2d2v),   dimension(:), allocatable  :: particles_list  ! name was: struct_markers_list

    ! replacing the flow grid, the markers grid and the deposition particles grid:
    type(sll_t_cartesian_mesh_4d), pointer                      :: particles_grid
    sll_real64                                                  :: particles_grid_h
    !> @}

    !> @name Sparse grid interpolation of the remapped density f, dimensions of remapping grid = particle grid
    !> @{
    sll_int32                                        :: remapped_f_interpolation_degree
    logical                                          :: domain_is_periodic(2)
    type(sll_t_sparse_grid_interpolator_4d)          :: sparse_grid_interpolator
    sll_int32,  dimension(4)                         :: sparse_grid_max_levels
    sll_real64, dimension(:), allocatable            :: remapped_f_sparse_grid_coefficients
    sll_real64, dimension(:), allocatable            :: tmp_f_values_on_remapping_sparse_grid !< used only during remapping
    !> @}

  contains

    !> @name Getters
    !> @{
    procedure :: get_x          => particle_lbf_4d_get_x      ! [[function%20particle_lbf_4d_get_x]]
    procedure :: get_v          => particle_lbf_4d_get_v      ! [[function%20particle_lbf_4d_get_v]]
    procedure :: get_charge     => particle_lbf_4d_get_charge ! [[function%20particle_lbf_4d_get_charge]]
    procedure :: get_mass       => particle_lbf_4d_get_mass   ! [[function%20particle_lbf_4d_get_mass]]
    !> @}
    
    !> @name Setters
    !> @{
    procedure :: set_x              => particle_lbf_4d_set_x
    procedure :: set_v              => particle_lbf_4d_set_v
    !> @}

    !> @name Sampling and resampling
    !> @{
    procedure :: resample
    !> @}

    procedure :: particle_lbf_4d_write_hat_density_on_remapping_grid        !> this evaluates an analytic f0
    procedure :: particle_lbf_4d_write_landau_density_on_remapping_grid     !> this evaluates an analytic f0

    !todo: fusionner
    procedure :: particle_lbf_4d_reset_markers_position
    procedure :: reset_deposition_particles_coordinates
    procedure :: reset_deposition_particles_weights_with_direct_interpolation

    procedure :: particle_lbf_4d_set_markers_connectivity

    procedure :: particle_lbf_4d_reconstruct_f
    procedure :: particle_lbf_4d_interpolate_value_of_remapped_f

    procedure :: get_ltp_deformation_matrix                           !> the local bwd flow using structured flow markers

    ! Initializer
    procedure :: init => initialize_particle_lbf_4d_group   !> Initialization function
    procedure :: free => delete_particle_lbf_4d_group       !> Destructor

  end type sll_t_particle_group_lbf_2d2v

contains

  ! gets the physical coordinates of a 'particle'.  Simpler version of
  ! [[file:~/sllrzg/src/particle_methods/pic_remapped/pic_lbfr/sll_m_pic_lbfr_4d_group.F90::function%20pic_lbfr_4d_get_x]]
  
  pure function particle_lbf_4d_get_x( self, i ) result( r )
    ! AAA_ALH_TODO
  end function particle_lbf_4d_get_x

  ! Simpler version of
  ! [[file:~/sllrzg/src/particle_methods/pic_remapped/pic_lbfr/sll_m_pic_lbfr_4d_group.F90::function%20pic_lbfr_4d_get_v]]
  
  ! AAA_ALH_TODO procedure :: get_v          => particle_lbf_4d_get_v

  ! Simpler version of
  ! [[file:~/sllrzg/src/particle_methods/pic_remapped/pic_lbfr/sll_m_pic_lbfr_4d_group.F90::function%20pic_lbfr_4d_get_charge]]
  
  ! AAA_ALH_TODO procedure :: get_charge     => particle_lbf_4d_get_charge
  
  ! Simpler version of
  ! [[file:~/sllrzg/src/particle_methods/pic_remapped/pic_lbfr/sll_m_pic_lbfr_4d_group.F90::function%20pic_lbfr_4d_get_mass]]
  
  ! AAA_ALH_TODO procedure :: get_mass       => particle_lbf_4d_get_mass
  
  ! Simpler version of
  ! [[file:~/sllrzg/src/particle_methods/pic_remapped/pic_lbfr/sll_m_pic_lbfr_4d_group.F90::subroutine%20pic_lbfr_4d_set_x]]
  
  ! AAA_ALH_TODO procedure :: set_x              => particle_lbf_4d_set_x
  
  ! Simpler version of
  ! [[file:~/sllrzg/src/particle_methods/pic_remapped/pic_lbfr/sll_m_pic_lbfr_4d_group.F90::subroutine%20pic_lbfr_4d_set_v]]
  
  ! AAA_ALH_TODO procedure :: set_v              => particle_lbf_4d_set_v
  
  ! AAA_ALH_TODO procedure :: resample
  ! AAA_ALH_TODO procedure :: particle_lbf_4d_write_hat_density_on_remapping_grid        !> this evaluates an analytic f0
  ! AAA_ALH_TODO procedure :: particle_lbf_4d_write_landau_density_on_remapping_grid     !> this evaluates an analytic f0
  
  ! AAA_ALH_TODO fusionner
  !procedure :: particle_lbf_4d_reset_markers_position
  !procedure :: reset_deposition_particles_coordinates
  !procedure :: reset_deposition_particles_weights_with_direct_interpolation

  ! AAA_ALH_TODO procedure :: particle_lbf_4d_set_markers_connectivity

  ! AAA_ALH_TODO procedure :: particle_lbf_4d_reconstruct_f
  ! AAA_ALH_TODO procedure :: particle_lbf_4d_interpolate_value_of_remapped_f

  ! AAA_ALH_TODO procedure :: get_ltp_deformation_matrix                           !> the local bwd flow using structured flow markers

  ! AAA_ALH_TODO procedure :: init => initialize_particle_lbf_4d_group   !> Initialization function
  ! AAA_ALH_TODO procedure :: free => delete_particle_lbf_4d_group       !> Destructor


end module sll_m_particle_group_lbf_2d2v
