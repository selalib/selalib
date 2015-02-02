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

!> @file
!> @brief Originally copied from [[file:particle_group_2d.F90]] ALH MCP

module sll_mlt_pic_2d_group_module

#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"
#include "particle_representation.h"

  use sll_cartesian_meshes
  use sll_mlt_pic_2d_particle_module ! [[file:mlt_pic_2d_particle.F90]]

  implicit none

  !> ALH MCP <<level_part_indices_2d>> multi-dimentional array of particles, indexed by (k_x, k_vx, etc) which contains
  !> the integer index of particles. This is used to find the neighbours of each particle. The integer index of a
  !> particle corresponds is the index of that particle in [[p_list]] which is of type [[level_part_list_2d]].

  type :: level_part_indices_2d
     sll_int32 :: level

     !> <<p_level_indices>> contains the  particles
     sll_int32,dimension(:,:),pointer :: p_level_indices
  end type level_part_indices_2d

  ! <<level_part_list_2d>> a <<list>> is a one-dimensional array of existing particles, numbered from 1. There is one
  ! separate particle list for each refinement level.

  type :: level_part_list_2d
     sll_int32 :: level
     sll_int32 :: number_active_particles

     !> <<p_level_list>> contains pointers to the actual particles from
     !> [[file:mlt_pic_2d_particle.F90::sll_mlt_pic_2d_particle]] for each refinement level. Particles pointers are
     !> only allocated on demand. The array of Fortran pointers is defined as in
     !> [[http://nf.nci.org.au/training/FortranAdvanced/slides/slides.035.html]].

     type(sll_mlt_pic_2d_particle_pointer),dimension(:),pointer :: p_level_list

  end type level_part_list_2d

  !> <<level_target_values_2d>> ALH MCP Target values that particles must have at the end of the remapping process

  type :: level_target_values_2d
     sll_int32 :: level
     sll_int32 :: number_parts_x
     sll_int32 :: number_parts_vx
     type(sll_cartesian_mesh_2d),pointer :: level_remapping_grid
     sll_real64,dimension(:,:),pointer :: level_target_values
  end type level_target_values_2d

  ! <<sll_mlt_pic_2d_group>>

  type :: sll_mlt_pic_2d_group

     ! ALH_MCP_06_2014
     sll_int32 :: spline_degree

     sll_int32 :: number_particles! peut etre a faire en SLL_PRIV
     sll_int32 :: active_particles! tout ça doit passer en 32

     ! "guard" means particles going out of the domain
     sll_int32 :: guard_list_size! tout ça doit passer en 32

     ! qoverm = q/m (charge/mass)
     sll_real64 :: qoverm 

     ! ALH_MCP_06_2014
     logical :: domain_is_x_periodic

     ! <<x_mesh>>
     type(sll_cartesian_mesh_1d), pointer :: x_mesh

     !> We need a separate list of particles for each level of [[level_part_list_2d]]

     type(level_part_list_2d), dimension(:), pointer :: p_list ! <<p_list>>
     type(sll_mlt_pic_2d_particle_guard), dimension(:), pointer :: p_guard

     !> ALH MCP One full matrix of indices for each refinement level of particles - p_ is for "particle_" (cf
     !> [[level_part_indices_2d]])

     type(level_part_indices_2d), dimension(:),pointer :: p_indices

     ! <<target_values>> uses [[level_target_values_2d]]

     type(level_target_values_2d), dimension(:),pointer :: target_values

     !> <<min_level>> <<max_level>> Minimum and maximum refinement level

     sll_int32 :: min_level
     sll_int32 :: max_level

  end type sll_mlt_pic_2d_group

  interface sll_delete
     module procedure sll_mlt_pic_2d_group_delete
  end interface sll_delete

contains

  function sll_mlt_pic_2d_group_new( &
       spline_degree,&
       remap_grid_vx_min,                       &
       remap_grid_vx_max,                       &
       particle_array_size,                     &
       guard_list_size,                         &
       qoverm,                                  &
       domain_is_x_periodic,                    &
       x_mesh,                                    &
       min_level,                               &
       max_level                                &
       ) result(p_group)

    type(sll_mlt_pic_2d_group), pointer :: p_group

    sll_int32,intent(in)::spline_degree
    sll_real64, intent(in) :: remap_grid_vx_min
    sll_real64, intent(in) :: remap_grid_vx_max
    sll_int32,  intent(in) :: particle_array_size
    sll_int32,  intent(in) :: guard_list_size
    sll_real64, intent(in) :: qoverm
    logical,    intent(in) :: domain_is_x_periodic

    ! cf [[file:../meshes/sll_cartesian_meshes.F90::sll_cartesian_mesh_2d]]

    type(sll_cartesian_mesh_1d), pointer :: x_mesh

    sll_int32,  intent(in) :: min_level
    sll_int32,  intent(in) :: max_level
    sll_int32 :: number_particles
    sll_int32 :: level
    sll_int32 :: ierr
    sll_int32 :: level_number_part_x
    sll_int32 :: level_number_part_vx

    sll_int32 :: remap_grid_number_cells_x
    sll_int32 :: remap_grid_number_cells_vx
    sll_real64 :: remap_grid_x_min
    sll_real64 :: remap_grid_x_max

    ! ALH MCP 06/2014 - <<Allocating_the_3D_structure>> (refinement level,x index,v index) of particle indices
    SLL_ALLOCATE(p_group,ierr)
    p_group%spline_degree=spline_degree
    number_particles = 0
    p_group%min_level = min_level
    p_group%max_level = max_level
    p_group%domain_is_x_periodic = domain_is_x_periodic
    SLL_ALLOCATE(p_group%p_indices(max_level),ierr)
    SLL_ALLOCATE(p_group%p_list(max_level),ierr)

    ! cf [[target_values]]
    SLL_ALLOCATE(p_group%target_values(max_level),ierr)
    do level = min_level,max_level

       ! cf [[number_parts_x]]
       if(domain_is_x_periodic)then
          level_number_part_x = 2**level
       else
          level_number_part_x = 2**level+1
       endif
       level_number_part_vx = 2**level+1

       p_group%target_values(level)%number_parts_x = level_number_part_x
       p_group%target_values(level)%number_parts_vx = level_number_part_vx

       number_particles = number_particles &
            + p_group%target_values(level)%number_parts_x*p_group%target_values(level)%number_parts_vx

       p_group%p_indices(level)%level = level

       ! <<allocating_p_level_indices>>
       SLL_ALLOCATE(p_group%p_indices(level)%p_level_indices(level_number_part_x,level_number_part_vx),ierr)
       p_group%p_indices(level)%p_level_indices = 0

       p_group%p_list(level)%level = level

       ! <allocating_p_level_list>> size may be optimized in [[p_level_list]]

       SLL_ALLOCATE(p_group%p_list(level)%p_level_list(level_number_part_x * level_number_part_vx),ierr)

       p_group%target_values(level)%level = level

       SLL_ALLOCATE(p_group%target_values(level)%level_remapping_grid,ierr)
       remap_grid_number_cells_x  = 2**level
       remap_grid_number_cells_vx = 2**level
       remap_grid_x_min           = x_mesh%eta_min ! cf [[file:../meshes/sll_cartesian_meshes.F90::eta_min]]
       remap_grid_x_max           = x_mesh%eta_max ! cf [[file:../meshes/sll_cartesian_meshes.F90::eta_max]]
       p_group%target_values(level)%level_remapping_grid => new_cartesian_mesh_2d( &
            remap_grid_number_cells_x,                                       &
            remap_grid_number_cells_vx,                                      &
            remap_grid_x_min,                                                &
            remap_grid_x_max,                                                &
            remap_grid_vx_min,                                               &
            remap_grid_vx_max)

       ! <<allocating_level_target_values>>
       SLL_ALLOCATE(p_group%target_values(level)%level_target_values(2**level,2**level),ierr)

       p_group%target_values(level)%level_target_values = 0
    enddo

    if( number_particles > particle_array_size ) then
       print *, 'sll_mlt_pic_2d_group_new(): ERROR, number_particles=',number_particles,' should not ', &
            'be greater than the memory size requested, particle_array_size=',particle_array_size,'.'
       STOP
    end if

    p_group%number_particles = number_particles
    p_group%active_particles = number_particles
    p_group%guard_list_size  = guard_list_size
    p_group%qoverm           = qoverm

    SLL_ALLOCATE(p_group%p_guard(guard_list_size),ierr)

    if (.not.associated(x_mesh) ) then
       print*, 'error: passed x_mesh not associated'
    endif

    ! for p_group%x_mesh see [[x_mesh]]
    p_group%x_mesh => x_mesh

  end function sll_mlt_pic_2d_group_new

  ! <<sll_mlt_pic_2d_group_delete>>

  subroutine sll_mlt_pic_2d_group_delete(p_group)
    type(sll_mlt_pic_2d_group), pointer :: p_group
    sll_int32 :: level
    sll_int32 :: ierr

    if(.not. associated(p_group) ) then
       print *, 'sll_mlt_pic_2d_group_delete(): ERROR, passed group was not associated.'
    end if

    ! ALH_MCP_06_2014 Deallocations corresponding to [[Allocating_the_3D_structure]]

    do level = p_group%min_level,p_group%max_level
       SLL_DEALLOCATE(p_group%target_values(level)%level_remapping_grid,ierr)

       ! cf [[allocating_level_target_values]]
       SLL_DEALLOCATE(p_group%target_values(level)%level_target_values,ierr)

       ! cf [[allocating_p_level_indices]]
       SLL_DEALLOCATE(p_group%p_indices(level)%p_level_indices,ierr)

       ! cf [[allocating_p_level_list]]
       SLL_DEALLOCATE(p_group%p_list(level)%p_level_list,ierr)
    enddo

    SLL_DEALLOCATE(p_group%target_values,ierr)
    SLL_DEALLOCATE(p_group%p_indices,ierr)

    SLL_DEALLOCATE(p_group%p_list,ierr)
    SLL_DEALLOCATE(p_group%p_guard,ierr)
    SLL_DEALLOCATE(p_group,ierr)

  end subroutine sll_mlt_pic_2d_group_delete
end module sll_mlt_pic_2d_group_module
