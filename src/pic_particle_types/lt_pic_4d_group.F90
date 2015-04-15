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


module sll_lt_pic_4d_group_module
!module sll_lt_particle_group_4d_module   ! old name

#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"
#include "particle_representation.h"

  use sll_lt_pic_particle_module

!  use sll_particle_representations
  use sll_cartesian_meshes
  implicit none

  ! <<sll_lt_pic_4d_group>>

  type :: sll_lt_pic_4d_group
!  type :: sll_lt_particle_group_4d  ! old name
     sll_int32 :: spline_degree
     sll_int32 :: number_parts_x
     sll_int32 :: number_parts_y
     sll_int32 :: number_parts_vx
     sll_int32 :: number_parts_vy
     sll_int32 :: number_particles
     sll_int32 :: active_particles
     sll_int32 :: guard_list_size
     sll_real64 :: qoverm 
     logical    :: domain_is_x_periodic
     logical    :: domain_is_y_periodic
     logical    :: track_markers_outside_domain   !! default value is true
     type(sll_cartesian_mesh_2d), pointer                     :: mesh

     ! <<sll_lt_pic_4d_group-p_list>> uses [[file:lt_pic_4d_particle.F90::sll_lt_pic_4d_particle]]
     type(sll_lt_pic_4d_particle), dimension(:), pointer        :: p_list
     type(sll_lt_pic_4d_particle_guard_ptr), dimension(:), pointer :: p_guard

     ! num_postprocess_particles: an array indexed by the thread number (if any),
     ! of the number of particles to post-process after the main loop
     sll_int32, dimension(:), pointer :: num_postprocess_particles

     ! <<sll_lt_pic_4d_group-remapping_grid>> uses [[file:../meshes/sll_cartesian_meshes.F90::sll_cartesian_mesh_4d]]
     type(sll_cartesian_mesh_4d),               pointer       :: remapping_grid
!     sll_real64, dimension(:,:,:,:),          pointer       :: remapped_values     ! old name
     sll_real64, dimension(:,:,:,:),          pointer       :: target_values

     ! MCP [DEBUG]
     sll_real64, dimension(:,:,:,:,:,:),          pointer       :: debug_bsl_remap

     sll_real64, dimension(:),                pointer       :: ltpic_interpolation_coefs

    ! temp fields (landau damping)
    sll_real64 :: thermal_speed
    sll_real64 :: alpha_landau
    sll_real64 :: k_landau


  end type sll_lt_pic_4d_group


  interface sll_delete
     module procedure sll_lt_pic_4d_group_delete
!     module procedure delete_lt_particle_4d_group  ! old name
  end interface sll_delete

contains


  function sll_lt_pic_4d_group_new( &
!  function new_lt_particle_4d_group( &   ! old name
        spline_degree,     &
        number_parts_x,    &
        number_parts_y,    &
        number_parts_vx,   &
        number_parts_vy,   &
        remap_grid_vx_min, &
        remap_grid_vx_max, &
        remap_grid_vy_min, &
        remap_grid_vy_max, &       
        !        particle_array_size, &
        !        guard_list_size,     &
        qoverm,              &
        domain_is_x_periodic,&
        domain_is_y_periodic,&
        mesh ) result(res)

    type(sll_lt_pic_4d_group), pointer :: res

    sll_int32, intent(in)   :: spline_degree
    sll_int32, intent(in)   :: number_parts_x
    sll_int32, intent(in)   :: number_parts_y
    sll_int32, intent(in)   :: number_parts_vx
    sll_int32, intent(in)   :: number_parts_vy
    sll_real64, intent(in)  :: remap_grid_vx_min   
    sll_real64, intent(in)  :: remap_grid_vx_max   
    sll_real64, intent(in)  :: remap_grid_vy_min   
    sll_real64, intent(in)  :: remap_grid_vy_max   
    sll_real64, intent(in)  :: qoverm
    logical, intent(in)      :: domain_is_x_periodic
    logical, intent(in)      :: domain_is_y_periodic
    
    !        sll_int32, intent(in)   :: particle_array_size
    !        sll_int32, intent(in)   :: guard_list_size
    sll_int32   :: particle_array_size
    sll_int32   :: guard_list_size

    type(sll_cartesian_mesh_2d), pointer :: mesh
    sll_int32   :: ierr
    sll_int32   :: number_particles
    sll_int32   :: remap_grid_number_cells_x
    sll_int32   :: remap_grid_number_cells_y
    sll_int32   :: remap_grid_number_cells_vx
    sll_int32   :: remap_grid_number_cells_vy
    sll_real64  :: remap_grid_x_min   
    sll_real64  :: remap_grid_x_max   
    sll_real64  :: remap_grid_y_min   
    sll_real64  :: remap_grid_y_max   
        
    number_particles = number_parts_x * number_parts_y * number_parts_vx * number_parts_vy
    particle_array_size = number_particles
    guard_list_size = number_particles
    !    if( number_particles > particle_array_size ) then
    !       print *, 'sll_lt_pic_4d_group_new(): ERROR (code=6454357),  number_particles should not ', &
    !            'be greater than the requested memory size, particle_array_size.'
    !       print *, 'note: number_particles, particle_array_size = ', number_particles, particle_array_size
    !       STOP
    !    end if

    SLL_ALLOCATE( res, ierr )
    res%spline_degree = spline_degree
    res%number_parts_x   = number_parts_x
    res%number_parts_y   = number_parts_y
    res%number_parts_vx  = number_parts_vx
    res%number_parts_vy  = number_parts_vy
    res%number_particles = number_particles
    res%active_particles = number_particles
    res%guard_list_size  = guard_list_size
    res%qoverm           = qoverm
    res%domain_is_x_periodic = domain_is_x_periodic
    res%domain_is_y_periodic = domain_is_y_periodic
    res%track_markers_outside_domain = .false.
    SLL_ALLOCATE( res%p_list(particle_array_size), ierr )
    SLL_ALLOCATE( res%p_guard(guard_list_size), ierr )
    SLL_ALLOCATE( res%target_values(number_parts_x, number_parts_y, number_parts_vx, number_parts_vy), ierr )

    ! MCP [DEBUG]
    SLL_ALLOCATE( res%debug_bsl_remap(number_parts_x, number_parts_y, number_parts_vx, number_parts_vy,4,3), ierr )

    ! physical mesh, used for Poisson solver
    if (.not.associated(mesh) ) then
       print*, 'sll_lt_pic_4d_group_new(): ERROR, passed mesh not associated'
    endif
    res%mesh => mesh

    ! phase space grid to initialize and remap the particles    
    remap_grid_x_min = mesh%eta1_min
    remap_grid_x_max = mesh%eta1_max
    remap_grid_y_min = mesh%eta2_min
    remap_grid_y_max = mesh%eta2_max
    
    if( domain_is_x_periodic )then
        remap_grid_number_cells_x = number_parts_x
    else
        remap_grid_number_cells_x = number_parts_x - 1
    end if
    if( domain_is_y_periodic )then
        remap_grid_number_cells_y = number_parts_y
    else
        remap_grid_number_cells_y = number_parts_y - 1
    end if
    remap_grid_number_cells_vx = number_parts_vx - 1
    remap_grid_number_cells_vy = number_parts_vy - 1
    res%remapping_grid => new_cartesian_mesh_4d(  remap_grid_number_cells_x,        &
                                                remap_grid_number_cells_y,        & 
                                                remap_grid_number_cells_vx,       &
                                                remap_grid_number_cells_vy,       &
                                                remap_grid_x_min,   &
                                                remap_grid_x_max,   &
                                                remap_grid_y_min,   &
                                                remap_grid_y_max,   &
                                                remap_grid_vx_min,  &
                                                remap_grid_vx_max,  &
                                                remap_grid_vy_min,  &
                                                remap_grid_vy_max   & 
                                              )
  
 

    if( spline_degree == 3) then  ! ( otherwise, no need... )   
        SLL_ALLOCATE( res%ltpic_interpolation_coefs(-1:1), ierr )
           res%ltpic_interpolation_coefs(-1) = -1./6.
           res%ltpic_interpolation_coefs(0)  =  8./6.
           res%ltpic_interpolation_coefs(1)  =  -1./6.
    end if

  
  
  end function sll_lt_pic_4d_group_new


  subroutine sll_lt_pic_4d_group_delete(p_group)
    type(sll_lt_pic_4d_group), pointer :: p_group
    sll_int32 :: ierr

    if(.not. associated(p_group) ) then
       print *, 'delete_lt_particle_group_4d(): ERROR, passed group was not ', &
            'associated.'
    end if
    SLL_DEALLOCATE(p_group%p_list, ierr)
    SLL_DEALLOCATE(p_group%p_guard, ierr)
    SLL_DEALLOCATE(p_group%target_values, ierr)

    ! MCP [DEBUG]
    SLL_DEALLOCATE( p_group%debug_bsl_remap, ierr )

    SLL_DEALLOCATE(p_group, ierr)
    
  end subroutine sll_lt_pic_4d_group_delete


end module sll_lt_pic_4d_group_module
