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
!> @namespace sll_multilevel_lt_particle_representations

! ALH_MCP_06_2014 - copied from [[file:sll_particle_init4D.F90]]

module sll_mlt_pic_2d_initializers

#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"
#include "particle_representation.h"

  use sll_constants, only: sll_pi
  use sll_mlt_pic_2d_group_module ! [[file:../pic_particle_types/mlt_pic_2d_group.F90::sll_mlt_pic_2d_group_module]]
  use sll_cartesian_meshes
  use deviators
  use sll_representation_conversion_module

  implicit none
  
!   private sll_init_spatial_particle2D, suite_hamm
! !  In the future, we have to build a directory called
! !  random_number_generators 


contains

  ! <<sll_initialize_multilevel_lt_particles_2d>> initialize the lt_particle group with the landau f0 distribution
  subroutine sll_initialize_multilevel_lt_particles_2d( &
              thermal_speed, alpha, k,                  &
              p_group )
    sll_real64, intent(in) :: thermal_speed, alpha, k

    ! ALH_MCP_06_2014 [[file:../pic_particle_types/mlt_pic_2d_group.F90::sll_mlt_pic_2d_group]]

    type(sll_mlt_pic_2d_group), pointer, intent(inout) :: p_group

    call write_landau_density_on_multilevel_remap_grid( thermal_speed, alpha, k_landau, p_group )
    !call plot_2d_slice_remapping_grid("init_values_on_rg.dat", p_group )
    call sll_compute_new_multilevel_lt_particles_2d( p_group ) ! [[sll_compute_new_multilevel_lt_particles_2d]]
    
  end subroutine sll_initialize_multilevel_lt_particles_2d

  ! ALH_MCP_06_2014 <<write_landau_density_on_multilevel_remap_grid>> Setting the Landau density values on each level of the remap grids

  subroutine write_landau_density_on_multilevel_remap_grid(    &
       thermal_speed, alpha, k_landau,                         &
       p_group                                                 &
       )

    sll_real64, intent(in)                                             :: thermal_speed, alpha, k_landau

    ! [[file:../pic_particle_types/mlt_pic_2d_group.F90::sll_mlt_pic_2d_group]]
    type(sll_mlt_pic_2d_group), pointer, intent(inout)  :: p_group

    sll_int64                                                          :: k_x
    sll_int64                                                          :: k_vx
    sll_int64                                                          :: number_particles
    sll_int64                                                          :: number_parts_x
    sll_int64                                                          :: number_parts_vx
    sll_real64                                                         :: h_parts_x    
    sll_real64                                                         :: h_parts_vx   
    sll_real64                                                         :: parts_x_min  
    sll_real64                                                         :: parts_vx_min 
    sll_real64                                                         :: one_over_thermal_velocity
    sll_real64                                                         :: one_over_two_pi
    sll_real64                                                         :: x_k
    sll_real64                                                         :: vx_k
    type(sll_cartesian_mesh_1d),      pointer                            :: m1d
    sll_real64                                                         :: f_x, f_vx
    sll_int32                                                          :: level

    number_particles = p_group%number_particles
    one_over_thermal_velocity = 1./thermal_speed   
    one_over_two_pi = 1./(2*sll_pi)
    
    ! Poisson mesh associated with the particles [[file:../pic_particle_types/mlt_pic_2d_group.F90::x_mesh]]
    m1d => p_group%x_mesh

    do level = p_group%min_level,p_group%max_level
       h_parts_x    = p_group%target_values(level)%level_remapping_grid%delta_eta1
       h_parts_vx   = p_group%target_values(level)%level_remapping_grid%delta_eta2

       parts_x_min    = p_group%target_values(level)%level_remapping_grid%eta1_min
       parts_vx_min   = p_group%target_values(level)%level_remapping_grid%eta2_min

       number_parts_x  = p_group%target_values(level)%number_parts_x
       number_parts_vx = p_group%target_values(level)%number_parts_vx

       ! compute the values of f0 on the (cartesian, phase-space) remapping grid at that level
       x_k = parts_x_min

       ! number_parts = number of particles
       ! cf [[file:../pic_particle_types/multilevel_lt_particle_group_2d.F90::number_parts_x]]
       do k_x = 1, number_parts_x
          f_x = eval_landau(alpha, k_landau, x_k) ! [[eval_landau]]
          vx_k = parts_vx_min
          do k_vx = 1, number_parts_vx
             f_vx = one_over_thermal_velocity * exp(-0.5*(vx_k*one_over_thermal_velocity)**2)

             ! [[file:../pic_particle_types/multilevel_lt_particle_group_2d.F90::level_target_values]]
             p_group%target_values(level)%level_target_values(k_x,k_vx) = one_over_two_pi * f_x * f_vx
             vx_k = vx_k + h_parts_vx
          end do
          x_k = x_k + h_parts_x
       end do
    end do
  end subroutine write_landau_density_on_multilevel_remap_grid
    
  ! <<sll_compute_new_multilevel_lt_particles_2d>> position the particle on the cartesian remapping grid and compute new
  ! weights in order to approximate the point values stored in the remapping grid

  subroutine sll_compute_new_multilevel_lt_particles_2d(p_group)

    ! [[file:../pic_particle_types/mlt_pic_2d_group.F90::sll_mlt_pic_2d_group]]
    type(sll_mlt_pic_2d_group), pointer, intent(inout) :: p_group
    sll_int64 :: k, k_ngb
    sll_int64 :: k_x
    sll_int64 :: j_y
    sll_int64 :: k_vx
    sll_int64 :: j_vy
    sll_int32 :: ierr
    sll_int64 :: number_parts_x
    sll_int64 :: number_parts_y
    sll_int64 :: number_parts_vx
    sll_int64 :: number_parts_vy
    sll_real64 :: h_parts_x    
    sll_real64 :: h_parts_y    
    sll_real64 :: h_parts_vx   
    sll_real64 :: h_parts_vy   
    sll_real64 :: parts_x_min  
    sll_real64 :: parts_y_min  
    sll_real64 :: parts_vx_min 
    sll_real64 :: parts_vy_min 
    sll_real32 :: w_k
    sll_real64 :: x_k
    sll_real64 :: y_j
     sll_real64 :: vx_k
    sll_real64 :: vy_j
    sll_real32 :: d_vol
    type(sll_cartesian_mesh_1d),      pointer  :: m1d
    sll_int32 :: level

    ! Poisson mesh associated with the particles [[file:../pic_particle_types/mlt_pic_2d_group.F90::x_mesh]]

    m1d => p_group%x_mesh
    
    ! [[file:../pic_particle_types/mlt_pic_2d_group.F90::min_level]]
    do level = p_group % min_level, p_group % max_level

       ! computing d_j = A_j(g_j-g_{j-1}) as in [[file:~/mcp/maltpic/algos-maltpic.tex::hierarchic_approx]]

       ! <<calling_substract_multilevel_lt_pic_density_from_remap_grid>> at
       ! [[substract_multilevel_lt_pic_density_from_remap_grid]]

       call substract_multilevel_lt_pic_density_from_remap_grid(p_group,level)

       !<<current>>
       aaa

!       ! compute the particle weights from the values of f0 on the (cartesian, phase-space) remapping grid
!       ! k = particle index
!       k = 0
!       x_k = parts_x_min
!       do k_x = 1, number_parts_x
!          vx_k = parts_vx_min
!          do k_vx = 1, number_parts_vx
!            
!             k = k+1
!             if( p_group%spline_degree == 1 )then
!                w_k = d_vol * real( p_group%target_values(level)%level_target_values(k_x,k_vx) ,f32)
!            else
!               print *, 'sll_compute_new_multilevel_lt_particles_2d(): ERROR, value of p_group%spline_degree ', &
!                        ' is invalid: ', p_group%spline_degree
!               STOP
!            end if
!!            print *, '-- -- -- sll_lt_pic_initialize_some4Dfunction(): ', &
!!                        'number_parts_x, k_x = ', number_parts_x, k_x, &
!!                        'number_parts_y, j_y = ', number_parts_y, j_y
!!            print *, 'm2d%eta1_max, p_group%remapping_grid%eta1_max , x_k = ', &
!!                        m2d%eta1_max, p_group%remapping_grid%eta1_max , x_k
!!            print *, 'm2d%eta2_max, p_group%remapping_grid%eta2_max , y_j = ', &
!!                        m2d%eta2_max, p_group%remapping_grid%eta2_max , y_j
!                        
!!            if( w_k > 0 )then
!            p_group%p_list(k)%q = w_k
!
!            ! each particle knows its indices back in (x,vx) form
!            ! [[file:../pic_particle_types/multilevel_lt_particle_group_2d.F90::level_sll_sll_particle_2d]]
!            p_group%p_list(k)%k_x = k_x
!            p_group%p_list(k)%k_vx = k_vx
!            p_group%p_indices(level)%p_level_indices(k_x,k_vx) = k
!!            end if
!            
!            call global_to_cell_offset (  &
!                 x_k, &
!                 m1d, &
!                 p_group%p_list(k)%ic, &
!                 p_group%p_list(k)%dx )
!            p_group%p_list(k)%vx = vx_k 
!
!            vx_k = vx_k + h_parts_vx
!         end do
!         x_k = x_k + h_parts_x
!      end do
!   end do

  end subroutine sll_compute_new_multilevel_lt_particles_2d

  ! <<eval_landau>>
    
  function eval_landau(alp, kx, x)
    sll_real64 :: alp, kx, x
    sll_real64 :: eval_landau
    eval_landau = 1._f64 + alp * cos(kx * x)
  end function eval_landau

  ! <<substract_multilevel_lt_pic_density_from_remap_grid>>

  subroutine substract_multilevel_lt_pic_density_from_remap_grid(p_group,level)
    type(sll_multilevel_lt_particle_group_2d), pointer, intent(inout) :: p_group
    sll_int32 :: level

    ! number of particles at this level 
    number_parts_x  = p_group % target_values(level) % number_parts_x
    number_parts_vx = p_group % target_values(level) % number_parts_vx

    ! remove d_j from target values on finer levels
    do j_aux = level+1, max_level 

       ! for every point of the current level, remove its value from all finer levels

       do k_x = 1, number_parts_x
          do k_vx = 1, number_parts_vx

             ! find the location of the same particle on the finer grid

             k_x_aux = (k_x - 1) * 2 ** (j_aux - level) + 1
             k_vx_aux = (k_vx - 1) * 2 ** (j_aux - level) + 1
       
             ! remove the value from current level to level j_aux

             ! do we need to remove the value from all the fine particles in the support of the bessel function? But
             ! with which weight?

             xxx

       h_parts_x    = p_group % target_values(level) % level_remapping_grid % delta_eta1
       h_parts_vx   = p_group % target_values(level) % level_remapping_grid % delta_eta2

       parts_x_min    = p_group % target_values(level) % level_remapping_grid % eta1_min
       parts_vx_min   = p_group % target_values(level) % level_remapping_grid % eta2_min


       d_vol = h_parts_x*h_parts_vx

       xxx
    enddo
    xxx
  end subroutine substract_multilevel_lt_pic_density_from_remap_grid

end module sll_mlt_pic_2d_initializers
