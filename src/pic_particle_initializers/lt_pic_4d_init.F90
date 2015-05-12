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

module sll_lt_pic_4d_init
!module sll_lt_particle_initializers   ! old name
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"
#include "particle_representation.h"

  use sll_constants, only: sll_pi
  use sll_lt_pic_4d_group_module
  use sll_pic_utilities
  use sll_lt_pic_4d_utilities
  use sll_cartesian_meshes
  use sll_representation_conversion_module

  implicit none
  
!   private sll_init_spatial_particle2D, suite_hamm
! !  In the future, we have to build a directory called
! !  random_number_generators 


contains

  
   
  ! initialize the lt_pic group with the landau f0 distribution
!  subroutine sll_init_lt_particles_4d_landau (       &  ! old name
  subroutine sll_lt_pic_4d_init_landau (       &
              thermal_speed, alpha, k_landau,   &
              p_group)

    sll_real64, intent(in)                                  :: thermal_speed, alpha, k_landau
    type(sll_lt_pic_4d_group), pointer, intent(inout)       :: p_group

    call write_landau_density_on_remap_grid( thermal_speed, alpha, k_landau, p_group)
!    call plot_2d_slice_remapping_grid("init_values_on_rg.dat", p_group )
    call sll_lt_pic_4d_compute_new_particles( p_group )
    
  end subroutine sll_lt_pic_4d_init_landau

  ! initialize the particle density with a tensor product hat function with max value f_max attained at (x0,y0,vx0,vy0)
!  subroutine sll_init_lt_particles_4d_hat_f (           &          ! old name
  subroutine sll_lt_pic_4d_init_hat_f (           &
             x0,y0,vx0,vy0,r_x,r_y,r_vx,r_vy, basis_height, hat_shift,    &
             p_group )

    sll_real64, intent(in)                                  :: x0,y0,vx0,vy0,r_x,r_y,r_vx,r_vy, basis_height, hat_shift
    type(sll_lt_pic_4d_group), pointer, intent(inout)  :: p_group

    ! [[write_hat_density_on_remap_grid]]
    call write_hat_density_on_remap_grid( x0,y0,vx0,vy0,r_x,r_y,r_vx,r_vy, basis_height, hat_shift, p_group )
!    call plot_2d_slice_remapping_grid("init_values_on_rg.dat", p_group )
    call sll_lt_pic_4d_compute_new_particles( p_group )
    
  end subroutine sll_lt_pic_4d_init_hat_f



  subroutine sll_lt_pic_4d_remap( p_group )

    type(sll_lt_pic_4d_group), pointer, intent(inout) :: p_group

    ! n_virtual arguments: here to define the virtual cells where the backward flow is linearized
    sll_int32 :: n_virtual_x_for_remapping
    sll_int32 :: n_virtual_y_for_remapping
    sll_int32 :: n_virtual_vx_for_remapping
    sll_int32 :: n_virtual_vy_for_remapping

    ! 1st try: take just 1...
    n_virtual_x_for_remapping = 1
    n_virtual_y_for_remapping = 1
    n_virtual_vx_for_remapping = 1
    n_virtual_vy_for_remapping = 1

    call sll_lt_pic_4d_write_f_on_remapping_grid( p_group,                        &
                                                  n_virtual_x_for_remapping,      &
                                                  n_virtual_y_for_remapping,      &
                                                  n_virtual_vx_for_remapping,     &
                                                  n_virtual_vy_for_remapping)
    !    call sll_lt_pic_4d_write_f_on_remap_grid( p_group )

    call sll_lt_pic_4d_compute_new_particles( p_group )
    
  end subroutine sll_lt_pic_4d_remap



  subroutine write_landau_density_on_remap_grid(    &
              thermal_speed, alpha, k_landau,       &
              p_group                               &
              )

    sll_real64, intent(in)                                  :: thermal_speed, alpha, k_landau
    type(sll_lt_pic_4d_group), pointer, intent(inout)       :: p_group
    sll_int32 :: j_x
    sll_int32 :: j_y
    sll_int32 :: j_vx
    sll_int32 :: j_vy
    sll_int32 :: number_particles
    sll_int32 :: number_parts_x
    sll_int32 :: number_parts_y
    sll_int32 :: number_parts_vx
    sll_int32 :: number_parts_vy
    sll_real64 :: h_parts_x    
    sll_real64 :: h_parts_y    
    sll_real64 :: h_parts_vx   
    sll_real64 :: h_parts_vy   
    sll_real64 :: parts_x_min  
    sll_real64 :: parts_y_min  
    sll_real64 :: parts_vx_min 
    sll_real64 :: parts_vy_min 
    sll_real64 :: one_over_thermal_velocity
    sll_real64 :: one_over_two_pi
    sll_real64 :: x_j
    sll_real64 :: y_j
    sll_real64 :: vx_j
    sll_real64 :: vy_j
    type(sll_cartesian_mesh_2d),      pointer  :: m2d
    sll_real64 :: f_x, f_y, f_vx, f_vy

    number_particles = p_group%number_particles
    one_over_thermal_velocity = 1./thermal_speed   
    one_over_two_pi = 1./(2*sll_pi)
    
    number_parts_x  = p_group%number_parts_x
    number_parts_y  = p_group%number_parts_y
    number_parts_vx = p_group%number_parts_vx
    number_parts_vy = p_group%number_parts_vy

    h_parts_x    = p_group%remapping_grid%delta_eta1
    h_parts_y    = p_group%remapping_grid%delta_eta2
    h_parts_vx   = p_group%remapping_grid%delta_eta3
    h_parts_vy   = p_group%remapping_grid%delta_eta4

    parts_x_min    = p_group%remapping_grid%eta1_min
    parts_y_min    = p_group%remapping_grid%eta2_min
    parts_vx_min   = p_group%remapping_grid%eta3_min
    parts_vy_min   = p_group%remapping_grid%eta4_min

    ! Poisson mesh associated to the particles
    m2d => p_group%mesh

    ! compute the values of f0 on the (cartesian, phase-space) remapping grid
    x_j = parts_x_min
    do j_x = 1, number_parts_x
      f_x = eval_landau_fx(alpha, k_landau, x_j)
      y_j = parts_y_min
      do j_y = 1, number_parts_y
        vx_j = parts_vx_min
        do j_vx = 1, number_parts_vx
          f_vx = one_over_thermal_velocity * exp(-0.5*(vx_j*one_over_thermal_velocity)**2)
          vy_j = parts_vy_min
          do j_vy = 1, number_parts_vy
            f_vy = one_over_thermal_velocity * exp(-0.5*(vy_j*one_over_thermal_velocity)**2)
            p_group%target_values(j_x,j_y,j_vx,j_vy) = one_over_two_pi * f_x * f_vx * f_vy
            vy_j = vy_j + h_parts_vy
          end do
          vx_j = vx_j + h_parts_vx
        end do
        y_j = y_j + h_parts_y
      end do
      x_j = x_j + h_parts_x
    end do

  end subroutine write_landau_density_on_remap_grid
    
  function eval_hat_function(x0,y0,vx0,vy0,r_x,r_y,r_vx,r_vy, basis_height, hat_shift, x, y, vx, vy)
    sll_real64 :: x0,y0,vx0,vy0             ! centers of the hat
    sll_real64 :: r_x,r_y,r_vx,r_vy         ! radii of the hat
    sll_real64 :: basis_height, hat_shift
    sll_real64 :: x, y, vx, vy
    sll_real64 :: eval_hat_function
    
    sll_real64 :: inv_r_x 
    sll_real64 :: inv_r_y 
    sll_real64 :: inv_r_vx
    sll_real64 :: inv_r_vy
    
    inv_r_x  = 1./r_x
    inv_r_y  = 1./r_y
    inv_r_vx = 1./r_vx
    inv_r_vy = 1./r_vy

    eval_hat_function = basis_height + hat_shift * max(0._f64, 1. - inv_r_x*abs(x-x0) )             &
                                                 * max(0._f64, 1. - inv_r_y*abs(y-y0) )             &
                                                 * max(0._f64, 1. - inv_r_vx*abs(vx-vx0) )          &
                                                 * max(0._f64, 1. - inv_r_vy*abs(vy-vy0) )
  end function eval_hat_function
  
  

  ! <<write_hat_density_on_remap_grid>>
  subroutine write_hat_density_on_remap_grid ( &
        x0,y0,vx0,vy0,          &
        r_x,r_y,r_vx,r_vy,      &
        basis_height, hat_shift,&
        p_group                 &
      )

    sll_real64, intent(in)                                  :: x0,y0,vx0,vy0
    sll_real64, intent(in)                                  :: r_x,r_y,r_vx,r_vy
    sll_real64, intent(in)                                  :: basis_height, hat_shift
    type(sll_lt_pic_4d_group), pointer, intent(inout)       :: p_group
    sll_int32 :: j_x
    sll_int32 :: j_y
    sll_int32 :: j_vx
    sll_int32 :: j_vy
    sll_int32 :: number_particles
    sll_int32 :: number_parts_x
    sll_int32 :: number_parts_y
    sll_int32 :: number_parts_vx
    sll_int32 :: number_parts_vy
    sll_real64 :: h_parts_x    
    sll_real64 :: h_parts_y    
    sll_real64 :: h_parts_vx   
    sll_real64 :: h_parts_vy   
    sll_real64 :: parts_x_min  
    sll_real64 :: parts_y_min  
    sll_real64 :: parts_vx_min 
    sll_real64 :: parts_vy_min 
    sll_real64 :: x_j
    sll_real64 :: y_j
    sll_real64 :: vx_j
    sll_real64 :: vy_j
!    sll_real64 :: inv_r_x
!    sll_real64 :: inv_r_y
!    sll_real64 :: inv_r_vx
!    sll_real64 :: inv_r_vy
    type(sll_cartesian_mesh_2d),      pointer  :: m2d
    sll_real64 :: f_x, f_y, f_vx, f_vy

    number_particles = p_group%number_particles
    
!    inv_r_x  = 1./r_x
!    inv_r_y  = 1./r_y
!    inv_r_vx = 1./r_vx
!    inv_r_vy = 1./r_vy
!    
    number_parts_x  = p_group%number_parts_x
    number_parts_y  = p_group%number_parts_y
    number_parts_vx = p_group%number_parts_vx
    number_parts_vy = p_group%number_parts_vy

    h_parts_x    = p_group%remapping_grid%delta_eta1
    h_parts_y    = p_group%remapping_grid%delta_eta2
    h_parts_vx   = p_group%remapping_grid%delta_eta3
    h_parts_vy   = p_group%remapping_grid%delta_eta4

    parts_x_min    = p_group%remapping_grid%eta1_min
    parts_y_min    = p_group%remapping_grid%eta2_min
    parts_vx_min   = p_group%remapping_grid%eta3_min
    parts_vy_min   = p_group%remapping_grid%eta4_min

    ! Poisson mesh associated to the particles
    m2d => p_group%mesh

    ! compute the values of f0 on the (cartesian, phase-space) remapping grid
    x_j = parts_x_min
    do j_x = 1, number_parts_x
      y_j = parts_y_min
      do j_y = 1, number_parts_y
        vx_j = parts_vx_min
        do j_vx = 1, number_parts_vx
          vy_j = parts_vy_min
          do j_vy = 1, number_parts_vy
            p_group%target_values(j_x,j_y,j_vx,j_vy) = eval_hat_function(x0,y0,vx0,vy0,r_x,r_y,r_vx,r_vy, &
                                                                         basis_height, hat_shift,         &
                                                                         x_j, y_j, vx_j, vy_j)
            vy_j = vy_j + h_parts_vy
          end do
          vx_j = vx_j + h_parts_vx
        end do
        y_j = y_j + h_parts_y
      end do
      x_j = x_j + h_parts_x
    end do

  end subroutine write_hat_density_on_remap_grid
        

  ! position the particle on the cartesian remapping grid
  ! and compute new weights in order to approximate the point values stored in the remapping grid
  !  subroutine sll_compute_new_lt_particles_4d( &        ! old name
  subroutine sll_lt_pic_4d_compute_new_particles( &
              p_group )

    type(sll_lt_pic_4d_group), pointer, intent(inout) :: p_group
    sll_int32 :: k, k_ngb
    sll_int32 :: k_temp_debug
    sll_int32 :: j_x
    sll_int32 :: j_y
    sll_int32 :: j_vx
    sll_int32 :: j_vy
    sll_int32 :: j_aux_x
    sll_int32 :: j_aux_y
    sll_int32 :: j_aux_vx
    sll_int32 :: j_aux_vy
    sll_int32 :: l_x
    sll_int32 :: l_y
    sll_int32 :: l_vx
    sll_int32 :: l_vy
    sll_int32 :: ierr
    sll_int32 :: number_parts_x
    sll_int32 :: number_parts_y
    sll_int32 :: number_parts_vx
    sll_int32 :: number_parts_vy
    sll_real64 :: h_parts_x    
    sll_real64 :: h_parts_y    
    sll_real64 :: h_parts_vx   
    sll_real64 :: h_parts_vy   
    sll_real64 :: parts_x_min  
    sll_real64 :: parts_y_min  
    sll_real64 :: parts_vx_min 
    sll_real64 :: parts_vy_min 
    sll_real32 :: w_k
    sll_real64 :: x_j
    sll_real64 :: y_j
    sll_real64 :: vx_j
    sll_real64 :: vy_j
    sll_real32 :: d_vol
    type(sll_cartesian_mesh_2d),      pointer  :: m2d
    sll_int32,  dimension(:,:,:,:), pointer  :: particle_indices    !  why pointer ?
    
    number_parts_x  = p_group%number_parts_x
    number_parts_y  = p_group%number_parts_y
    number_parts_vx = p_group%number_parts_vx
    number_parts_vy = p_group%number_parts_vy

    h_parts_x    = p_group%remapping_grid%delta_eta1
    h_parts_y    = p_group%remapping_grid%delta_eta2
    h_parts_vx   = p_group%remapping_grid%delta_eta3
    h_parts_vy   = p_group%remapping_grid%delta_eta4
    
    d_vol = real( h_parts_x*h_parts_y*h_parts_vx*h_parts_vy,f32)
    
    parts_x_min    = p_group%remapping_grid%eta1_min
    parts_y_min    = p_group%remapping_grid%eta2_min
    parts_vx_min   = p_group%remapping_grid%eta3_min
    parts_vy_min   = p_group%remapping_grid%eta4_min

    ! Poisson mesh associated to the particles
    m2d => p_group%mesh
    
    SLL_ALLOCATE( particle_indices(number_parts_x, number_parts_y, number_parts_vx, number_parts_vy), ierr )
    particle_indices(:,:,:,:) = 0

    ! compute the particle weights from the values of f0 on the (cartesian, phase-space) remapping grid
    k_temp_debug = 0
    x_j = parts_x_min
    do j_x = 1, number_parts_x
      y_j = parts_y_min
      do j_y = 1, number_parts_y
        vx_j = parts_vx_min
        do j_vx = 1, number_parts_vx
          vy_j = parts_vy_min
          do j_vy = 1, number_parts_vy
            
            k_temp_debug = k_temp_debug + 1
            call get_particle_index_from_initial_position_on_cartesian_grid(            &
                j_x, j_y, j_vx, j_vy,                                                   &
                number_parts_x, number_parts_y, number_parts_vx, number_parts_vy,       &
                k                                                                       &
            )

            SLL_ASSERT(k == k_temp_debug)

            if( p_group%spline_degree == 1 )then
                w_k = d_vol * real( p_group%target_values(j_x,j_y,j_vx,j_vy) ,f32)
            else if( p_group%spline_degree == 3 )then
                w_k = 0
                do l_x = -1, 1
                    j_aux_x = j_x + l_x
                    if( p_group%domain_is_x_periodic )then
                        if( j_aux_x < 1 ) j_aux_x = j_aux_x + number_parts_x
                        if( j_aux_x > number_parts_x ) j_aux_x = j_aux_x - number_parts_x
                    end if
                    if( j_aux_x >= 1 .and. j_aux_x <= number_parts_x )then
                        do l_y = -1, 1
                            j_aux_y = j_y + l_y
                            if( p_group%domain_is_y_periodic )then
                                if( j_aux_y < 1 ) j_aux_y = j_aux_y + number_parts_y
                                if( j_aux_y > number_parts_y ) j_aux_y = j_aux_y - number_parts_y
                            end if
                            if( j_aux_y >= 1 .and. j_aux_y <= number_parts_y )then
                                do l_vx = -1, 1
                                    j_aux_vx = j_vx + l_vx
                                    if( j_aux_vx >= 1 .and. j_aux_vx <= number_parts_vx )then
                                        do l_vy = -1, 1
                                            j_aux_vy = j_vy + l_vy
                                            if( j_aux_vy >= 1 .and. j_aux_vy <= number_parts_vy )then
                                                ! MCP: here by discarding outside loop instances we assume 
                                                ! that non-periodic bc = zero bc. We should instead
                                                ! keep those loop instances and use the specified bounddary condition when 
                                                ! the node is in some "fat boundary zone" )
                                                w_k = w_k +         &
                                                    real(           &
                                                    p_group%ltpic_interpolation_coefs(l_x )  *  &
                                                    p_group%ltpic_interpolation_coefs(l_y )  *  &
                                                    p_group%ltpic_interpolation_coefs(l_vx)  *  &
                                                    p_group%ltpic_interpolation_coefs(l_vy)  *  &
                                                    p_group%target_values(j_aux_x,j_aux_y,j_aux_vx,j_aux_vy) ,f32)
                                            end if
                                        end do
                                    end if                        
                                end do
                            end if
                        end do
                    end if
                end do
                w_k = d_vol * w_k
            else
               print *, 'sll_lt_pic_initialize_some4Dfunction(): ERROR, value of p_group%spline_degree ', &
                        ' is invalid: ', p_group%spline_degree
               STOP
            end if

            p_group%p_list(k)%q = w_k
            particle_indices(j_x,j_y,j_vx,j_vy) = k

            call global_to_cell_offset_extended(    &
                    x_j, y_j, &
                    m2d,      &
                    p_group%p_list(k)%ic_x, &
                    p_group%p_list(k)%ic_y, &
                    p_group%p_list(k)%dx, &
                    p_group%p_list(k)%dy )

            !            call global_to_cell_offset (  &
            !                 x_j, y_j, &
            !                 m2d, &
            !                 p_group%p_list(k)%ic, &
            !                 p_group%p_list(k)%dx, &
            !                 p_group%p_list(k)%dy )
            p_group%p_list(k)%vx = vx_j 
            p_group%p_list(k)%vy = vy_j

            ! set the particle connectivity
            if(j_x == 1)then
                ! [neighbor index = own index] means: no neighbor
                ! in the x-periodic case this will be changed when dealing the last particle in the x dimension
                p_group%p_list(k)%ngb_xleft_index = k
            else
                ! set the connectivity (in both directions) with left neighbor
                k_ngb = particle_indices(j_x-1,j_y,j_vx,j_vy)
                p_group%p_list(k)%ngb_xleft_index = k_ngb
                p_group%p_list(k_ngb)%ngb_xright_index = k
                if(j_x == number_parts_x)then
                    if( p_group%domain_is_x_periodic )then
                        ! set the connectivity (in both directions) with right neighbor
                        k_ngb = particle_indices(1,j_y,j_vx,j_vy)
                        p_group%p_list(k)%ngb_xright_index = k_ngb
                        p_group%p_list(k_ngb)%ngb_xleft_index = k
                    else
                        ! [neighbor index = own index] means: no neighbor
                        p_group%p_list(k)%ngb_xright_index = k
                    end if
                end if
            end if
            if(j_y == 1)then
                ! [neighbor index = own index] means: no neighbor
                ! in the y-periodic case this will be changed when dealing the last particle in the y dimension
                p_group%p_list(k)%ngb_yleft_index = k
            else
                ! set the connectivity (in both directions) with left neighbor
                k_ngb = particle_indices(j_x,j_y-1,j_vx,j_vy)
                p_group%p_list(k)%ngb_yleft_index = k_ngb
                p_group%p_list(k_ngb)%ngb_yright_index = k
                if(j_y == number_parts_y)then
                    if( p_group%domain_is_y_periodic )then
                        ! set the connectivity (in both directions) with right neighbor
                        k_ngb = particle_indices(j_x,1,j_vx,j_vy)
                        p_group%p_list(k)%ngb_yright_index = k_ngb
                        p_group%p_list(k_ngb)%ngb_yleft_index = k
                    else
                        ! [neighbor index = own index] means: no neighbor
                        p_group%p_list(k)%ngb_yright_index = k
                    end if
                end if
            end if
            if(j_vx == 1)then
                ! [neighbor index = own index] means: no neighbor
                p_group%p_list(k)%ngb_vxleft_index = k
            else
                ! set the connectivity (in both directions) with left neighbor
                k_ngb = particle_indices(j_x,j_y,j_vx-1,j_vy)
                p_group%p_list(k)%ngb_vxleft_index = k_ngb
                p_group%p_list(k_ngb)%ngb_vxright_index = k
                if(j_vx == number_parts_vx)then
                    ! [neighbor index = own index] means: no neighbor
                    p_group%p_list(k)%ngb_vxright_index = k
                end if
            end if
            if(j_vy == 1)then
                ! [neighbor index = own index] means: no neighbor
                p_group%p_list(k)%ngb_vyleft_index = k
            else
                ! set the connectivity (in both directions) with left neighbor
                k_ngb = particle_indices(j_x,j_y,j_vx,j_vy-1)
                p_group%p_list(k)%ngb_vyleft_index = k_ngb
                p_group%p_list(k_ngb)%ngb_vyright_index = k
                if(j_vy == number_parts_vy)then
                    ! [neighbor index = own index] means: no neighbor
                    p_group%p_list(k)%ngb_vyright_index = k
                end if
            end if
            
            vy_j = vy_j + h_parts_vy          
          end do
          vx_j = vx_j + h_parts_vx
        end do
        y_j = y_j + h_parts_y
      end do
      x_j = x_j + h_parts_x
    end do

  end subroutine sll_lt_pic_4d_compute_new_particles
    

  function eval_landau_fx(alpha, kx, x)
    sll_real64 :: alpha, kx, x
    sll_real64 :: eval_landau_fx
    eval_landau_fx = 1._f64 + alpha * cos(kx * x)
  end function eval_landau_fx

end module sll_lt_pic_4d_init
