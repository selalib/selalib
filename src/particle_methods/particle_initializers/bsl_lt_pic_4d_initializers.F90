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


!!!!!!!!!!!!!!!!!!!           THIS MODULE DOES NOT SEEM TO BE USED...

module sll_bsl_lt_pic_4d_initializers
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"
!#include "particle_representation.h"

  use sll_constants, only: sll_pi
  use sll_cartesian_meshes
!  use sll_bsl_lt_pic_4d_group_module
  use gaussian

  use sll_bsl_lt_pic_4d_utilities_module
  use sll_cartesian_meshes
!  use sll_representation_conversion_module -- better not use this

  implicit none


contains




  ! initialize the bsl_lt_pic group with the landau f0 distribution
!  subroutine sll_bsl_lt_pic_4d_initializer_landau_f0 (       &
!              thermal_speed, alpha, k_landau,   &
!              p_group)
!
!    sll_real64, intent(in)                                  :: thermal_speed, alpha, k_landau
!    type(sll_bsl_lt_pic_4d_group), pointer, intent(inout)       :: p_group
!
!    call write_landau_density_on_remap_grid( thermal_speed, alpha, k_landau, p_group)
!    call sll_bsl_lt_pic_4d_compute_new_particles( p_group )
!
!  end subroutine sll_bsl_lt_pic_4d_initializer_landau_f0

  ! initialize the particle density with a tensor product hat function with max value f_max attained at (x0,y0,vx0,vy0)
  subroutine bsl_lt_pic_4d_init_hat_f (           &
             x0,y0,vx0,vy0,r_x,r_y,r_vx,r_vy, basis_height, hat_shift,    &
              )

    sll_real64, intent(in)                                  :: x0,y0,vx0,vy0,r_x,r_y,r_vx,r_vy, basis_height, hat_shift
    type(sll_bsl_lt_pic_4d_group), pointer, intent(inout)  :: p_group

    ! [[write_hat_density_on_remap_grid]]
    call write_hat_density_on_remap_grid( x0,y0,vx0,vy0,r_x,r_y,r_vx,r_vy, basis_height, hat_shift, p_group )
!    call plot_2d_slice_remapping_grid("init_values_on_rg.dat", p_group )
    call sll_bsl_lt_pic_4d_compute_new_particles( p_group )

  end subroutine bsl_lt_pic_4d_init_hat_f

  ! todo:  this has nothing to do with the initial f0, put it somewhere else (utilities?)
  subroutine sll_bsl_lt_pic_4d_remap( p_group )

    type(sll_bsl_lt_pic_4d_group), pointer, intent(inout) :: p_group

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

    call sll_bsl_lt_pic_4d_write_f_on_remapping_grid( p_group,                        &
                                                  n_virtual_x_for_remapping,      &
                                                  n_virtual_y_for_remapping,      &
                                                  n_virtual_vx_for_remapping,     &
                                                  n_virtual_vy_for_remapping)

    call sll_bsl_lt_pic_4d_compute_new_particles( p_group )
    
  end subroutine sll_bsl_lt_pic_4d_remap



    
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
    type(sll_bsl_lt_pic_4d_group), pointer, intent(inout)       :: p_group
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
    m2d => p_group%space_mesh_2d

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
        


end module sll_bsl_lt_pic_4d_initializers
