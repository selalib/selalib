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

!module sll_lt_pic_utilities        ! old name
module sll_lt_pic_4d_utilities
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"
#include "sll_accumulators.h" 
#include "particle_representation.h"
  use sll_lt_pic_4d_group_module
  use sll_representation_conversion_module
  use sll_timer

implicit none 

contains

  subroutine sll_first_lt_pic_charge_accumulation_4d( p_group, q_accum )
    type(sll_lt_pic_4d_group), pointer      :: p_group
    type(sll_charge_accumulator_2d), pointer     :: q_accum
    type(sll_lt_pic_4d_particle), dimension(:), pointer :: p
    sll_int64  :: i
    sll_int64  :: number_particles
    sll_real64 :: tmp1
    sll_real64 :: tmp2

    SLL_ASSERT( associated(p_group) .and. associated(q_accum))
    number_particles =  p_group%number_particles
    p             => p_group%p_list

    do i=1,number_particles
       SLL_ACCUMULATE_PARTICLE_CHARGE(q_accum,p(i),tmp1,tmp2)
    end do
  end subroutine sll_first_lt_pic_charge_accumulation_4d


  ! <<get_ltp_deformation_matrix>>
  !> Compute the coefficients of the particle 'deformation' matrix
  !! which approximates the Jacobian matrix of the exact backward flow
  !! (defined between the current time and that of the particle initialization -- or the last particle remapping)
  !! at the current particle position.
  !! It is computed in two steps:
  !!    * first we compute a finite-difference approximation of the forward Jacobian matrix by using the fact that the
  !!      initial (or remapped) particles were located on a cartesian grid in phase space. Specifically, the entries of
  !!      the forward Jacobian matrix (say, (J^n_k)_xy = d(F^n_x)/dy (z^0_k) -- with z is the phase space coordinate)
  !!      are approximated by finite differences: here, by  DF / 2*h_y
  !!      with DF = F^n_x(z^0_k + (0,h_y,0,0)) - F^n_x(z^0_k - (0,h_y,0,0))
  !!    * then we invert that approximated forward Jacobian matrix
  !!
  !! Note: when computing finite differences in a periodic dimension (say x), one must be careful since two values of F_x
  !!    can be close in the periodic interval but distant by almost L_x in the (stored) [0,L_x[ representation.
  !!    To account for this (frequent) phenomenon we do the following:
  !!    when the difference DF (see example above) is larger than L_x/2, we make it smaller by adding +/- L_x to it.
  !!    Note that here we could very well be making the slope smaller than what it should actually be: indeed if the function
  !!    F^n_x is having a steep slope at z^0_k which adds one (or several) periods L_x to DF then our modification will
  !!    artificially lower the slope. But this is impossible to prevent in full generality (indeed: a steep slope crossing the
  !!    0 or L_x value will be lowered anyway in the [0,L_x[ representation) and should not be frequent (indeed it only happens
  !!    when F^n_x has high derivatives and the particle grid is coarse, which will lead to bad approximations anyhow).

  subroutine get_ltp_deformation_matrix (       &
                        k,                      &    
                        particles_m2d,          &
                        p_list,                 &
                        domain_is_x_periodic,   &
                        domain_is_y_periodic,   &
                        mesh_period_x,          &
                        mesh_period_y,          &
                        h_parts_x,              &    
                        h_parts_y,              &    
                        h_parts_vx,             &   
                        h_parts_vy,             &   
                        inv_h_parts_x,          &    
                        inv_h_parts_y,          &    
                        inv_h_parts_vx,         &   
                        inv_h_parts_vy,         &   
                        ref_radius,             &
                        x_k, y_k,               &
                        vx_k, vy_k,             &
                        d11,d12,d13,d14,        &
                        d21,d22,d23,d24,        &
                        d31,d32,d33,d34,        &
                        d41,d42,d43,d44,        &
                        part_radius_x,          &
                        part_radius_y,          &
                        part_radius_vx,         &
                        part_radius_vy          &
                        )                            
    
        sll_int64, intent(in) :: k
!        sll_int32, intent(in) :: part_degree
        type(sll_cartesian_mesh_2d),               pointer,  intent(in) :: particles_m2d  ! Poisson mesh associated with the particles
        type(sll_lt_pic_4d_particle), dimension(:),  pointer,  intent(in) :: p_list
    
        LOGICAL, intent(in) :: domain_is_x_periodic
        LOGICAL, intent(in) :: domain_is_y_periodic    
        sll_real64, intent(in)  :: mesh_period_x
        sll_real64, intent(in)  :: mesh_period_y
        
        sll_real64, intent(in)  :: h_parts_x    
        sll_real64, intent(in)  :: h_parts_y    
        sll_real64, intent(in)  :: h_parts_vx   
        sll_real64, intent(in)  :: h_parts_vy   
        sll_real64, intent(in)  :: inv_h_parts_x    
        sll_real64, intent(in)  :: inv_h_parts_y    
        sll_real64, intent(in)  :: inv_h_parts_vx   
        sll_real64, intent(in)  :: inv_h_parts_vy   
        sll_real64, intent(in)  :: ref_radius      

        sll_real64, intent(out) :: x_k, y_k         ! particle center in physical space
        sll_real64, intent(out) :: vx_k, vy_k       ! particle center in velocity space
        sll_real64, intent(out) :: d11,d12,d13,d14  ! coefs of matrix D (backward Jacobian)
        sll_real64, intent(out) :: d21,d22,d23,d24
        sll_real64, intent(out) :: d31,d32,d33,d34
        sll_real64, intent(out) :: d41,d42,d43,d44
        sll_real64, intent(out) :: part_radius_x    ! radius of particle support, in x dimension
        sll_real64, intent(out) :: part_radius_y    ! radius of particle support, in y dimension
        sll_real64, intent(out) :: part_radius_vx   ! radius of particle support, in vx dimension
        sll_real64, intent(out) :: part_radius_vy   ! radius of particle support, in vy dimension
    
        sll_int64   :: k_ngb
        sll_real64  :: x_k_left,  x_k_right
        sll_real64  :: y_k_left,  y_k_right
        sll_real64  :: vx_k_left, vx_k_right
        sll_real64  :: vy_k_left, vy_k_right
    
        sll_real64  :: j11,j12,j13,j14   ! coefs of matrix J = D^-1 (forward Jacobian)
        sll_real64  :: j21,j22,j23,j24
        sll_real64  :: j31,j32,j33,j34
        sll_real64  :: j41,j42,j43,j44
        sll_real64  :: factor, det_J, inv_det_J
            
        
        call cell_offset_to_global( p_list(k)%dx, &
                                    p_list(k)%dy, &
                                    p_list(k)%ic, &
                                    particles_m2d, x_k, y_k )
        vx_k   = p_list(k)%vx
        vy_k   = p_list(k)%vy
!        w_k    = p_list(k)%q           

        ! Compute the forward Jacobian matrix J_k        

        ! ------   d/d_x terms 
        factor = 0.5*inv_h_parts_x        

        k_ngb  = p_list(k)%ngb_xright_index
        if( k_ngb == k )then           
           ! no right neighbor is available, use a non-centered finite difference
           factor = 2*factor
           x_k_right   = x_k
           y_k_right   = y_k
           vx_k_right  = vx_k
           vy_k_right  = vy_k
        else
            call cell_offset_to_global( p_list(k_ngb)%dx, &
                                    p_list(k_ngb)%dy, &
                                    p_list(k_ngb)%ic, &
                                    particles_m2d, x_k_right, y_k_right )
            vx_k_right = p_list(k_ngb)%vx
            vy_k_right = p_list(k_ngb)%vy        
            if( domain_is_x_periodic .and. x_k_right < x_k - 0.5*mesh_period_x ) x_k_right = x_k_right + mesh_period_x
            if( domain_is_x_periodic .and. x_k_right > x_k + 0.5*mesh_period_x ) x_k_right = x_k_right - mesh_period_x
            if( domain_is_y_periodic .and. y_k_right < y_k - 0.5*mesh_period_y ) y_k_right = y_k_right + mesh_period_y
            if( domain_is_y_periodic .and. y_k_right > y_k + 0.5*mesh_period_y ) y_k_right = y_k_right - mesh_period_y
        end if
        

        k_ngb  = p_list(k)%ngb_xleft_index
        if( k_ngb == k )then           
           ! no left neighbor is available, use a non-centered finite difference
           factor = 2*factor
           x_k_left   = x_k
           y_k_left   = y_k
           vx_k_left  = vx_k
           vy_k_left  = vy_k
        else
            call cell_offset_to_global( p_list(k_ngb)%dx, &
                                        p_list(k_ngb)%dy, &
                                        p_list(k_ngb)%ic, &
                                        particles_m2d, x_k_left, y_k_left )
            vx_k_left  = p_list(k_ngb)%vx
            vy_k_left  = p_list(k_ngb)%vy
            if( domain_is_x_periodic .and. x_k_left < x_k - 0.5*mesh_period_x ) x_k_left = x_k_left + mesh_period_x
            if( domain_is_x_periodic .and. x_k_left > x_k + 0.5*mesh_period_x ) x_k_left = x_k_left - mesh_period_x
            if( domain_is_y_periodic .and. y_k_left < y_k - 0.5*mesh_period_y ) y_k_left = y_k_left + mesh_period_y
            if( domain_is_y_periodic .and. y_k_left > y_k + 0.5*mesh_period_y ) y_k_left = y_k_left - mesh_period_y
        end if
        
           j11 = factor * ( x_k_right  - x_k_left  )
           j21 = factor * ( y_k_right  - y_k_left  )
           j31 = factor * ( vx_k_right - vx_k_left )
           j41 = factor * ( vy_k_right - vy_k_left )

        ! ------   d/d_y terms 
        factor = 0.5*inv_h_parts_y

        k_ngb  = p_list(k)%ngb_yright_index
        if( k_ngb == k )then           
           ! no right neighbor is available, use a non-centered finite difference
           factor = 2*factor
           x_k_right   = x_k
           y_k_right   = y_k
           vx_k_right  = vx_k
           vy_k_right  = vy_k
        else
            call cell_offset_to_global(  p_list(k_ngb)%dx, &
                                        p_list(k_ngb)%dy, &
                                        p_list(k_ngb)%ic, &
                                        particles_m2d, x_k_right, y_k_right )
            vx_k_right = p_list(k_ngb)%vx
            vy_k_right = p_list(k_ngb)%vy
            if( domain_is_x_periodic .and. x_k_right < x_k - 0.5*mesh_period_x ) x_k_right = x_k_right + mesh_period_x
            if( domain_is_x_periodic .and. x_k_right > x_k + 0.5*mesh_period_x ) x_k_right = x_k_right - mesh_period_x
            if( domain_is_y_periodic .and. y_k_right < y_k - 0.5*mesh_period_y ) y_k_right = y_k_right + mesh_period_y
            if( domain_is_y_periodic .and. y_k_right > y_k + 0.5*mesh_period_y ) y_k_right = y_k_right - mesh_period_y
        end if

        k_ngb  = p_list(k)%ngb_yleft_index
        if( k_ngb == k )then           
           ! no left neighbor is available, use a non-centered finite difference
           factor = 2*factor
           x_k_left   = x_k
           y_k_left   = y_k
           vx_k_left  = vx_k
           vy_k_left  = vy_k
        else
            call cell_offset_to_global( p_list(k_ngb)%dx, &
                                        p_list(k_ngb)%dy, &
                                        p_list(k_ngb)%ic, &
                                        particles_m2d, x_k_left, y_k_left )
            vx_k_left  = p_list(k_ngb)%vx
            vy_k_left  = p_list(k_ngb)%vy
            if( domain_is_x_periodic .and. x_k_left < x_k - 0.5*mesh_period_x ) x_k_left = x_k_left + mesh_period_x
            if( domain_is_x_periodic .and. x_k_left > x_k + 0.5*mesh_period_x ) x_k_left = x_k_left - mesh_period_x
            if( domain_is_y_periodic .and. y_k_left < y_k - 0.5*mesh_period_y ) y_k_left = y_k_left + mesh_period_y
            if( domain_is_y_periodic .and. y_k_left > y_k + 0.5*mesh_period_y ) y_k_left = y_k_left - mesh_period_y
        end if
        
           j12 = factor * ( x_k_right  - x_k_left  )
           j22 = factor * ( y_k_right  - y_k_left  )
           j32 = factor * ( vx_k_right - vx_k_left )
           j42 = factor * ( vy_k_right - vy_k_left )


        ! ------   d/d_vx terms 
        factor = 0.5*inv_h_parts_vx

        k_ngb  = p_list(k)%ngb_vxright_index
        if( k_ngb == k )then           
           ! no right neighbor is available, use a non-centered finite difference
           factor = 2*factor
           x_k_right   = x_k
           y_k_right   = y_k
           vx_k_right  = vx_k
           vy_k_right  = vy_k
        else
           call cell_offset_to_global(  p_list(k_ngb)%dx, &
                                        p_list(k_ngb)%dy, &
                                        p_list(k_ngb)%ic, &
                                        particles_m2d, x_k_right, y_k_right )
            vx_k_right = p_list(k_ngb)%vx
            vy_k_right = p_list(k_ngb)%vy
            if( domain_is_x_periodic .and. x_k_right < x_k - 0.5*mesh_period_x ) x_k_right = x_k_right + mesh_period_x
            if( domain_is_x_periodic .and. x_k_right > x_k + 0.5*mesh_period_x ) x_k_right = x_k_right - mesh_period_x
            if( domain_is_y_periodic .and. y_k_right < y_k - 0.5*mesh_period_y ) y_k_right = y_k_right + mesh_period_y
            if( domain_is_y_periodic .and. y_k_right > y_k + 0.5*mesh_period_y ) y_k_right = y_k_right - mesh_period_y
        end if

        k_ngb  = p_list(k)%ngb_vxleft_index
        if( k_ngb == k )then           
            ! no left neighbor is available, use a non-centered finite difference
            factor = 2*factor
            x_k_left   = x_k
            y_k_left   = y_k
            vx_k_left  = vx_k
            vy_k_left  = vy_k
        else
            call cell_offset_to_global( p_list(k_ngb)%dx, &
                                        p_list(k_ngb)%dy, &
                                        p_list(k_ngb)%ic, &
                                        particles_m2d, x_k_left, y_k_left )
            vx_k_left  = p_list(k_ngb)%vx
            vy_k_left  = p_list(k_ngb)%vy
            if( domain_is_x_periodic .and. x_k_left < x_k - 0.5*mesh_period_x ) x_k_left = x_k_left + mesh_period_x
            if( domain_is_x_periodic .and. x_k_left > x_k + 0.5*mesh_period_x ) x_k_left = x_k_left - mesh_period_x
            if( domain_is_y_periodic .and. y_k_left < y_k - 0.5*mesh_period_y ) y_k_left = y_k_left + mesh_period_y
            if( domain_is_y_periodic .and. y_k_left > y_k + 0.5*mesh_period_y ) y_k_left = y_k_left - mesh_period_y
        end if
        
           j13 = factor * ( x_k_right  - x_k_left  )
           j23 = factor * ( y_k_right  - y_k_left  )
           j33 = factor * ( vx_k_right - vx_k_left )
           j43 = factor * ( vy_k_right - vy_k_left )


        ! ------   d/d_vy terms 
        factor = 0.5*inv_h_parts_vy

        k_ngb  = p_list(k)%ngb_vyright_index
        if( k_ngb == k )then           
           ! no right neighbor is available, use a non-centered finite difference
           factor = 2*factor
           x_k_right   = x_k
           y_k_right   = y_k
           vx_k_right  = vx_k
           vy_k_right  = vy_k
        else
           call cell_offset_to_global(  p_list(k_ngb)%dx, &
                                        p_list(k_ngb)%dy, &
                                        p_list(k_ngb)%ic, &
                                        particles_m2d, x_k_right, y_k_right )
            vx_k_right = p_list(k_ngb)%vx
            vy_k_right = p_list(k_ngb)%vy
            if( domain_is_x_periodic .and. x_k_right < x_k - 0.5*mesh_period_x ) x_k_right = x_k_right + mesh_period_x
            if( domain_is_x_periodic .and. x_k_right > x_k + 0.5*mesh_period_x ) x_k_right = x_k_right - mesh_period_x
            if( domain_is_y_periodic .and. y_k_right < y_k - 0.5*mesh_period_y ) y_k_right = y_k_right + mesh_period_y
            if( domain_is_y_periodic .and. y_k_right > y_k + 0.5*mesh_period_y ) y_k_right = y_k_right - mesh_period_y
        end if

        k_ngb  = p_list(k)%ngb_vyleft_index
        if( k_ngb == k )then           
            ! no left neighbor is available, use a non-centered finite difference
            factor = 2*factor
            x_k_left   = x_k
            y_k_left   = y_k
            vx_k_left  = vx_k
            vy_k_left  = vy_k
        else
            call cell_offset_to_global( p_list(k_ngb)%dx, &
                                        p_list(k_ngb)%dy, &
                                        p_list(k_ngb)%ic, &
                                        particles_m2d, x_k_left, y_k_left )
            vx_k_left  = p_list(k_ngb)%vx
            vy_k_left  = p_list(k_ngb)%vy
            if( domain_is_x_periodic .and. x_k_left < x_k - 0.5*mesh_period_x ) x_k_left = x_k_left + mesh_period_x
            if( domain_is_x_periodic .and. x_k_left > x_k + 0.5*mesh_period_x ) x_k_left = x_k_left - mesh_period_x
            if( domain_is_y_periodic .and. y_k_left < y_k - 0.5*mesh_period_y ) y_k_left = y_k_left + mesh_period_y
            if( domain_is_y_periodic .and. y_k_left > y_k + 0.5*mesh_period_y ) y_k_left = y_k_left - mesh_period_y
        end if
        
           j14 = factor * ( x_k_right  - x_k_left  )
           j24 = factor * ( y_k_right  - y_k_left  )
           j34 = factor * ( vx_k_right - vx_k_left )
           j44 = factor * ( vy_k_right - vy_k_left )
           
           
           ! Compute D_k the inverse of J_k
           
        det_J = j11 * j22 * j33 * j44 &
              + j11 * j23 * j34 * j42 &
               + j11 * j24 * j32 * j43 &
              + j12 * j21 * j34 * j43 &
              + j12 * j23 * j31 * j44 &
              + j12 * j24 * j33 * j41 &
              + j13 * j21 * j32 * j44 &
              + j13 * j22 * j34 * j41 &
              + j13 * j24 * j31 * j42 &
              + j14 * j21 * j33 * j42 &
              + j14 * j22 * j31 * j43 &
              + j14 * j23 * j32 * j41 &
              - j11 * j22 * j34 * j43 &
              - j11 * j23 * j32 * j44 &
              - j11 * j24 * j33 * j42 &
              - j12 * j21 * j33 * j44 &
              - j12 * j23 * j34 * j41 &
              - j12 * j24 * j31 * j43 &
              - j13 * j21 * j34 * j42 &
              - j13 * j22 * j31 * j44 &
              - j13 * j24 * j32 * j41 &
              - j14 * j21 * j32 * j43 &
              - j14 * j22 * j33 * j41 &
              - j14 * j23 * j31 * j42 
        !print*,'det  =',det    
        inv_det_J = 1./det_J
        
        d11 = j22 * j33 * j44 &
            + j23 * j34 * j42 &
            + j24 * j32 * j43 &
            - j22 * j34 * j43 &
            - j23 * j32 * j44 &
            - j24 * j33 * j42
        d11 = inv_det_J * d11
                         
        d12 = j12 * j34 * j43 &
            + j13 * j32 * j44 &
            + j14 * j33 * j42 &
            - j12 * j33 * j44 &
            - j13 * j34 * j42 &
            - j14 * j32 * j43
        d12 = inv_det_J * d12
                             
        d13 = j12 * j23 * j44 &
            + j13 * j24 * j42 &
            + j14 * j22 * j43 &
            - j12 * j24 * j43 &
            - j13 * j22 * j44 &
            - j14 * j23 * j42 
        d13 = inv_det_J * d13    
                     
        d14 = j12 * j24 * j33 &
            + j13 * j22 * j34 &
            + j14 * j23 * j32 &
            - j12 * j23 * j34 &
            - j13 * j24 * j32 &
            - j14 * j22 * j33 
        d14 = inv_det_J * d14
                         
        d21 = j21 * j34 * j43 &
            + j23 * j31 * j44 &
            + j24 * j33 * j41 &
            - j21 * j33 * j44 &
            - j23 * j34 * j41 &
            - j24 * j31 * j43 
        d21 = inv_det_J * d21
                         
        d22 = j11 * j33 * j44 &
            + j13 * j34 * j41 &
            + j14 * j31 * j43 &
            - j11 * j34 * j43 &
            - j13 * j31 * j44 &
            - j14 * j33 * j41
        d22 = inv_det_J * d22
                         
        d23 = j11 * j24 * j43 &
            + j13 * j21 * j44 &
            + j14 * j23 * j41 &
            - j11 * j23 * j44 &
            - j13 * j24 * j41 &
            - j14 * j21 * j43
        d23 = inv_det_J * d23
                     
        d24 = j11 * j23 * j34 &
            + j13 * j24 * j31 &
            + j14 * j21 * j33 &
            - j11 * j24 * j33 &
            - j13 * j21 * j34 &
            - j14 * j23 * j31
        d24 = inv_det_J * d24
                         
        d31 = j21 * j32 * j44 &
            + j22 * j34 * j41 &
            + j24 * j31 * j42 &
            - j21 * j34 * j42 &
            - j22 * j31 * j44 &
            - j24 * j32 * j41
        d31 = inv_det_J * d31
                     
        d32 = j11 * j34 * j42 &
            + j12 * j31 * j44 &
            + j14 * j32 * j41 &
            - j11 * j32 * j44 &
            - j12 * j34 * j41 &
            - j14 * j31 * j42
        d32 = inv_det_J * d32
                     
        d33 = j11 * j22 * j44 &
            + j12 * j24 * j41 &
            + j14 * j21 * j42 &
            - j11 * j24 * j42 &
            - j12 * j21 * j44 &
            - j14 * j22 * j41
        d33 = inv_det_J * d33
                     
        d34 = j11 * j24 * j32 &
            + j12 * j21 * j34 &
            + j14 * j22 * j31 &
            - j11 * j22 * j34 &
            - j12 * j24 * j31 &
            - j14 * j21 * j32
        d34 = inv_det_J * d34
                     
        d41 = j21 * j33 * j42 &
            + j22 * j31 * j43 &
            + j23 * j32 * j41 &
            - j21 * j32 * j43 &
            - j22 * j33 * j41 &
            - j23 * j31 * j42
        d41 = inv_det_J * d41
                     
        d42 = j11 * j32 * j43 &
            + j12 * j33 * j41 &
            + j13 * j31 * j42 &
            - j11 * j33 * j42 &
            - j12 * j31 * j43 &
            - j13 * j32 * j41
        d42 = inv_det_J * d42
                     
        d43 = j11 * j23 * j42 &
            + j12 * j21 * j43 &
            + j13 * j22 * j41 &
            - j11 * j22 * j43 &
            - j12 * j23 * j41 &
            - j13 * j21 * j42
        d43 = inv_det_J * d43
                     
        d44 = j11 * j22 * j33 &
            + j12 * j23 * j31 &
            + j13 * j21 * j32 &
            - j11 * j23 * j32 &
            - j12 * j21 * j33 &
            - j13 * j22 * j31
        d44 = inv_det_J * d44
        
        !Compute an approximation of the support    
!        ref_radius = 0.5*(part_degree+1)
        part_radius_x  =  ref_radius *( abs(j11) *h_parts_x + abs(j12) * h_parts_y + abs(j13) * h_parts_vx + abs(j14) * h_parts_vy )
        part_radius_y  =  ref_radius *( abs(j21) *h_parts_x + abs(j22) * h_parts_y + abs(j23) * h_parts_vx + abs(j24) * h_parts_vy )
        part_radius_vx =  ref_radius *( abs(j31) *h_parts_x + abs(j32) * h_parts_y + abs(j33) * h_parts_vx + abs(j34) * h_parts_vy )
        part_radius_vy =  ref_radius *( abs(j41) *h_parts_x + abs(j42) * h_parts_y + abs(j43) * h_parts_vx + abs(j44) * h_parts_vy )        
      
      
!        print *, "*    -- det_J = ", det_J
!        print *, "**** -- part_radius_x = ", part_radius_x
        
end subroutine get_ltp_deformation_matrix

  ! <<in_bounds_periodic>> extracted from [[file:../simulation/simulation_4d_vp_lt_pic_cartesian.F90::function
  ! in_bounds]]

  function in_bounds_periodic( x, y, mesh, x_periodic, y_periodic ) result(res)

    use sll_cartesian_meshes

    ! [[file:../working_precision/sll_working_precision.h]]
    use sll_working_precision

    logical :: res
    sll_real64, intent(in) :: x
    sll_real64, intent(in) :: y
    logical, intent(in) :: x_periodic
    logical, intent(in) :: y_periodic
    type(sll_cartesian_mesh_2d), pointer :: mesh
    if((x >= mesh%eta1_min) &
         .and. ((x < mesh%eta1_max .and. x_periodic) .or. (x <= mesh%eta1_max .and. .not. x_periodic)) &
         .and. (y >= mesh%eta2_min) &
         .and. ((y < mesh%eta2_max .and. y_periodic) .or. (y <= mesh%eta2_max .and. .not. y_periodic))) then
       res = .true.
    else
       res = .false.
    end if
  end function in_bounds_periodic

  ! <<apply_periodic_bc>> extracted from [[file:../simulation/simulation_4d_vp_lt_pic_cartesian.F90::subroutine
  ! apply_periodic_bc]]

  subroutine apply_periodic_bc( mesh, x, y )

    use sll_cartesian_meshes

    ! [[file:../working_precision/sll_working_precision.h]]
    use sll_working_precision

    type(sll_cartesian_mesh_2d), pointer :: mesh
    sll_real64, intent(inout) :: x
    sll_real64, intent(inout) :: y
    sll_real64 :: xmin
    sll_real64 :: xmax
    sll_real64 :: ymin
    sll_real64 :: ymax

    xmin = mesh%eta1_min
    xmax = mesh%eta1_max
    ymin = mesh%eta2_min
    ymax = mesh%eta2_max
!    if( x < xmin ) x = x + xmax-xmin
!    if( x > xmax ) x = x - xmax-xmin
!    if( y < ymin ) y = y + ymax-ymin
!    if( y > ymax ) y = y - ymax-ymin
    do while( x < xmin ) 
      x = x + (xmax-xmin)
    end do
    do while( x >= xmax ) 
      x = x - (xmax-xmin)
    end do
    do while( y < ymin ) 
      y = y + (ymax-ymin)
    end do
    do while( y >= ymax )
      y = y - (ymax-ymin)
    end do

    ! and the condition that the particle is in-bounds should trigger some
    ! alarm as this would not be supposed to happen here!
  end subroutine apply_periodic_bc

  ! <<onestep>> <<ALH>> utility function for finding the neighbours of a particle, used by [[ONESTEPMACRO]]. "dim"
  ! corresponds to one of x,y,vx,vy.

  subroutine onestep(dim,dim_t0,kprime,p_list,h_parts_dim)

    sll_int :: dim
    sll_real64 :: dim_t0
    sll_int64 :: neighbour

    ! [[file:~/mcp/selalib/src/pic_particle_types/lt_pic_4d_group.F90::sll_lt_pic_4d_group-p_list]]
    type(sll_lt_pic_4d_particle), dimension(:), pointer,intent(in) :: p_list
    
    sll_int64 :: ngb_dim_right_index
    sll_int64 :: ngb_dim_left_index
    sll_real64 :: h_parts_dim
    sll_int64 :: kprime
    sll_int32 :: j,jumps

    ! <<up>> means that kprime needs to go up ie increase in coordinate
    
    logical :: up

    ! Move to a closer neighbour only if dim_t0 is not located in a cell of size h_parts_dim and with a left bound of
    ! dim_t0

    if(kprime == 0)return

    ! How many jumps do we need to do in that direction to reduce the distance 'dim_t0' to a minimum?

    jumps = int(abs(dim_t0/h_parts_dim))
    
    ! dim_t0 < 0 means that the virtual particle is at the left of kprime (dim_t0 is a relative coordinate). If dim_t0
    ! is negative, add one step to move to a positive relative signed distance between dim_t0 and kprime.
    
    up = .true.
    if(dim_t0 < 0) then
       jumps = jumps + 1
       up = .false.
    end if

    ! resulting signed distance between marker at t0 and kprime
    
    if(up) then
       dim_t0 = dim_t0 - jumps * h_parts_dim
    else
       dim_t0 = dim_t0 + jumps * h_parts_dim
    endif

    ! do as many jumps as required through the neighbour pointers in the given dimension (1:x, 2:y, 3:vx, 4:vy). kprime
    ! can become zero (ie fall out of the domain) in non-periodic dimensions.
    
    j = 1
    do while(j<=jumps .and. kprime/=0)
       
       ! going through neighbours
       ! [[file:~/mcp/selalib/src/pic_particle_types/lt_pic_4d_particle.F90::neighbour_pointers]]
       
       select case (dim)
#define ALONG_X 1
       case(ALONG_X)
          if(up) then
             neighbour = p_list(kprime)%ngb_xright_index
          else
             neighbour = p_list(kprime)%ngb_xleft_index
          endif
#define ALONG_Y 2
       case(ALONG_Y)
          if(up) then
             neighbour = p_list(kprime)%ngb_yright_index
          else
             neighbour = p_list(kprime)%ngb_yleft_index
          endif
#define ALONG_VX 3
       case(ALONG_VX)
          if(up) then
             neighbour = p_list(kprime)%ngb_vxright_index
          else
             neighbour = p_list(kprime)%ngb_vxleft_index
          endif
#define ALONG_VY 4
       case(ALONG_VY)
          if(up) then
             neighbour = p_list(kprime)%ngb_vyright_index
          else
             neighbour = p_list(kprime)%ngb_vyleft_index
          endif
       case default
          SLL_ASSERT(.false.)
       end select

       ! The convention in
       ! [[file:~/mcp/selalib/src/pic_particle_initializers/lt_pic_4d_init.F90::sll_lt_pic_4d_compute_new_particles]] is
       ! that if there is no neighbour then a neighbour index is equal to the particle index
          
       if(neighbour/=kprime) then
          kprime = neighbour
       else
          kprime = 0
       endif
       j = j + 1
    end do
  end subroutine onestep

  ! <<sll_lt_pic_4d_write_bsl_f_on_remap_grid>> <<ALH>> write the density on the (phase-space) remapping
  ! grid, using the method described in the "BSL-remapping" notes (version of december 2, 2014) cf
  ! [[file:~/mcp/maltpic/ltpic-bsl.tex::BSL_remapping]] and more precisely
  ! [[file:~/mcp/maltpic/ltpic-bsl.tex::BSL_remapping_step_1]].  Algorithm from
  ! [[file:~/mcp/maltpic/ltpic-bsl.tex::algo:pic-vr]] (but without the deposition step)
  !
  ! -- this function should be a faster alternative to [[sll_lt_pic_4d_write_f_on_remap_grid]] --
  ! 
  ! Note: the (x,y)-projection of the remapping grid may be larger than the "Poisson" 2d mesh associated with the
  ! particle group (in particular if the (x,y) domain is not periodic)

  subroutine sll_lt_pic_4d_write_bsl_f_on_remap_grid (p_group)

    ! [[file:../pic_particle_types/lt_pic_4d_group.F90::sll_lt_pic_4d_group]] p_group contains both the existing
    ! particles and the virtual remapping grid

    type(sll_lt_pic_4d_group),pointer,intent(inout) :: p_group
    
    ! cf [[file:~/mcp/maltpic/ltpic-bsl.tex::N*]]

    sll_int32 :: n_virtual = 2 ! <<n_virtual>>
    sll_int32 :: num_virtual_cells_x
    sll_int32 :: num_virtual_cells_y
    sll_int32 :: num_virtual_cells_vx
    sll_int32 :: num_virtual_cells_vy

    ! [[file:~/mcp/maltpic/ltpic-bsl.tex::h_parts_x]] and h_parts_y, h_parts_vx, h_parts_vy

    sll_real64 :: h_parts_x
    sll_real64 :: h_parts_y
    sll_real64 :: h_parts_vx
    sll_real64 :: h_parts_vy

    sll_real64 :: inv_h_parts_x    
    sll_real64 :: inv_h_parts_y    
    sll_real64 :: inv_h_parts_vx   
    sll_real64 :: inv_h_parts_vy   

    sll_real64 :: parts_x_min  
    sll_real64 :: parts_y_min  
    sll_real64 :: parts_vx_min 
    sll_real64 :: parts_vy_min 

    ! same as \delta{x,y,vx,vy} in [[file:~/mcp/maltpic/ltpic-bsl.tex::h_parts_x]]
    sll_real64 :: h_virtual_cell_x
    sll_real64 :: h_virtual_cell_y
    sll_real64 :: h_virtual_cell_vx
    sll_real64 :: h_virtual_cell_vy

    sll_real64 :: x
    sll_real64 :: y
    sll_real64 :: vx
    sll_real64 :: vy

    ! working values (offsets in virtual cell)

    sll_real32 :: dx
    sll_real32 :: dy
    sll_real32 :: dvx
    sll_real32 :: dvy

    ! working space

    sll_real64 :: tmp

    ! index of particle closest to the center of each virtual cell. Array dimensions defined by the contents of
    ! [[file:../pic_particle_types/lt_pic_4d_group.F90::sll_lt_pic_4d_group-remapping_grid]]. If [[n_virtual]] is
    ! greater than 1, the size of this array is smaller than the number of real remapping_grid cells.

    sll_int64,dimension(:,:,:,:),allocatable :: closest_particle
    sll_real64,dimension(:,:,:,:),allocatable :: closest_particle_distance

    sll_int32 :: i ! x dimension
    sll_int32 :: j ! y dimension
    sll_int64 :: k,kprime ! particle index
    sll_int64 :: neighbour ! particle index for local use
    sll_int32 :: l ! vx dimension
    sll_int32 :: m ! vy dimension

    ! indices in a virtual cell (go from 1 to [[n_virtual]])

    sll_int :: ivirt ! x dimension
    sll_int :: jvirt ! y dimension
    sll_int :: lvirt ! vx dimension
    sll_int :: mvirt ! vy dimension

    sll_int :: i_x,i_y,i_vx,i_vy

    ! particle pointer [[file:../pic_particle_types/lt_pic_4d_particle.F90::sll_lt_pic_4d_particle]]
    type(sll_lt_pic_4d_particle),pointer :: p

    ! <<g>> cartesian grid pointer to
    ! [[file:../pic_particle_types/lt_pic_4d_group.F90::sll_lt_pic_4d_group-remapping_grid]]

    type(sll_cartesian_mesh_4d),pointer :: g

    ! periodicity

    LOGICAL :: domain_is_x_periodic
    LOGICAL :: domain_is_y_periodic    
    sll_real64 :: mesh_period_x
    sll_real64 :: mesh_period_y
    sll_real64 :: inv_period_x
    sll_real64 :: inv_period_y

    ! results from [[get_ltp_deformation_matrix]]

    sll_real64 :: d11,d12,d13,d14 ! coefs of matrix D (backward Jacobian)
    sll_real64 :: d21,d22,d23,d24
    sll_real64 :: d31,d32,d33,d34
    sll_real64 :: d41,d42,d43,d44

    ! coordinates of particle k at time n and time 0
    sll_real64 :: x_k,y_k,vx_k,vy_k
    sll_real64 :: x_k_t0,y_k_t0,vx_k_t0,vy_k_t0

    sll_real64 :: part_radius_x 
    sll_real64 :: part_radius_y 
    sll_real64 :: part_radius_vx    
    sll_real64 :: part_radius_vy

    ! coordinates of a virtual particle at time 0 relative to the coordinates of one real particle

    sll_real64 :: x_t0,y_t0,vx_t0,vy_t0

    sll_int32 :: part_degree

    sll_int32 :: ierr

    ! temporary workspace
    sll_real64 :: x_aux
    sll_real64 :: y_aux
    sll_real64 :: vx_aux
    sll_real64 :: vy_aux

    sll_real64 :: length

    !aaa
    sll_real64 :: x_kprime_t0
    sll_real64 :: y_kprime_t0
    sll_real64 :: vx_kprime_t0
    sll_real64 :: vy_kprime_t0

    ! value 1 or 2 points to each side of an hypercube in direction x,y,vx or vy
    sll_int :: side_x,side_y,side_vx,side_vy
    sll_int64,dimension(2,2,2,2) :: hcube

    sll_int64 :: j_x,j_y,j_vx,j_vy

    sll_int64 :: number_parts_x
    sll_int64 :: number_parts_y
    sll_int64 :: number_parts_vx
    sll_int64 :: number_parts_vy
    
    ! --- end of declarations

    g => p_group%remapping_grid
    part_degree = p_group%spline_degree

    ! Preparatory work: find out the particle which is closest to each cell center by looping over all particles and
    ! noting which virtual cell contains it. The leftmost virtual cell in each dimension may not be complete.

    num_virtual_cells_x =  int(g%num_cells1/n_virtual)+1
    num_virtual_cells_y =  int(g%num_cells2/n_virtual)+1
    num_virtual_cells_vx = int(g%num_cells3/n_virtual)+1
    num_virtual_cells_vy = int(g%num_cells4/n_virtual)+1

    SLL_ALLOCATE(closest_particle(num_virtual_cells_x,num_virtual_cells_y,num_virtual_cells_vx,num_virtual_cells_vy),ierr)
    closest_particle(:,:,:,:) = 0

    SLL_ALLOCATE(closest_particle_distance(num_virtual_cells_x,num_virtual_cells_y,num_virtual_cells_vx,num_virtual_cells_vy),ierr)
    closest_particle_distance(:,:,:,:) = 0

    ! remapping grid cell size - same as in [[write_f_on_remap_grid-h_parts_x]]

    h_parts_x    = g%delta_eta1
    h_parts_y    = g%delta_eta2
    h_parts_vx   = g%delta_eta3
    h_parts_vy   = g%delta_eta4

    inv_h_parts_x  = 1./h_parts_x
    inv_h_parts_y  = 1./h_parts_y
    inv_h_parts_vx = 1./h_parts_vx
    inv_h_parts_vy = 1./h_parts_vy

    parts_x_min    = p_group%remapping_grid%eta1_min
    parts_y_min    = p_group%remapping_grid%eta2_min
    parts_vx_min   = p_group%remapping_grid%eta3_min
    parts_vy_min   = p_group%remapping_grid%eta4_min

    number_parts_x = p_group%number_parts_x
    number_parts_y = p_group%number_parts_y
    number_parts_vx = p_group%number_parts_vx
    number_parts_vy = p_group%number_parts_vy
    
    ! virtual cell size

    h_virtual_cell_x  = n_virtual * h_parts_x
    h_virtual_cell_y  = n_virtual * h_parts_y
    h_virtual_cell_vx = n_virtual * h_parts_vx
    h_virtual_cell_vy = n_virtual * h_parts_vy
    
    ! preparatory loop to fill the [[closest_particle]] array containing the particle closest to the center of each
    ! virtual cell

    do k=1,p_group%number_particles ! [[file:../pic_particle_types/lt_pic_4d_group.F90::number_particles]]
       p => p_group%p_list(k)

       ! find absolute (x,y,vx,vy) coordinates for this particle. Uses
       ! [[file:sll_representation_conversion.F90::cell_offset_to_global]]

       call cell_offset_to_global(p%dx,p%dy,p%ic,p_group%mesh,x,y)
       vx = p%vx
       vy = p%vy

       ! which _virtual_ cell is this particle in? uses
       ! [[file:sll_representation_conversion.F90::compute_cell_and_offset]] and [[g]]

       call compute_cell_and_offset(x,g%eta1_min,1./h_virtual_cell_x,i,dx)
       i=i+1
       SLL_ASSERT(i>0)
       SLL_ASSERT(dx>=0)
       SLL_ASSERT(dx<=1)
       call compute_cell_and_offset(y,g%eta2_min,1./h_virtual_cell_y,j,dy)
       j=j+1
       SLL_ASSERT(j>0)
       SLL_ASSERT(dy>=0)
       SLL_ASSERT(dy<=1)
       call compute_cell_and_offset(vx,g%eta3_min,1./h_virtual_cell_vx,l,dvx)
       l=l+1
       SLL_ASSERT(l>0)
       SLL_ASSERT(dvx>=0)
       SLL_ASSERT(dvx<=1)
       call compute_cell_and_offset(vy,g%eta4_min,1./h_virtual_cell_vy,m,dvy)
       m=m+1
       SLL_ASSERT(m>0)
       SLL_ASSERT(dvy>=0)
       SLL_ASSERT(dvy<=1)

       ! what is the distance from this particle to the virtual cell center? Speed things up a bit by skipping the
       ! square root calculation that will not change the final comparison of distances. Use adimensional values because
       ! adding values from different dimensions with different units makes little sense.

       tmp =  (dx - 0.5)**2.    &
            + (dy - 0.5)**2.    &
            + (dvx - 0.5)**2.    &
            + (dvy - 0.5)**2.

       ! if new particle is closer to center, keep the new one

       if(closest_particle(i,j,l,m) == 0 .or. tmp < closest_particle_distance(i,j,l,m)) then
          closest_particle(i,j,l,m) = k
          closest_particle_distance(i,j,l,m) = tmp
       end if
    end do

    ! Periodicity treatments copied from [[sll_lt_pic_4d_write_f_on_remap_grid-periodicity]]

    domain_is_x_periodic = .true.   ! temp
    domain_is_y_periodic = .true.   ! temp
    
    if(domain_is_x_periodic) then
      ! here the domain corresponds to the Poisson mesh
      mesh_period_x = p_group%mesh%eta1_max - p_group%mesh%eta1_min
      inv_period_x = 1./mesh_period_x
    else
      mesh_period_x = 0
      inv_period_x = 0
    end if
    
    if(domain_is_y_periodic) then
      ! here the domain corresponds to the Poisson mesh
      mesh_period_y = p_group%mesh%eta2_max - p_group%mesh%eta2_min
      inv_period_y = 1./mesh_period_y
    else
      mesh_period_y = 0
      inv_period_y = 0
    end if  

    ! initialize [[file:../pic_particle_types/lt_pic_4d_group.F90::target_values]]
    
    p_group%target_values(:,:,:,:) = 0

    ! MCP: [DEBUG] store the (computed) absolute initial position of the virtual particle
    p_group%debug_bsl_remap = -100


    ! <<loop_on_virtual_cells>> [[file:~/mcp/maltpic/ltpic-bsl.tex::algo:pic-vr:loop_over_all_cells]]
    ! Loop over all cells of indices i,j,l,m which contain at least one particle

    do i = 1,num_virtual_cells_x
       do j = 1,num_virtual_cells_y
          do l = 1,num_virtual_cells_vx
             do m = 1,num_virtual_cells_vy

                ! [[file:~/mcp/maltpic/ltpic-bsl.tex::algo:pic-vr:create_virtual_particles]] Create a temporary set of
                ! virtual particles inside the cell.  Note: in our case the virtual particles coincide with the existing
                ! remapping_grid [[file:../pic_particle_types/lt_pic_4d_group.F90::sll_lt_pic_4d_group-remapping_grid]]
                ! defined in p_group [[file:../pic_particle_types/lt_pic_4d_group.F90::sll_lt_pic_4d_group]]. So nothing
                ! more to do.

                ! [[file:~/mcp/maltpic/ltpic-bsl.tex::algo:pic-vr:find_closest_real_particle]] Find the real particle
                ! which is closest to the cell center.  Note: speed-wise, it may be necessary to find a way not to scan
                ! all the particles for every cell.  We avoid scanning all the particles for each cell by using the
                ! precomputed array [[closest_particle]]. Virtual cells which do not contain any particle are skipped.

                k = closest_particle(i,j,l,m)
                if(k /= 0) then

                   ! [[file:~/mcp/maltpic/ltpic-bsl.tex::hat-bz*]] Compute backward image of l-th virtual node by the
                   ! k-th backward flow. MCP -> oui, avec la matrice de deformation calculée avec la fonction
                   ! [[get_ltp_deformation_matrix]] pour la particule k. Calling [[get_ltp_deformation_matrix]]
                   ! with parameters inspired from [[sll_lt_pic_4d_write_f_on_remap_grid-get_ltp_deformation_matrix]]
    
                   call get_ltp_deformation_matrix (               &
                        k,                                         &    
                        p_group%mesh,                              &
                        p_group%p_list,                            &
                        domain_is_x_periodic,                      &
                        domain_is_y_periodic,                      &
                        mesh_period_x,                             &
                        mesh_period_y,                             &
                        h_parts_x,                                 &    
                        h_parts_y,                                 &    
                        h_parts_vx,                                &   
                        h_parts_vy,                                &   
                        1./h_parts_x,                              &    
                        1./h_parts_y,                              &    
                        1./h_parts_vx,                             &   
                        1./h_parts_vy,                             &   
                        0.5_f64*(part_degree+1),                   &
                        x_k,y_k,vx_k,vy_k,                         &
                        d11,d12,d13,d14,                           &
                        d21,d22,d23,d24,                           &
                        d31,d32,d33,d34,                           &
                        d41,d42,d43,d44,                           &
                        part_radius_x,                             &
                        part_radius_y,                             &
                        part_radius_vx,                            &
                        part_radius_vy                             &
                        )

                   ! Find position of particle k at time 0
                   ! [[get_initial_position_on_cartesian_grid_from_particle_index]]

                   call get_initial_position_on_cartesian_grid_from_particle_index(k, &
                        number_parts_x,number_parts_y,number_parts_vx,number_parts_vy, &
                        j_x,j_y,j_vx,j_vy)
                   x_k_t0 =  parts_x_min  + (j_x-1)  * h_parts_x
                   y_k_t0 =  parts_y_min  + (j_y-1)  * h_parts_y
                   vx_k_t0 = parts_vx_min + (j_vx-1) * h_parts_vx
                   vy_k_t0 = parts_vy_min + (j_vy-1) * h_parts_vy

                   ! <<loop_on_virtual_particles_in_one_virtual_cell>>
                   ! [[file:~/mcp/maltpic/ltpic-bsl.tex::algo:pic-vr:find_f0_for_each_virtual_particle]] Loop over all
                   ! virtual particles in the cell to compute the value of f0 at that point (Following
                   ! [[file:~/mcp/maltpic/ltpic-bsl.tex::BSL_remapping_algo]])

                   do ivirt = 1,n_virtual
                      do jvirt = 1,n_virtual
                         do lvirt = 1,n_virtual
                            do mvirt = 1,n_virtual

                               ! real index of the virtual particle in
                               ! [[file:../pic_particle_types/lt_pic_4d_group.F90::target_values]]

                               i_x =  (i-1)*n_virtual + ivirt
                               SLL_ASSERT(i_x>0)
                               i_y =  (j-1)*n_virtual + jvirt
                               SLL_ASSERT(i_y>0)
                               i_vx = (l-1)*n_virtual + lvirt
                               SLL_ASSERT(i_vx>0)
                               i_vy = (m-1)*n_virtual + mvirt
                               SLL_ASSERT(i_vy>0)

                               ! The index may go out of the domain for higher values of x,y,vx,vy in each dimension
                               ! (because the corners of the corresponding virtual cell do not correspond to existing
                               ! real particles). In that case, just ignore that value.

                               if(i_x<=p_group%number_parts_x            &
                                    .and. i_y<=p_group%number_parts_y    &
                                    .and. i_vx<=p_group%number_parts_vx  &
                                    .and. i_vy<=p_group%number_parts_vy) then

                                  ! Location of virtual particle (ivirt,jvirt,lvirt,mvirt) at time n
                                  
                                  x =  parts_x_min  + (i-1)*h_virtual_cell_x  + (ivirt-1)*h_parts_x
                                  y =  parts_y_min  + (j-1)*h_virtual_cell_y  + (jvirt-1)*h_parts_y
                                  vx = parts_vx_min + (l-1)*h_virtual_cell_vx + (lvirt-1)*h_parts_vx
                                  vy = parts_vy_min + (m-1)*h_virtual_cell_vy + (mvirt-1)*h_parts_vy

                                  ! particle k has to be inside the current virtual cell

                                  SLL_ASSERT(abs(x-x_k) < h_virtual_cell_x)
                                  SLL_ASSERT(abs(y-y_k) < h_virtual_cell_y)
                                  SLL_ASSERT(abs(vx-vx_k) < h_virtual_cell_vx)
                                  SLL_ASSERT(abs(vy-vy_k) < h_virtual_cell_vy)
                                  
                                  ! Location of virtual particle (ivirt,jvirt,lvirt,mvirt) at time 0 _relative_ to the
                                  ! position of particle k at time 0 (x_k_t0,y_k_t0,vx_k_t0,vy_k_t0) according to flow
                                  ! deformation

                                  x_t0  = d11 * (x - x_k) + d12 * (y - y_k) + d13 * (vx - vx_k) + d14 * (vy - vy_k)
                                  y_t0  = d21 * (x - x_k) + d22 * (y - y_k) + d23 * (vx - vx_k) + d24 * (vy - vy_k)
                                  vx_t0 = d31 * (x - x_k) + d32 * (y - y_k) + d33 * (vx - vx_k) + d34 * (vy - vy_k)
                                  vy_t0 = d41 * (x - x_k) + d42 * (y - y_k) + d43 * (vx - vx_k) + d44 * (vy - vy_k)

                                  ! MCP: [DEBUG] store the (computed) absolute initial position of the virtual particle
                                  p_group%debug_bsl_remap(i_x,i_y,i_vx,i_vy,1,1) = x_k_t0 + x_t0
                                  p_group%debug_bsl_remap(i_x,i_y,i_vx,i_vy,2,1) = y_k_t0 + y_t0
                                  p_group%debug_bsl_remap(i_x,i_y,i_vx,i_vy,3,1) = vx_k_t0 + vx_t0
                                  p_group%debug_bsl_remap(i_x,i_y,i_vx,i_vy,4,1) = vy_k_t0 + vy_t0

                                  ! Minimize the path to take along periodic boundaries

                                  if(domain_is_x_periodic) then
                                     if(x_t0 > mesh_period_x/2) then
                                        x_t0 = x_t0 - mesh_period_x
                                     else
                                        if(x_t0 < -mesh_period_x/2) then
                                           x_t0 = x_t0 + mesh_period_x
                                        end if
                                     end if
                                  endif
                                  
                                  if(domain_is_y_periodic) then
                                     if(y_t0 > mesh_period_y/2) then
                                        y_t0 = y_t0 - mesh_period_y
                                     else
                                        if(y_t0 < -mesh_period_y/2) then
                                           y_t0 = y_t0 + mesh_period_y
                                        end if
                                     end if
                                  endif
                                  
                                  ! MCP: [DEBUG] store the (computed) absolute initial position of the virtual particle
                                  p_group%debug_bsl_remap(i_x,i_y,i_vx,i_vy,1,2) = x_k_t0 + x_t0
                                  p_group%debug_bsl_remap(i_x,i_y,i_vx,i_vy,2,2) = y_k_t0 + y_t0
                                  p_group%debug_bsl_remap(i_x,i_y,i_vx,i_vy,3,2) = vx_k_t0 + vx_t0
                                  p_group%debug_bsl_remap(i_x,i_y,i_vx,i_vy,4,2) = vy_k_t0 + vy_t0

                                  ! [[file:~/mcp/maltpic/ltpic-bsl.tex::neighbors-grid-0]] find the neighbours of the
                                  ! virtual particle (ivirt,jvirt,lvirt,mvirt) at time 0 through the "logical
                                  ! neighbours" pointers of particle k. To reduce the amount of code, start with finding
                                  ! the closest neighbour which has lower coordinates in all directions. The particle
                                  ! located at (x_t0,y_t0,vx_t0,vy_t0) (coordinates relative to particle k to start
                                  ! with) gets progressively closer to kprime step by step (ie from neighbour to
                                  ! neighbour).

                                  kprime = k
                                  
                                  ! Calls [[onestep]]. "dim" can be x,y,vx,vy. cf
                                  ! [[file:../pic_particle_types/lt_pic_4d_particle.F90::sll_lt_pic_4d_particle]] for
                                  ! pointers to neighbours.

#define ONESTEPMACRO(dimpos,dimname) call onestep(dimpos,dimname/**/_t0,kprime,p_group%p_list,h_parts_/**/dimname)

                                  ONESTEPMACRO(ALONG_X,x)
                                  ONESTEPMACRO(ALONG_Y,y)
                                  ONESTEPMACRO(ALONG_VX,vx)
                                  ONESTEPMACRO(ALONG_VY,vy)

                                  ! If we end up with kprime == 0, it means that we have not found a cell that contains
                                  ! the particle so we just set that particle value to zero

                                  if (kprime /= 0) then

                                     ! kprime is the left-most vertex of the hypercube. find all the other vertices
                                     ! through the neighbour pointers in
                                     ! [[file:../pic_particle_types/lt_pic_4d_particle.F90::neighbour_pointers]]

                                     hcube(1,1,1,1) = kprime

                                     hcube(2,1,1,1) = p_group%p_list(kprime)%ngb_xright_index ! 1 step
                                     hcube(1,2,1,1) = p_group%p_list(kprime)%ngb_yright_index
                                     hcube(1,1,2,1) = p_group%p_list(kprime)%ngb_vxright_index
                                     hcube(1,1,1,2) = p_group%p_list(kprime)%ngb_vyright_index

                                     ! if any of the first four vertices is undefined (the convention in
                                     ! [[file:~/mcp/selalib/src/pic_particle_initializers/lt_pic_4d_init.F90::sll_lt_pic_4d_compute_new_particles]]
                                     ! is that the neighbour index is then equal to the particle index), it means that
                                     ! we reached the mesh border. just set the value of f for that particle as zero as
                                     ! before.

                                     if (hcube(2,1,1,1) /= kprime        &
                                          .and. hcube(1,2,1,1) /= kprime &
                                          .and. hcube(1,1,2,1) /= kprime &
                                          .and. hcube(1,1,1,2) /= kprime) then

                                        ! remaining vertices of the hypercube. they should all exist now that the first
                                        ! 4 vertices are checked.

                                        ! 1 step in x + 1 other step
                                        hcube(2,2,1,1) = p_group%p_list(hcube(2,1,1,1))%ngb_yright_index
                                        hcube(2,1,2,1) = p_group%p_list(hcube(2,1,1,1))%ngb_vxright_index
                                        hcube(2,1,1,2) = p_group%p_list(hcube(2,1,1,1))%ngb_vyright_index

                                        ! 1 step in y + 1 other step
                                        hcube(1,2,2,1) = p_group%p_list(hcube(1,2,1,1))%ngb_vxright_index
                                        hcube(1,2,1,2) = p_group%p_list(hcube(1,2,1,1))%ngb_vyright_index

                                        ! 1 step in vx + 1 other step
                                        hcube(1,1,2,2) = p_group%p_list(hcube(1,1,2,1))%ngb_vyright_index

                                        ! all combinations of 3 steps
                                        hcube(1,2,2,2) = p_group%p_list(hcube(1,2,2,1))%ngb_vyright_index
                                        hcube(2,1,2,2) = p_group%p_list(hcube(2,1,2,1))%ngb_vyright_index
                                        hcube(2,2,1,2) = p_group%p_list(hcube(2,2,1,1))%ngb_vyright_index
                                        hcube(2,2,2,1) = p_group%p_list(hcube(2,2,1,1))%ngb_vxright_index

                                        ! 4 steps
                                        hcube(2,2,2,2) = p_group%p_list(hcube(2,2,2,1))%ngb_vyright_index

                                        ! [[file:~/mcp/maltpic/ltpic-bsl.tex::affine-fn*]] use the values of f0 at these
                                        ! neighbours to interpolate the value of f0 at
                                        ! [[file:~/mcp/maltpic/ltpic-bsl.tex::hat-bz*]]. MCP -> oui. Ici si tu utilises
                                        ! des particules affines (part_deg = 1) la valeur de f0 se déduit de celle du
                                        ! poids de la particule.  En fait tu peux utiliser une formule semblable à celle
                                        ! qui est utilisée dans la fonction sll_lt_pic_4d_write_f_on_remap_grid, mais
                                        ! sans faire intervenir la matrice de déformation à l'intérieur des splines.

                                        ! place the resulting value of f on the virtual particle in
                                        ! p_group%target_values

                                          ! MCP: [BEGIN-DEBUG] store the (computed) absolute initial position of the virtual particle
                                        ! [[get_initial_position_on_cartesian_grid_from_particle_index]]
                                           call get_initial_position_on_cartesian_grid_from_particle_index(kprime, &
                                                number_parts_x,number_parts_y,number_parts_vx,number_parts_vy, &
                                                j_x,j_y,j_vx,j_vy)
                                           x_kprime_t0 =  parts_x_min  + (j_x-1)  * h_parts_x
                                           y_kprime_t0 =  parts_y_min  + (j_y-1)  * h_parts_y
                                           vx_kprime_t0 = parts_vx_min + (j_vx-1) * h_parts_vx
                                           vy_kprime_t0 = parts_vy_min + (j_vy-1) * h_parts_vy

                                           p_group%debug_bsl_remap(i_x,i_y,i_vx,i_vy,1,3) = x_kprime_t0 + x_t0
                                           p_group%debug_bsl_remap(i_x,i_y,i_vx,i_vy,2,3) = y_kprime_t0 + y_t0
                                           p_group%debug_bsl_remap(i_x,i_y,i_vx,i_vy,3,3) = vx_kprime_t0 + vx_t0
                                           p_group%debug_bsl_remap(i_x,i_y,i_vx,i_vy,4,3) = vy_kprime_t0 + vy_t0

                                          ! MCP [END-DEBUG]


                                        do side_x = 1,2
                                           if(side_x == 1)then
                                              x_aux = x_t0
                                           else
                                              x_aux = h_parts_x - x_t0
                                           end if
                                           do side_y = 1,2
                                              if(side_y == 1)then
                                                 y_aux = y_t0
                                              else
                                                 y_aux = h_parts_y - y_t0
                                              end if
                                              do side_vx = 1,2
                                                 if(side_vx == 1)then
                                                    vx_aux = vx_t0
                                                 else
                                                    vx_aux = h_parts_vx - vx_t0
                                                 end if
                                                 do side_vy = 1,2
                                                    if(side_vy == 1)then
                                                       vy_aux = vy_t0
                                                    else
                                                       vy_aux = h_parts_vy - vy_t0
                                                    end if

                                                    ! uses [[sll_pic_shape]]
                                                    p_group%target_values(i_x,i_y,i_vx,i_vy) =                    &
                                                         p_group%target_values(i_x,i_y,i_vx,i_vy)                 &
                                                         + p_group%p_list(hcube(side_x,side_y,side_vx,side_vy))%q &
                                                         * sll_pic_shape(part_degree,x_aux,y_aux,vx_aux,vy_aux,   &
                                                         inv_h_parts_x,inv_h_parts_y,inv_h_parts_vx,inv_h_parts_vy)
                                                 end do
                                              end do
                                           end do
                                        end do
                                     end if
                                  end if
                               end if
                            end do
                         end do
                      end do
                   end do
                end if
             end do
          end do
       end do
    end do
  end subroutine sll_lt_pic_4d_write_bsl_f_on_remap_grid

  ! <<sll_lt_pic_4d_write_f_on_remap_grid>>
  ! write the density on the (phase-space) remapping grid
  ! Note: the (x,y)-projection of the remapping grid may be larger than the "Poisson" 2d mesh associated with the particle group
  ! (in particular if the (x,y) domain is not periodic)
!  subroutine write_lt_pic_density_on_remap_grid ( &            ! old name
  subroutine sll_lt_pic_4d_write_f_on_remap_grid ( &
                                p_group          &
                                )
    type(sll_lt_pic_4d_group), pointer,  intent(in)  :: p_group
    
    type(sll_cartesian_mesh_2d),      pointer  :: particles_m2d
    type(sll_lt_pic_4d_particle), dimension(:),  pointer :: p

    sll_int32 :: part_degree
    
    sll_int32 :: mx_min
    sll_int32 :: mx_max
    sll_int32 :: my_min
    sll_int32 :: my_max
    sll_int32 :: mx, my

    sll_int32 :: number_parts_x
    sll_int32 :: number_parts_y
    sll_int32 :: number_parts_vx
    sll_int32 :: number_parts_vy
    sll_int64 :: i_x_min,  i_x_max,  i_x
    sll_int64 :: i_y_min,  i_y_max,  i_y
    sll_int64 :: i_vx_min, i_vx_max, i_vx
    sll_int64 :: i_vy_min, i_vy_max, i_vy
    sll_int64 :: k, k_ngb
    
    sll_real64 :: h_parts_x    
    sll_real64 :: h_parts_y    
    sll_real64 :: h_parts_vx   
    sll_real64 :: h_parts_vy   
    sll_real64 :: inv_h_parts_x    
    sll_real64 :: inv_h_parts_y    
    sll_real64 :: inv_h_parts_vx   
    sll_real64 :: inv_h_parts_vy   
    sll_real64 :: parts_x_min  
    sll_real64 :: parts_y_min  
    sll_real64 :: parts_vx_min 
    sll_real64 :: parts_vy_min 
    sll_real64 :: mesh_period_x
    sll_real64 :: mesh_period_y
    sll_real64 :: inv_period_x
    sll_real64 :: inv_period_y

    sll_real64 :: x_k_left,  x_k_right
    sll_real64 :: y_k_left,  y_k_right
    sll_real64 :: vx_k_left, vx_k_right
    sll_real64 :: vy_k_left, vy_k_right

    sll_real64 :: d11,d12,d13,d14 ! coefs of matrix D (backward Jacobian)
    sll_real64 :: d21,d22,d23,d24
    sll_real64 :: d31,d32,d33,d34
    sll_real64 :: d41,d42,d43,d44
    sll_real64 :: j11,j12,j13,j14 ! coefs of matrix J = D^-1 (forward Jacobian)
    sll_real64 :: j21,j22,j23,j24
    sll_real64 :: j31,j32,j33,j34
    sll_real64 :: j41,j42,j43,j44
    sll_real64 :: factor, det_J, inv_det_J
    
    sll_real64 :: x_k, x_k_per          !! CHECK ----  x_k_per not used ??   
    sll_real64 :: y_k, y_k_per
    sll_real64 :: vx_k, vy_k, w_k
    
    sll_real64 :: ref_radius
    sll_real64 :: part_radius_x 
    sll_real64 :: part_radius_y 
    sll_real64 :: part_radius_vx    
    sll_real64 :: part_radius_vy

    sll_real64 :: x,  x_aux
    sll_real64 :: y,  y_aux
    sll_real64 :: vx, vx_aux
    sll_real64 :: vy, vy_aux
    
    LOGICAL :: domain_is_x_periodic
    LOGICAL :: domain_is_y_periodic    

    ! Poisson mesh associated to the particles
    particles_m2d => p_group%mesh
    
    p            => p_group%p_list
    part_degree  = p_group%spline_degree
    ref_radius = 0.5*(part_degree+1)        

    number_parts_x  = p_group%number_parts_x
    number_parts_y  = p_group%number_parts_y
    number_parts_vx = p_group%number_parts_vx
    number_parts_vy = p_group%number_parts_vy

    ! array of density values on the remapping grid
    p_group%target_values(:,:,:,:) = 0
        
    h_parts_x    = p_group%remapping_grid%delta_eta1 ! <<write_f_on_remap_grid-h_parts_x>>
    h_parts_y    = p_group%remapping_grid%delta_eta2
    h_parts_vx   = p_group%remapping_grid%delta_eta3
    h_parts_vy   = p_group%remapping_grid%delta_eta4

    inv_h_parts_x  = 1./h_parts_x
    inv_h_parts_y  = 1./h_parts_y
    inv_h_parts_vx = 1./h_parts_vx
    inv_h_parts_vy = 1./h_parts_vy

    parts_x_min    = p_group%remapping_grid%eta1_min
    parts_y_min    = p_group%remapping_grid%eta2_min
    parts_vx_min   = p_group%remapping_grid%eta3_min
    parts_vy_min   = p_group%remapping_grid%eta4_min
        

    ! <<sll_lt_pic_4d_write_f_on_remap_grid-periodicity>>

    domain_is_x_periodic = .true.   ! temp
    domain_is_y_periodic = .true.   ! temp
    
    if(domain_is_x_periodic) then
      ! here the domain corresponds to the Poisson mesh
      mesh_period_x = particles_m2d%eta1_max - particles_m2d%eta1_min
      inv_period_x = 1./mesh_period_x
    else
      mesh_period_x = 0
      inv_period_x = 0
    end if
    
    if(domain_is_y_periodic) then
      ! here the domain corresponds to the Poisson mesh
      mesh_period_y = particles_m2d%eta2_max - particles_m2d%eta2_min
      inv_period_y = 1./mesh_period_y
    else
      mesh_period_y = 0
      inv_period_y = 0
    end if  
    
    ! loop over the particles first, then over the relevant grid points 
    do k = 1, p_group%number_particles

!        call cell_offset_to_global( p(k)%dx, &
!                                    p(k)%dy, &
!                                    p(k)%ic, &
!                                    particles_m2d, x_k, y_k )
!        vx_k   = p(k)%vx
!        vy_k   = p(k)%vy
        w_k    = p(k)%q           

        ! <<sll_lt_pic_4d_write_f_on_remap_grid-get_ltp_deformation_matrix>>

        call get_ltp_deformation_matrix(                    &
                                    k,                      &
                                    particles_m2d,          &
                                    p,                      &
                                    domain_is_x_periodic,   &
                                    domain_is_y_periodic,   &
                                    mesh_period_x,          &
                                    mesh_period_y,          &
                                    h_parts_x,              &     
                                    h_parts_y,              &    
                                    h_parts_vx,             &   
                                    h_parts_vy,             &   
                                    inv_h_parts_x,          &    
                                    inv_h_parts_y,          &    
                                    inv_h_parts_vx,         &      
                                    inv_h_parts_vy,         &   
                                    ref_radius,             &
                                    x_k, y_k,               &
                                    vx_k, vy_k,             &
                                    d11,d12,d13,d14,        &
                                    d21,d22,d23,d24,        &
                                    d31,d32,d33,d34,        &
                                    d41,d42,d43,d44,        &
                                    part_radius_x,          &
                                    part_radius_y,          &
                                    part_radius_vx,         &
                                    part_radius_vy          &
                                   )

        if(domain_is_x_periodic) then
           mx_min = int(floor(   inv_period_x * (particles_m2d%eta1_min - x_k - part_radius_x) ) ) + 1
           mx_max = int(ceiling( inv_period_x * (particles_m2d%eta1_max - x_k + part_radius_x) ) ) - 1
        else
           mx_min = 0
           mx_max = 0
        end if
        
        if(domain_is_y_periodic) then
           my_min = int(floor(   inv_period_y * (particles_m2d%eta2_min - y_k - part_radius_y) ) ) + 1
           my_max = int(ceiling( inv_period_y * (particles_m2d%eta2_max - y_k + part_radius_y) ) ) - 1
        else
           my_min = 0
           my_max = 0
        end if
               
        ! Now determine the grid overlapping the support of the particle (see doc sheet)
        do mx = mx_min, mx_max
            do my = my_min, my_max
            
                i_x_min  = max( 1, &
                               int(floor(   inv_h_parts_x  * (x_k + real(mx*mesh_period_x,f64)      &
                                                                - parts_x_min  - part_radius_x ) ) ) +2 )
                i_x_max  = min( number_parts_x, &
                               int(ceiling( inv_h_parts_x  * (x_k + real(mx*mesh_period_x,f64)      &
                                                                - parts_x_min  + part_radius_x ) ), i32 )  )
                i_y_min  = max( 1, &
                               int(floor(   inv_h_parts_y  * (y_k + real(my*mesh_period_y,f64)      &
                                                                - parts_y_min  - part_radius_y ) ) ) +2 )
                i_y_max  = min( number_parts_y, &
                               int(ceiling( inv_h_parts_y  * (y_k + real(my*mesh_period_y,f64)      &
                                                                - parts_y_min  + part_radius_y ) ), i32 ) )
                i_vx_min = max( 1, &
                               int(floor(   inv_h_parts_vx * (vx_k - parts_vx_min - part_radius_vx ) ) ) +2 )
                i_vx_max = min( number_parts_vx, &
                               int(ceiling( inv_h_parts_vx * (vx_k - parts_vx_min + part_radius_vx ) ), i32 ) )
                i_vy_min = max( 1, &
                               int(floor(   inv_h_parts_vy * (vy_k - parts_vy_min - part_radius_vy ) ) )  +2 )
                i_vy_max = min( number_parts_vy, &
                               int(ceiling( inv_h_parts_vy * (vy_k - parts_vy_min + part_radius_vy ) ), i32 )  )
                ! loop over the grid points 
                do i_x = i_x_min, i_x_max
                    x = parts_x_min + (i_x-1) * h_parts_x - mx * mesh_period_x
                    do i_y = i_y_min, i_y_max
                        y = parts_y_min + (i_y-1) * h_parts_y - my * mesh_period_y
                        do i_vx = i_vx_min, i_vx_max
                            vx = parts_vx_min + (i_vx-1) * h_parts_vx
                            do i_vy = i_vy_min, i_vy_max
                               vy = parts_vy_min + (i_vy-1) * h_parts_vy

                                x_aux =   d11 * (x - x_k)  &
                                         + d12 * (y - y_k)  &
                                         + d13 * (vx - vx_k)&
                                         + d14 * (vy - vy_k)
                                 
                                y_aux =   d21 * (x - x_k)  &
                                         + d22 * (y - y_k)  &
                                         + d23 * (vx - vx_k)&
                                         + d24 * (vy - vy_k)
                         
                                vx_aux =  d31 * (x - x_k)  &
                                         + d32 * (y - y_k)  &
                                         + d33 * (vx - vx_k)&
                                         + d34 * (vy - vy_k)
                                 
                                vy_aux =  d41 * (x - x_k)  &
                                         + d42 * (y - y_k)  &
                                         + d43 * (vx - vx_k)&
                                         + d44 * (vy - vy_k)

                                p_group%target_values(i_x, i_y, i_vx, i_vy)                 &
                                      = p_group%target_values(i_x, i_y, i_vx, i_vy)         &
                                                    + w_k * sll_pic_shape( part_degree,     &
                                                                            x_aux,          &
                                                                            y_aux,          &
                                                                            vx_aux,         &
                                                                            vy_aux,         &
                                                                            inv_h_parts_x,  &
                                                                            inv_h_parts_y,  &
                                                                            inv_h_parts_vx, &
                                                                            inv_h_parts_vy  &
                                                                          )
                            end do
                        end do
                    end do
                end do
            end do
        end do
    end do ! end loop on particles
    
  end subroutine sll_lt_pic_4d_write_f_on_remap_grid


!  subroutine plot_lt_pic_density_on_2d_grid( &              ! old name
  subroutine sll_lt_pic_4d_plot_f_on_2d_grid( &
      file_name,        &
      plot_dim1,        &
      plot_dim2,        &
      plot_cst_dim3,    &
      plot_cst_dim4,    &
      x3_plot_cst,      &
      x4_plot_cst,      &      
      plot_m2d,         &
      p_group           &
      )
    character(len=*),                         intent(in)    :: file_name
    type(sll_cartesian_mesh_2d),      pointer,  intent(in)    :: plot_m2d
    type(sll_lt_pic_4d_group), pointer,  intent(in)    :: p_group
    sll_int32, intent(in)  ::  plot_dim1
    sll_int32, intent(in)  :: plot_dim2
    sll_int32, intent(in)  :: plot_cst_dim3
    sll_int32, intent(in)  :: plot_cst_dim4
    sll_real64, intent(in) :: x3_plot_cst
    sll_real64, intent(in) :: x4_plot_cst    

    type(sll_cartesian_mesh_2d),      pointer                 :: particles_m2d
    sll_real64, dimension(:,:),     pointer                 :: plotted_values    
    type(sll_lt_pic_4d_particle),       dimension(:),  pointer  :: p

    sll_int32 :: ierr
    sll_int32 :: part_degree

    sll_int32 :: mx_min
    sll_int32 :: mx_max
    sll_int32 :: my_min
    sll_int32 :: my_max
    sll_int32 :: mx, my
    sll_int32 :: m1, m2
        
    sll_int64 :: num_points_1
    sll_int64 :: num_points_2
!    sll_int64 :: i_x, i_y, i_vx, i_vy
    sll_int64 :: i1_min, i1_max, i1
    sll_int64 :: i2_min, i2_max, i2
    sll_int64 :: k, k_ngb

    sll_real64 :: h_parts_x    
    sll_real64 :: h_parts_y    
    sll_real64 :: h_parts_vx   
    sll_real64 :: h_parts_vy   
    sll_real64 :: inv_h_parts_x    
    sll_real64 :: inv_h_parts_y    
    sll_real64 :: inv_h_parts_vx   
    sll_real64 :: inv_h_parts_vy   
    sll_real64 :: parts_x_min  
    sll_real64 :: parts_y_min  
    sll_real64 :: parts_vx_min 
    sll_real64 :: parts_vy_min 

    sll_real64 :: h_plot_1    
    sll_real64 :: h_plot_2
    sll_real64 :: inv_h_plot_1    
    sll_real64 :: inv_h_plot_2
    sll_real64 :: plot_x1_min  
    sll_real64 :: plot_x2_min  

    sll_real64 :: mesh_period_x
    sll_real64 :: mesh_period_y
    sll_real64 :: inv_period_x
    sll_real64 :: inv_period_y

    sll_real64 :: x_k_left,  x_k_right
    sll_real64 :: y_k_left,  y_k_right
    sll_real64 :: vx_k_left, vx_k_right
    sll_real64 :: vy_k_left, vy_k_right

    sll_real64 :: mesh_period_1
    sll_real64 :: mesh_period_2

    sll_real64 :: d11,d12,d13,d14 ! coefs of matrix D (backward Jacobian)
    sll_real64 :: d21,d22,d23,d24
    sll_real64 :: d31,d32,d33,d34
    sll_real64 :: d41,d42,d43,d44
    sll_real64 :: j11,j12,j13,j14 ! coefs of matrix J = D^-1 (forward Jacobian)
    sll_real64 :: j21,j22,j23,j24
    sll_real64 :: j31,j32,j33,j34
    sll_real64 :: j41,j42,j43,j44
    sll_real64 :: factor, det_J, inv_det_J
        
    sll_real64 :: ref_radius
    sll_real64 :: part_radius_x 
    sll_real64 :: part_radius_y 
    sll_real64 :: part_radius_vx    
    sll_real64 :: part_radius_vy

    sll_real64 :: part_radius_1 
    sll_real64 :: part_radius_2 

    sll_real64 :: x_k, x_k_per
    sll_real64 :: y_k, y_k_per
    sll_real64 :: vx_k, vy_k, w_k
    sll_real64 :: x1_k, x2_k
    sll_real64 :: x_aux
    sll_real64 :: y_aux
    sll_real64 :: vx_aux
    sll_real64 :: vy_aux
    sll_real64 :: x_plot
    sll_real64 :: y_plot
    sll_real64 :: vx_plot
    sll_real64 :: vy_plot
    sll_real64 :: x1_plot
    sll_real64 :: x2_plot
        
    LOGICAL :: domain_is_x_periodic
    LOGICAL :: domain_is_y_periodic    

    sll_real64 :: sum_part_radius_x
    sll_real64 :: sum_part_radius_y
    sll_real64 :: sum_part_radius_vx
    sll_real64 :: sum_part_radius_vy
    sll_int32 :: mx_min_sum
    sll_int32 :: mx_max_sum
    sll_int32 :: my_min_sum
    sll_int32 :: my_max_sum
    sll_int64 :: n_debug

    sll_real64 :: time
    type(sll_time_mark)  :: t0


    print*, "[sll_lt_pic_4d_plot_f_on_2d_grid] plotting density in file :", file_name
    print*, "[sll_lt_pic_4d_plot_f_on_2d_grid] plot_m2d = ", plot_m2d
    print*, "[sll_lt_pic_4d_plot_f_on_2d_grid] plot_m2d%num_cells1 = ", plot_m2d%num_cells1
    print*, "[sll_lt_pic_4d_plot_f_on_2d_grid] plot_m2d%num_cells2 = ", plot_m2d%num_cells2
    call sll_set_time_mark(t0)
    sum_part_radius_x = 0
    sum_part_radius_y = 0
    sum_part_radius_vx = 0
    sum_part_radius_vy = 0
    mx_min_sum = 0
    mx_max_sum = 0
    my_min_sum = 0
    my_max_sum = 0
    n_debug = int(0.8 * p_group%number_particles)

    p             => p_group%p_list
    part_degree   = p_group%spline_degree
    ref_radius = 0.5*(part_degree+1)
    
    num_points_1  = plot_m2d%num_cells1+1
    num_points_2  = plot_m2d%num_cells2+1

    ! Poisson mesh associated to the particles
    particles_m2d => p_group%mesh

    ! array of density values on the remapping grid
    SLL_ALLOCATE( plotted_values(num_points_1, num_points_2), ierr )
    plotted_values(:,:) = 0
        
    ! parameters of the particle (remapping) grid
    h_parts_x    = p_group%remapping_grid%delta_eta1
    h_parts_y    = p_group%remapping_grid%delta_eta2
    h_parts_vx   = p_group%remapping_grid%delta_eta3
    h_parts_vy   = p_group%remapping_grid%delta_eta4

    inv_h_parts_x  = 1./h_parts_x
    inv_h_parts_y  = 1./h_parts_y
    inv_h_parts_vx = 1./h_parts_vx
    inv_h_parts_vy = 1./h_parts_vy

    parts_x_min    = p_group%remapping_grid%eta1_min
    parts_y_min    = p_group%remapping_grid%eta2_min
    parts_vx_min   = p_group%remapping_grid%eta3_min
    parts_vy_min   = p_group%remapping_grid%eta4_min

    ! parameters of the plotting grid
    h_plot_1    = plot_m2d%delta_eta1
    h_plot_2    = plot_m2d%delta_eta2

    inv_h_plot_1  = 1./h_plot_1
    inv_h_plot_2  = 1./h_plot_2

    plot_x1_min    = plot_m2d%eta1_min
    plot_x2_min    = plot_m2d%eta2_min
        
    domain_is_x_periodic = .true.   ! temp
    domain_is_y_periodic = .true.   ! temp
    
    if(domain_is_x_periodic) then
      ! here the domain corresponds to the Poisson mesh
      mesh_period_x = particles_m2d%eta1_max - particles_m2d%eta1_min
      inv_period_x = 1./mesh_period_x
    else
      mesh_period_x = 0
      inv_period_x = 0
    end if
    
    if(domain_is_y_periodic) then
      ! here the domain corresponds to the Poisson mesh
      mesh_period_y = particles_m2d%eta2_max - particles_m2d%eta2_min
      inv_period_y = 1./mesh_period_y
    else
      mesh_period_y = 0
      inv_period_y = 0
    end if  
    
    ! loop over the particles first, then over the relevant grid points 
    do k = 1, p_group%number_particles


!        call cell_offset_to_global( p(k)%dx, &
!                                    p(k)%dy, &
!                                    p(k)%ic, &
!                                    particles_m2d, x_k, y_k )
!        vx_k   = p(k)%vx
!        vy_k   = p(k)%vy
        
        w_k    = p(k)%q           

        call get_ltp_deformation_matrix(                    &
                                    k,                      &
                                    particles_m2d,          &
                                    p,                      &
                                    domain_is_x_periodic,   &
                                    domain_is_y_periodic,   &
                                    mesh_period_x,          &
                                    mesh_period_y,          &
                                    h_parts_x,              &     
                                    h_parts_y,              &    
                                    h_parts_vx,             &   
                                    h_parts_vy,             &   
                                    inv_h_parts_x,          &    
                                    inv_h_parts_y,          &    
                                    inv_h_parts_vx,         &      
                                    inv_h_parts_vy,         &   
                                    ref_radius,             &
                                    x_k, y_k,               &
                                    vx_k, vy_k,             &
                                    d11,d12,d13,d14,        &
                                    d21,d22,d23,d24,        &
                                    d31,d32,d33,d34,        &
                                    d41,d42,d43,d44,        &
                                    part_radius_x,          &
                                    part_radius_y,          &
                                    part_radius_vx,         &
                                    part_radius_vy          &
                                   )

        select case( plot_dim1 )
            case( 1 )
                x1_k = x_k
                part_radius_1 = part_radius_x
            case( 2 )
                x1_k = y_k
                part_radius_1 = part_radius_y
            case( 3 )
                x1_k = vx_k
                part_radius_1 = part_radius_vx
            case( 4 )
                x1_k = vy_k
                part_radius_1 = part_radius_vy
            case default
                print*, "WARNING (code=569824526), invalid value for plot_dim1: ", plot_dim1, &
                        " taking 1 instead"
                x1_k = x_k
                part_radius_1 = part_radius_x
        end select

        select case( plot_dim2 )
            case( 1 )
                x2_k = x_k
                part_radius_2 = part_radius_x
            case( 2 )
                x2_k = y_k
                part_radius_2 = part_radius_y
            case( 3 )
                x2_k = vx_k
                part_radius_2 = part_radius_vx
            case( 4 )
                x2_k = vy_k
                part_radius_2 = part_radius_vy
            case default
                print*, "WARNING (code=569824527), invalid value for plot_dim2: ", plot_dim2, &
                        " taking 2 instead"
                x2_k = y_k
                part_radius_2 = part_radius_y
        end select

        select case( plot_cst_dim3 )
            case( 1 )
!                x3_k    = x_k
!                part_radius_3 = part_radius_x
                x_plot  = x3_plot_cst
            case( 2 )
!                x3_k    = y_k
!                part_radius_3 = part_radius_y
                y_plot  = x3_plot_cst
            case( 3 )
!                x3_k    = vx_k
!                part_radius_3 = part_radius_vx
                vx_plot = x3_plot_cst
            case( 4 )
!                x3_k    = vy_k
!                part_radius_3 = part_radius_vy
                vy_plot = x3_plot_cst
            case default
                print*, "WARNING (code=569824528), invalid value for plot_cst_dim3: ", plot_cst_dim3, &
                        " taking 3 instead"
!                x3_k    = vx_k
!                part_radius_3 = part_radius_vx
                vx_plot = x3_plot_cst
        end select

        select case( plot_cst_dim4 )
            case( 1 )
!                x4_k    = x_k
!                part_radius_4 = part_radius_x
                x_plot  = x4_plot_cst
            case( 2 )
!                x4_k    = y_k
!                part_radius_4 = part_radius_y
                y_plot  = x4_plot_cst
            case( 3 )
!                x4_k    = vx_k
!                part_radius_4 = part_radius_vx
                vx_plot = x4_plot_cst
            case( 4 )
!                x4_k    = vy_k
!                part_radius_4 = part_radius_vy
                vy_plot = x4_plot_cst
            case default
                print*, "WARNING (code=569824529), invalid value for plot_cst_dim4: ", plot_cst_dim4, &
                        " taking 4 instead"
!                x4_k    = vy_k
!                part_radius_4 = part_radius_vy
                vy_plot = x4_plot_cst
        end select


        if(domain_is_x_periodic) then
           mx_min = int(floor(   inv_period_x * (particles_m2d%eta1_min - x_k - part_radius_x) ) ) + 1
           mx_max = int(ceiling( inv_period_x * (particles_m2d%eta1_max - x_k + part_radius_x) ) ) - 1
        else
           mx_min = 0
           mx_max = 0
        end if
        
        if(domain_is_y_periodic) then
           my_min = int(floor(   inv_period_y * (particles_m2d%eta2_min - y_k - part_radius_y) ) ) + 1
           my_max = int(ceiling( inv_period_y * (particles_m2d%eta2_max - y_k + part_radius_y) ) ) - 1
        else
           my_min = 0
           my_max = 0
        end if
        ! debug
        sum_part_radius_x = sum_part_radius_x + part_radius_x
        sum_part_radius_y = sum_part_radius_y + part_radius_y
        sum_part_radius_vx = sum_part_radius_vx + part_radius_vx
        sum_part_radius_vy = sum_part_radius_vy + part_radius_vy
        mx_min_sum = mx_min_sum + mx_min
        mx_max_sum = mx_max_sum + mx_max
        my_min_sum = my_min_sum + my_min
        my_max_sum = my_max_sum + my_max
        if (mod(k,n_debug)==0) then 
          print*, "[sll_lt_pic_4d_plot_f_on_2d_grid] plotting part k = ", k
          print*, " -- average mx_min, mx_max = ", mx_min_sum * 1./n_debug, mx_max_sum * 1./n_debug
          print*, " -- average my_min, my_max = ", my_min_sum * 1./n_debug, my_max_sum * 1./n_debug
          print*, " -- -- -- -- -- -- "

          print*, " -- (average part_radius_x) / h_parts_x = ", (sum_part_radius_x * 1./n_debug) / h_parts_x
          print*, " -- (average part_radius_y) / h_parts_y = ", (sum_part_radius_y * 1./n_debug) / h_parts_y
          print*, " -- (average part_radius_vx) / h_parts_vx = ", (sum_part_radius_vx * 1./n_debug) / h_parts_vx
          print*, " -- (average part_radius_vy) / h_parts_vy = ", (sum_part_radius_vy * 1./n_debug) / h_parts_vy
          print*, " -- -- -- -- -- -- "

          time = sll_time_elapsed_since(t0)
          print*, 'TIME=', time, ' -- ', n_debug*1./time, 'average particles treated/sec'
          call sll_set_time_mark(t0)
          sum_part_radius_x = 0
          sum_part_radius_y = 0
          sum_part_radius_vx = 0
          sum_part_radius_vy = 0
          mx_min_sum = 0
          mx_max_sum = 0
          my_min_sum = 0
          my_max_sum = 0
          
        end if
        
        ! Now determine the grid points overlapping the support of the particle
        do mx = mx_min, mx_max            
            do my = my_min, my_max

                x_k_per = x_k + mx * mesh_period_x
                y_k_per = y_k + my * mesh_period_y

                m1 = 0
                m2 = 0
                mesh_period_1 = 0
                mesh_period_2 = 0
                if( plot_dim1 == 1 )then
                    m1 = mx
                    mesh_period_1 = mesh_period_x
                else if( plot_dim1 == 2 )then
                    m1 = my
                    mesh_period_1 = mesh_period_y
                end if
                if( plot_dim2 == 1 )then
                    m2 = mx
                    mesh_period_2 = mesh_period_x
                else if( plot_dim2 == 2 )then
                    m2 = my
                    mesh_period_2 = mesh_period_y
                end if
!                if( plot_dim3 == 1 )then
!                    m3 = mx
!                    mesh_period_3 = mesh_period_x
!                else if( plot_dim3 == 2 )then
!                    m3 = my
!                    mesh_period_3 = mesh_period_y
!                end if
!                if( plot_dim4 == 1 )then
!                    m4 = mx
!                    mesh_period_4 = mesh_period_x
!                else if( plot_dim4 == 2 )then
!                    m4 = my
!                    mesh_period_4 = mesh_period_y
!                end if
            
                i1_min  = max( 1, &
                               int(floor(   inv_h_plot_1  * ( x1_k + real( m1*mesh_period_1, f64)        &
                                                                - plot_x1_min  - part_radius_1      ) ) ) + 2 )
                i1_max  = min( num_points_1, &
                               int(ceiling( inv_h_plot_1  * ( x1_k + real( m1*mesh_period_1, f64)        &
                                                                - plot_x1_min  + part_radius_1      ) ), i64 ) )
                i2_min  = max( 1, &
                               int(floor(   inv_h_plot_2  * ( x2_k + real( m2*mesh_period_2, f64)        &
                                                                - plot_x2_min  - part_radius_2      ) ) ) + 2 )
                i2_max  = min( num_points_2, &
                               int(ceiling( inv_h_plot_2  * ( x2_k + real( m2*mesh_period_2, f64)        &
                                                                - plot_x2_min  + part_radius_2      ) ), i64 ) )

                ! loop over the grid points 
                do i1 = i1_min, i1_max
                    x1_plot = plot_x1_min + (i1-1) * h_plot_1
                    do i2 = i2_min, i2_max
                        x2_plot = plot_x2_min + (i2-1) * h_plot_2

                        select case( plot_dim1 )
                            case( 1 )
                                x_plot  = x1_plot
                            case( 2 )
                                y_plot  = x1_plot
                            case( 3 )
                                vx_plot = x1_plot
                            case( 4 )
                                vy_plot = x1_plot
                            case default
                                print*, "WARNING (code=569824526), invalid value for plot_dim1: ", plot_dim1, &
                                        " taking 1 instead"
                                x_plot  = x1_plot
                        end select
                
                        select case( plot_dim2 )
                            case( 1 )
                                x_plot  = x2_plot
                            case( 2 )
                                y_plot  = x2_plot
                            case( 3 )
                                vx_plot = x2_plot
                            case( 4 )
                                vy_plot = x2_plot
                            case default
                                print*, "WARNING (code=569824527), invalid value for plot_dim2: ", plot_dim2, &
                                        " taking 2 instead"
                                y_plot  = x2_plot
                        end select

    
                        x_aux =   d11 * (x_plot  - x_k_per)  &
                                 + d12 * (y_plot  - y_k_per)  &
                                 + d13 * (vx_plot - vx_k)&
                                 + d14 * (vy_plot - vy_k)
                         
                        y_aux =   d21 * (x_plot  - x_k_per)  &
                                 + d22 * (y_plot  - y_k_per)  &
                                 + d23 * (vx_plot - vx_k)&
                                 + d24 * (vy_plot - vy_k)
                 
                        vx_aux =  d31 * (x_plot  - x_k_per)  &
                                 + d32 * (y_plot  - y_k_per)  &
                                 + d33 * (vx_plot - vx_k)&
                                 + d34 * (vy_plot - vy_k)
                         
                        vy_aux =  d41 * (x_plot  - x_k_per)  &
                                 + d42 * (y_plot  - y_k_per)  &
                                 + d43 * (vx_plot - vx_k)&
                                 + d44 * (vy_plot - vy_k)
                                         
                                plotted_values(i1, i2)               &
                                      = plotted_values(i1, i2)       &
                                        + w_k * sll_pic_shape( part_degree,     &
                                                                x_aux,          &
                                                                y_aux,          &
                                                                vx_aux,         &
                                                                vy_aux,         &
                                                                inv_h_parts_x,      &
                                                                inv_h_parts_y,      &
                                                                inv_h_parts_vx,     &
                                                                inv_h_parts_vy      &
                                                                )
!                        if (mod(k,100_i64)==0) then 
!                            print*, "[sll_lt_pic_4d_plot_f_on_2d_grid] plotted_values(i1, i2) = ", plotted_values(i1, i2)
!                        end if
                                
                    end do
                end do
            end do
        end do
    end do ! end loop on particles

!    open(90,file='lt_pic_density.dat')    
    print*, '#  plot_dim1 = ',  plot_dim1
    print*, '#  plot_dim2 = ',  plot_dim2

    open(90,file=file_name)

    write(90,*) '#  plot_dim1 = ',  plot_dim1
    write(90,*) '#  plot_x1_min = ',  plot_x1_min
    write(90,*) '#  plot_x1_max = num_points_1 * h_plot_1 + plot_x1_min = ',  num_points_1 * h_plot_1 + plot_x1_min
    

    write(90,*) '#  plot_dim2 = ',  plot_dim2
    write(90,*) '#  plot_x2_min = ',  plot_x2_min
    write(90,*) '#  plot_x2_max = num_points_2 * h_plot_2 + plot_x2_min = ',  num_points_2 * h_plot_2 + plot_x2_min

    write(90,*) '#  plot_cst_dim3 = ',  plot_cst_dim3
    write(90,*) '#  plot_cst_dim4 = ',  plot_cst_dim4
    write(90,*) '#  x3_plot_cst = ', x3_plot_cst, '  x4_plot_cst = ',x4_plot_cst
    write(90,*) '#  '
    write(90,*) '#  x1 (variable)        x2 (variable)        f( x1, x2, x3_plot_cst, x4_plot_cst) '

    do i1 = 1, num_points_1
        x1_plot = plot_x1_min + (i1-1) * h_plot_1
        do i2 = 1, num_points_2 
            x2_plot = plot_x2_min + (i2-1) * h_plot_2
            write(90,*) x1_plot, x2_plot, plotted_values(i1, i2)
        end do
        write(90,*) 
    end do
    write(90,*) 
    close(90)

        
  end subroutine sll_lt_pic_4d_plot_f_on_2d_grid





  ! mostly for debug 
  subroutine plot_2d_slice_remapping_grid( &
      file_name,        &
      p_group ) 

    character(len=*),                         intent(in)    :: file_name
    type(sll_lt_pic_4d_group), pointer,  intent(in)    :: p_group

    sll_int32 :: j_x, j_y, j_vx, j_vy
    sll_real64 :: x_j, y_j, vx_j, vy_j

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

    print*, "[plot_2d_slice_remapping_grid] plotting slice in file :", file_name
 
 
    open(90,file=file_name)


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

    j_y = int(0.5*number_parts_y)
    y_j = parts_y_min + (j_y-1) * h_parts_y
    j_vy = int(0.5*number_parts_vy)
    vy_j = parts_vy_min + (j_vy-1) * h_parts_vy
    write(90,*) "# j_y = ", j_y, " *** y_j = ", y_j, " *** j_vy = ", j_vy, " ***  vy_j = ", vy_j
    x_j = parts_x_min
    do j_x = 1, number_parts_x
        vx_j = parts_vx_min
        do j_vx = 1, number_parts_vx
            write(90,*) x_j, vx_j, p_group%target_values(j_x,j_y,j_vx,j_vy)
            vx_j = vx_j + h_parts_vx
        end do
        write(90,*) 
        x_j = x_j + h_parts_x
    end do
    write(90,*) 
    close(90)

        
  end subroutine plot_2d_slice_remapping_grid



  ! mostly for debug 
  subroutine plot_lt_particles( &
      file_name,        &
      p_group ) 

    character(len=*),                         intent(in)    :: file_name
    type(sll_lt_pic_4d_group), pointer,  intent(in)    :: p_group

    sll_int32 :: k
    sll_real64 :: x_k    
    sll_real64 :: y_k    
    type(sll_lt_pic_4d_particle),      pointer    :: particle


    print*, "[plot_lt_particles] plotting slice in file :", file_name
 
 
    open(90,file=file_name)

    do k = 1, p_group%number_particles
        particle => p_group%p_list(k)
        call cell_offset_to_global( particle%dx, particle%dy, particle%ic, &
                                      p_group%mesh, x_k, y_k )
        write(90,*) x_k, "  ", y_k, "  ", particle%vx , "  ", particle%vy,  "  ", particle%q

    end do
    write(90,*) 
    close(90)

        
  end subroutine plot_lt_particles



  ! <<sll_pic_shape>>
  !> sll_pic_shape(degree, x, y, vx, vy, inv_hx, inv_hy, inv_hvx, inv_hvy)
  !! computes the value of the 4d B-spline particle shape (no particle transformation here)
  !!
  !! note:
  !!  - inv_hx, inv_hy, inv_hvx, inv_hvy are the inverse
  !!    of the inter-particle distances (hx, hy, hvx, hvy) on the initial (or remapping) grid
  !!
  !!  - x, y, vx, vy are phase-space coordinates relative to the particle center
  !!    and they are not scaled: for instance if degree = 1 (resp. if degree = 3),
  !!    and |x| >= hx (resp. |x| >= 2*hx), then sll_pic_shape returns 0
  !!    (and similarly for y, vx, vy)

  function sll_pic_shape( &
    degree,             &
    x, y, vx, vy,       &
    inv_hx, inv_hy, inv_hvx, inv_hvy   &
    ) &
    result(shape)
      
    sll_int32 :: degree
    sll_real64 :: x, y, vx, vy
    sll_real64 :: inv_hx, inv_hy, inv_hvx, inv_hvy
    sll_real64 :: shape
!    sll_int32 :: nc_eta2
!    sll_int32 :: ierr

    ! uses [[sll_b_spline]]
    shape =   inv_hx  * sll_b_spline(degree,inv_hx  * x ) &
            * inv_hy  * sll_b_spline(degree,inv_hy  * y ) &
            * inv_hvx * sll_b_spline(degree,inv_hvx * vx) &
            * inv_hvy * sll_b_spline(degree,inv_hvy * vy) 
  end function sll_pic_shape


  !> added by MCP
  !! 4d-reference B-spline shape function (independent of the grid resolution), function support is [-part_cp, part_cp]^4
  !! with part_cp = (degree +1)/2
  function sll_ref_pic_shape( &
    degree,     &
    x,y,vx,vy &
    ) &
    result(ref_shape)
      
    sll_int32 :: degree
    sll_real64 :: x, y, vx, vy
    sll_real64 :: ref_shape
!    sll_int32 :: nc_eta2
!    sll_int32 :: ierr

    ref_shape = sll_b_spline(degree,x) * sll_b_spline(degree,y) * sll_b_spline(degree,vx) * sll_b_spline(degree, vy)        
  end function sll_ref_pic_shape
  
  
  ! added by MCP
  ! <<sll_b_spline>>
  ! univariate centered (and reference, ie independent of the grid resolution) B-splines. Support is ( -(degree+1)/2, (degree+1)/2 )
  function sll_b_spline( &
    degree, &
    x       &
    ) &
    result(res)
    sll_int32 :: degree
    sll_real64 :: x
    sll_real64 :: x_aux
    sll_real64 :: res
    
    res = 0
    if ( degree == 1 ) then
        ! pw affine spline
       if ( x <= -1 .or. x >= 1 ) then
            return 
        end if
        if ( x < 0 ) then
            res = x+1
            return
        end if
        res = 1-x
        return
    end if
    if ( degree == 3 ) then
        x_aux = x+2
        if ( x_aux <= 0 .or. x_aux >= 4 ) then
            return 
        end if
        if ( x_aux < 1 ) then
            res = 0.16666666666666666*(x_aux**3)
            return 
        end if
        if ( x_aux < 2 ) then
            res = 0.6666666666666666 - 2.*x_aux + 2.*(x_aux**2) - 0.5*(x_aux**3)
            return 
        end if
        if ( x_aux < 3 ) then
            res = -7.333333333333333 + 10.*x_aux - 4.*(x_aux**2) + 0.5*(x_aux**3)
            return 
        end if
        res = 10.66666666666666 - 8.*x_aux + 2.*(x_aux**2) - 0.16666666666666666*(x_aux**3)
        return        
    end if
    if ( degree == 5 ) then
        x_aux = x+3
        if ( x_aux <= 0 .or. x_aux >= 6 ) then
            res = 0
            return 
        end if
        if ( x_aux < 1 ) then
            res = 0.00833333333333333*(x_aux**5)
            return 
        end if
        if ( x_aux < 2 ) then
            res = -0.0416666666667 *(x_aux**5) + 0.25 *(x_aux**4) - 0.5 *(x_aux**3) + 0.5 *(x_aux**2) - 0.25 *(x_aux) + 0.05
            return 
        end if
        if ( x_aux < 3 ) then
            res = 0.0833333333333 *(x_aux**5) - 1.0 *(x_aux**4) + 4.5 *(x_aux**3) - 9.5 *(x_aux**2) + 9.75 *(x_aux) - 3.95
            return 
        end if
        if ( x_aux < 4 ) then
            res = -0.0833333333333 *(x_aux**5) + 1.5 *(x_aux**4) - 10.5 *(x_aux**3) + 35.5 *(x_aux**2) - 57.75 *(x_aux) + 36.55
            return 
        end if        
        if ( x_aux < 5 ) then
            res = 0.0416666666667 *(x_aux**5) - 1.0 *(x_aux**4) + 9.5 *(x_aux**3) - 44.5 *(x_aux**2) + 102.25 *(x_aux) - 91.45
            return 
        end if
        res = -0.008333333333333333 * ((x_aux-6.)**5)
        return
    end if
    print *, 'sll_b_spline(): ERROR, invalid value of argument degree ', degree
    STOP
    return 
  end function sll_b_spline


  ! <<get_initial_position_on_cartesian_grid_from_particle_index>>
  
subroutine get_initial_position_on_cartesian_grid_from_particle_index (k,                                               &
                                                                       n_parts_x, n_parts_y, n_parts_vx, n_parts_vy,    &
                                                                       j_x, j_y, j_vx, j_vy                             &
                                                                       )
    sll_int64, intent(in) :: k
    sll_int64, intent(in) :: n_parts_x
    sll_int64, intent(in) :: n_parts_y
    sll_int64, intent(in) :: n_parts_vx
    sll_int64, intent(in) :: n_parts_vy
    sll_int64, intent(out) :: j_x
    sll_int64, intent(out) :: j_y
    sll_int64, intent(out) :: j_vx
    sll_int64, intent(out) :: j_vy
    sll_int64              :: k_aux

    k_aux = k-1
    ! here, k_aux = (j_vy-1) + (j_vx-1) * n_parts_vy + (j_y-1) * n_parts_vy * n_parts_vx + (j_x-1) * n_parts_vy * n_parts_vx * n_parts_y
    j_vy = mod(k_aux, n_parts_vy) + 1

    k_aux = (k_aux - (j_vy-1)) / n_parts_vy
    ! here, k_aux = (j_vx-1) + (j_y-1) * n_parts_vx + (j_x-1) * n_parts_vx * n_parts_y
    j_vx = mod(k_aux, n_parts_vx) + 1

    k_aux = (k_aux - (j_vx-1)) / n_parts_vx
    ! here, k_aux = (j_y-1) + (j_x-1) * n_parts_y
    j_y = mod(k_aux, n_parts_y) + 1

    k_aux = (k_aux - (j_y-1)) / n_parts_y
    ! here, k_aux = (j_x-1)
    j_x = k_aux + 1

end subroutine

! <<get_particle_index_from_initial_position_on_cartesian_grid>>

subroutine get_particle_index_from_initial_position_on_cartesian_grid (j_x, j_y, j_vx, j_vy,                            &
                                                                       n_parts_x, n_parts_y, n_parts_vx, n_parts_vy,    &
                                                                       k                                                &
                                                                       )
    sll_int64, intent(in) :: j_x
    sll_int64, intent(in) :: j_y
    sll_int64, intent(in) :: j_vx
    sll_int64, intent(in) :: j_vy
    sll_int64, intent(in) :: n_parts_x
    sll_int64, intent(in) :: n_parts_y
    sll_int64, intent(in) :: n_parts_vx
    sll_int64, intent(in) :: n_parts_vy
    sll_int64, intent(out) :: k

!    k = 1+ (j_x-1) + (j_y-1) * n_parts_x + (j_vx-1) * n_parts_x * n_parts_y + (j_vy-1) * n_parts_x * n_parts_y * n_parts_vx
    k = 1+ (j_vy-1) + (j_vx-1) * n_parts_vy + (j_y-1) * n_parts_vy * n_parts_vx + (j_x-1) * n_parts_vy * n_parts_vx * n_parts_y

end subroutine


end module sll_lt_pic_4d_utilities
