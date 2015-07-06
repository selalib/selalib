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

  use sll_constants, only: sll_pi
  use sll_lt_pic_4d_group_module
  use sll_representation_conversion_module
  use sll_timer
  use sll_utilities, only: sll_new_file_id, int2string
  use sll_gnuplot

implicit none 

contains

    !  subroutine sll_first_lt_pic_charge_accumulation_4d( p_group, q_accum )
    !    type(sll_lt_pic_4d_group), pointer      :: p_group
    !    type(sll_charge_accumulator_2d), pointer     :: q_accum
    !    type(sll_lt_pic_4d_particle), dimension(:), pointer :: p
    !    sll_int32  :: i
    !    sll_int32  :: number_particles
    !    sll_real64 :: tmp1
    !    sll_real64 :: tmp2
    !
    !    SLL_ASSERT( associated(p_group) .and. associated(q_accum))
    !    number_particles =  p_group%number_particles
    !    p             => p_group%p_list
    !
    !    do i=1,number_particles
    !       SLL_ACCUMULATE_PARTICLE_CHARGE(q_accum,p(i),tmp1,tmp2)
    !    end do
    !  end subroutine sll_first_lt_pic_charge_accumulation_4d


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
                        track_markers_outside_domain, &
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
    
        sll_int32, intent(in) :: k
!        sll_int32, intent(in) :: part_degree
        type(sll_cartesian_mesh_2d),               pointer,  intent(in) :: particles_m2d  ! Poisson mesh associated with the particles
        type(sll_lt_pic_4d_particle), dimension(:),  pointer,  intent(in) :: p_list
    
        LOGICAL, intent(in) :: domain_is_x_periodic
        LOGICAL, intent(in) :: domain_is_y_periodic    
        LOGICAL, intent(in) :: track_markers_outside_domain
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
    
        sll_int32   :: k_ngb
        sll_real64  :: x_k_left,  x_k_right
        sll_real64  :: y_k_left,  y_k_right
        sll_real64  :: vx_k_left, vx_k_right
        sll_real64  :: vy_k_left, vy_k_right
    
        sll_real64  :: j11,j12,j13,j14   ! coefs of matrix J = D^-1 (forward Jacobian)
        sll_real64  :: j21,j22,j23,j24
        sll_real64  :: j31,j32,j33,j34
        sll_real64  :: j41,j42,j43,j44
        sll_real64  :: factor, det_J, inv_det_J
            
        
        call cell_offset_to_global_extended(p_list(k)%dx,   &
                                            p_list(k)%dy,   &
                                            p_list(k)%ic_x,   &
                                            p_list(k)%ic_y,   &
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
            call cell_offset_to_global_extended(p_list(k_ngb)%dx,   &
                                                p_list(k_ngb)%dy,   &
                                                p_list(k_ngb)%ic_x,   &
                                                p_list(k_ngb)%ic_y,   &
                                                particles_m2d, x_k_right, y_k_right )
            vx_k_right = p_list(k_ngb)%vx
            vy_k_right = p_list(k_ngb)%vy        
            if( .not. track_markers_outside_domain )then
                if( domain_is_x_periodic .and. x_k_right < x_k - 0.5*mesh_period_x ) x_k_right = x_k_right + mesh_period_x
                if( domain_is_x_periodic .and. x_k_right > x_k + 0.5*mesh_period_x ) x_k_right = x_k_right - mesh_period_x
                if( domain_is_y_periodic .and. y_k_right < y_k - 0.5*mesh_period_y ) y_k_right = y_k_right + mesh_period_y
                if( domain_is_y_periodic .and. y_k_right > y_k + 0.5*mesh_period_y ) y_k_right = y_k_right - mesh_period_y
            end if
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
            call cell_offset_to_global_extended(p_list(k_ngb)%dx,   &
                                                p_list(k_ngb)%dy,   &
                                                p_list(k_ngb)%ic_x,   &
                                                p_list(k_ngb)%ic_y,   &
                                                particles_m2d, x_k_left, y_k_left )
            vx_k_left  = p_list(k_ngb)%vx
            vy_k_left  = p_list(k_ngb)%vy
            if( .not. track_markers_outside_domain )then
                if( domain_is_x_periodic .and. x_k_left < x_k - 0.5*mesh_period_x ) x_k_left = x_k_left + mesh_period_x
                if( domain_is_x_periodic .and. x_k_left > x_k + 0.5*mesh_period_x ) x_k_left = x_k_left - mesh_period_x
                if( domain_is_y_periodic .and. y_k_left < y_k - 0.5*mesh_period_y ) y_k_left = y_k_left + mesh_period_y
                if( domain_is_y_periodic .and. y_k_left > y_k + 0.5*mesh_period_y ) y_k_left = y_k_left - mesh_period_y
            end if
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
            call cell_offset_to_global_extended(p_list(k_ngb)%dx,   &
                                                p_list(k_ngb)%dy,   &
                                                p_list(k_ngb)%ic_x,   &
                                                p_list(k_ngb)%ic_y,   &
                                                particles_m2d, x_k_right, y_k_right )
            vx_k_right = p_list(k_ngb)%vx
            vy_k_right = p_list(k_ngb)%vy
            if( .not. track_markers_outside_domain )then
                if( domain_is_x_periodic .and. x_k_right < x_k - 0.5*mesh_period_x ) x_k_right = x_k_right + mesh_period_x
                if( domain_is_x_periodic .and. x_k_right > x_k + 0.5*mesh_period_x ) x_k_right = x_k_right - mesh_period_x
                if( domain_is_y_periodic .and. y_k_right < y_k - 0.5*mesh_period_y ) y_k_right = y_k_right + mesh_period_y
                if( domain_is_y_periodic .and. y_k_right > y_k + 0.5*mesh_period_y ) y_k_right = y_k_right - mesh_period_y
            end if
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
            call cell_offset_to_global_extended(p_list(k_ngb)%dx,   &
                                                p_list(k_ngb)%dy,   &
                                                p_list(k_ngb)%ic_x,   &
                                                p_list(k_ngb)%ic_y,   &
                                                particles_m2d, x_k_left, y_k_left )
            vx_k_left  = p_list(k_ngb)%vx
            vy_k_left  = p_list(k_ngb)%vy
            if( .not. track_markers_outside_domain )then
                if( domain_is_x_periodic .and. x_k_left < x_k - 0.5*mesh_period_x ) x_k_left = x_k_left + mesh_period_x
                if( domain_is_x_periodic .and. x_k_left > x_k + 0.5*mesh_period_x ) x_k_left = x_k_left - mesh_period_x
                if( domain_is_y_periodic .and. y_k_left < y_k - 0.5*mesh_period_y ) y_k_left = y_k_left + mesh_period_y
                if( domain_is_y_periodic .and. y_k_left > y_k + 0.5*mesh_period_y ) y_k_left = y_k_left - mesh_period_y
            end if
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
            call cell_offset_to_global_extended(p_list(k_ngb)%dx,   &
                                                p_list(k_ngb)%dy,   &
                                                p_list(k_ngb)%ic_x,   &
                                                p_list(k_ngb)%ic_y,   &
                                                particles_m2d, x_k_right, y_k_right )
            vx_k_right = p_list(k_ngb)%vx
            vy_k_right = p_list(k_ngb)%vy
            if( .not. track_markers_outside_domain )then
                if( domain_is_x_periodic .and. x_k_right < x_k - 0.5*mesh_period_x ) x_k_right = x_k_right + mesh_period_x
                if( domain_is_x_periodic .and. x_k_right > x_k + 0.5*mesh_period_x ) x_k_right = x_k_right - mesh_period_x
                if( domain_is_y_periodic .and. y_k_right < y_k - 0.5*mesh_period_y ) y_k_right = y_k_right + mesh_period_y
                if( domain_is_y_periodic .and. y_k_right > y_k + 0.5*mesh_period_y ) y_k_right = y_k_right - mesh_period_y
            end if
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
            call cell_offset_to_global_extended(p_list(k_ngb)%dx,   &
                                                p_list(k_ngb)%dy,   &
                                                p_list(k_ngb)%ic_x,   &
                                                p_list(k_ngb)%ic_y,   &
                                                particles_m2d, x_k_left, y_k_left )
            vx_k_left  = p_list(k_ngb)%vx
            vy_k_left  = p_list(k_ngb)%vy
            if( .not. track_markers_outside_domain )then
                if( domain_is_x_periodic .and. x_k_left < x_k - 0.5*mesh_period_x ) x_k_left = x_k_left + mesh_period_x
                if( domain_is_x_periodic .and. x_k_left > x_k + 0.5*mesh_period_x ) x_k_left = x_k_left - mesh_period_x
                if( domain_is_y_periodic .and. y_k_left < y_k - 0.5*mesh_period_y ) y_k_left = y_k_left + mesh_period_y
                if( domain_is_y_periodic .and. y_k_left > y_k + 0.5*mesh_period_y ) y_k_left = y_k_left - mesh_period_y
            end if
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
            call cell_offset_to_global_extended(p_list(k_ngb)%dx,   &
                                                p_list(k_ngb)%dy,   &
                                                p_list(k_ngb)%ic_x,   &
                                                p_list(k_ngb)%ic_y,   &
                                                particles_m2d, x_k_right, y_k_right )
            vx_k_right = p_list(k_ngb)%vx
            vy_k_right = p_list(k_ngb)%vy
            if( .not. track_markers_outside_domain )then
                if( domain_is_x_periodic .and. x_k_right < x_k - 0.5*mesh_period_x ) x_k_right = x_k_right + mesh_period_x
                if( domain_is_x_periodic .and. x_k_right > x_k + 0.5*mesh_period_x ) x_k_right = x_k_right - mesh_period_x
                if( domain_is_y_periodic .and. y_k_right < y_k - 0.5*mesh_period_y ) y_k_right = y_k_right + mesh_period_y
                if( domain_is_y_periodic .and. y_k_right > y_k + 0.5*mesh_period_y ) y_k_right = y_k_right - mesh_period_y
            end if
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
            call cell_offset_to_global_extended(p_list(k_ngb)%dx,   &
                                                p_list(k_ngb)%dy,   &
                                                p_list(k_ngb)%ic_x,   &
                                                p_list(k_ngb)%ic_y,   &
                                                particles_m2d, x_k_left, y_k_left )
            vx_k_left  = p_list(k_ngb)%vx
            vy_k_left  = p_list(k_ngb)%vy
            if( .not. track_markers_outside_domain )then
                if( domain_is_x_periodic .and. x_k_left < x_k - 0.5*mesh_period_x ) x_k_left = x_k_left + mesh_period_x
                if( domain_is_x_periodic .and. x_k_left > x_k + 0.5*mesh_period_x ) x_k_left = x_k_left - mesh_period_x
                if( domain_is_y_periodic .and. y_k_left < y_k - 0.5*mesh_period_y ) y_k_left = y_k_left + mesh_period_y
                if( domain_is_y_periodic .and. y_k_left > y_k + 0.5*mesh_period_y ) y_k_left = y_k_left - mesh_period_y
            end if
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


  function in_bounds_periodic_x( x, mesh, x_periodic ) result(res)

    use sll_cartesian_meshes

    ! [[file:../working_precision/sll_working_precision.h]]
    use sll_working_precision

    logical :: res
    sll_real64, intent(in) :: x
    logical, intent(in) :: x_periodic
    type(sll_cartesian_mesh_2d), pointer :: mesh
    if( (x >= mesh%eta1_min)                                                                    &
        .and.                                                                                   &
        ((x < mesh%eta1_max .and. x_periodic) .or. (x <= mesh%eta1_max .and. .not. x_periodic)) &
      ) then
       res = .true.
    else
       res = .false.
    end if
  end function in_bounds_periodic_x


  function in_bounds_periodic_y( y, mesh, y_periodic ) result(res)

    use sll_cartesian_meshes

    ! [[file:../working_precision/sll_working_precision.h]]
    use sll_working_precision

    logical :: res
    sll_real64, intent(in) :: y
    logical, intent(in) :: y_periodic
    type(sll_cartesian_mesh_2d), pointer :: mesh
    if( (y >= mesh%eta2_min)                                                                    &
        .and.                                                                                   &
        ((y < mesh%eta2_max .and. y_periodic) .or. (y <= mesh%eta2_max .and. .not. y_periodic)) &
      ) then
       res = .true.
    else
       res = .false.
    end if
  end function in_bounds_periodic_y

  ! <<apply_periodic_bc>> extracted from [[file:../simulation/simulation_4d_vp_lt_pic_cartesian.F90::subroutine
  ! apply_periodic_bc]]

  subroutine apply_periodic_bc_old( mesh, x, y )

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
  end subroutine apply_periodic_bc_old

  subroutine apply_periodic_bc( mesh, x, y )

    use sll_cartesian_meshes
    ! [[file:../working_precision/sll_working_precision.h]]
    use sll_working_precision

    type(sll_cartesian_mesh_2d), pointer :: mesh
    sll_real64, intent(inout) :: x
    sll_real64, intent(inout) :: y

    x = mesh%eta1_min + modulo(x - mesh%eta1_min, mesh%eta1_max - mesh%eta1_min)
    y = mesh%eta2_min + modulo(y - mesh%eta2_min, mesh%eta2_max - mesh%eta2_min)
  end subroutine apply_periodic_bc


  subroutine apply_periodic_bc_x( mesh, x)

    use sll_cartesian_meshes
    ! [[file:../working_precision/sll_working_precision.h]]
    use sll_working_precision

    type(sll_cartesian_mesh_2d), pointer :: mesh
    sll_real64, intent(inout) :: x

    x = mesh%eta1_min + modulo(x - mesh%eta1_min, mesh%eta1_max - mesh%eta1_min)
  end subroutine apply_periodic_bc_x


  subroutine apply_periodic_bc_y( mesh, y)

    use sll_cartesian_meshes
    ! [[file:../working_precision/sll_working_precision.h]]
    use sll_working_precision

    type(sll_cartesian_mesh_2d), pointer, intent(in) :: mesh
    sll_real64, intent(inout) :: y

    y = mesh%eta2_min + modulo(y - mesh%eta2_min, mesh%eta2_max - mesh%eta2_min)
  end subroutine apply_periodic_bc_y


  ! puts the point (x,y) back into the computational domain if periodic in x or y (or both)
  ! otherwise, does nothing
  subroutine periodic_correction(p_group, x, y)
    type(sll_lt_pic_4d_group), pointer,  intent(in)    :: p_group
    sll_real64, intent(inout) :: x
    sll_real64, intent(inout) :: y
    type(sll_cartesian_mesh_2d), pointer :: mesh

    mesh => p_group%mesh

    if( p_group%domain_is_x_periodic                        &
        .and.                                               &
        ( (x < mesh%eta1_min) .or. (x >= mesh%eta1_max) )   &
      ) then
          call apply_periodic_bc_x( mesh, x)
    end if
    if( p_group%domain_is_y_periodic                        &
        .and.                                               &
        ( (y < mesh%eta2_min) .or. (y >= mesh%eta2_max) )   &
      ) then
          call apply_periodic_bc_y( mesh, y)
    end if
  end subroutine periodic_correction


  ! <<onestep>> <<ALH>> utility function for finding the neighbours of a particle, used by [[ONESTEPMACRO]]. "dim"
  ! corresponds to one of x,y,vx,vy.

  subroutine onestep(dim,dim_t0,kprime,p_list,h_parts_dim)

    sll_int :: dim
    sll_real64 :: dim_t0
    sll_int32 :: neighbour

    ! [[file:~/mcp/selalib/src/pic_particle_types/lt_pic_4d_group.F90::sll_lt_pic_4d_group-p_list]]
    type(sll_lt_pic_4d_particle), dimension(:), pointer,intent(in) :: p_list
    
    sll_int32 :: ngb_dim_right_index
    sll_int32 :: ngb_dim_left_index
    sll_real64 :: h_parts_dim
    sll_int32 :: kprime
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

  ! note (March 25): a new version of this routine is being written, that also allows to deposit the charge
  !       when finished the new routine (called sll_lt_pic_4d_remap_or_deposit_f) should replace this one.

!  subroutine sll_lt_pic_4d_write_bsl_f_on_remap_grid (p_group)
!
!    ! [[file:../pic_particle_types/lt_pic_4d_group.F90::sll_lt_pic_4d_group]] p_group contains both the existing
!    ! particles and the virtual remapping grid
!
!    type(sll_lt_pic_4d_group),pointer,intent(inout) :: p_group
!
!    ! cf [[file:~/mcp/maltpic/ltpic-bsl.tex::N*]]
!
!    sll_int32 :: n_virtual = 2 ! <<n_virtual>>
!    sll_int32 :: num_virtual_cells_x
!    sll_int32 :: num_virtual_cells_y
!    sll_int32 :: num_virtual_cells_vx
!    sll_int32 :: num_virtual_cells_vy
!
!    ! [[file:~/mcp/maltpic/ltpic-bsl.tex::h_parts_x]] and h_parts_y, h_parts_vx, h_parts_vy
!
!    sll_real64 :: h_parts_x
!    sll_real64 :: h_parts_y
!    sll_real64 :: h_parts_vx
!    sll_real64 :: h_parts_vy
!
!    sll_real64 :: inv_h_parts_x
!    sll_real64 :: inv_h_parts_y
!    sll_real64 :: inv_h_parts_vx
!    sll_real64 :: inv_h_parts_vy
!
!    sll_real64 :: parts_x_min
!    sll_real64 :: parts_y_min
!    sll_real64 :: parts_vx_min
!    sll_real64 :: parts_vy_min
!
!    ! same as \delta{x,y,vx,vy} in [[file:~/mcp/maltpic/ltpic-bsl.tex::h_parts_x]]
!    sll_real64 :: h_virtual_cell_x
!    sll_real64 :: h_virtual_cell_y
!    sll_real64 :: h_virtual_cell_vx
!    sll_real64 :: h_virtual_cell_vy
!
!    sll_real64 :: x
!    sll_real64 :: y
!    sll_real64 :: vx
!    sll_real64 :: vy
!
!    ! working values (offsets in virtual cell)
!
!    sll_real32 :: dx
!    sll_real32 :: dy
!    sll_real32 :: dvx
!    sll_real32 :: dvy
!
!    ! working space
!
!    sll_real64 :: tmp
!
!    ! index of particle closest to the center of each virtual cell. Array dimensions defined by the contents of
!    ! [[file:../pic_particle_types/lt_pic_4d_group.F90::sll_lt_pic_4d_group-remapping_grid]]. If [[n_virtual]] is
!    ! greater than 1, the size of this array is smaller than the number of real remapping_grid cells.
!
!    sll_int64,dimension(:,:,:,:),allocatable :: closest_particle
!    sll_real64,dimension(:,:,:,:),allocatable :: closest_particle_distance
!
!    sll_int32 :: i ! x dimension
!    sll_int32 :: j ! y dimension
!    sll_int64 :: k,kprime ! particle index
!    sll_int64 :: neighbour ! particle index for local use
!    sll_int32 :: l ! vx dimension
!    sll_int32 :: m ! vy dimension
!
!    ! indices in a virtual cell (go from 1 to [[n_virtual]])
!
!    sll_int :: ivirt ! x dimension
!    sll_int :: jvirt ! y dimension
!    sll_int :: lvirt ! vx dimension
!    sll_int :: mvirt ! vy dimension
!
!    sll_int :: i_x,i_y,i_vx,i_vy
!
!    ! particle pointer [[file:../pic_particle_types/lt_pic_4d_particle.F90::sll_lt_pic_4d_particle]]
!    type(sll_lt_pic_4d_particle),pointer :: p
!
!    ! <<g>> cartesian grid pointer to
!    ! [[file:../pic_particle_types/lt_pic_4d_group.F90::sll_lt_pic_4d_group-remapping_grid]]
!
!    type(sll_cartesian_mesh_4d),pointer :: g
!
!    ! periodicity
!
!    LOGICAL :: domain_is_x_periodic
!    LOGICAL :: domain_is_y_periodic
!    sll_real64 :: mesh_period_x
!    sll_real64 :: mesh_period_y
!    sll_real64 :: inv_period_x
!    sll_real64 :: inv_period_y
!
!    ! results from [[get_ltp_deformation_matrix]]
!
!    sll_real64 :: d11,d12,d13,d14 ! coefs of matrix D (backward Jacobian)
!    sll_real64 :: d21,d22,d23,d24
!    sll_real64 :: d31,d32,d33,d34
!    sll_real64 :: d41,d42,d43,d44
!
!    ! coordinates of particle k at time n and time 0
!    sll_real64 :: x_k,y_k,vx_k,vy_k
!    sll_real64 :: x_k_t0,y_k_t0,vx_k_t0,vy_k_t0
!
!    sll_real64 :: part_radius_x
!    sll_real64 :: part_radius_y
!    sll_real64 :: part_radius_vx
!    sll_real64 :: part_radius_vy
!
!    ! coordinates of a virtual particle at time 0 relative to the coordinates of one real particle
!
!    sll_real64 :: x_t0,y_t0,vx_t0,vy_t0
!
!    sll_int32 :: part_degree
!
!    sll_int32 :: ierr
!
!    ! temporary workspace
!    sll_real64 :: x_aux
!    sll_real64 :: y_aux
!    sll_real64 :: vx_aux
!    sll_real64 :: vy_aux
!
!    sll_real64 :: length
!
!    sll_real64 :: x_kprime_t0
!    sll_real64 :: y_kprime_t0
!    sll_real64 :: vx_kprime_t0
!    sll_real64 :: vy_kprime_t0
!
!    ! value 1 or 2 points to each side of an hypercube in direction x,y,vx or vy
!    sll_int :: side_x,side_y,side_vx,side_vy
!    sll_int64,dimension(2,2,2,2) :: hcube
!
!    sll_int64 :: j_x,j_y,j_vx,j_vy
!
!    sll_int64 :: number_parts_x
!    sll_int64 :: number_parts_y
!    sll_int64 :: number_parts_vx
!    sll_int64 :: number_parts_vy
!
!    ! --- end of declarations
!
!    g => p_group%remapping_grid
!    part_degree = p_group%spline_degree
!
!    ! Preparatory work: find out the particle which is closest to each cell center by looping over all particles and
!    ! noting which virtual cell contains it. The leftmost virtual cell in each dimension may not be complete.
!
!    num_virtual_cells_x =  int(g%num_cells1/n_virtual)+1
!    num_virtual_cells_y =  int(g%num_cells2/n_virtual)+1
!    num_virtual_cells_vx = int(g%num_cells3/n_virtual)+1
!    num_virtual_cells_vy = int(g%num_cells4/n_virtual)+1
!
!    SLL_ALLOCATE(closest_particle(num_virtual_cells_x,num_virtual_cells_y,num_virtual_cells_vx,num_virtual_cells_vy),ierr)
!    closest_particle(:,:,:,:) = 0
!
!    SLL_ALLOCATE(closest_particle_distance(num_virtual_cells_x,num_virtual_cells_y,num_virtual_cells_vx,num_virtual_cells_vy),ierr)
!    closest_particle_distance(:,:,:,:) = 0
!
!    ! remapping grid cell size - same as in [[write_f_on_remap_grid-h_parts_x]]
!
!    h_parts_x    = g%delta_eta1
!    h_parts_y    = g%delta_eta2
!    h_parts_vx   = g%delta_eta3
!    h_parts_vy   = g%delta_eta4
!
!    inv_h_parts_x  = 1./h_parts_x
!    inv_h_parts_y  = 1./h_parts_y
!    inv_h_parts_vx = 1./h_parts_vx
!    inv_h_parts_vy = 1./h_parts_vy
!
!    parts_x_min    = p_group%remapping_grid%eta1_min
!    parts_y_min    = p_group%remapping_grid%eta2_min
!    parts_vx_min   = p_group%remapping_grid%eta3_min
!    parts_vy_min   = p_group%remapping_grid%eta4_min
!
!    number_parts_x = p_group%number_parts_x
!    number_parts_y = p_group%number_parts_y
!    number_parts_vx = p_group%number_parts_vx
!    number_parts_vy = p_group%number_parts_vy
!
!    ! virtual cell size
!
!    h_virtual_cell_x  = n_virtual * h_parts_x
!    h_virtual_cell_y  = n_virtual * h_parts_y
!    h_virtual_cell_vx = n_virtual * h_parts_vx
!    h_virtual_cell_vy = n_virtual * h_parts_vy
!
!    ! preparatory loop to fill the [[closest_particle]] array containing the particle closest to the center of each
!    ! virtual cell
!
!    do k=1,p_group%number_particles ! [[file:../pic_particle_types/lt_pic_4d_group.F90::number_particles]]
!       p => p_group%p_list(k)
!
!       ! find absolute (x,y,vx,vy) coordinates for this particle. Uses
!       ! [[file:sll_representation_conversion.F90::cell_offset_to_global]]
!
!       call cell_offset_to_global(p%dx,p%dy,p%ic,p_group%mesh,x,y)
!       vx = p%vx
!       vy = p%vy
!
!       ! which _virtual_ cell is this particle in? uses
!       ! [[file:sll_representation_conversion.F90::compute_cell_and_offset]] and [[g]]
!
!       x_aux = x - g%eta1_min
!       i = int( x_aux / h_virtual_cell_x ) + 1
!
!       y_aux = y - g%eta2_min
!       j = int( y_aux / h_virtual_cell_y ) + 1
!
!       vx_aux = vx - g%eta3_min
!       l = int( vx_aux / h_virtual_cell_vx ) + 1
!
!       vy_aux = vy - g%eta4_min
!       m = int( vy_aux / h_virtual_cell_vy ) + 1
!
!       ! discard particles in virtual cells off-bounds
!       if(  i >= 1 .and. i <= num_virtual_cells_x .and. &
!            j >= 1 .and. j <= num_virtual_cells_y .and. &
!            l >= 1 .and. l <= num_virtual_cells_vx .and. &
!            m >= 1 .and. m <= num_virtual_cells_vy  )then
!
!          call update_closest_particle_arrays(k,                         &
!                                              x_aux, y_aux, vx_aux, vy_aux,   &
!                                              i, j, l, m,                     &
!                                              h_virtual_cell_x, h_virtual_cell_y, h_virtual_cell_vx, h_virtual_cell_vy,   &
!                                              closest_particle,               &
!                                              closest_particle_distance)
!
!        end if
!
!    end do
!
!    ! Periodicity treatments copied from [[sll_lt_pic_4d_write_f_on_remap_grid-periodicity]]
!
!    domain_is_x_periodic = .true.   ! temp
!    domain_is_y_periodic = .true.   ! temp
!
!    if(domain_is_x_periodic) then
!      ! here the domain corresponds to the Poisson mesh
!      mesh_period_x = p_group%mesh%eta1_max - p_group%mesh%eta1_min
!      inv_period_x = 1./mesh_period_x
!    else
!      mesh_period_x = 0
!      inv_period_x = 0
!    end if
!
!    if(domain_is_y_periodic) then
!      ! here the domain corresponds to the Poisson mesh
!      mesh_period_y = p_group%mesh%eta2_max - p_group%mesh%eta2_min
!      inv_period_y = 1./mesh_period_y
!    else
!      mesh_period_y = 0
!      inv_period_y = 0
!    end if
!
!    ! initialize [[file:../pic_particle_types/lt_pic_4d_group.F90::target_values]]
!
!    p_group%target_values(:,:,:,:) = 0
!
!    ! MCP: [DEBUG] store the (computed) absolute initial position of the virtual particle
!    p_group%debug_bsl_remap = -100
!
!
!    ! <<loop_on_virtual_cells>> [[file:~/mcp/maltpic/ltpic-bsl.tex::algo:pic-vr:loop_over_all_cells]]
!    ! Loop over all cells of indices i,j,l,m which contain at least one particle
!
!    do i = 1,num_virtual_cells_x
!       do j = 1,num_virtual_cells_y
!          do l = 1,num_virtual_cells_vx
!             do m = 1,num_virtual_cells_vy
!
!                ! [[file:~/mcp/maltpic/ltpic-bsl.tex::algo:pic-vr:create_virtual_particles]] Create a temporary set of
!                ! virtual particles inside the cell.  Note: in our case the virtual particles coincide with the existing
!                ! remapping_grid [[file:../pic_particle_types/lt_pic_4d_group.F90::sll_lt_pic_4d_group-remapping_grid]]
!                ! defined in p_group [[file:../pic_particle_types/lt_pic_4d_group.F90::sll_lt_pic_4d_group]]. So nothing
!                ! more to do.
!
!                ! [[file:~/mcp/maltpic/ltpic-bsl.tex::algo:pic-vr:find_closest_real_particle]] Find the real particle
!                ! which is closest to the cell center.  Note: speed-wise, it may be necessary to find a way not to scan
!                ! all the particles for every cell.  We avoid scanning all the particles for each cell by using the
!                ! precomputed array [[closest_particle]]. Virtual cells which do not contain any particle are skipped.
!
!                k = closest_particle(i,j,l,m)
!                if(k /= 0) then
!
!                   ! [[file:~/mcp/maltpic/ltpic-bsl.tex::hat-bz*]] Compute backward image of l-th virtual node by the
!                   ! k-th backward flow. MCP -> oui, avec la matrice de deformation calculée avec la fonction
!                   ! [[get_ltp_deformation_matrix]] pour la particule k. Calling [[get_ltp_deformation_matrix]]
!                   ! with parameters inspired from [[sll_lt_pic_4d_write_f_on_remap_grid-get_ltp_deformation_matrix]]
!
!                   call get_ltp_deformation_matrix (               &
!                        k,                                         &
!                        p_group%mesh,                              &
!                        p_group%p_list,                            &
!                        domain_is_x_periodic,                      &
!                        domain_is_y_periodic,                      &
!                        mesh_period_x,                             &
!                        mesh_period_y,                             &
!                        h_parts_x,                                 &
!                        h_parts_y,                                 &
!                        h_parts_vx,                                &
!                        h_parts_vy,                                &
!                        1./h_parts_x,                              &
!                        1./h_parts_y,                              &
!                        1./h_parts_vx,                             &
!                        1./h_parts_vy,                             &
!                        0.5_f64*(part_degree+1),                   &
!                        x_k,y_k,vx_k,vy_k,                         &
!                        d11,d12,d13,d14,                           &
!                        d21,d22,d23,d24,                           &
!                        d31,d32,d33,d34,                           &
!                        d41,d42,d43,d44,                           &
!                        part_radius_x,                             &
!                        part_radius_y,                             &
!                        part_radius_vx,                            &
!                        part_radius_vy                             &
!                        )
!
!                   ! Find position of particle k at time 0
!                   ! [[get_initial_position_on_cartesian_grid_from_particle_index]]
!
!                   call get_initial_position_on_cartesian_grid_from_particle_index(k, &
!                        number_parts_x,number_parts_y,number_parts_vx,number_parts_vy, &
!                        j_x,j_y,j_vx,j_vy)
!                   x_k_t0 =  parts_x_min  + (j_x-1)  * h_parts_x
!                   y_k_t0 =  parts_y_min  + (j_y-1)  * h_parts_y
!                   vx_k_t0 = parts_vx_min + (j_vx-1) * h_parts_vx
!                   vy_k_t0 = parts_vy_min + (j_vy-1) * h_parts_vy
!
!                   ! <<loop_on_virtual_particles_in_one_virtual_cell>>
!                   ! [[file:~/mcp/maltpic/ltpic-bsl.tex::algo:pic-vr:find_f0_for_each_virtual_particle]] Loop over all
!                   ! virtual particles in the cell to compute the value of f0 at that point (Following
!                   ! [[file:~/mcp/maltpic/ltpic-bsl.tex::BSL_remapping_algo]])
!
!                   do ivirt = 1,n_virtual
!                      do jvirt = 1,n_virtual
!                         do lvirt = 1,n_virtual
!                            do mvirt = 1,n_virtual
!
!                               ! real index of the virtual particle in
!                               ! [[file:../pic_particle_types/lt_pic_4d_group.F90::target_values]]
!
!                               i_x =  (i-1)*n_virtual + ivirt
!                               SLL_ASSERT(i_x>0)
!                               i_y =  (j-1)*n_virtual + jvirt
!                               SLL_ASSERT(i_y>0)
!                               i_vx = (l-1)*n_virtual + lvirt
!                               SLL_ASSERT(i_vx>0)
!                               i_vy = (m-1)*n_virtual + mvirt
!                               SLL_ASSERT(i_vy>0)
!
!                               ! The index may go out of the domain for higher values of x,y,vx,vy in each dimension
!                               ! (because the corners of the corresponding virtual cell do not correspond to existing
!                               ! real particles). In that case, just ignore that value.
!
!                               if(i_x<=p_group%number_parts_x            &
!                                    .and. i_y<=p_group%number_parts_y    &
!                                    .and. i_vx<=p_group%number_parts_vx  &
!                                    .and. i_vy<=p_group%number_parts_vy) then
!
!                                  ! Location of virtual particle (ivirt,jvirt,lvirt,mvirt) at time n
!
!                                  x =  parts_x_min  + (i-1)*h_virtual_cell_x  + (ivirt-1)*h_parts_x
!                                  y =  parts_y_min  + (j-1)*h_virtual_cell_y  + (jvirt-1)*h_parts_y
!                                  vx = parts_vx_min + (l-1)*h_virtual_cell_vx + (lvirt-1)*h_parts_vx
!                                  vy = parts_vy_min + (m-1)*h_virtual_cell_vy + (mvirt-1)*h_parts_vy
!
!                                  ! particle k has to be inside the current virtual cell
!
!                                  SLL_ASSERT(abs(x-x_k) < h_virtual_cell_x)
!                                  SLL_ASSERT(abs(y-y_k) < h_virtual_cell_y)
!                                  SLL_ASSERT(abs(vx-vx_k) < h_virtual_cell_vx)
!                                  SLL_ASSERT(abs(vy-vy_k) < h_virtual_cell_vy)
!
!                                  ! Location of virtual particle (ivirt,jvirt,lvirt,mvirt) at time 0 _relative_ to the
!                                  ! position of particle k at time 0 (x_k_t0,y_k_t0,vx_k_t0,vy_k_t0) according to flow
!                                  ! deformation
!
!                                  x_t0  = d11 * (x - x_k) + d12 * (y - y_k) + d13 * (vx - vx_k) + d14 * (vy - vy_k)
!                                  y_t0  = d21 * (x - x_k) + d22 * (y - y_k) + d23 * (vx - vx_k) + d24 * (vy - vy_k)
!                                  vx_t0 = d31 * (x - x_k) + d32 * (y - y_k) + d33 * (vx - vx_k) + d34 * (vy - vy_k)
!                                  vy_t0 = d41 * (x - x_k) + d42 * (y - y_k) + d43 * (vx - vx_k) + d44 * (vy - vy_k)
!
!                                  ! MCP: [DEBUG] store the (computed) absolute initial position of the virtual particle
!                                  p_group%debug_bsl_remap(i_x,i_y,i_vx,i_vy,1,1) = x_k_t0 + x_t0
!                                  p_group%debug_bsl_remap(i_x,i_y,i_vx,i_vy,2,1) = y_k_t0 + y_t0
!                                  p_group%debug_bsl_remap(i_x,i_y,i_vx,i_vy,3,1) = vx_k_t0 + vx_t0
!                                  p_group%debug_bsl_remap(i_x,i_y,i_vx,i_vy,4,1) = vy_k_t0 + vy_t0
!
!                                  ! Minimize the path to take along periodic boundaries
!
!                                  if(domain_is_x_periodic) then
!                                     if(x_t0 > mesh_period_x/2) then
!                                        x_t0 = x_t0 - mesh_period_x
!                                     else
!                                        if(x_t0 < -mesh_period_x/2) then
!                                           x_t0 = x_t0 + mesh_period_x
!                                        end if
!                                     end if
!                                  endif
!
!                                  if(domain_is_y_periodic) then
!                                     if(y_t0 > mesh_period_y/2) then
!                                        y_t0 = y_t0 - mesh_period_y
!                                     else
!                                        if(y_t0 < -mesh_period_y/2) then
!                                           y_t0 = y_t0 + mesh_period_y
!                                        end if
!                                     end if
!                                  endif
!
!                                  ! MCP: [DEBUG] store the (computed) absolute initial position of the virtual particle
!                                  p_group%debug_bsl_remap(i_x,i_y,i_vx,i_vy,1,2) = x_k_t0 + x_t0
!                                  p_group%debug_bsl_remap(i_x,i_y,i_vx,i_vy,2,2) = y_k_t0 + y_t0
!                                  p_group%debug_bsl_remap(i_x,i_y,i_vx,i_vy,3,2) = vx_k_t0 + vx_t0
!                                  p_group%debug_bsl_remap(i_x,i_y,i_vx,i_vy,4,2) = vy_k_t0 + vy_t0
!
!                                  ! [[file:~/mcp/maltpic/ltpic-bsl.tex::neighbors-grid-0]] find the neighbours of the
!                                  ! virtual particle (ivirt,jvirt,lvirt,mvirt) at time 0 through the "logical
!                                  ! neighbours" pointers of particle k. To reduce the amount of code, start with finding
!                                  ! the closest neighbour which has lower coordinates in all directions. The particle
!                                  ! located at (x_t0,y_t0,vx_t0,vy_t0) (coordinates relative to particle k to start
!                                  ! with) gets progressively closer to kprime step by step (ie from neighbour to
!                                  ! neighbour).
!
!                                  kprime = k
!
!                                  ! Calls [[onestep]]. "dim" can be x,y,vx,vy. cf
!                                  ! [[file:../pic_particle_types/lt_pic_4d_particle.F90::sll_lt_pic_4d_particle]] for
!                                  ! pointers to neighbours.
!
!#define ONESTEPMACRO(dimpos,dimname) call onestep(dimpos,dimname/**/_t0,kprime,p_group%p_list,h_parts_/**/dimname)
!
!                                  ONESTEPMACRO(ALONG_X,x)
!                                  ONESTEPMACRO(ALONG_Y,y)
!                                  ONESTEPMACRO(ALONG_VX,vx)
!                                  ONESTEPMACRO(ALONG_VY,vy)
!
!                                  ! If we end up with kprime == 0, it means that we have not found a cell that contains
!                                  ! the particle so we just set that particle value to zero
!
!                                  if (kprime /= 0) then
!
!                                     ! kprime is the left-most vertex of the hypercube. find all the other vertices
!                                     ! through the neighbour pointers in
!                                     ! [[file:../pic_particle_types/lt_pic_4d_particle.F90::neighbour_pointers]]
!
!                                     hcube(1,1,1,1) = kprime
!
!                                     hcube(2,1,1,1) = p_group%p_list(kprime)%ngb_xright_index ! 1 step
!                                     hcube(1,2,1,1) = p_group%p_list(kprime)%ngb_yright_index
!                                     hcube(1,1,2,1) = p_group%p_list(kprime)%ngb_vxright_index
!                                     hcube(1,1,1,2) = p_group%p_list(kprime)%ngb_vyright_index
!
!                                     ! if any of the first four vertices is undefined (the convention in
!                                     ! [[file:~/mcp/selalib/src/pic_particle_initializers/lt_pic_4d_init.F90::sll_lt_pic_4d_compute_new_particles]]
!                                     ! is that the neighbour index is then equal to the particle index), it means that
!                                     ! we reached the mesh border. just set the value of f for that particle as zero as
!                                     ! before.
!
!                                     if (hcube(2,1,1,1) /= kprime        &
!                                          .and. hcube(1,2,1,1) /= kprime &
!                                          .and. hcube(1,1,2,1) /= kprime &
!                                          .and. hcube(1,1,1,2) /= kprime) then
!
!                                        ! remaining vertices of the hypercube. they should all exist now that the first
!                                        ! 4 vertices are checked.
!
!                                        ! 1 step in x + 1 other step
!                                        hcube(2,2,1,1) = p_group%p_list(hcube(2,1,1,1))%ngb_yright_index
!                                        hcube(2,1,2,1) = p_group%p_list(hcube(2,1,1,1))%ngb_vxright_index
!                                        hcube(2,1,1,2) = p_group%p_list(hcube(2,1,1,1))%ngb_vyright_index
!
!                                        ! 1 step in y + 1 other step
!                                        hcube(1,2,2,1) = p_group%p_list(hcube(1,2,1,1))%ngb_vxright_index
!                                        hcube(1,2,1,2) = p_group%p_list(hcube(1,2,1,1))%ngb_vyright_index
!
!                                        ! 1 step in vx + 1 other step
!                                        hcube(1,1,2,2) = p_group%p_list(hcube(1,1,2,1))%ngb_vyright_index
!
!                                        ! all combinations of 3 steps
!                                        hcube(1,2,2,2) = p_group%p_list(hcube(1,2,2,1))%ngb_vyright_index
!                                        hcube(2,1,2,2) = p_group%p_list(hcube(2,1,2,1))%ngb_vyright_index
!                                        hcube(2,2,1,2) = p_group%p_list(hcube(2,2,1,1))%ngb_vyright_index
!                                        hcube(2,2,2,1) = p_group%p_list(hcube(2,2,1,1))%ngb_vxright_index
!
!                                        ! 4 steps
!                                        hcube(2,2,2,2) = p_group%p_list(hcube(2,2,2,1))%ngb_vyright_index
!
!                                        ! [[file:~/mcp/maltpic/ltpic-bsl.tex::affine-fn*]] use the values of f0 at these
!                                        ! neighbours to interpolate the value of f0 at
!                                        ! [[file:~/mcp/maltpic/ltpic-bsl.tex::hat-bz*]]. MCP -> oui. Ici si tu utilises
!                                        ! des particules affines (part_deg = 1) la valeur de f0 se déduit de celle du
!                                        ! poids de la particule.  En fait tu peux utiliser une formule semblable à celle
!                                        ! qui est utilisée dans la fonction sll_lt_pic_4d_write_f_on_remap_grid, mais
!                                        ! sans faire intervenir la matrice de déformation à l'intérieur des splines.
!
!                                        ! place the resulting value of f on the virtual particle in
!                                        ! p_group%target_values
!
!                                          ! MCP: [BEGIN-DEBUG] store the (computed) absolute initial position of the virtual particle
!                                        ! [[get_initial_position_on_cartesian_grid_from_particle_index]]
!                                           call get_initial_position_on_cartesian_grid_from_particle_index(kprime, &
!                                                number_parts_x,number_parts_y,number_parts_vx,number_parts_vy, &
!                                                j_x,j_y,j_vx,j_vy)
!                                           x_kprime_t0 =  parts_x_min  + (j_x-1)  * h_parts_x
!                                           y_kprime_t0 =  parts_y_min  + (j_y-1)  * h_parts_y
!                                           vx_kprime_t0 = parts_vx_min + (j_vx-1) * h_parts_vx
!                                           vy_kprime_t0 = parts_vy_min + (j_vy-1) * h_parts_vy
!
!                                           p_group%debug_bsl_remap(i_x,i_y,i_vx,i_vy,1,3) = x_kprime_t0 + x_t0
!                                           p_group%debug_bsl_remap(i_x,i_y,i_vx,i_vy,2,3) = y_kprime_t0 + y_t0
!                                           p_group%debug_bsl_remap(i_x,i_y,i_vx,i_vy,3,3) = vx_kprime_t0 + vx_t0
!                                           p_group%debug_bsl_remap(i_x,i_y,i_vx,i_vy,4,3) = vy_kprime_t0 + vy_t0
!
!                                          ! MCP [END-DEBUG]
!
!
!                                        do side_x = 1,2
!                                           if(side_x == 1)then
!                                              x_aux = x_t0
!                                           else
!                                              x_aux = h_parts_x - x_t0
!                                           end if
!                                           do side_y = 1,2
!                                              if(side_y == 1)then
!                                                 y_aux = y_t0
!                                              else
!                                                 y_aux = h_parts_y - y_t0
!                                              end if
!                                              do side_vx = 1,2
!                                                 if(side_vx == 1)then
!                                                    vx_aux = vx_t0
!                                                 else
!                                                    vx_aux = h_parts_vx - vx_t0
!                                                 end if
!                                                 do side_vy = 1,2
!                                                    if(side_vy == 1)then
!                                                       vy_aux = vy_t0
!                                                    else
!                                                       vy_aux = h_parts_vy - vy_t0
!                                                    end if
!
!                                                    ! uses [[sll_pic_shape]]
!                                                    p_group%target_values(i_x,i_y,i_vx,i_vy) =                    &
!                                                         p_group%target_values(i_x,i_y,i_vx,i_vy)                 &
!                                                         + p_group%p_list(hcube(side_x,side_y,side_vx,side_vy))%q &
!                                                         * sll_pic_shape(part_degree,x_aux,y_aux,vx_aux,vy_aux,   &
!                                                         inv_h_parts_x,inv_h_parts_y,inv_h_parts_vx,inv_h_parts_vy)
!                                                 end do
!                                              end do
!                                           end do
!                                        end do
!                                     end if
!                                  end if
!                               end if
!                            end do
!                         end do
!                      end do
!                   end do
!                end if
!             end do
!          end do
!       end do
!    end do
!  end subroutine sll_lt_pic_4d_write_bsl_f_on_remap_grid

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
    sll_int32 :: i_x_min,  i_x_max,  i_x
    sll_int32 :: i_y_min,  i_y_max,  i_y
    sll_int32 :: i_vx_min, i_vx_max, i_vx
    sll_int32 :: i_vy_min, i_vy_max, i_vy
    sll_int32 :: k, k_ngb
    
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
    LOGICAL :: track_markers_outside_domain

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
        
    track_markers_outside_domain = p_group%track_markers_outside_domain

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
                                    track_markers_outside_domain, &
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

  ! WARNING: this function uses the classical LTP representation of f -- hence NOT the virtual-bsl one --
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
        
    sll_int32 :: num_points_1
    sll_int32 :: num_points_2
    sll_int32 :: i1_min, i1_max, i1
    sll_int32 :: i2_min, i2_max, i2
    sll_int32 :: k, k_ngb

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

!    sll_real64 :: x_k_left,  x_k_right
!    sll_real64 :: y_k_left,  y_k_right
!    sll_real64 :: vx_k_left, vx_k_right
!    sll_real64 :: vy_k_left, vy_k_right

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
    LOGICAL :: track_markers_outside_domain

    sll_real64 :: sum_part_radius_x
    sll_real64 :: sum_part_radius_y
    sll_real64 :: sum_part_radius_vx
    sll_real64 :: sum_part_radius_vy
    sll_int32 :: mx_min_sum
    sll_int32 :: mx_max_sum
    sll_int32 :: my_min_sum
    sll_int32 :: my_max_sum
    sll_int32 :: n_debug

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
    track_markers_outside_domain = p_group%track_markers_outside_domain

    num_points_1  = int(plot_m2d%num_cells1+1, i32)
    num_points_2  = int(plot_m2d%num_cells2+1, i32)

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
                                    track_markers_outside_domain,   &
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
                                                                - plot_x1_min  - part_radius_1      ) ), i32 ) + 2 )
                i1_max  = min( num_points_1, &
                               int(ceiling( inv_h_plot_1  * ( x1_k + real( m1*mesh_period_1, f64)        &
                                                                - plot_x1_min  + part_radius_1      ) ), i32 ) )
                i2_min  = max( 1, &
                               int(floor(   inv_h_plot_2  * ( x2_k + real( m2*mesh_period_2, f64)        &
                                                                - plot_x2_min  - part_radius_2      ) ), i32 ) + 2 )
                i2_max  = min( num_points_2, &
                               int(ceiling( inv_h_plot_2  * ( x2_k + real( m2*mesh_period_2, f64)        &
                                                                - plot_x2_min  + part_radius_2      ) ), i32 ) )

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

    sll_int32                             :: file_id    !< file unit number
    sll_int32                             :: error

    print*, "[plot_2d_slice_remapping_grid] plotting slice in file :", file_name

    call sll_new_file_id(file_id, error)
    open(file_id, file=file_name, iostat=error)


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
    write(file_id,*) "# j_y = ", j_y, " *** y_j = ", y_j, " *** j_vy = ", j_vy, " ***  vy_j = ", vy_j
    x_j = parts_x_min
    do j_x = 1, number_parts_x
        vx_j = parts_vx_min
        do j_vx = 1, number_parts_vx
            write(file_id,*) x_j, vx_j, p_group%target_values(j_x,j_y,j_vx,j_vy)
            vx_j = vx_j + h_parts_vx
        end do
        write(file_id,*)
        x_j = x_j + h_parts_x
    end do
    write(file_id,*)
    close(file_id)

        
  end subroutine plot_2d_slice_remapping_grid



  ! mostly for debug 
  subroutine plot_lt_particles( &
      file_name,        &
      p_group ) 

    character(len=*),                    intent(in)    :: file_name
    type(sll_lt_pic_4d_group), pointer,  intent(in)    :: p_group
    sll_int32                             :: file_id    !< file unit number
    sll_int32                             :: error
    sll_int32 :: k
    sll_real64 :: x_k    
    sll_real64 :: y_k    
    type(sll_lt_pic_4d_particle),      pointer    :: particle


    print*, "[plot_lt_particles] plotting slice in file :", file_name
 
    call sll_new_file_id(file_id, error)
    open(file_id, file=file_name, iostat=error)

    do k = 1, p_group%number_particles
        particle => p_group%p_list(k)

        call cell_offset_to_global_extended(particle%dx, particle%dy, particle%ic_x, particle%ic_y,   &
                                            p_group%mesh, x_k, y_k )

        call periodic_correction(p_group, x_k, y_k)

        write(file_id,*) x_k, "  ", y_k, "  ", particle%vx , "  ", particle%vy,  "  ", particle%q

    end do
    write(file_id,*)
    close(file_id)

  end subroutine plot_lt_particles


!! plot_f_slice_x_vx  plots an approximation of f_x_vx = \int \int f(x,y,v_x,v_y) d y d v_y
!  using a 4d grid:
!   - the nodes in x and v_x are used to plot the values of f_x_vx
!     and
!     the nodes in y and v_y are used to approximate the integrals along y and v_y
!
!  - Grid dimensions: here we give
!       -)  n_virtual_cells_D: the number of "virtual cells" (where the backward flow is linearized) per dimension D
!           ("virtual" corresponds to the terminology used in the function sll_lt_pic_4d_write_f_on_grid_or_deposit)
!       -)  n_virtual_D: the number of virtual nodes (or virtual particles) per virtual cells, in the dimension D.
!    Then the number of plotted ("virtual") nodes is
!       n_virtual_D * n_virtual_cells_D per dimension D

subroutine plot_f_slice_x_vx(p_group,           &
                             plot_grid_x_min,   &
                             plot_grid_x_max,   &
                             plot_grid_y_min,   &
                             plot_grid_y_max,   &
                             plot_grid_vx_min,  &
                             plot_grid_vx_max,  &
                             plot_grid_vy_min,  &
                             plot_grid_vy_max,  &
                             n_virtual_cells_x, n_virtual_cells_y, n_virtual_cells_vx, n_virtual_cells_vy,       &
                             n_virtual_x, n_virtual_y, n_virtual_vx, n_virtual_vy,                   &
                             array_name, iplot)

  type(sll_lt_pic_4d_group), pointer,  intent(inout)  :: p_group
  character(len=*),         intent(in)  :: array_name !< field name
  sll_real64, intent(in)                :: plot_grid_x_min
  sll_real64, intent(in)                :: plot_grid_x_max
  sll_real64, intent(in)                :: plot_grid_y_min
  sll_real64, intent(in)                :: plot_grid_y_max
  sll_real64, intent(in)                :: plot_grid_vx_min
  sll_real64, intent(in)                :: plot_grid_vx_max
  sll_real64, intent(in)                :: plot_grid_vy_min
  sll_real64, intent(in)                :: plot_grid_vy_max
  sll_int32, intent(in)                 :: n_virtual_cells_x
  sll_int32, intent(in)                 :: n_virtual_cells_y
  sll_int32, intent(in)                 :: n_virtual_cells_vx
  sll_int32, intent(in)                 :: n_virtual_cells_vy
  sll_int32, intent(in)                 :: n_virtual_x
  sll_int32, intent(in)                 :: n_virtual_y
  sll_int32, intent(in)                 :: n_virtual_vx
  sll_int32, intent(in)                 :: n_virtual_vy
  sll_int32, intent(in)                 :: iplot      !< plot counter
  character(len=4)                      :: cplot
  sll_int32 :: ierr
  sll_int32 :: num_virtual_parts_x
  sll_int32 :: num_virtual_parts_y
  sll_int32 :: num_virtual_parts_vx
  sll_int32 :: num_virtual_parts_vy

  logical       :: scenario_is_deposition
  logical       :: use_remapping_grid
  sll_real64, dimension(:,:),       pointer :: x_vx_grid_values
  type(sll_cartesian_mesh_4d),      pointer :: plotting_grid_4d
  type(sll_charge_accumulator_2d),  pointer :: dummy_q_accumulator

  nullify(dummy_q_accumulator)

  ! number of points in the grid
  num_virtual_parts_x =  n_virtual_cells_x *  n_virtual_x
  num_virtual_parts_y =  n_virtual_cells_y *  n_virtual_y
  num_virtual_parts_vx = n_virtual_cells_vx * n_virtual_vx
  num_virtual_parts_vy = n_virtual_cells_vy * n_virtual_vy

  SLL_ALLOCATE( x_vx_grid_values(num_virtual_parts_x,num_virtual_parts_vx),ierr)

  plotting_grid_4d => new_cartesian_mesh_4d( num_virtual_parts_x-1,     &
                                             num_virtual_parts_y-1,     &
                                             num_virtual_parts_vx-1,    &
                                             num_virtual_parts_vy-1,    &
                                             plot_grid_x_min,           &
                                             plot_grid_x_max,           &
                                             plot_grid_y_min,           &
                                             plot_grid_y_max,           &
                                             plot_grid_vx_min,          &
                                             plot_grid_vx_max,          &
                                             plot_grid_vy_min,          &
                                             plot_grid_vy_max           &
                                            )

  scenario_is_deposition = .false.
  use_remapping_grid = .false.

  call sll_lt_pic_4d_write_f_on_grid_or_deposit (p_group, dummy_q_accumulator,      &
                                                       scenario_is_deposition,      &
                                                       use_remapping_grid,          &
                                                       plotting_grid_4d,            &
                                                       x_vx_grid_values,            &
                                                       n_virtual_x,                 &
                                                       n_virtual_y,                 &
                                                       n_virtual_vx,                &
                                                       n_virtual_vy)

  call sll_gnuplot_2d(  plot_grid_x_min,  plot_grid_x_max,  num_virtual_parts_x,    &
                        plot_grid_vx_min, plot_grid_vx_max, num_virtual_parts_vx,   &
                        x_vx_grid_values, array_name, iplot, ierr )

end subroutine plot_f_slice_x_vx


subroutine plot_f_slice_1st_try(p_group,        &
                                array_name, iplot)
  type(sll_lt_pic_4d_group), pointer,  intent(inout)  :: p_group
  character(len=*),         intent(in)  :: array_name !< field name
  sll_int32, intent(in)                 :: iplot      !< plot counter
!  sll_int32                             :: npoints
!  sll_int32                             :: ipoints
  character(len=4)                      :: cplot
  sll_int32                             :: n_virtual_x_for_remapping
  sll_int32                             :: n_virtual_y_for_remapping
  sll_int32                             :: n_virtual_vx_for_remapping
  sll_int32                             :: n_virtual_vy_for_remapping

  n_virtual_x_for_remapping = 2
  n_virtual_y_for_remapping = 2
  n_virtual_vx_for_remapping = 2
  n_virtual_vy_for_remapping = 2

  call sll_lt_pic_4d_write_f_on_remapping_grid( p_group,                        &
                                                n_virtual_x_for_remapping,      &
                                                n_virtual_y_for_remapping,      &
                                                n_virtual_vx_for_remapping,     &
                                                n_virtual_vy_for_remapping)

  call int2string(iplot, cplot)

  call plot_2d_slice_remapping_grid( array_name//"_"//cplot//".dat",        &
                                     p_group )

end subroutine plot_f_slice_1st_try


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
    sll_int32, intent(in) :: k
    sll_int32, intent(in) :: n_parts_x
    sll_int32, intent(in) :: n_parts_y
    sll_int32, intent(in) :: n_parts_vx
    sll_int32, intent(in) :: n_parts_vy
    sll_int32, intent(out) :: j_x
    sll_int32, intent(out) :: j_y
    sll_int32, intent(out) :: j_vx
    sll_int32, intent(out) :: j_vy
    sll_int32              :: k_aux

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
    sll_int32, intent(in) :: j_x
    sll_int32, intent(in) :: j_y
    sll_int32, intent(in) :: j_vx
    sll_int32, intent(in) :: j_vy
    sll_int32, intent(in) :: n_parts_x
    sll_int32, intent(in) :: n_parts_y
    sll_int32, intent(in) :: n_parts_vx
    sll_int32, intent(in) :: n_parts_vy
    sll_int32, intent(out) :: k

    if(         j_x <= 0  .or. j_x > n_parts_x      &
        .or.    j_y <= 0  .or. j_y > n_parts_y      &
        .or.    j_vx <= 0 .or. j_vx > n_parts_vx    &
        .or.    j_vy <= 0 .or. j_vy > n_parts_vy  )then
        k = 0
    else
        k = 1+ (j_vy-1) + (j_vx-1) * n_parts_vy + (j_y-1) * n_parts_vy * n_parts_vx + (j_x-1) * n_parts_vy * n_parts_vx * n_parts_y
    end if

    SLL_ASSERT( k <= n_parts_x * n_parts_y * n_parts_vx * n_parts_vy )

end subroutine


  ! <<sll_lt_pic_4d_deposit_charge_on_2d_mesh>>
  subroutine sll_lt_pic_4d_deposit_charge_on_2d_mesh( p_group, q_accumulator,           &
                                                      n_virtual_x_for_deposition,       &
                                                      n_virtual_y_for_deposition,       &
                                                      n_virtual_vx_for_deposition,      &
                                                      n_virtual_vy_for_deposition,      &
                                                      given_total_density)

    type(sll_lt_pic_4d_group),pointer,intent(inout) :: p_group
    type(sll_charge_accumulator_2d), pointer, intent(inout) :: q_accumulator
    sll_int32, intent(in) :: n_virtual_x_for_deposition
    sll_int32, intent(in) :: n_virtual_y_for_deposition
    sll_int32, intent(in) :: n_virtual_vx_for_deposition
    sll_int32, intent(in) :: n_virtual_vy_for_deposition
    type(sll_cartesian_mesh_4d),    pointer :: dummy_grid_4d
    sll_real64, dimension(:,:),     pointer :: dummy_array_2d

    sll_real64, intent(in), optional :: given_total_density
    logical       :: scenario_is_deposition
    logical       :: use_remapping_grid

    scenario_is_deposition = .true.
    use_remapping_grid = .false.
    nullify(dummy_grid_4d)
    nullify(dummy_array_2d)

    if( present(given_total_density) )then
        call sll_lt_pic_4d_write_f_on_grid_or_deposit (p_group, q_accumulator,          &
                                                       scenario_is_deposition,          &
                                                       use_remapping_grid,              &
                                                       dummy_grid_4d,                   &
                                                       dummy_array_2d,                  &
                                                       n_virtual_x_for_deposition,      &
                                                       n_virtual_y_for_deposition,      &
                                                       n_virtual_vx_for_deposition,     &
                                                       n_virtual_vy_for_deposition,     &
                                                       given_total_density)
    else
        call sll_lt_pic_4d_write_f_on_grid_or_deposit (p_group, q_accumulator,          &
                                                       scenario_is_deposition,          &
                                                       use_remapping_grid,              &
                                                       dummy_grid_4d,                   &
                                                       dummy_array_2d,                  &
                                                       n_virtual_x_for_deposition,      &
                                                       n_virtual_y_for_deposition,      &
                                                       n_virtual_vx_for_deposition,     &
                                                       n_virtual_vy_for_deposition)
    end if

  end subroutine sll_lt_pic_4d_deposit_charge_on_2d_mesh


  subroutine sll_lt_pic_4d_write_f_on_remapping_grid( p_group,                         &
                                                      n_virtual_x_for_remapping,       &
                                                      n_virtual_y_for_remapping,       &
                                                      n_virtual_vx_for_remapping,      &
                                                      n_virtual_vy_for_remapping )

    type(sll_lt_pic_4d_group),pointer,intent(inout) :: p_group
    sll_int32, intent(in) :: n_virtual_x_for_remapping
    sll_int32, intent(in) :: n_virtual_y_for_remapping
    sll_int32, intent(in) :: n_virtual_vx_for_remapping
    sll_int32, intent(in) :: n_virtual_vy_for_remapping
    type(sll_charge_accumulator_2d),pointer :: dummy_q_accumulator
    type(sll_cartesian_mesh_4d),    pointer :: dummy_grid_4d
    sll_real64, dimension(:,:),     pointer :: dummy_array_2d
    logical       :: scenario_is_deposition
    logical       :: use_remapping_grid

    nullify(dummy_q_accumulator)
    nullify(dummy_grid_4d)
    nullify(dummy_array_2d)

    scenario_is_deposition = .false.
    use_remapping_grid = .true.

    call sll_lt_pic_4d_write_f_on_grid_or_deposit (p_group, dummy_q_accumulator,    &
                                                   scenario_is_deposition,          &
                                                   use_remapping_grid,              &
                                                   dummy_grid_4d,                   &
                                                   dummy_array_2d,                  &
                                                   n_virtual_x_for_remapping,       &
                                                   n_virtual_y_for_remapping,       &
                                                   n_virtual_vx_for_remapping,      &
                                                   n_virtual_vy_for_remapping)

  end subroutine sll_lt_pic_4d_write_f_on_remapping_grid


  ! <<sll_lt_pic_4d_write_f_on_grid_or_deposit>> <<ALH>> has two scenarios:
  !  - 1.  the "write f" scenario:
  !        write the density on the (phase-space) remapping grid, using the method described
  !        in the "BSL-remapping" notes (version of december 2, 2014) cf
  !        [[file:~/mcp/maltpic/ltpic-bsl.tex::BSL_remapping]] and more precisely
  !        [[file:~/mcp/maltpic/ltpic-bsl.tex::BSL_remapping_step_1]].  Algorithm from
  !        [[file:~/mcp/maltpic/ltpic-bsl.tex::algo:pic-vr]] (but without the deposition step)
  !
  !        -- this function should be a faster alternative to [[sll_lt_pic_4d_write_f_on_remap_grid]] --
  !
  !        Note: the (x,y)-projection of the remapping grid may be larger than the "Poisson" 2d mesh associated with the
  !        particle group (in particular if the (x,y) domain is not periodic)
  !
  !  - 2.  the "deposition" scenario:
  !        deposit the charge on the Poisson cells given in the given charge_accumulator_2d object
  !
  !  In every case, this routine computes (approximated) values of the density f(t_n) on 'virtual' nodes (sometimes called
  !  'virtual particles') located on a cartesian grid of the 4d phase-space.
  !  For different reasons, these virtual nodes are gathered in 'virtual' cells, and the given arguments
  !  n_virtual_x, n_virtual_y, n_virtual_vx, n_virtual_vy is the number of virtual nodes per virtual cell, in every dimension
  !
  !  -> in the "write f" scenario, these virtual nodes are in fact the nodes of a given grid -- either the remapping one,
  !     or some given grid -- hence they are not really virtual.
  !     On the contrary the virtual cells are really virtual, in the sense that they are just a way to gather the computations
  !
  !  -> in the "deposition" scenario, the virtual nodes are really virtual: they correspond to temporary particles which
  !     are deposited with a standard PIC procedure. And the virtual cells are "half-virtual" in the sense that their (x,y)
  !     projection coincides with the cells of the Poisson mesh, whereas in the velocity dimensions they are created to gather
  !     the computations, just as in the "remapping" scenario
  !
  !     In particular, taking a larger value for n_virtual has the following effect:
  !     -> in the  "write f" scenario, larger virtual cells will be used to compute the approximated values of f(t_n) on the
  !        nodes of the (remapping or given) grid. This will speed-up the code and is morally ok if the characteristic flow
  !        is smooth
  !     -> in the "deposition" scenario, finer grids of virtual point particles will be (temporarily) created and
  !        deposited. This will slow down the code and is morally required if the density f(t_n) is not locally smooth
  !
  !  - Note: in the "write f" scenario, the grid is known through:
  !        the number of grid points (not cells) in every dimension
  !        and the max and min coordinates of these points, in every dimension
  !
  !  - given_total_density is an optional argument that is given to make the deposition method conservative
  !    (note: it may be used also in the 'write_f' scenario, but one has to define what conservative means in this case)
  !
  !  Note: This routine is an evolution from sll_lt_pic_4d_write_bsl_f_on_remap_grid (which will be eventually discarded)

  ! todo: Treat the non-periodic case. In this case we can place the virtual particles slightly off the boundaries of the
  ! todo: virtual cells, so that we do not need a special treatment for the particles on the right (x and y) domain boundaries

  subroutine sll_lt_pic_4d_write_f_on_grid_or_deposit (p_group, q_accumulator,      &
                                                       scenario_is_deposition,      &
                                                       use_remapping_grid,          &
                                                       given_grid_4d,               &
                                                       given_array_2d,              &
                                                       n_virtual_x,                 &
                                                       n_virtual_y,                 &
                                                       n_virtual_vx,                &
                                                       n_virtual_vy,                &
                                                       given_total_density)

    ! [[file:../pic_particle_types/lt_pic_4d_group.F90::sll_lt_pic_4d_group]] p_group contains both the existing
    ! particles and the virtual remapping grid
    type(sll_lt_pic_4d_group),pointer,intent(inout) :: p_group
    type(sll_charge_accumulator_2d), pointer, intent(inout) :: q_accumulator
    logical, intent(in) :: scenario_is_deposition            ! if false, then scenario is "write on grid"
    logical, intent(in) :: use_remapping_grid                ! if false, then grid must be given
    type(sll_cartesian_mesh_4d), pointer, intent(in)        :: given_grid_4d
    sll_real64, dimension(:,:),     pointer, intent(inout)  :: given_array_2d   ! assumed in x, vx for now
    ! <<n_virtual>>      ! see comments above for the meaning
    sll_int32, intent(in) :: n_virtual_x
    sll_int32, intent(in) :: n_virtual_y
    sll_int32, intent(in) :: n_virtual_vx
    sll_int32, intent(in) :: n_virtual_vy

    sll_real64, intent(in), optional :: given_total_density

    type(charge_accumulator_cell_2d), pointer :: charge_accumulator_cell

    sll_real64 :: deposited_density
    sll_real64 :: charge_correction_factor

    ! cf [[file:~/mcp/maltpic/ltpic-bsl.tex::N*]]

    sll_int32 :: num_virtual_cells_x
    sll_int32 :: num_virtual_cells_y
    sll_int32 :: num_virtual_cells_vx
    sll_int32 :: num_virtual_cells_vy

    sll_int32 :: number_virtual_particles_x
    sll_int32 :: number_virtual_particles_y
    sll_int32 :: number_virtual_particles_vx
    sll_int32 :: number_virtual_particles_vy

    ! [[file:~/mcp/maltpic/ltpic-bsl.tex::h_parts_x]] and h_parts_y, h_parts_vx, h_parts_vy

    sll_real64 :: h_parts_x
    sll_real64 :: h_parts_y
    sll_real64 :: h_parts_vx
    sll_real64 :: h_parts_vy

    sll_real64 :: inv_h_parts_x
    sll_real64 :: inv_h_parts_y
    sll_real64 :: inv_h_parts_vx
    sll_real64 :: inv_h_parts_vy

    sll_real64 :: h_virtual_parts_x
    sll_real64 :: h_virtual_parts_y
    sll_real64 :: h_virtual_parts_vx
    sll_real64 :: h_virtual_parts_vy

    sll_real64 :: inv_h_virtual_parts_x
    sll_real64 :: inv_h_virtual_parts_y
    sll_real64 :: inv_h_virtual_parts_vx
    sll_real64 :: inv_h_virtual_parts_vy

    sll_real64 :: phase_space_virtual_dvol

    sll_real64 :: parts_x_min
    sll_real64 :: parts_y_min
    sll_real64 :: parts_vx_min
    sll_real64 :: parts_vy_min

    sll_real64 :: virtual_parts_x_min
    sll_real64 :: virtual_parts_y_min
    sll_real64 :: virtual_parts_vx_min
    sll_real64 :: virtual_parts_vy_min

    sll_real64 :: virtual_cells_x_min
    sll_real64 :: virtual_cells_y_min
    sll_real64 :: virtual_cells_vx_min
    sll_real64 :: virtual_cells_vy_min

    sll_real64 :: virtual_grid_x_min  ! do we need this? should be = virtual_parts_x_min...
    sll_real64 :: virtual_grid_x_max
    sll_real64 :: virtual_grid_y_min
    sll_real64 :: virtual_grid_y_max
    sll_real64 :: virtual_grid_vx_min
    sll_real64 :: virtual_grid_vx_max
    sll_real64 :: virtual_grid_vy_min
    sll_real64 :: virtual_grid_vy_max

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

    !    sll_real32 :: dx
    !    sll_real32 :: dy
    !    sll_real32 :: dvx
    !    sll_real32 :: dvy

    sll_real64 :: closest_particle_distance_to_first_corner
    sll_real64 :: particle_distance_to_first_corner

    ! working space

    sll_real64 :: tmp, tmp1, tmp2
    sll_real32 :: tmp_dx, tmp_dy


    ! index of particle closest to the center of each virtual cell. Array dimensions defined by the contents of
    ! [[file:../pic_particle_types/lt_pic_4d_group.F90::sll_lt_pic_4d_group-remapping_grid]]. If [[n_virtual]] is
    ! greater than 1, the size of this array is smaller than the number of real remapping_grid cells.

    sll_int32,dimension(:,:,:,:),allocatable :: closest_particle
    sll_real64,dimension(:,:,:,:),allocatable :: closest_particle_distance

    sll_int32 :: i ! x dimension
    sll_int32 :: j ! y dimension
    sll_int32 :: k,kprime ! particle index
    sll_int32 :: neighbour ! particle index for local use
    sll_int32 :: l ! vx dimension
    sll_int32 :: m ! vy dimension

    sll_int32 :: k_neighbor
    sll_int32 :: k_particle_closest_to_first_corner

    ! indices in a virtual cell (go from 1 to [[n_virtual]])

    sll_int :: ivirt ! x dimension
    sll_int :: jvirt ! y dimension
    sll_int :: lvirt ! vx dimension
    sll_int :: mvirt ! vy dimension

    sll_int :: i_x,i_y,i_vx,i_vy

    ! particle pointer [[file:../pic_particle_types/lt_pic_4d_particle.F90::sll_lt_pic_4d_particle]]
    type(sll_lt_pic_4d_particle),pointer :: particle

    ! <<g>> cartesian grid pointer to
    ! [[file:../pic_particle_types/lt_pic_4d_group.F90::sll_lt_pic_4d_group-remapping_grid]]

    type(sll_cartesian_mesh_4d),pointer :: g

    ! periodicity

    LOGICAL :: domain_is_x_periodic
    LOGICAL :: domain_is_y_periodic
    LOGICAL :: track_markers_outside_domain

    LOGICAL :: find_k_prime_step_by_step

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

    sll_real64 :: x_to_xk, y_to_yk, vx_to_vxk, vy_to_vyk

    sll_real64 :: x_t0_to_xkprime_t0
    sll_real64 :: y_t0_to_ykprime_t0
    sll_real64 :: vx_t0_to_vxkprime_t0
    sll_real64 :: vy_t0_to_vykprime_t0

    sll_real64 :: d1_x, d1_y, d1_vx, d1_vy
    sll_real64 :: d2_x, d2_y, d2_vx, d2_vy
    sll_real64 :: d3_x, d3_y, d3_vx, d3_vy
    sll_real64 :: d4_x, d4_y, d4_vx, d4_vy

    sll_real64 :: part_radius_x
    sll_real64 :: part_radius_y
    sll_real64 :: part_radius_vx
    sll_real64 :: part_radius_vy

    sll_real64 :: dx_in_virtual_cell
    sll_real64 :: dy_in_virtual_cell

    sll_real64 :: x_center_virtual_cell
    sll_real64 :: y_center_virtual_cell

    sll_real64 :: f_value_on_virtual_particle
    sll_real64 :: virtual_charge

    ! coordinates of a virtual particle at time 0 relative to the coordinates of one real particle

    sll_real64 :: x_t0,y_t0,vx_t0,vy_t0

    sll_real64 :: x_t0_to_xk_t0
    sll_real64 :: y_t0_to_yk_t0
    sll_real64 :: vx_t0_to_vxk_t0
    sll_real64 :: vy_t0_to_vyk_t0

    sll_int32 :: part_degree

    sll_int32 :: ierr
    sll_int32 :: i_cell

    ! temporary workspace
    sll_real64 :: x_aux
    sll_real64 :: y_aux
    sll_real64 :: vx_aux
    sll_real64 :: vy_aux

    sll_real64 :: length

    ! value 1 or 2 points to each side of an hypercube in direction x,y,vx or vy
    sll_int :: side_x,side_y,side_vx,side_vy
    sll_int32,dimension(2,2,2,2) :: hcube

    sll_int32 :: j_x,j_y,j_vx,j_vy

    sll_real64 :: one_over_two_pi
    sll_real64 :: one_over_thermal_velocity_squared
    sll_real64 :: alpha_landau
    sll_real64 :: k_landau

    ! pw-affine approximations of exp and cos for fast (?) evaluations
    sll_int32 :: i_table, ncells_table
    sll_real64, dimension(:),allocatable :: cos_table, exp_table
    sll_real64 :: s, s_aux
    sll_real64 :: hs_cos_table, hs_exp_table
    sll_real64 :: ds_table
    sll_real64 :: smin_cos_table, smax_cos_table
    sll_real64 :: smin_exp_table, smax_exp_table
    sll_real64 :: si_cos, si_exp

    sll_real64 :: cos_approx, exp_approx

    ! --- end of declarations

    ! -- creating g the virtual grid [begin] --

    if( scenario_is_deposition )then


        ! todo: modify the code with offset virtual particles to deposit the charge on the Poisson cells

        if( p_group%domain_is_x_periodic )then
            num_virtual_cells_x = p_group%mesh%num_cells1
            virtual_grid_x_min = p_group%mesh%eta1_min
            virtual_grid_x_max = p_group%mesh%eta1_max
        else
            print *, "error (87585758769753486576676543): change code here, place the virtual nodes inside virtual (Poisson) cells"
            print *, "error (87585758769753486576676543): so that the virtual cells can be just the Poisson cells -- "
            stop


            ! an extra cell is needed outside (in every direction) so that the approximation of f(t_n) by regular
            ! splines located at the virtual nodes is accurate close to the domain boundaries
            num_virtual_cells_x = p_group%mesh%num_cells1 + 2
            virtual_grid_x_min = p_group%mesh%eta1_min - p_group%mesh%delta_eta1
            virtual_grid_x_max = p_group%mesh%eta1_max + p_group%mesh%delta_eta1
        end if

        if( p_group%domain_is_y_periodic )then
            num_virtual_cells_y = p_group%mesh%num_cells2
            virtual_grid_y_min = p_group%mesh%eta2_min
            virtual_grid_y_max = p_group%mesh%eta2_max
        else
            ! same reason than for num_virtual_cells_x
            num_virtual_cells_y = p_group%mesh%num_cells2 + 2
            virtual_grid_y_min = p_group%mesh%eta2_min - p_group%mesh%delta_eta2
            virtual_grid_y_max = p_group%mesh%eta2_max + p_group%mesh%delta_eta2
        end if

        ! Because the Poisson mesh does not prescribe any resolution in velocity
        ! the resolution of the 'virtual' cells in the velocity dimensions is inferred from the remapping (or initial) grid
        num_virtual_cells_vx = p_group%number_parts_vx
        virtual_grid_vx_min = p_group%remapping_grid%eta3_min
        virtual_grid_vx_max = p_group%remapping_grid%eta3_max

        num_virtual_cells_vy = p_group%number_parts_vy
        virtual_grid_vy_min = p_group%remapping_grid%eta4_min
        virtual_grid_vy_max = p_group%remapping_grid%eta4_max

        number_virtual_particles_x =  n_virtual_x  * num_virtual_cells_x
        number_virtual_particles_y =  n_virtual_y  * num_virtual_cells_y
        number_virtual_particles_vx = n_virtual_vx * num_virtual_cells_vx
        number_virtual_particles_vy = n_virtual_vy * num_virtual_cells_vy

        g => new_cartesian_mesh_4d( number_virtual_particles_x,        &
                                    number_virtual_particles_y,        &
                                    number_virtual_particles_vx,       &
                                    number_virtual_particles_vy,       &
                                    virtual_grid_x_min,   &
                                    virtual_grid_x_max,   &
                                    virtual_grid_y_min,   &
                                    virtual_grid_y_max,   &
                                    virtual_grid_vx_min,  &
                                    virtual_grid_vx_max,  &
                                    virtual_grid_vy_min,  &
                                    virtual_grid_vy_max   &
                                   )

        if( present(given_total_density) )then
            deposited_density = 0
        end if

    else

        if( use_remapping_grid )then

            g => p_group%remapping_grid

            number_virtual_particles_x = p_group%number_parts_x
            number_virtual_particles_y = p_group%number_parts_y
            number_virtual_particles_vx = p_group%number_parts_vx
            number_virtual_particles_vy = p_group%number_parts_vy

            num_virtual_cells_x =  int(ceiling(number_virtual_particles_x * 1. / n_virtual_x) )
            num_virtual_cells_y =  int(ceiling(number_virtual_particles_y * 1. / n_virtual_y) )
            num_virtual_cells_vx = int(ceiling(number_virtual_particles_vx * 1. / n_virtual_vx))
            num_virtual_cells_vy = int(ceiling(number_virtual_particles_vy * 1. / n_virtual_vy))


            ! initialize [[file:../pic_particle_types/lt_pic_4d_group.F90::target_values]]
            p_group%target_values(:,:,:,:) = 0

            !print *, "6453 before remap -> DEBUG: ", p_group%number_parts_x/2,    &
            !                               p_group%number_parts_y/2,    &
            !                               p_group%number_parts_vx/2,   &
            !                               p_group%number_parts_vy/2
            !print *, "6454 before remap -> DEBUG: ", p_group%target_values(p_group%number_parts_x/2,  p_group%number_parts_y/2, &
            !                                                     p_group%number_parts_vx/2, p_group%number_parts_vy/2)



        else

            ! then use the given 4d grid and write values in given (x, vx for now) array given_array_2d
            g => given_grid_4d

            number_virtual_particles_x = given_grid_4d%num_cells1 + 1
            number_virtual_particles_y = given_grid_4d%num_cells2 + 1
            number_virtual_particles_vx = given_grid_4d%num_cells3 + 1
            number_virtual_particles_vy = given_grid_4d%num_cells4 + 1

            SLL_ASSERT( mod(number_virtual_particles_x,  n_virtual_x)  == 0 )
            SLL_ASSERT( mod(number_virtual_particles_y,  n_virtual_y)  == 0 )
            SLL_ASSERT( mod(number_virtual_particles_vx, n_virtual_vx) == 0 )
            SLL_ASSERT( mod(number_virtual_particles_vy, n_virtual_vy) == 0 )

            num_virtual_cells_x = number_virtual_particles_x / n_virtual_x
            num_virtual_cells_y = number_virtual_particles_y / n_virtual_y
            num_virtual_cells_vx = number_virtual_particles_vx / n_virtual_vx
            num_virtual_cells_vy = number_virtual_particles_vy / n_virtual_vy

            ! for now we assume that given_array_2d is in (x, vx) space
            SLL_ASSERT(size(given_array_2d,1) == number_virtual_particles_x)
            SLL_ASSERT(size(given_array_2d,2) == number_virtual_particles_vx)
            given_array_2d(:,:) = 0

        end if

    end if

    ! -- creating g the virtual grid [end] --

    part_degree = p_group%spline_degree

    track_markers_outside_domain = p_group%track_markers_outside_domain

    ! Preparatory work: find out the particle which is closest to each cell center by looping over all particles and
    ! noting which virtual cell contains it. The leftmost virtual cell in each dimension may not be complete.

    SLL_ALLOCATE(closest_particle(num_virtual_cells_x,num_virtual_cells_y,num_virtual_cells_vx,num_virtual_cells_vy),ierr)
    closest_particle(:,:,:,:) = 0

    SLL_ALLOCATE(closest_particle_distance(num_virtual_cells_x,num_virtual_cells_y,num_virtual_cells_vx,num_virtual_cells_vy),ierr)
    closest_particle_distance(:,:,:,:) = 0

    ! remapping grid cell size - same as in [[write_f_on_remap_grid-h_parts_x]]

    h_parts_x    = p_group%remapping_grid%delta_eta1
    h_parts_y    = p_group%remapping_grid%delta_eta2
    h_parts_vx   = p_group%remapping_grid%delta_eta3
    h_parts_vy   = p_group%remapping_grid%delta_eta4

    inv_h_parts_x  = 1./h_parts_x
    inv_h_parts_y  = 1./h_parts_y
    inv_h_parts_vx = 1./h_parts_vx
    inv_h_parts_vy = 1./h_parts_vy

    h_virtual_parts_x    = g%delta_eta1
    h_virtual_parts_y    = g%delta_eta2
    h_virtual_parts_vx   = g%delta_eta3
    h_virtual_parts_vy   = g%delta_eta4

    inv_h_virtual_parts_x  = 1./h_virtual_parts_x
    inv_h_virtual_parts_y  = 1./h_virtual_parts_y
    inv_h_virtual_parts_vx = 1./h_virtual_parts_vx
    inv_h_virtual_parts_vy = 1./h_virtual_parts_vy

    phase_space_virtual_dvol = h_virtual_parts_x * h_virtual_parts_y * h_virtual_parts_vx * h_virtual_parts_vy

    parts_x_min    = p_group%remapping_grid%eta1_min
    parts_y_min    = p_group%remapping_grid%eta2_min
    parts_vx_min   = p_group%remapping_grid%eta3_min
    parts_vy_min   = p_group%remapping_grid%eta4_min

    virtual_parts_x_min    = g%eta1_min
    virtual_parts_y_min    = g%eta2_min
    virtual_parts_vx_min   = g%eta3_min
    virtual_parts_vy_min   = g%eta4_min

    ! offset: the first virtual particles are not on the left boundary of the first virtual cells
    virtual_cells_x_min    = g%eta1_min - 0.5 * g%delta_eta1
    virtual_cells_y_min    = g%eta2_min - 0.5 * g%delta_eta2
    virtual_cells_vx_min   = g%eta3_min - 0.5 * g%delta_eta3
    virtual_cells_vy_min   = g%eta4_min - 0.5 * g%delta_eta4

    ! virtual cell size
    h_virtual_cell_x  = n_virtual_x * h_virtual_parts_x
    h_virtual_cell_y  = n_virtual_y * h_virtual_parts_y
    h_virtual_cell_vx = n_virtual_vx * h_virtual_parts_vx
    h_virtual_cell_vy = n_virtual_vy * h_virtual_parts_vy

    if( scenario_is_deposition )then
        SLL_ASSERT( h_virtual_cell_x == p_group%mesh%delta_eta1)
        SLL_ASSERT( h_virtual_cell_y == p_group%mesh%delta_eta2)
    end if

    !! this is for the Landau damping case:
    one_over_two_pi = 1./(2*sll_pi)
    one_over_thermal_velocity_squared = (1./p_group%thermal_speed)**2
    alpha_landau = p_group%alpha_landau
    k_landau = p_group%k_landau


    ! preparatory loop to fill the [[closest_particle]] array containing the particle closest to the center of each
    ! virtual cell

    closest_particle_distance_to_first_corner = 1e30
    k_particle_closest_to_first_corner = 0

    do k=1, p_group%number_particles ! [[file:../pic_particle_types/lt_pic_4d_group.F90::number_particles]]
       particle => p_group%p_list(k)

       ! find absolute (x,y,vx,vy) coordinates for this particle. Uses
       ! [[file:sll_representation_conversion.F90::cell_offset_to_global]]

       call cell_offset_to_global_extended(particle%dx,particle%dy,particle%ic_x,particle%ic_y,p_group%mesh,x,y)
       call periodic_correction(p_group, x,y)
       vx = particle%vx
       vy = particle%vy

       ! which _virtual_ cell is this particle in? uses
       ! [[file:sll_representation_conversion.F90::compute_cell_and_offset]] and [[g]]

       x_aux = x - virtual_cells_x_min
       i = int( x_aux / h_virtual_cell_x ) + 1

       y_aux = y - virtual_cells_y_min
       j = int( y_aux / h_virtual_cell_y ) + 1

       vx_aux = vx - virtual_cells_vx_min
       l = int( vx_aux / h_virtual_cell_vx ) + 1

       vy_aux = vy - virtual_cells_vy_min
       m = int( vy_aux / h_virtual_cell_vy ) + 1

       ! discard particles in virtual cells off-bounds
       if(  i >= 1 .and. i <= num_virtual_cells_x .and. &
            j >= 1 .and. j <= num_virtual_cells_y .and. &
            l >= 1 .and. l <= num_virtual_cells_vx .and. &
            m >= 1 .and. m <= num_virtual_cells_vy  )then

          call update_closest_particle_arrays(k,                         &
                                              x_aux, y_aux, vx_aux, vy_aux,   &
                                              i, j, l, m,                     &
                                              h_virtual_cell_x, h_virtual_cell_y, h_virtual_cell_vx, h_virtual_cell_vy,   &
                                              closest_particle,               &
                                              closest_particle_distance)

       end if

       particle_distance_to_first_corner = abs(x_aux) + abs(y_aux) + abs(vx_aux) + abs(vy_aux)      !  (why not L1 after all)
       if( particle_distance_to_first_corner < closest_particle_distance_to_first_corner )then
            closest_particle_distance_to_first_corner = particle_distance_to_first_corner
            k_particle_closest_to_first_corner = k
       end if
    end do

    closest_particle(1,1,1,1) = k_particle_closest_to_first_corner

    ! Periodicity treatments copied from [[sll_lt_pic_4d_write_f_on_remap_grid-periodicity]]

    domain_is_x_periodic = .true.   ! temp --- todo: treat the general case
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

    ! MCP: [DEBUG] store the (computed) absolute initial position of the virtual particle
    !p_group%debug_bsl_remap = -100

    if( p_group%use_exact_f0 )then

        ! create two tables for fast evaluation of cosine and exponential
        ! -- THIS IS USEFUL FOR THE LANDAU f0 AND MAYBE NOT FOR OTHER INITIAL DENSITIES --

        ncells_table = 1000
        SLL_ALLOCATE(cos_table(0:ncells_table),ierr)
        SLL_ALLOCATE(exp_table(0:ncells_table),ierr)

        !  cos table range
        smin_cos_table = 0.
        smax_cos_table = 2 * sll_pi
        hs_cos_table = (smax_cos_table - smin_cos_table) / ncells_table

        !  exp table range
        smin_exp_table = -50.
        smax_exp_table =  50.
        hs_exp_table = (smax_exp_table - smin_exp_table) / ncells_table

        si_cos = smin_cos_table
        si_exp = smin_exp_table
        do i_table = 0, ncells_table
            cos_table(i_table) = cos(si_cos)
            exp_table(i_table) = exp(si_exp)
            si_cos = si_cos + hs_cos_table
            si_exp = si_exp + hs_exp_table
        end do

    end if

    ! <<loop_on_virtual_cells>> [[file:~/mcp/maltpic/ltpic-bsl.tex::algo:pic-vr:loop_over_all_cells]]
    ! Loop over all cells of indices i,j,l,m which contain at least one particle

    do i = 1, num_virtual_cells_x
       do j = 1, num_virtual_cells_y

          if( scenario_is_deposition )then

              ! index of the Poisson cell from i and j (see global_to_cell_offset)
              i_cell = i + (j-1) * p_group%mesh%num_cells1

              charge_accumulator_cell => q_accumulator%q_acc(i_cell)

          end if

          do l = 1,num_virtual_cells_vx
             do m = 1,num_virtual_cells_vy

                ! [[file:~/mcp/maltpic/ltpic-bsl.tex::algo:pic-vr:create_virtual_particles]] Create a temporary set of
                ! virtual particles inside the cell.
                ! Note: as written above in the remapping scenario the virtual particles coincide with the existing
                ! remapping_grid [[file:../pic_particle_types/lt_pic_4d_group.F90::sll_lt_pic_4d_group-remapping_grid]]
                ! defined in p_group [[file:../pic_particle_types/lt_pic_4d_group.F90::sll_lt_pic_4d_group]].
                ! In the deposition scenario the virtual particles are used to deposit the charge and they are not stored.
                ! So nothing more to do.

                ! [[file:~/mcp/maltpic/ltpic-bsl.tex::algo:pic-vr:find_closest_real_particle]] Find the real particle
                ! which is closest to the cell center.  Note: speed-wise, it may be necessary to find a way not to scan
                ! all the particles for every cell.  We avoid scanning all the particles for each cell by using the
                ! precomputed array [[closest_particle]]. Virtual cells which do not contain any particle are skipped.

                k = closest_particle(i,j,l,m)

        ! before, we used           call cell_offset_to_global(particle%dx,particle%dy,particle%ic,p_group%mesh,x,y) ;\

#define UPDATE_CLOSEST_PARTICLE_ARRAYS_USING_NEIGHBOR_CELLS(di,dj,dl,dm)                                                \
    do;                                                                                                                 \
        k_neighbor = closest_particle(i+(di), j+(dj), l+(dl), m+(dm));                                                  \
;                                                                                                                       \
        if(k_neighbor /= 0) then;  do          ;                                                                        \
            particle => p_group%p_list(k_neighbor);                                                                     \
            call cell_offset_to_global_extended(particle%dx,particle%dy,particle%ic_x,particle%ic_y,p_group%mesh,x,y) ; \
            call periodic_correction(p_group,x,y) ;                                                                     \
            vx = particle%vx;                                                                                           \
            vy = particle%vy;                                                                                           \
            x_aux = x - g%eta1_min;                                                                                     \
            y_aux = y - g%eta2_min;                                                                                     \
            vx_aux = vx - g%eta3_min;                                                                                   \
            vy_aux = vy - g%eta4_min;                                                                                   \
            call update_closest_particle_arrays(k_neighbor,                                                             \
                                                x_aux, y_aux, vx_aux, vy_aux,                                           \
                                                i, j, l, m,                                                             \
                                                h_virtual_cell_x, h_virtual_cell_y,                                     \
                                                h_virtual_cell_vx, h_virtual_cell_vy,                                   \
                                                closest_particle,                                                       \
                                                closest_particle_distance);                                             \
        exit;                                                                                                           \
        end do;                                                                                                         \
        end if;                                                                                                         \
    exit;                                                                                                               \
    end do



                if(k == 0) then

                    if( i > 1 )then
                        UPDATE_CLOSEST_PARTICLE_ARRAYS_USING_NEIGHBOR_CELLS(-1,0,0,0)
                    end if
                    if( i < num_virtual_cells_x )then
                        UPDATE_CLOSEST_PARTICLE_ARRAYS_USING_NEIGHBOR_CELLS( 1,0,0,0)
                    end if

                    if( j > 1 )then
                        UPDATE_CLOSEST_PARTICLE_ARRAYS_USING_NEIGHBOR_CELLS(0,-1,0,0)
                    end if
                    if( j < num_virtual_cells_y )then
                        UPDATE_CLOSEST_PARTICLE_ARRAYS_USING_NEIGHBOR_CELLS(0, 1,0,0)
                    end if

                    if( l > 1 )then
                        UPDATE_CLOSEST_PARTICLE_ARRAYS_USING_NEIGHBOR_CELLS(0,0,-1,0)
                    end if
                    if( l < num_virtual_cells_vx )then
                        UPDATE_CLOSEST_PARTICLE_ARRAYS_USING_NEIGHBOR_CELLS(0,0, 1,0)
                    end if

                    if( m > 1 )then
                        UPDATE_CLOSEST_PARTICLE_ARRAYS_USING_NEIGHBOR_CELLS(0,0,0,-1)
                    end if
                    if( m < num_virtual_cells_vy )then
                        UPDATE_CLOSEST_PARTICLE_ARRAYS_USING_NEIGHBOR_CELLS(0,0,0, 1)
                    end if

                end if

                k = closest_particle(i,j,l,m)
                SLL_ASSERT(k /= 0)

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
                    track_markers_outside_domain,              &
                    mesh_period_x,                             &
                    mesh_period_y,                             &
                    h_parts_x,                                 &
                    h_parts_y,                                 &
                    h_parts_vx,                                &
                    h_parts_vy,                                &
                    inv_h_parts_x,                             &
                    inv_h_parts_y,                             &
                    inv_h_parts_vx,                            &
                    inv_h_parts_vy,                            &
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

               call get_initial_position_on_cartesian_grid_from_particle_index(k,   &
                    p_group%number_parts_x, p_group%number_parts_y,                 &
                    p_group%number_parts_vx, p_group%number_parts_vy,               &
                    j_x,j_y,j_vx,j_vy)
               x_k_t0 =  parts_x_min  + (j_x-1)  * h_parts_x
               y_k_t0 =  parts_y_min  + (j_y-1)  * h_parts_y
               vx_k_t0 = parts_vx_min + (j_vx-1) * h_parts_vx
               vy_k_t0 = parts_vy_min + (j_vy-1) * h_parts_vy

               ! <<loop_on_virtual_particles_in_one_virtual_cell>>
               ! [[file:~/mcp/maltpic/ltpic-bsl.tex::algo:pic-vr:find_f0_for_each_virtual_particle]] Loop over all
               ! virtual particles in the cell to compute the value of f0 at that point (Following
               ! [[file:~/mcp/maltpic/ltpic-bsl.tex::BSL_remapping_algo]])



               ! i_x, i_y, i_vx, i_vy: real index of the virtual particle in
               ! [[file:../pic_particle_types/lt_pic_4d_group.F90::target_values]]

               ! x, y, vx, vy = will be the location of the virtual particle at time n
               ! x =  virtual_parts_x_min  + (i-1)*h_virtual_cell_x  + (ivirt-1)*h_virtual_parts_x
               ! y =  virtual_parts_y_min  + (j-1)*h_virtual_cell_y  + (jvirt-1)*h_virtual_parts_y
               ! vx = virtual_parts_vx_min + (l-1)*h_virtual_cell_vx + (lvirt-1)*h_virtual_parts_vx
               ! vy = virtual_parts_vy_min + (m-1)*h_virtual_cell_vy + (mvirt-1)*h_virtual_parts_vy

               i_x = (i-1) * n_virtual_x    ! this index is needed in the "write f on grid" scenario
               x =  virtual_parts_x_min  + (i_x-1) * h_virtual_parts_x
               dx_in_virtual_cell = - 0.5 * h_virtual_parts_x
               x_to_xk = x - x_k

               do ivirt = 1, n_virtual_x

                  i_x = i_x + 1
                  SLL_ASSERT( i_x == (i-1)*n_virtual_x + ivirt )

                  !if( (.not. scenario_is_deposition) .and. use_remapping_grid) then
                  !  if( (i > p_group%number_parts_x * 95/100) &
                  !      .and. (j == p_group%number_parts_y /2) &
                  !      .and. (l == p_group%number_parts_vx /2) &
                  !      .and. (m == p_group%number_parts_vy /2) )then
                  !      print *, "675465437545 -- i, ivirt, i_x = ", i, ivirt, i_x
                  !  end if
                  !end if

                  x =       x       + h_virtual_parts_x
                  x_to_xk = x_to_xk + h_virtual_parts_x

                  dx_in_virtual_cell = dx_in_virtual_cell + h_virtual_parts_x       ! for the deposition scenario

                  d1_x = d11 * x_to_xk
                  d2_x = d21 * x_to_xk
                  d3_x = d31 * x_to_xk
                  d4_x = d41 * x_to_xk

                  i_y = (j-1) * n_virtual_y     ! this index is needed in the "write f on grid" scenario
                  y =  virtual_parts_y_min + (i_y-1)*h_virtual_parts_y
                  dy_in_virtual_cell = - 0.5 * h_virtual_parts_y
                  y_to_yk = y - y_k
                  do jvirt = 1, n_virtual_y
                     i_y = i_y + 1
                     SLL_ASSERT( i_y == (j-1)*n_virtual_y + jvirt )
                     y =       y       + h_virtual_parts_y
                     y_to_yk = y_to_yk + h_virtual_parts_y

                     dy_in_virtual_cell = dy_in_virtual_cell + h_virtual_parts_y        ! for the deposition scenario

                     d1_y = d12 * y_to_yk
                     d2_y = d22 * y_to_yk
                     d3_y = d32 * y_to_yk
                     d4_y = d42 * y_to_yk

                     i_vx = (l-1) * n_virtual_vx     ! this index is needed in the "write f on grid" scenario
                     vx = virtual_parts_vx_min + (i_vx-1)*h_virtual_parts_vx
                     vx_to_vxk = vx - vx_k
                     do lvirt = 1, n_virtual_vx
                        i_vx = i_vx + 1
                        SLL_ASSERT( i_vx == (l-1)*n_virtual_vx + lvirt )

                        vx =        vx        + h_virtual_parts_vx
                        vx_to_vxk = vx_to_vxk + h_virtual_parts_vx

                        d1_vx = d13 * vx_to_vxk
                        d2_vx = d23 * vx_to_vxk
                        d3_vx = d33 * vx_to_vxk
                        d4_vx = d43 * vx_to_vxk

                        i_vy = (m-1) * n_virtual_vy      ! this index is needed in the "write f on grid" scenario
                        vy = virtual_parts_vy_min + (i_vy-1)*h_virtual_parts_vy
                        vy_to_vyk = vy - vy_k
                        do mvirt = 1, n_virtual_vy
                           i_vy = i_vy + 1
                           SLL_ASSERT( i_vy == (m-1)*n_virtual_vy + mvirt )

                           vy =        vy        + h_virtual_parts_vy
                           vy_to_vyk = vy_to_vyk + h_virtual_parts_vy

                           d1_vy = d14 * vy_to_vyk
                           d2_vy = d24 * vy_to_vyk
                           d3_vy = d34 * vy_to_vyk
                           d4_vy = d44 * vy_to_vyk


                           ! The index may go out of the domain for higher values of x,y,vx,vy in each dimension
                           ! (because the corners of the corresponding virtual cell do not correspond to existing
                           ! real particles). In that case, just ignore that value.

                           if( scenario_is_deposition                                    &
                                .or. (      i_x  <= number_virtual_particles_x     &
                                      .and. i_y  <= number_virtual_particles_y     &
                                      .and. i_vx <= number_virtual_particles_vx    &
                                      .and. i_vy <= number_virtual_particles_vy) ) then

                              ! Location of virtual particle (ivirt,jvirt,lvirt,mvirt) at time 0 _relative_ to the
                              ! position of particle k at time 0 (x_k_t0,y_k_t0,vx_k_t0,vy_k_t0) according to flow
                              ! deformation


!                              x_t0  = d11 * (x - x_k) + d12 * (y - y_k) + d13 * (vx - vx_k) + d14 * (vy - vy_k)
!                              y_t0  = d21 * (x - x_k) + d22 * (y - y_k) + d23 * (vx - vx_k) + d24 * (vy - vy_k)
!                              vx_t0 = d31 * (x - x_k) + d32 * (y - y_k) + d33 * (vx - vx_k) + d34 * (vy - vy_k)
!                              vy_t0 = d41 * (x - x_k) + d42 * (y - y_k) + d43 * (vx - vx_k) + d44 * (vy - vy_k)


                              if( p_group%use_exact_f0 )then

                                ! here (x_t0, y_t0, vx_t0, vy_t0) is the (approx) position of the virtual particle at time t=0

                                x_t0  = d1_x + d1_y + d1_vx + d1_vy + x_k_t0
                                y_t0  = d2_x + d2_y + d2_vx + d2_vy + y_k_t0
                                vx_t0 = d3_x + d3_y + d3_vx + d3_vy + vx_k_t0
                                vy_t0 = d4_x + d4_y + d4_vx + d4_vy + vy_k_t0

                                ! WARNING -- this is only valid for the landau damping case...

                                ! -- fast (?) approximation of cos(x_t0)

                                ! s = smin_table + hs_table * (i_table + ds_table) with integer i_table >= 0 and 0 <= ds_table <1
                                s = modulo( k_landau * x_t0, 2*sll_pi)
                                s_aux =  s / hs_cos_table       ! because xmin_cos_table = 0
                                i_table = int(s_aux)
                                SLL_ASSERT(i_table >= 0)
                                SLL_ASSERT(i_table < ncells_table)
                                ds_table = s_aux - i_table
                                cos_approx = (1-ds_table) * cos_table(i_table) + ds_table * cos_table(i_table+1)

                                ! -- fast (?) approximation of exp(vx_t0**2 + vy_t0**2)
                                s = -0.5 * one_over_thermal_velocity_squared * (vx_t0**2 + vy_t0**2)
                                s_aux = (s - smin_exp_table) / hs_exp_table
                                i_table = int(s_aux)
                                if(i_table >= 0 .and. i_table < ncells_table)then
                                    ds_table = s_aux - i_table
                                    exp_approx = (1-ds_table) * exp_table(i_table) + ds_table * exp_table(i_table+1)
                                else
                                    exp_approx = 0
                                end if
                                f_value_on_virtual_particle = one_over_two_pi                   &
                                    * (1._f64 + alpha_landau * cos_approx)                      &
                                    * one_over_thermal_velocity_squared                         &
                                    * exp_approx

                                !f_value_on_virtual_particle = one_over_two_pi                   &
                                !    * (1._f64 + alpha_landau * cos(k_landau * (x_t0))) &
                                !    * one_over_thermal_velocity_squared                         &
                                !    * exp(-0.5 * one_over_thermal_velocity_squared              &
                                !          * (vx_t0**2 + vy_t0**2)       &
                                !         )

                              else

                                  find_k_prime_step_by_step = .false.

                                  if( find_k_prime_step_by_step )then

                                      ! here (x_t0, y_t0, vx_t0, vy_t0) is the (approx) position of the virtual particle at time t=0,
                                      ! RELATIVE to the k-th particle (marker) position at time t=0

                                      print *,"ERROR 876549 -- here x_t0 and the other variables should be recalled -- as RELATIVE"
                                      stop
                                      x_t0  = d1_x + d1_y + d1_vx + d1_vy
                                      y_t0  = d1_x + d1_y + d1_vx + d1_vy
                                      vx_t0 = d1_x + d1_y + d1_vx + d1_vy
                                      vy_t0 = d1_x + d1_y + d1_vx + d1_vy

                                        !! MCP: [DEBUG] store the (computed) absolute initial position of the virtual particle
                                        !if((.not. scenario_is_deposition) .and. use_remapping_grid)then
                                        ! p_group%debug_bsl_remap(i_x,i_y,i_vx,i_vy,1,1) = x_k_t0 + x_t0
                                        ! p_group%debug_bsl_remap(i_x,i_y,i_vx,i_vy,2,1) = y_k_t0 + y_t0
                                        ! p_group%debug_bsl_remap(i_x,i_y,i_vx,i_vy,3,1) = vx_k_t0 + vx_t0
                                        ! p_group%debug_bsl_remap(i_x,i_y,i_vx,i_vy,4,1) = vy_k_t0 + vy_t0
                                        !endif

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

                                        !! MCP: [DEBUG] store the (computed) absolute initial position of the virtual particle
                                        !if((.not. scenario_is_deposition) .and. use_remapping_grid)then
                                        !p_group%debug_bsl_remap(i_x,i_y,i_vx,i_vy,1,2) = x_k_t0 + x_t0
                                        !p_group%debug_bsl_remap(i_x,i_y,i_vx,i_vy,2,2) = y_k_t0 + y_t0
                                        !p_group%debug_bsl_remap(i_x,i_y,i_vx,i_vy,3,2) = vx_k_t0 + vx_t0
                                        !p_group%debug_bsl_remap(i_x,i_y,i_vx,i_vy,4,2) = vy_k_t0 + vy_t0
                                        !endif

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

                                      !if( (.not. scenario_is_deposition) .and. use_remapping_grid) then
                                      !  if( (i > p_group%number_parts_x * 95/100) &
                                      !      .and. (j == p_group%number_parts_y /2) &
                                      !      .and. (l == p_group%number_parts_vx /2) &
                                      !      .and. (m == p_group%number_parts_vy /2) )then
                                      !      print *, "77654864 -- i, ivirt, i_x = ", i, ivirt, i_x
                                      !      print *, "77654865 -- found step by step: kprime = ", kprime
                                      !  end if
                                      !end if

                                  else

                                    ! find directly the index kprime of the lower left particle on the remapping grid
                                    ! that contains z_t0, the position of the virtual particle at time = 0

                                    ! here (x_t0, y_t0, vx_t0, vy_t0) is the (approx) position of the virtual particle at time t=0

                                    x_t0_to_xk_t0   = d1_x + d1_y + d1_vx + d1_vy
                                    y_t0_to_yk_t0   = d2_x + d2_y + d2_vx + d2_vy
                                    vx_t0_to_vxk_t0 = d3_x + d3_y + d3_vx + d3_vy
                                    vy_t0_to_vyk_t0 = d4_x + d4_y + d4_vx + d4_vy

                                    x_t0 = x_t0_to_xk_t0 + x_k_t0
                                    y_t0 = y_t0_to_yk_t0 + y_k_t0
                                    vx_t0 = vx_t0_to_vxk_t0 + vx_k_t0
                                    vy_t0 = vy_t0_to_vyk_t0 + vy_k_t0

                                    call periodic_correction(p_group,x_t0,y_t0)

                                    tmp = (x_t0 - parts_x_min) / h_parts_x
                                    j_x  = 1 + int(floor(tmp))
                                    x_t0_to_xkprime_t0 = x_t0 - (parts_x_min + (j_x-1)*h_parts_x)

                                    tmp = (y_t0 - parts_y_min) / h_parts_y
                                    j_y  = 1 + int(floor(tmp))
                                    y_t0_to_ykprime_t0 = y_t0 - (parts_y_min + (j_y-1) * h_parts_y)

                                    tmp = (vx_t0 - parts_vx_min) / h_parts_vx
                                    j_vx  = 1 + int(floor(tmp))
                                    vx_t0_to_vxkprime_t0 = vx_t0 - (parts_vx_min + (j_vx-1) * h_parts_vx)

                                    tmp = (vy_t0 - parts_vy_min) / h_parts_vy
                                    j_vy  = 1 + int(floor(tmp))
                                    vy_t0_to_vykprime_t0 = vy_t0 - (parts_vy_min + (j_vy-1) * h_parts_vy)

                                    call get_particle_index_from_initial_position_on_cartesian_grid (               &
                                                            j_x, j_y, j_vx, j_vy,                                   &
                                                            p_group%number_parts_x, p_group%number_parts_y,         &
                                                            p_group%number_parts_vx, p_group%number_parts_vy,       &
                                                            kprime )


                                      !if( (.not. scenario_is_deposition) .and. use_remapping_grid) then
                                      !  if( (i > p_group%number_parts_x * 95/100) &
                                      !      .and. (j == p_group%number_parts_y /2) &
                                      !      .and. (l == p_group%number_parts_vx /2) &
                                      !      .and. (m == p_group%number_parts_vy /2) )then
                                      !      print *, "77654878 -- i, ivirt, i_x = ", i, ivirt, i_x
                                      !      print *, "77654879 -- found directly: kprime = ", kprime
                                      !      print *, "77654879 ---- "
                                      !      tmp = parts_x_min + (j_x-1)* h_parts_x
                                      !      print *, "10 ---- xmin_part_cell  = ", tmp
                                      !      print *, "10 ---- x_t0            = ", x_t0
                                      !      print *, "10 ---- xmax_part_cell  = ", tmp + h_parts_vx
                                      !      tmp = parts_y_min + (j_y-1)* h_parts_y
                                      !      print *, "11 ---- ymin_part_cell  = ", tmp
                                      !      print *, "11 ---- y_t0            = ", y_t0
                                      !      print *, "11 ---- ymax_part_cell  = ", tmp + h_parts_y
                                      !      tmp = parts_vx_min + (j_vx-1)* h_parts_vx
                                      !      print *, "12 ---- vxmin_part_cell  = ", tmp
                                      !      print *, "12 ---- vx_t0            = ", vx_t0
                                      !      print *, "12 ---- vxmax_part_cell  = ", tmp + h_parts_vx
                                      !      tmp = parts_vy_min + (j_vy-1)* h_parts_vy
                                      !      print *, "13 ---- vymin_part_cell  = ", tmp
                                      !      print *, "13 ---- vy_t0            = ", vy_t0
                                      !      print *, "13 ---- vymax_part_cell  = ", tmp + h_parts_vy
                                      !  end if
                                      !end if



                                    !    if( kprime <= 0 .or. kprime > p_group%number_particles )then
                                    !    ! this means we are outside of the domain defined by the initial particle grid
                                    !    kprime = 0
                                    !    end if

                                  end if

                                  SLL_ASSERT(kprime >= 0)

                                  ! If we end up with kprime == 0, it means that we have not found a cell that contains
                                  ! the particle so we just set that (virtual) particle value to zero

                                  f_value_on_virtual_particle = 0

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


                                          !if( (.not. scenario_is_deposition) .and. use_remapping_grid) then
                                          !  if( (i > p_group%number_parts_x * 85/100) &
                                          !      .and. (j == p_group%number_parts_y /2) &
                                          !      .and. (l == p_group%number_parts_vx /2) &
                                          !      .and. (m == p_group%number_parts_vy /2) )then
                                          !      print *, "77654889 -- i, ivirt, i_x = ", i, ivirt, i_x
                                          !      print *, "77654890 -- hcube ok -- --  "
                                          !  end if
                                          !end if

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

                                          ! MCP: [BEGIN-DEBUG] store the (computed) absolute initial position of the virtual particle
                                          ! [[get_initial_position_on_cartesian_grid_from_particle_index]]
                                            !   call get_initial_position_on_cartesian_grid_from_particle_index(kprime, &
                                            !        p_group%number_parts_x, p_group%number_parts_y,     &
                                            !        p_group%number_parts_vx, p_group%number_parts_vy,   &
                                            !        j_x,j_y,j_vx,j_vy)
                                            !   x_kprime_t0 =  parts_x_min  + (j_x-1)  * h_parts_x
                                            !   y_kprime_t0 =  parts_y_min  + (j_y-1)  * h_parts_y
                                            !   vx_kprime_t0 = parts_vx_min + (j_vx-1) * h_parts_vx
                                            !   vy_kprime_t0 = parts_vy_min + (j_vy-1) * h_parts_vy
                                            !
                                            !   if((.not. scenario_is_deposition) .and. use_remapping_grid)then
                                            !      p_group%debug_bsl_remap(i_x,i_y,i_vx,i_vy,1,3) = x_kprime_t0 + x_t0
                                            !      p_group%debug_bsl_remap(i_x,i_y,i_vx,i_vy,2,3) = y_kprime_t0 + y_t0
                                            !      p_group%debug_bsl_remap(i_x,i_y,i_vx,i_vy,3,3) = vx_kprime_t0 + vx_t0
                                            !      p_group%debug_bsl_remap(i_x,i_y,i_vx,i_vy,4,3) = vy_kprime_t0 + vy_t0
                                            !   endif

                                          ! MCP [END-DEBUG]

                                        ! [[file:~/mcp/maltpic/ltpic-bsl.tex::affine-fn*]] use the values of f0 at these
                                        ! neighbours to interpolate the value of f0 at
                                        ! [[file:~/mcp/maltpic/ltpic-bsl.tex::hat-bz*]]. MCP -> oui. Ici si tu utilises
                                        ! des particules affines (part_deg = 1) la valeur de f0 se déduit de celle du
                                        ! poids de la particule.  En fait tu peux utiliser une formule semblable à celle
                                        ! qui est utilisée dans la fonction sll_lt_pic_4d_write_f_on_remap_grid, mais
                                        ! sans faire intervenir la matrice de déformation à l'intérieur des splines.

                                        ! place the resulting value of f on the virtual particle in
                                        ! p_group%target_values

                                        do side_x = 1,2
                                           if(side_x == 1)then
                                              x_aux = x_t0_to_xkprime_t0
                                           else
                                              x_aux = h_parts_x - x_t0_to_xkprime_t0
                                           end if
                                           do side_y = 1,2
                                              if(side_y == 1)then
                                                 y_aux = y_t0_to_ykprime_t0
                                              else
                                                 y_aux = h_parts_y - y_t0_to_ykprime_t0
                                              end if
                                              do side_vx = 1,2
                                                 if(side_vx == 1)then
                                                    vx_aux = vx_t0_to_vxkprime_t0
                                                 else
                                                    vx_aux = h_parts_vx - vx_t0_to_vxkprime_t0
                                                 end if
                                                 do side_vy = 1,2
                                                    if(side_vy == 1)then
                                                       vy_aux = vy_t0_to_vykprime_t0
                                                    else
                                                       vy_aux = h_parts_vy - vy_t0_to_vykprime_t0
                                                    end if

                                                    ! uses [[sll_pic_shape]]
                                                    f_value_on_virtual_particle =                                                  &
                                                        f_value_on_virtual_particle                                                &
                                                          + p_group%p_list(hcube(side_x,side_y,side_vx,side_vy))%q                 &
                                                         * sll_pic_shape(part_degree,x_aux,y_aux,vx_aux,vy_aux,                    &
                                                                         inv_h_parts_x,inv_h_parts_y,inv_h_parts_vx,inv_h_parts_vy)

                                                 end do
                                              end do
                                           end do
                                        end do

                                          !if( (.not. scenario_is_deposition) .and. use_remapping_grid) then
                                          !  if( (i > p_group%number_parts_x * 95/100) &
                                          !      .and. (j == p_group%number_parts_y /2) &
                                          !      .and. (l == p_group%number_parts_vx /2) &
                                          !      .and. (m == p_group%number_parts_vy /2) )then
                                          !      print *, "6754654 -- i, ivirt, i_x = ", i, ivirt, i_x
                                          !     print *, "6754654 -- f_value_on_virtual_particle = ", f_value_on_virtual_particle
                                          ! end if
                                          !end if




                                     end if     ! test to see whether the hyper cube was inside initial particle grid
                                  end if     ! test on (k_prime \=0 )
                              end if    ! test on (use_exact_f0)

                              ! now f_value_on_virtual_particle has been computed we can use it

                              if( f_value_on_virtual_particle /= 0 )then

                                     if( scenario_is_deposition )then

                                        virtual_charge = f_value_on_virtual_particle * phase_space_virtual_dvol

                                        tmp1 = (1.0_f64 - dx_in_virtual_cell)
                                        tmp2 = (1.0_f64 - dy_in_virtual_cell)

                                        charge_accumulator_cell%q_sw = charge_accumulator_cell%q_sw             &
                                                + virtual_charge * tmp1 * tmp2

                                        charge_accumulator_cell%q_se = charge_accumulator_cell%q_se             &
                                                + virtual_charge *  dx_in_virtual_cell * tmp2

                                        charge_accumulator_cell%q_nw = charge_accumulator_cell%q_nw             &
                                                + virtual_charge * tmp1 *  dy_in_virtual_cell

                                        charge_accumulator_cell%q_ne = charge_accumulator_cell%q_ne             &
                                                + virtual_charge *  dx_in_virtual_cell *  dy_in_virtual_cell

                                        if( present(given_total_density) )then
                                            deposited_density = deposited_density + virtual_charge
                                        end if


                                    else

                                        if( use_remapping_grid )then
                                            p_group%target_values(i_x,i_y,i_vx,i_vy) = f_value_on_virtual_particle
                                        else
                                            given_array_2d(i_x,i_vx) = given_array_2d(i_x,i_vx)         &
                                                    + f_value_on_virtual_particle * h_virtual_parts_y * h_virtual_parts_vy
                                        end if

                                    end if

                              end if

                           end if
                        end do
                     end do
                  end do
               end do
             end do
          end do
       end do
    end do

    !print *, "654 end remap -> DEBUG: ", p_group%number_parts_x,    &
    !                               p_group%number_parts_y/2,    &
    !                               p_group%number_parts_vx/2,   &
    !                               p_group%number_parts_vy/2
    !
    !print *, "655 end remap  , p_group%number_parts_x = ", p_group%number_parts_x
    !print *, "655 end remap  , number_virtual_particles_x = ", number_virtual_particles_x
    !
    !print *, "655 end remap -> DEBUG: p_group%target_values(p_group%number_parts_x - 2,  p_group%number_parts_y/2, "
    !print *,                                                "p_group%number_parts_vx/2, p_group%number_parts_vy/2) = "
    !print *, "               =  ", p_group%target_values(p_group%number_parts_x - 2,  p_group%number_parts_y/2, &
    !                                                     p_group%number_parts_vx/2, p_group%number_parts_vy/2)
    !
    !print *, "655 end remap -> DEBUG: p_group%target_values(p_group%number_parts_x - 1,  p_group%number_parts_y/2, "
    !print *,                                                "p_group%number_parts_vx/2, p_group%number_parts_vy/2) = "
    !print *, "               =  ", p_group%target_values(p_group%number_parts_x - 1,  p_group%number_parts_y/2, &
    !                                                     p_group%number_parts_vx/2, p_group%number_parts_vy/2)
    !
    !print *, "655 end remap -> DEBUG: p_group%target_values(p_group%number_parts_x,  p_group%number_parts_y/2, "
    !print *,                                                "p_group%number_parts_vx/2, p_group%number_parts_vy/2) = "
    !print *, "               =  ", p_group%target_values(p_group%number_parts_x,  p_group%number_parts_y/2, &
    !                                                     p_group%number_parts_vx/2, p_group%number_parts_vy/2)

    if( scenario_is_deposition .and. present(given_total_density) )then

        if( deposited_density == 0 )then
            print *, "WARNING (76576537475) -- total deposited charge is zero, which is strange..."
            print *, "                      -- (no charge correction in this case) "

        else
            charge_correction_factor = given_total_density / deposited_density

            do i = 1,num_virtual_cells_x
               do j = 1,num_virtual_cells_y

                  ! determining the index of the Poisson cell from i and j
                  i_cell = i + (j-1) * p_group%mesh%num_cells1

                  charge_accumulator_cell => q_accumulator%q_acc(i_cell)

                  charge_accumulator_cell%q_sw = charge_accumulator_cell%q_sw * charge_correction_factor
                  charge_accumulator_cell%q_se = charge_accumulator_cell%q_se * charge_correction_factor
                  charge_accumulator_cell%q_nw = charge_accumulator_cell%q_nw * charge_correction_factor
                  charge_accumulator_cell%q_ne = charge_accumulator_cell%q_ne * charge_correction_factor
               end do
            end do
        end if
    end if

  end subroutine sll_lt_pic_4d_write_f_on_grid_or_deposit



  ! update the arrays closest_particle and closest_particle_distance with the index of the given particle
  ! if closer to what had been stored up to now.

  ! x_aux : x_particle - x_min_virtual_mesh   and  similarly for y, vx, vy
  subroutine update_closest_particle_arrays(k_part,                         &
                                            x_aux, y_aux, vx_aux, vy_aux,   &
                                            i, j, l, m,                     &
                                            h_virtual_cell_x, h_virtual_cell_y, h_virtual_cell_vx, h_virtual_cell_vy,   &
                                            closest_particle,               &
                                            closest_particle_distance)

      sll_int32, intent(in) :: k_part
      sll_int32, intent(in) :: i, j, l, m
      sll_real64, intent(in) :: x_aux, y_aux, vx_aux, vy_aux
      sll_real64, intent(in) :: h_virtual_cell_x, h_virtual_cell_y, h_virtual_cell_vx, h_virtual_cell_vy

      sll_int32, dimension(:,:,:,:)  :: closest_particle
      sll_real64, dimension(:,:,:,:) :: closest_particle_distance

      sll_real64 :: square_dist_to_cell_center

      square_dist_to_cell_center = (x_aux  - (i + 0.5) * h_virtual_cell_x )**2.    &
                                 + (y_aux  - (j + 0.5) * h_virtual_cell_y )**2.    &
                                 + (vx_aux - (l + 0.5) * h_virtual_cell_vx )**2.    &
                                 + (vy_aux - (m + 0.5) * h_virtual_cell_vy )**2.

      ! if new particle is closer to center, keep the new one

      if(closest_particle(i,j,l,m) == 0 .or. square_dist_to_cell_center < closest_particle_distance(i,j,l,m)) then
         closest_particle(i,j,l,m) = k_part
         closest_particle_distance(i,j,l,m) = square_dist_to_cell_center
      end if

  end subroutine update_closest_particle_arrays


end module sll_lt_pic_4d_utilities
