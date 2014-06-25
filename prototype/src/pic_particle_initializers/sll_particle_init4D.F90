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

module sll_particle_initializers
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"
#include "particle_representation.h"

  use sll_constants, only: sll_pi
  use sll_particle_group_2d_module
  use sll_logical_meshes
  use gaussian
  use hammersley
  use sll_representation_conversion_module

  implicit none
  
!   private sll_init_spatial_particle2D, suite_hamm
! !  In the future, we have to build a directory called
! !  random_number_generators 


contains
  subroutine sll_initial_particles_4d( &
              thermal_speed, alpha, k, &
              m2d, &
              num_particles, &
              p_group )
    sll_real64, intent(in) :: thermal_speed, alpha, k
    type(sll_logical_mesh_2d), intent(in) :: m2d
    sll_int32, intent(in)  :: num_particles
    type(sll_particle_group_2d), pointer, intent(inout) :: p_group
    sll_int32  :: j
    sll_int32  :: ierr, ncx, ic_x,ic_y
    sll_real64 :: x, y, vx, vy, nu, xmin, ymin, rdx, rdy
    sll_real32 :: weight, off_x,off_y
    sll_real64 :: tmp1, tmp2

    weight = (m2d%eta1_max - m2d%eta1_min) * &
           (m2d%eta2_max - m2d%eta2_min)/real(num_particles,f32)
 
    rdx = 1._f64/m2d%delta_eta1
    rdy = 1._f64/m2d%delta_eta2
    xmin = m2d%eta1_min
    ymin = m2d%eta2_min
    ncx  = m2d%num_cells1

    open(90,file='initialparticles.dat')
    write(90,*) '#  POSITIONS in 2d    |||    VELOCITIES in 2d'
    j=1
    do while ( j <= num_particles )
       call random_number(x)
       x = (m2d%eta1_max - xmin)*x + xmin
       call random_number(y)
       y = 2._f64 * y! (2._f64*alpha)*y + 1._f64 - alpha
       if (eval_landau(alpha, k, x) >= y ) then
          y = (m2d%eta2_max - ymin)*suite_hamm(j,3) + ymin
          !
          nu = thermal_speed*sqrt( -2.0_f64*log(1.0_f64 - &
               (real(j,f64)-0.5_f64)/real(num_particles,f64)) )
          vx = nu * cos(suite_hamm(j,2)*2.0_f64*sll_pi)
          vy = nu * sin(suite_hamm(j,2)*2.0_f64*sll_pi)
          write(90,*) x, y, vx, vy 
!          call set_group_particle_values( p_group, j, x, y, vx, vy, weight)
          SET_PARTICLE_VALUES(p_group%p_list(j),x,y,vx,vy,weight,xmin,ymin,ncx,ic_x,ic_y,off_x,off_y,rdx,rdy,tmp1,tmp2)
          j = j + 1          
       endif
    end do
    print*, 'nb d essais', j
    close(90)
    
  end subroutine sll_initial_particles_4d
  
!!$  subroutine sll_initialize_some4Dfunction( &
!!$              thermal_speed, alpha, k, &
!!$              m2d, &
!!$              num_particles, &
!!$              p_group )
!!$
!!$    sll_real64, intent(in) :: thermal_speed, alpha, k
!!$    type(sll_logical_mesh_2d), intent(in) :: m2d
!!$    sll_int64, intent(in)  :: num_particles
!!$    type(sll_particle_group_2d), pointer, intent(inout) :: p_group
!!$    sll_int64 :: j
!!$    sll_int32 :: ierr
!!$    sll_real64 :: xn, yn, interm_y, nu
!!$
!!$!    p_group%p_list(:)%q = 1._f32/real(num_particles,f32) !  CHANGE 8/07
!!$    p_group%p_list(:)%q = (m2d%eta1_max - m2d%eta1_min) * &
!!$         (m2d%eta2_max - m2d%eta2_min)/real(num_particles,f32) !  CHANGE 11/07
!!$
!!$    do j=1,num_particles! standard Gaussian function in velocity 2D
!!$       nu = thermal_speed*sqrt( -2.0_f64*log(1.0_f64 - &
!!$            (real(j,f64)-0.5_f64)/real(num_particles,f64)) )
!!$       p_group%p_list(j)%vx = nu * cos(suite_hamm(j,2)*2.0_f64*sll_pi)
!!$       p_group%p_list(j)%vy = nu * sin(suite_hamm(j,2)*2.0_f64*sll_pi)
!!$    enddo
!!$! --- Here I should introduce particles by using 
!!$! --- call initialize_particle_2d()
!!$! --- in particle_representation.F90. 
!!$
!!$    !Rejection sampling for the function x --> 1+alpha*cos(k*x)
!!$    open(90,file='initial_random_parts.dat')
!!$    write(90,*) '#  POSITIONS         VITESSES'
!!$    j=1
!!$    do while (j<=num_particles)
!!$       call random_number(xn)
!!$       xn = (m2d%eta1_max - m2d%eta1_min)*xn + m2d%eta1_min
!!$       call random_number(yn)
!!$       yn = 2._f64 * yn! (2._f64*alpha)*yn + 1._f64 - alpha
!!$       if (eval_landau(alpha, k, xn) >= yn ) then
!!$          interm_y = (m2d%eta2_max - m2d%eta2_min)*suite_hamm(j,3) & 
!!$                     + m2d%eta2_min
!!$          write(90,*) xn, interm_y, p_group%p_list(j)%vx, p_group%p_list(j)%vy 
!!$          call global_to_cell_offset (  &
!!$               xn, interm_y, &
!!$               m2d, &
!!$               p_group%p_list(j)%ic, &
!!$               p_group%p_list(j)%dx, &
!!$               p_group%p_list(j)%dy )
!!$          j = j + 1          
!!$       endif
!!$    enddo
!!$    print*, 'nb d essais', j
!!$    close(90)
!!$    
!!$  end subroutine sll_initialize_some4Dfunction

  function eval_landau(alp, kx, x)
    sll_real64 :: alp, kx, x
    sll_real64 :: eval_landau
    eval_landau = 1._f64 + alp * cos(kx * x)
  end function eval_landau

end module sll_particle_initializers
