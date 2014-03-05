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

  use sll_constants, only: sll_pi
  use sll_particle_group_2d_module
  use sll_logical_meshes
  use sll_visu_pic

  implicit none
  
  private suite_hamm, sll_init_spatial_particle2D
! !  In the future, we have to build a directory called
! !  random_number_generators 


contains

  subroutine sll_initialize_some4Dfunction( &
              thermal_speed, alpha, k, &
              m2d, &
              num_particles, &
              p_group )

    sll_real64, intent(in) :: thermal_speed, alpha, k
    type(sll_logical_mesh_2d), intent(in) :: m2d
    sll_int64, intent(in)  :: num_particles
    type(sll_particle_group_2d), pointer, intent(inout) :: p_group
    sll_int64 :: j
    sll_int32 :: ierr
    sll_real64 :: xn, yn, interm_y, nu
    sll_real64, dimension(1:num_particles)  :: particles_X, particles_Y, poids
!!$    sll_real64, dimension(m2d%num_cells1,m2d%num_cells2) :: dens
!!!    character(len=5) :: nom


    do j=1,num_particles! standard Gaussian function in velocity 2D
       nu = thermal_speed*sqrt( -2.0_f64*log(1.0_f64 - &
            (real(j,f64)-0.5_f64)/real(num_particles,f64)) )
       p_group%p_list(j)%vx = nu * cos(suite_hamm(j,2)*2.0_f64*sll_pi)
       p_group%p_list(j)%vy = nu * sin(suite_hamm(j,2)*2.0_f64*sll_pi)
    enddo

    !Rejection sampling for the function x --> 1+alpha*cos(k*x)
    j=1
    open(90,file='initial_random_parts.dat')
    write(90,*) '#  POSITIONS         VITESSES'
    do while (j<=num_particles)
       call random_number(xn)
       xn = (m2d%eta1_max - m2d%eta1_min)*xn + m2d%eta1_min
       call random_number(yn)
       yn = (2._f64*alpha)*yn + 1._f64 - alpha
       if (eval_landau(alpha, k, xn) >= yn ) then
          interm_y = (m2d%eta2_max - m2d%eta2_min)*suite_hamm(j,3) & 
                     + m2d%eta2_min
          write(90,*) xn, interm_y
          particles_X(j) = xn
          particles_Y(j) = interm_y
          poids(j) = 1._f64
          call sll_init_spatial_particle2D(  &
               xn, interm_y, &
               m2d, &
               p_group%p_list(j)%ic, &
               p_group%p_list(j)%dx, &
               p_group%p_list(j)%dy     )
          j = j + 1          
       endif
    enddo
    print*, 'nb d essais', j
    close(90)
    
!!$    call compute_df_cic(particles_X, particles_Y, poids, &
!!$         m2d%eta1_min, m2d%eta1_max, m2d%num_cells1, &
!!$         m2d%eta2_min, m2d%eta2_max, m2d%num_cells2, dens)
    
    call distribution_xv_gnuplot( 'TEST', particles_X, particles_Y, &
         m2d%eta1_min, m2d%eta1_max, m2d%num_cells1, &
         m2d%eta2_min, m2d%eta2_max, m2d%num_cells2, 1, 0._f64)
    print*, 'OK for the writing'

  end subroutine sll_initialize_some4Dfunction


  subroutine compute_cell_and_offset( &
     x, &
     xmin, &
     dx, &
     i_cell, &
     offset )
    sll_real64, intent(in)  ::  x, xmin, dx
    sll_int32,  intent(out) ::  i_cell
    sll_real32, intent(out) ::  offset
    sll_real64  :: temp

    temp = (x - xmin)/dx
    i_cell = int(temp)
    offset = temp - real(i_cell,f64)

!!$    i_cell = int( (x - xmin)/dx )
!!$    offset = mod( x - xmin, dx )
!!$    offset = offset/dx! the cell for a charge accumulator is [0,1]x[0,1]
!                           and not [0,delta_x]x[0,delta_y]  !!!
  end subroutine compute_cell_and_offset

  subroutine sll_init_spatial_particle2D( x, y, &
                      m2d,   &
                      icell, &
                      offset_x, offset_y )
! transforms a particle position (x,y) in our type (icell, dx, dy)
    sll_real64, intent(in)  :: x, y
    type(sll_logical_mesh_2d), intent(in) :: m2d
    sll_int32,  intent(out) :: icell
    sll_real32, intent(out) :: offset_x, offset_y
    sll_int32               :: icell_x, icell_y

    call compute_cell_and_offset(x, m2d%eta1_min, m2d%delta_eta1, icell_x, offset_x)
    if ( (icell_x<0).or.(icell_x>m2d%num_cells1) ) print*,'ERROR: bad icell_x', icell_x
    if ( (offset_x<0).or.(offset_x.ge.1) ) print*, 'ERROR: bad offset_x', offset_x

    call compute_cell_and_offset(y, m2d%eta2_min, m2d%delta_eta2, icell_y, offset_y)
    if ( (icell_y<0).or.(icell_y>m2d%num_cells2) ) print*,'ERROR: bad icell_y', icell_y
    if ( (offset_y<0).or.(offset_y.ge.1) ) print*, 'ERROR: bad offset_y', offset_y

    icell = icell_x + 1 + icell_y * m2d%num_cells1
    if ( (icell<1).or.(icell>(m2d%num_cells1*m2d%num_cells2)) ) print*,'ERROR: bad icell', icell
  end subroutine sll_init_spatial_particle2D

!! !> @brief Returns
!! !> @param[in]
  function suite_hamm (n,b)
    sll_real64 :: suite_hamm
    sll_int64  :: n,m
    sll_int32  :: b
    sll_int64  :: u
    sll_real64 :: h,k
    
    k=0
    h=0
    m=n
    do while (m>0)
       k = k+1
       u = m / b
       h = h +  b**(-k) * (m - u*b)
       m = u
    enddo
    suite_hamm = h 
  end function suite_hamm

  function eval_landau(alp, kx, x)
    sll_real64 :: alp, kx, x
    sll_real64 :: eval_landau
    eval_landau = 1._f64 + alp * cos(kx * x)
  end function eval_landau

end module sll_particle_initializers
