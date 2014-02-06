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

  implicit none
  
  private suite_hamm! ! In the future, we have to build a directory called
!!!  ---                random_number_generators 


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
    sll_real64, dimension(:), allocatable :: nu
    sll_int32 :: ierr
    sll_real64 :: xn, yn, interm_value
    
    SLL_ALLOCATE( nu(1:num_particles), ierr )
    do j=1,num_particles
       nu(j) = thermal_speed*sqrt( -2.0_f64*log(1.0_f64-(real(j)-0.5_f64) &
               /num_particles) )
       p_group%p_list(j)%vx = nu(j)*cos(suite_hamm(j,2)*2.0_f64*sll_pi)
       p_group%p_list(j)%vy = nu(j)*sin(suite_hamm(j,2)*2.0_f64*sll_pi)
!
       interm_value = (m2d%eta2_max - m2d%eta2_min)*suite_hamm(j,3) & 
                     + m2d%eta2_min
       p_group%p_list(j)%ic = int( (interm_value-m2d%eta2_min)/m2d%delta_eta2 )
       p_group%p_list(j)%dy = mod( interm_value-m2d%eta2_min, m2d%delta_eta2 )
    enddo
    SLL_DEALLOCATE_ARRAY(nu, ierr)

    !Methode du rejet pour la fonction 1+alpha*cos(k*x)
    j=1
    do while (j<=num_particles)
       call random_number(xn)
       xn = (m2d%eta1_max - m2d%eta1_min)*xn + m2d%eta1_min
       call random_number(yn)
       yn = (2._f64*alpha)*yn + 1._f64 - alpha
       if (eval(alpha, k, xn) >= yn ) then
          p_group%p_list(j)%dx = mod( xn-m2d%eta1_min, m2d%delta_eta1 )
          p_group%p_list(j)%ic = p_group%p_list(j)%ic * m2d%num_cells1 &
               + int( (xn-m2d%eta1_min)/m2d%delta_eta1 ) + 1
          j = j + 1          
       endif
    enddo
    print*, 'nb d essais', j
    
  end subroutine sll_initialize_some4Dfunction


!!$  subroutine compute_cell_and_offset( &
!!$     x, &
!!$     nc, &
!!$     xmin, &
!!$     i_cell, &
!!$     offset )
!!$    intent(in) :: x, nc, xmin
!!$    intent(out) :: i_cell, offset
!!$  end subroutine compute_cell_and_offset
!!$
!!$  subroutine sll_initialize_particle(ncx, ncy, mesh_xmin, mesh_ymin, &
!!$                                     mesh_dx )
!!$    call compute_cell_and_offset(x, ncx, index_x)
!!$    call compute_cell_and_offset(y, ncy, index_y)
!!$    i_c= index_x
!!$  end subroutine sll_initialize_particle

!! !> @brief Returns
!! !> @param[in]
  function suite_hamm (n,b)
    sll_real64 :: suite_hamm
    sll_int64 :: n,m
    sll_int32 :: b
    sll_int64 :: u
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

  function eval(alp, kx, x)
    sll_real64 :: alp, kx, x
    sll_real64 :: eval
    eval = 1._f64 + alp * cos(kx * x)
  end function eval

end module sll_particle_initializers
