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
  implicit none

 

contains

  subroutine sll_initialize_maxwellian(thermal_speed,p_group)
    type(sll_particle_group_2d), pointer :: p_group
    sll_int64 :: j
    sll_real64, dimension(:), allocatable :: nu
    sll_int32 :: ierr

    SLL_ALLOCATE( nu(1:NUM8PARTICLES) )
    do j=1,NUM_PARTICLES
       nu(j) = spatial_variance*sqrt( -2.0_f64*log(1.0_f64-(real(j)-0.5_f64)/NUM_PARTICLES) )
       p_group%p_list(j)%vx = nu(j)*cos(suite_hamm(j,2)*2.0_f64*sll_pi)
       p_group%p_list(j)%vy = nu(j)*sin(suite_hamm(j,2)*2.0_f64*sll_pi)
    enddo
    SLL_DEALLOCATE(nu, ierr)

  end subroutine sll_initialize_maxwellian


  subroutine compute_cell_and_offset( &
     x, &
     nc, &
     xmin, &
     i_cell, &
     offset )
    intent(in) :: x, nc, xmin
    intent(out) :: i_cell, offset
  end subroutine compute_cell_and_offset

  subroutine sll_initialize_particle(ncx, ncy, mesh_xmin, mesh_ymin, &
                                     mesh_dx )
    call compute_cell_and_offset(x, ncx, index_x)
    call compute_cell_and_offset(y, ncy, index_y)
    i_c= index_x
  end subroutine sll_initialize_particle

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

end module sll_particle_initializers
