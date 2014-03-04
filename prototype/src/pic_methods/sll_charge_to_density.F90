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


module sll_charge_to_density
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"

  use sll_accumulators

  use sll_particle_representations

  implicit none
  
  type charge_accumulator_cell! for particles deposition on the grid
     sll_real64 :: q_sw
     sll_real64 :: q_se
     sll_real64 :: q_nw
     sll_real64 :: q_ne
  end type charge_accumulator_cell



contains

  subroutine sll_charge_meshdensity( all_charge, m2d, density )

    type(sll_logical_mesh_2d), intent(in) :: m2d
    type(charge_accumulator_cell), dimension(1:m2d%num_cells2*m2d%num_cells1), intent(in) :: all_charge
    sll_real64, dimension(1:1+m2d%num_cells1, 1:1+m2d%num_cells2) , intent(inout) :: density
    sll_int32 :: k, i, j
    sll_int32 :: ncx, ncy
    sll_real64 :: delta_x, delta_y

    ncx = m2d%num_cells1
    ncy = m2d%num_cells2
    delta_x = m2d%delta_eta1
    delta_y = m2d%delta_eta2

    do k = 1, ncx * ncy
       i = mod( k, ncx )
       j = int( k/ ncy  )
       density(i  ,j  ) = density(i,j)    + all_charge(k)%q_sw
       density(i+1,j  ) = density(i+1,j)  + all_charge(k)%q_se
       density(i  ,j+1) = density(i,j+1)  + all_charge(k)%q_nw
       density(i+1,j+1) = density(i+1,j+1)+ all_charge(k)%q_ne
    enddo
    density = density/(delta_x*delta_y)
  end subroutine sll_charge_meshdensity


  subroutine sll_accumulate_charge( &
                particle, &
                charge )

    type(sll_particle_2d), intent(in) :: particle
    type(charge_accumulator_cell), intent(inout) :: charge
    sll_int64    ::  j
    
    charge%q_sw = charge%q_sw + &
         particle%q * (1._f64 - particle%dx) * (1._f64 - particle%dy)

    charge%q_se = charge%q_se + &
         particle%q * particle%dx * (1._f64 - particle%dy)

    charge%q_nw = charge%q_nw + &
         particle%q * (1._f64 - particle%dx) * particle%dy

    charge%q_ne = charge%q_ne + &
         particle%q * particle%dx * particle%dy

  end subroutine sll_accumulate_charge


  function new_accumulate_charge( &
       cells_number ) result(cha)

    sll_int32, intent(in) :: cells_number
    type(charge_accumulator_cell), dimension(:), pointer :: cha
    sll_int32  :: ierr

    SLL_ALLOCATE( cha(1:cells_number), ierr)
    cha(:)%q_sw = 0._f64
    cha(:)%q_se = 0._f64
    cha(:)%q_nw = 0._f64
    cha(:)%q_ne = 0._f64

  end function new_accumulate_charge
  

 end module sll_charge_to_density
