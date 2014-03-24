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


module sll_accumulators
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"

  use sll_particle_representations
  implicit none
  
  type charge_accumulator_cell! for particles deposition on the grid
     sll_real64 :: q_sw
     sll_real64 :: q_se
     sll_real64 :: q_nw
     sll_real64 :: q_ne
  end type charge_accumulator_cell

  type field_accumulator_cell! for interpolating the elec field from
                     ! its grid values(in this type) in the particles
     sll_real64 :: Ex_sw
     sll_real64 :: Ex_se
     sll_real64 :: Ex_nw
     sll_real64 :: Ex_ne
     sll_real64 :: Ey_sw
     sll_real64 :: Ey_se
     sll_real64 :: Ey_nw
     sll_real64 :: Ey_ne
  end type field_accumulator_cell

  interface sll_delete
     module procedure delete_accumulate_charge
  end interface sll_delete

contains
  
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
  
  subroutine delete_accumulate_charge( allcharge )
     type(charge_accumulator_cell), dimension(:), pointer :: allcharge
     sll_int32  :: ierr
     
     if ( .not.associated(allcharge) ) then
        print*, 'delete_accumulate_charge: ERROR, charge was not associated'
     endif
     SLL_DEALLOCATE(allcharge, ierr)
   end subroutine delete_accumulate_charge

end module sll_accumulators
