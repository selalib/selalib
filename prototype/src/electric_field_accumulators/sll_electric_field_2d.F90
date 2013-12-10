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

module sll_electric_field_2d_accumulator
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"
  implicit none

  ! The idea of having this data structure is to precompute the values of
  ! the electric field, store them redundantly on a cell-based structure and
  ! reduce the memory thrashing that would inevitably occur if we were to
  ! compute the values of the electric field from the potential repeatedly.
  ! These costs ought to be amortized if the number of points to be advected,
  ! or particles to be advanced (per cell) is relatively large compared with
  ! the cost of filling out this data structure.
  !
  ! The values of the electric field within the structure in 2D are named 
  ! like this:
  !                NW (north-west)       NE (north-east)
  !                  +--------------------+
  !                  |                    |
  !                  |                    |
  !                  |                    |
  !                  |                    |
  !                  |      cell(i,j)     |
  !                  |                    |
  !                  |                    |
  !                  |                    |
  !                  |                    |
  !                  +--------------------+
  !                SW (south-west)       SE (south-east)
  !
  ! Each corner represents a point with two electric field components defined:
  ! ex and ey (x-component and y-component of the field respectively).

  type :: electric_field_2d_accumulator
     sll_real64 :: ex_sw
     sll_real64 :: ey_sw
     sll_real64 :: ex_se
     sll_real64 :: ey_se
     sll_real64 :: ex_nw
     sll_real64 :: ey_nw
     sll_real64 :: ex_ne
     sll_real64 :: ey_ne
  end type electric_field_2d_accumulator

  ! The following is a pseudo-accumulator. It is written for use while some
  ! problems with the parallel loading of the previous accumulator are resolved.
  ! This one is only useful to precompute the values of the electric field in
  ! a structure that can have the same layout as the potential array.
  type :: efield_2d_point
     sll_real64 :: ex
     sll_real64 :: ey
  end type efield_2d_point

end module sll_electric_field_2d_accumulator
