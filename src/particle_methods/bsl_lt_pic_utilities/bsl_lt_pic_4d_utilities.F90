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

module sll_bsl_lt_pic_4d_utilities_module

#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"
!#include "sll_accumulators.h"
!#include "particle_representation.h"

  use sll_constants, only: sll_pi
!  use sll_bsl_lt_pic_4d_group_module
!  use sll_bsl_lt_pic_4d_particle_module

!  use sll_representation_conversion_module
  use sll_timer
  use sll_utilities, only: sll_new_file_id, int2string
  use sll_gnuplot

implicit none 

contains

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

end subroutine get_particle_index_from_initial_position_on_cartesian_grid

end module  sll_bsl_lt_pic_4d_utilities_module