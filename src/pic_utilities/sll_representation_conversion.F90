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


module sll_representation_conversion_module
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"

  use sll_cartesian_meshes

  implicit none
  
contains

  !> given x, x_min and inverse_delta_x = 1/delta_x
  !! this computes i_cell_x (integer) and 0 <= offset_x < 1 such that
  !! x = x_min + (i_cell_x + offset_x) * delta_x

  subroutine compute_cell_and_offset( &
                                     x, &
                                     x_min, &
                                     inverse_delta_x, &
                                     i_cell_x, &
                                     offset_x )

    sll_real64, intent(in)  ::  x
    sll_real64, intent(in)  ::  x_min
    sll_real64, intent(in)  ::  inverse_delta_x
    sll_int32,  intent(out) ::  i_cell_x
    sll_real32, intent(out) ::  offset_x
    sll_real64 :: temp

    temp = (x - x_min) * inverse_delta_x
    i_cell_x  = int(temp)
    offset_x = temp - real(i_cell_x,f64)

  end subroutine compute_cell_and_offset


  !> given x, y and a 2d mesh m2d,
  !! this computes i_cell (integer) and 0 <= offset_x, offset_y < 1 such that
  !! x = x_min + (i_cell_x + offset_x) * delta_x
  !! and
  !! y = y_min + (i_cell_y + offset_y) * delta_y
  !! where i_cell = i_cell_x + 1 + i_cell_y * m2d%num_cells1
  !! and x_min,  y_min, delta_x, delta_y are derived from the 2d mesh
  !!
  !! WARNING -----  here i_cell_x goes from 0 to m2d%num_cells1 - 1
  !!                (and i_cell_y goes from 0 to m2d%num_cells2 - 1)
  !! BUT:
  !!                i_cell goes from 1 to m2d%num_cells1 * m2d%num_cells2  !!
  !! (beurk)

  subroutine global_to_cell_offset( x, y, &
                      m2d,   &
                      i_cell, &
                      offset_x, offset_y )
  ! transforms a standard particle position (x,y) in our type (i_cell, dx, dy)

    sll_real64, intent(in)  :: x, y
    type(sll_cartesian_mesh_2d), intent(in) :: m2d
    sll_int32,  intent(out) :: i_cell
    sll_real32, intent(out) :: offset_x, offset_y!    sll_real64, intent(out) :: offset_x, offset_y!
    sll_int32               :: i_cell_x, i_cell_y
    sll_real64              :: inverse_delta_x, inverse_delta_y

    inverse_delta_x = 1._f64/m2d%delta_eta1
    call compute_cell_and_offset( x, m2d%eta1_min, &
                                  inverse_delta_x, i_cell_x, &
                                  offset_x )

    if ( (i_cell_x<0).or.(i_cell_x.ge.m2d%num_cells1) ) &
         print*,'ERROR: bad i_cell_x', i_cell_x

    inverse_delta_y = 1._f64/m2d%delta_eta2
    call compute_cell_and_offset( y, m2d%eta2_min, &
                                  inverse_delta_y, i_cell_y, &
                                  offset_y )

    if ( (i_cell_y<0).or.(i_cell_y.ge.m2d%num_cells2) ) &
         print*,'ERROR: bad i_cell_y', i_cell_y

    i_cell = i_cell_x + 1 + i_cell_y * m2d%num_cells1

    if ( (i_cell<1).or.(i_cell>(m2d%num_cells1*m2d%num_cells2)) ) &
         print*,'ERROR: bad i_cell', i_cell
  end subroutine global_to_cell_offset


  !> performs the inverse transformation of global_to_cell_offset above
  !!

  subroutine cell_offset_to_global ( offset_x, offset_y, &
                                   i_cell, m2d, &
                                   x, y )
  ! transforms sll_type of a particle (i_cell, dx, dy) into the standard
  ! particle position (x,y)
    sll_real64, intent(out)  :: x, y
    type(sll_cartesian_mesh_2d), intent(in) :: m2d
    sll_int32,  intent(in) :: i_cell
    sll_real32, intent(in) :: offset_x, offset_y

    x = m2d%delta_eta1*( offset_x + real( mod(i_cell-1,m2d%num_cells1), f64) )
    y = m2d%delta_eta2*( offset_y + real( int( (i_cell-1)/m2d%num_cells1 ), f64))

  end subroutine cell_offset_to_global
	
end module sll_representation_conversion_module
