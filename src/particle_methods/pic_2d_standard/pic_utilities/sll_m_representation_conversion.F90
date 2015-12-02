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


module sll_m_representation_conversion
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"

  use sll_m_cartesian_meshes

  implicit none
  
contains

  !> given x, x_min and inverse_delta_x = 1/delta_x
  !! this computes i_cell_x_from_0 (integer) and 0 <= offset_x < 1 such that
  !! x = x_min + (i_cell_x_from_0 + offset_x) * delta_x

  subroutine compute_cell_and_offset( &
                                     x, &
                                     x_min, &
                                     inverse_delta_x, &
                                     i_cell_x_from_0, &
                                     offset_x )

    sll_real64, intent(in)  ::  x
    sll_real64, intent(in)  ::  x_min
    sll_real64, intent(in)  ::  inverse_delta_x
    sll_int32,  intent(out) ::  i_cell_x_from_0
    sll_real32, intent(out) ::  offset_x
    sll_real64 :: temp

    temp = (x - x_min) * inverse_delta_x
    i_cell_x_from_0  = int(temp)
    offset_x = real(temp - real(i_cell_x_from_0,f64), f32)

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

  subroutine global_to_cell_offset( x, y, &
                      m2d,   &
                      i_cell, &
                      offset_x, offset_y )
  ! transforms a standard particle position (x,y) in our type (i_cell, dx, dy)

    sll_real64, intent(in)  :: x, y
    type(sll_cartesian_mesh_2d), intent(in) :: m2d
    sll_int32,  intent(out) :: i_cell       !! watch out: from 1 to m2d%num_cells1 * m2d%num_cells2
    sll_real32, intent(out) :: offset_x, offset_y
    sll_int32               :: i_cell_x_from_0     !! watch out: from 0 to m2d%num_cells1 - 1
    sll_int32               :: i_cell_y_from_0     !! watch out: from 0 to m2d%num_cells2 - 1
    sll_real64              :: inverse_delta_x, inverse_delta_y

    inverse_delta_x = 1._f64/m2d%delta_eta1
    call compute_cell_and_offset( x, m2d%eta1_min, &
                                  inverse_delta_x, i_cell_x_from_0, &
                                  offset_x )

    if ( (i_cell_x_from_0 < 0).or.(i_cell_x_from_0 .ge. m2d%num_cells1) ) &
         print*,'ERROR: bad i_cell_x_from_0', i_cell_x_from_0

    inverse_delta_y = 1._f64/m2d%delta_eta2
    call compute_cell_and_offset( y, m2d%eta2_min, &
                                  inverse_delta_y, i_cell_y_from_0, &
                                  offset_y )

    if ( (i_cell_y_from_0 < 0).or.(i_cell_y_from_0 .ge. m2d%num_cells2) ) &
         print*,'ERROR: bad i_cell_y_from_0', i_cell_y_from_0

    i_cell = i_cell_x_from_0 + i_cell_y_from_0 * m2d%num_cells1 + 1

    if ( (i_cell <= 0) .or. (i_cell > (m2d%num_cells1*m2d%num_cells2)) ) &
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

    x = m2d%eta1_min+m2d%delta_eta1*(real(offset_x,f64)+real(mod(i_cell-1,m2d%num_cells1), f64) )
    y = m2d%eta2_min+m2d%delta_eta2*(real(offset_y,f64)+real(int( (i_cell-1)/m2d%num_cells1 ), f64))

  end subroutine cell_offset_to_global

  !! transforms a standard particle position (x,y) in (i_cell_x, i_cell_y, dx, dy)
  !! -> here the indices i_cell_x and i_cell_y do not need to be within [1, m2d%num_cells1] or [1, m2d%num_cells2]
  !!    so that: - in periodic domains, the flows are better represented (no information is lost using modulo)
  !!             - in non-periodic domains we can track outside particles (markers)
  subroutine global_to_cell_offset_extended( x, y, &
                      m2d,      &
                      i_cell_x, &
                      i_cell_y, &
                      offset_x, offset_y )

    sll_real64, intent(in)  :: x, y
    type(sll_cartesian_mesh_2d), intent(in) :: m2d
    sll_int32,  intent(out) :: i_cell_x       !! not necessarily in [1, m2d%num_cells1], see comments above
    sll_int32,  intent(out) :: i_cell_y       !! not necessarily in [1, m2d%num_cells2], see comments above
    sll_real32, intent(out) :: offset_x, offset_y
    sll_real64              :: temp

    temp = (x - m2d%eta1_min) / m2d%delta_eta1
    i_cell_x  = 1 + int(floor(temp))
    offset_x = real(temp - real(i_cell_x - 1,f64), f32)

    temp = (y - m2d%eta2_min) / m2d%delta_eta2
    i_cell_y  = 1 + int(floor(temp))
    offset_y = real(temp - real(i_cell_y - 1,f64), f32)
    SLL_ASSERT(offset_x >= 0)
    SLL_ASSERT(offset_x <= 1 )
    SLL_ASSERT(offset_y >= 0)
    SLL_ASSERT(offset_y <= 1 )

    !! note: the (integer) index of the Poisson cell (within space computational domain) is then obtained with get_poisson_cell_index

  end subroutine global_to_cell_offset_extended

  !> <<cell_offset_to_global_extended>> performs the inverse transformation of global_to_cell_offset_extended above
  !!
  subroutine cell_offset_to_global_extended ( offset_x, offset_y, &
                                   i_cell_x, i_cell_y, m2d, &
                                   x, y )
  ! transforms sll_type of a particle (i_cell, dx, dy) into the standard
  ! particle position (x,y)
    sll_real64, intent(out)  :: x, y
    type(sll_cartesian_mesh_2d), intent(in) :: m2d
    sll_int32,  intent(in) :: i_cell_x, i_cell_y
    sll_real32, intent(in) :: offset_x, offset_y

    x = m2d%eta1_min + m2d%delta_eta1*( real(offset_x + i_cell_x-1, f64) )
    y = m2d%eta2_min + m2d%delta_eta2*( real(offset_y + i_cell_y-1, f64) )

  end subroutine cell_offset_to_global_extended


  subroutine get_poisson_cell_index( m2d, i_cell_x, i_cell_y, i_cell )
    type(sll_cartesian_mesh_2d), intent(in) :: m2d
    sll_int32,  intent(in) :: i_cell_x       !! not necessarily in [1, m2d%num_cells1]
    sll_int32,  intent(in) :: i_cell_y      !! not necessarily in [1, m2d%num_cells2]
    sll_int32,  intent(out) :: i_cell        !! in [1, m2d%num_cells1 * m2d%num_cells2]

    i_cell = 1 + modulo(i_cell_x-1,  m2d%num_cells1) + modulo(i_cell_y-1,  m2d%num_cells2) * m2d%num_cells1

    SLL_ASSERT( i_cell >= 1)
    SLL_ASSERT( i_cell <= m2d%num_cells1 * m2d%num_cells2 )

  end subroutine get_poisson_cell_index


end module sll_m_representation_conversion
