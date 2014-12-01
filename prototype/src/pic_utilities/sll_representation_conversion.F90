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

  use sll_logical_meshes

  implicit none
  
contains
  subroutine compute_cell_and_offset( &
                                     x, &
                                     xmin, &
                                     rdelta_x, &! rdx =1/delta_x
                                     i_cell, &
                                     offset )
    sll_real64, intent(in)  ::  x, xmin, rdelta_x
    sll_int32,  intent(out) ::  i_cell
    sll_real32, intent(out) ::  offset!sll_real64, intent(out) ::  offset!     
    sll_real64 :: temp
!!$    sll_int64  :: ioio

    temp = (x - xmin)*rdelta_x
!!$    ioio = int(temp)
!!$    i_cell = ioio
    i_cell  = int(temp)
    offset = temp - real(i_cell,f64)
!!$    i_cell = int( (x - xmin)/dx )
!!$    offset = mod( x - xmin, dx )
!!$    offset = offset/dx! the cell for a charge accumulator is [0,1) x [0,1)
!                           and not [0,delta_x) x [0,delta_y)  !!!
  end subroutine compute_cell_and_offset

  subroutine global_to_cell_offset( x, y, &
                      m2d,   &
                      icell, &
                      offset_x, offset_y )
! transforms a standard particle position (x,y) in our type (icell, dx, dy)

    sll_real64, intent(in)  :: x, y
    type(sll_logical_mesh_2d), intent(in) :: m2d
    sll_int32,  intent(out) :: icell
    sll_real32, intent(out) :: offset_x, offset_y!    sll_real64, intent(out) :: offset_x, offset_y!
    sll_int32               :: icell_x, icell_y
    sll_real64              :: rdeltax, rdeltay

    rdeltax = 1._f64/m2d%delta_eta1
    call compute_cell_and_offset( x, m2d%eta1_min, &
                                  rdeltax, icell_x, &
                                  offset_x )

!    SLL_ASSERT( (icell_x>=0).and.(icell_x<m2d%num_cells1) )
    if ( (icell_x<0).or.(icell_x.ge.m2d%num_cells1) ) &
         print*,'ERROR: bad icell_x', icell_x

!!$    if ( (offset_x<0).or.(offset_x.ge.1) ) then
!!$       print*, 'ERROR: bad offset_x', offset_x
!!$       print*, m2d%eta1_min, '<', x, '<',m2d%eta1_max, 'Xcell is',icell_x
!!$       print*, 'Delta_x is', m2d%delta_eta1, 'ic*Delta_x is', icell_x*m2d%delta_eta1
!!$       print*, 'THUS, ic*Delta_x + Delta_x',  icell_x*m2d%delta_eta1 + m2d%delta_eta1
!!$    endif
    rdeltay = 1._f64/m2d%delta_eta2
    call compute_cell_and_offset( y, m2d%eta2_min, &
                                  rdeltay, icell_y, &
                                  offset_y )

    if ( (icell_y<0).or.(icell_y.ge.m2d%num_cells2) ) &
         print*,'ERROR: bad icell_y', icell_y

!!$    if ( (offset_y<0).or.(offset_y.ge.1) ) then
!!$       print*, 'ERROR: bad offset_y', offset_y
!!$       print*, m2d%eta2_min, '<', y, '<',m2d%eta2_max, 'Ycell is',icell_y
!!$       print*, 'Delta_y is', m2d%delta_eta2, 'ic*Delta_y is', icell_y*m2d%delta_eta2
!!$       print*, 'THUS, ic*Delta_y + Delta_y',  icell_y*m2d%delta_eta2 + m2d%delta_eta2
!!$    endif
    icell = icell_x + 1 + icell_y * m2d%num_cells1

    if ( (icell<1).or.(icell>(m2d%num_cells1*m2d%num_cells2)) ) &
         print*,'ERROR: bad icell', icell
  end subroutine global_to_cell_offset
  
  subroutine cell_offset_to_global ( offset_x, offset_y, &
                                   icell, m2d, &
                                   x, y )
! transforms sll_type of a particle (icell, dx, dy) into the standard
! particle position (x,y)
    sll_real64, intent(out)  :: x, y
    type(sll_logical_mesh_2d), intent(in) :: m2d
    sll_int32,  intent(in) :: icell
    sll_real32, intent(in) :: offset_x, offset_y!sll_real64, intent(in) :: offset_x, offset_y!

    x = m2d%delta_eta1*( offset_x + real( mod(icell-1,m2d%num_cells1), f64) )
    y = m2d%delta_eta2*( offset_y + real( int( (icell-1)/m2d%num_cells1 ), f64))

!    if (x.le.m2d%eta1_min .or. x.ge.m2d%eta1_max) print*, 'STOP', x

!    if (y.le.m2d%eta2_min .or. y.ge.m2d%eta2_max) print*, 'STOP', y

  end subroutine cell_offset_to_global
	
end module sll_representation_conversion_module
