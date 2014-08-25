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


module sll_particle_representations
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"

  use sll_logical_meshes

  implicit none

  type :: sll_particle_2d
     sll_int32  :: ic   ! cell index, linearly arranged
     sll_real32 :: dx!     sll_real64 :: dx
     sll_real32 :: dy!     sll_real64 :: dy!
     sll_real64 :: vx
     sll_real64 :: vy
     sll_real32 :: q
  end type sll_particle_2d

  type :: sll_particle_2d_guard
     type(sll_particle_2d), pointer :: p
  end type sll_particle_2d_guard

  type :: sll_particle_2d_guard_ptr
     type(sll_particle_2d_guard), dimension(:), pointer :: g_list
  end type sll_particle_2d_guard_ptr
   type :: sll_particle_1d1v
     sll_real64 :: x      ! This is not necessarily an offset
     sll_real64 :: vx
     sll_real64 :: w            !\gamma_k , weight
  end type sll_particle_1d1v


  type :: sll_particle_1d1v_idxx, extends(sll_particle_1d1v)
    sll_int32  :: icx !cell index
    sll_real32 :: dx
  end type  sll_particle_1d1v_idxx

  type :: sll_particle_1d1v_twoweight, extends(sll_particle_1d1v)
     sll_real64 :: wc      !c_k , constant weight
  endtype

  type :: sll_particle_1d1v_guard
     class(sll_particle_1d1v), pointer :: p

  end type sll_particle_1d1v_guard




!contains

!!$  subroutine initialize_particle_2d( &
!!$       p,  &
!!$       ic, &
!!$       dx, &
!!$       dy, &
!!$       vx, &
!!$       vy, &
!!$       q )
!!$
!!$    type(sll_particle_2d), intent(out) :: p
!!$    sll_int32,  intent(in) :: ic   ! cell index, linearly arranged
!!$    sll_real32, intent(in) :: dx
!!$    sll_real32, intent(in) :: dy
!!$    sll_real64, intent(in) :: vx
!!$    sll_real64, intent(in) :: vy
!!$    sll_real32, intent(in) :: q
!!$
!!$    p%ic = ic
!!$    p%dx = dx
!!$    p%dy = dy
!!$    p%vx = vx
!!$    p%vy = vy
!!$    p%q  = q
!!$  end subroutine initialize_particle_2d

!!$  subroutine compute_cell_and_offset( &
!!$                                     x, &
!!$                                     xmin, &
!!$                                     rdelta_x, &! rdx =1/delta_x
!!$                                     i_cell, &
!!$                                     offset )
!!$    sll_real64, intent(in)  ::  x, xmin, rdelta_x
!!$    sll_int32,  intent(out) ::  i_cell
!!$    sll_real32, intent(out) ::  offset
!!$    sll_real64 :: temp
!!$
!!$    temp = (x - xmin)*rdelta_x
!!$    i_cell  = int(temp)
!!$    offset = temp - real(i_cell,f32)
!!$! the cell for a charge accumulator is [0,1) x [0,1)
!!$!                  and not [0,delta_x) x [0,delta_y)  !
!!$  end subroutine compute_cell_and_offset
!!$
!!$
!!$  subroutine set_particle_values( particle,  &
!!$                                  x,  y,     &
!!$                                  vx, vy, q, &
!!$                                  mesh       )
!!$
!!$    type(sll_particle_2d), intent(out) :: particle
!!$    sll_real64, intent(in) :: x, y, vx, vy
!!$    sll_real32, intent(in) :: q
!!$    type(sll_logical_mesh_2d), pointer :: mesh
!!$    sll_real64             :: rdeltax, rdeltay
!!$    sll_int32              :: icell_x, icell_y
!!$    sll_int32              :: icell
!!$    sll_real32             :: offset_x, offset_y
!!$
!!$    rdeltax = 1._f64/mesh%delta_eta1
!!$    call compute_cell_and_offset( x, mesh%eta1_min, &
!!$                                  rdeltax, icell_x, &
!!$                                  offset_x )
!!$    rdeltay = 1._f64/mesh%delta_eta2
!!$    call compute_cell_and_offset( y, mesh%eta2_min, &
!!$                                  rdeltay, icell_y, &
!!$                                  offset_y )
!!$    icell = icell_x + 1 + icell_y * mesh%num_cells1
!!$
!!$    SLL_ASSERT( (icell_x>=0).and.(icell_x<mesh%num_cells1) )
!!$    SLL_ASSERT( (icell_y>=0).and.(icell_y<mesh%num_cells2) )
!!$    SLL_ASSERT( (icell>=1).and.(icell<=mesh%num_cells1*mesh%num_cells2) )
!!$
!!$    particle%ic = icell
!!$    particle%dx = offset_x
!!$    particle%dy = offset_y
!!$    particle%vx = vx
!!$    particle%vy = vy
!!$    particle%q  = q
!!$  end subroutine set_particle_values


end module sll_particle_representations
