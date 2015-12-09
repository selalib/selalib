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

module sll_m_remapped_pic_utilities

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_memory.h"
#include "sll_working_precision.h"

  use sll_m_cartesian_meshes, only: &
    sll_cartesian_mesh_2d

  implicit none

  public :: &
    apply_periodic_bc_on_cartesian_mesh_2d, &
    x_is_in_domain_2d

  private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

contains



  ! todo: put this in the right module (with the meshes?)
  ! tells whether the given point is in the given domain, with boolean arguments for the domain periodicity
  ! (taken from previous function in_bounds_periodic)
  function x_is_in_domain_2d( x, y, mesh, x_periodic, y_periodic ) result(res)

!    use sll_m_cartesian_meshes
    sll_real64,                     intent( in )            :: x, y
    type(sll_cartesian_mesh_2d),    intent( in ), pointer   :: mesh
    logical,                        intent( in )            :: x_periodic
    logical,                        intent( in )            :: y_periodic
    logical     :: res

    res = ( x >= mesh%eta1_min )                                                                                    &
          .and.                                                                                                     &
          ( ( x < mesh%eta1_max .and. x_periodic ) .or. ( x <= mesh%eta1_max .and. .not. x_periodic ) )             &
          .and.                                                                                                     &
          ( y >= mesh%eta2_min )                                                                                    &
          .and.                                                                                                     &
          ( ( y < mesh%eta2_max .and. y_periodic ) .or. ( y <= mesh%eta2_max .and. .not. y_periodic) )

  end function x_is_in_domain_2d

  ! <<apply_periodic_bc_on_cartesian_mesh_2d>>

  ! todo: put this in the right module (with the meshes?)
  subroutine apply_periodic_bc_on_cartesian_mesh_2d( mesh, x, y )

!    use sll_m_cartesian_meshes
    ! [[file:../working_precision/sll_m_working_precision.h]]
!    use sll_m_working_precision

    type(sll_cartesian_mesh_2d), pointer :: mesh
    sll_real64, intent(inout) :: x
    sll_real64, intent(inout) :: y

    x = mesh%eta1_min + modulo(x - mesh%eta1_min, mesh%eta1_max - mesh%eta1_min)
    y = mesh%eta2_min + modulo(y - mesh%eta2_min, mesh%eta2_max - mesh%eta2_min)
  end subroutine apply_periodic_bc_on_cartesian_mesh_2d

end module sll_m_remapped_pic_utilities
