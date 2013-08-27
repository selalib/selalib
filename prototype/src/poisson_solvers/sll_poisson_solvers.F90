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

!> @namespace sll_poisson_solvers
!> @author Pierre Navaro
!> @brief 
!> Library to solve Poisson equation in 2D and 3D
!>
!> - Add  :
!> \code
!> #include "sll_poisson_solvers.h"
!> \endcode

module sll_poisson_solvers

#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"
#include "sll_constants.h"

implicit none

type, public  :: poisson_2d
  sll_int32   :: nc_x, nc_y
  sll_real64  :: dx, dy
  sll_real64  :: x_min
  sll_real64  :: x_max
  sll_real64  :: y_min
  sll_real64  :: y_max
end type poisson_2d

type, public  :: poisson_3d
  sll_int32   :: nc_x, nc_y, nc_z
  sll_real64  :: dx, dy, dz
  sll_real64  :: x_min
  sll_real64  :: x_max
  sll_real64  :: y_min
  sll_real64  :: y_max
  sll_real64  :: z_min
  sll_real64  :: z_max
end type poisson_3d

end module sll_poisson_solvers
