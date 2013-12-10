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
!> @brief
!> Library to solve Poisson equation in 2D and 3D
!> 
!> - Modules available
!>   + sll_fishpack
!>   + sll_mudpack_cartesian
!>   + sll_mudpack_colella
!>   + sll_mudpack_polar
!>   + sll_poisson_1d_periodic
!>   + sll_poisson_2d_fem
!>   + sll_poisson_2d_periodic_fem
!>   + sll_poisson_2d_periodic
!>   + sll_poisson_2d_polar
!>   + sll_poisson_3d_periodic_seq
!>
!> - Parallel solvers
!>   + sll_poisson_2d_periodic_cartesian_par
!>   + sll_poisson_polar_parallel
!>
!> - Import module with  :
!> \code
!> #include "sll_poisson_solvers.h"
!> \endcode
module sll_poisson_solvers

#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"
#include "sll_constants.h"

implicit none

!> Poisson solver base class on 2d cartesian mesh
!> with periodic boundary conditions
!> We use FFT to solve Potential or Electric fields
type, public  :: poisson_2d
  sll_int32   :: nc_x       !< number of cells along x
  sll_int32   :: nc_y       !< number of cells along y
  sll_real64  :: dx         !< step size along x
  sll_real64  :: dy         !< step size along y
  sll_real64  :: x_min      !< left corner of x dimension
  sll_real64  :: x_max      !< right corner of y dimension
  sll_real64  :: y_min      !< left corner of y dimension
  sll_real64  :: y_max      !< right corner of y dimension
end type poisson_2d

!> Poisson solver base class on 3d cartesian mesh
!> with periodic boundary conditions
!> We use FFT to solve Potential or Electric fields
type, public  :: poisson_3d
  sll_int32   :: nc_x       !< number of cells along x
  sll_int32   :: nc_y       !< number of cells along y
  sll_int32   :: nc_z       !< number of cells along z
  sll_real64  :: dx         !< step size along x
  sll_real64  :: dy         !< step size along y
  sll_real64  :: dz         !< step size along z
  sll_real64  :: x_min      !< left corner of x dimension
  sll_real64  :: x_max      !< right corner of y dimension
  sll_real64  :: y_min      !< left corner of y dimension
  sll_real64  :: y_max      !< right corner of y dimension
  sll_real64  :: z_min      !< left corner of z dimension
  sll_real64  :: z_max      !< right corner of z dimension
end type poisson_3d

end module sll_poisson_solvers
