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
!> @author We will see
!> @brief 
!> Simulation class to solve slab drift kinetic equation in polar coordinates
!> (3d space (r,\theta,z) 1d velocity (v))
!> translation of vp4d_dk in simulation class
!> @details
!> Example of use in test program 
!> 
!> \code
!>
!>  use sll_simulation_4d_drift_kinetic_polar
!>  type(sll_simulation_4d_vp_polar)    :: simulation
!>  call simulation%init_from_file(trim(filename))
!>  call simulation%run()
!>  call delete(simulation)
!> \endcode


module sll_simulation_4d_drift_kinetic_polar
#include "sll_working_precision.h"
#include "sll_assert.h"
#include "sll_memory.h"
#include "sll_field_2d.h"
#include "sll_utilities.h"
  use sll_collective
  use sll_remapper
  use sll_constants
  use sll_cubic_spline_interpolator_1d
  use sll_test_4d_initializer
  use sll_poisson_2d_periodic_cartesian_par
  use sll_cubic_spline_interpolator_1d
  use sll_simulation_base
  implicit none

  
  type, extends(sll_simulation_base_class) :: &
    sll_simulation_4d_dk_polar
   contains
     procedure, pass(sim) :: run => run_dk4d_polar
     procedure, pass(sim) :: init_from_file => init_dk4d_polar
  end type sll_simulation_4d_dk_polar

  interface delete
     module procedure delete_dk4d_polar
  end interface delete

contains

  subroutine init_dk4d_polar( sim, filename )
    intrinsic :: trim
    class(sll_simulation_4d_dk_polar), intent(inout) :: sim
    character(len=*), intent(in)                                :: filename
  end subroutine init_dk4d_polar

  subroutine run_dk4d_polar(sim)
    class(sll_simulation_4d_dk_polar), intent(inout) :: sim
  end subroutine run_dk4d_polar


  subroutine delete_dk4d_polar( sim )
    class(sll_simulation_4d_dk_polar) :: sim
    sll_int32 :: ierr
  end subroutine delete_dk4d_polar

  

end module sll_simulation_4d_drift_kinetic_polar



