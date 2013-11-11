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
!> @author
!> Michel Mehrenberger (mehrenbe@math.unistra.fr)
!> Nicolas Crouseilles
!> Edwin Chacon Golcher (for checking code writing/interfaces/reusability)
!> @brief 
!> Simulation class to solve 4D vlasov poisson equation
!> parallelization follows simulation_4d_qns_general.F90
!> differs from simulation_4d_vlasov_poisson_cartesian.F90
!> as Poisson solvers is not intended to be parallel for the moment
!> should include new advectors and input files with explicit names
!> implementation and test of high order splitting is today's application goal  


module sll_simulation_4d_vlasov_parallel_poisson_sequential_cartesian

#include "sll_working_precision.h"
#include "sll_assert.h"
#include "sll_memory.h"
#include "sll_field_2d.h"
#include "sll_utilities.h"
#include "sll_poisson_solvers.h"
  use sll_collective
  use sll_remapper
  use sll_buffer_loader_utilities_module
  use sll_constants
  !for mesh
  use sll_logical_meshes
  
  use sll_gnuplot_parallel
  use sll_coordinate_transformation_2d_base_module
  use sll_module_coordinate_transformations_2d
  use sll_common_coordinate_transformations
  use sll_common_array_initializers_module
  use sll_parallel_array_initializer_module
  
  use sll_module_advection_1d_periodic

  use sll_poisson_1d_periodic
  
  use sll_fft

  use sll_simulation_base
  implicit none

end module sll_simulation_4d_vlasov_parallel_poisson_sequential_cartesian
