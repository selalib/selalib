
!***************************************************************************
!
! Selalib 2014     
! Module: unit_test_qn_solver_3d_polar_parallel_x1.F90
!
!> @brief 
!> test the quasi neutral solver 3D polar parallel in r
!> whhich is for the moment essentially a copy paste of sll_poisson_polar_parallel.F90
!> with adding of the third z direction
!> Start date: April 17, 2014
!> Last modification:: April 17, 2014
!   
!> @authors                    
!> Michel MEHRENBERGER (mehrebe@math.unistra.fr), 
!                                  
!***************************************************************************
program unit_test_qn_solver_3d_polar_parallel_x1
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"
  use sll_remapper
  use sll_constants
  use sll_collective
  use sll_qn_solver_3d_polar_parallel_x1_module
  use sll_boundary_condition_descriptors

implicit none

  print *,'#PASSED'

end program unit_test_qn_solver_3d_polar_parallel_x1
