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

!> @file
!> @brief Unit test for sll_mlt_pic_2d data structures

! Goal of this test: create then delete one mlt_pic particle group. This will in turn call all the defined entities that
! we want to check. If it works, then all the data type names are coherent. Tests data structures from
! [[file:mlt_pic_2d_particle.F90]] and [[file:mlt_pic_2d_group.F90]]. To compile check
! [[file:CMakeLists.txt::unit_test_mlt_pic_2d]]

program mlt_pic_2d_unit_test

  ! [[file:../assert/sll_assert.h]]
#include "sll_assert.h"

  ! [[file:../memory/sll_memory.h]]
#include "sll_memory.h"

  use sll_constants ! [[file:../constants/sll_constants.F90]] for [[sll_pi]]
  use sll_cartesian_meshes
  use sll_mlt_pic_2d_group_module ! [[file:mlt_pic_2d_group.F90::sll_mlt_pic_2d_group_module]]

  type(sll_mlt_pic_2d_group), pointer :: p_group

  ! values inspired from [[file:../simulation/simulation_4d_vp_lt_pic_cartesian.F90]]

#define KX   0.5_f64
#define NC_X 256_i32
#define XMIN 0._f64
#define XMAX (2._f64*sll_pi/KX)

  ! new mesh built with [[file:../meshes/sll_cartesian_meshes.F90::new_cartesian_mesh_1d]]
  
  type(sll_cartesian_mesh_1d), pointer :: m
  m => new_cartesian_mesh_1d(NC_X,XMIN,XMAX)

#define NUM_PARTICLES 1000000_i32
#define GUARD_SIZE    1000000_i32
#define PARTICLE_ARRAY_SIZE 2050000_i32
#define QOVERM 1._f64
#define REMAP_GRID_VX_MIN -6._f64
#define REMAP_GRID_VX_MAX  6._f64

  ! [[file:mlt_pic_2d_group.F90::sll_mlt_pic_2d_group_new]]

  p_group => sll_mlt_pic_2d_group_new(&
       1,&
       REMAP_GRID_VX_MIN,&
       REMAP_GRID_VX_MAX,&
       NUM_PARTICLES,&
       GUARD_SIZE,&
       QOVERM,&

       ! domain_is_x_periodic
       .true.,&

       m,&

       ! min_level defined in [[file:~/mcp/maltpic/algos-maltpic.tex::j0]]. Cannot be 0 since it is used as an index
       ! in array [[file:mlt_pic_2d_group.F90::target_values]] and this is Fortran
       1,&

       ! max_level defined in [[file:~/mcp/maltpic/algos-maltpic.tex::jmax]]
       5)

  ! cf [[file:../assert/sll_assert.h::SLL_ASSERT]]

  SLL_ASSERT(associated(p_group))
  SLL_ASSERT(p_group%min_level == 1)
  SLL_ASSERT(p_group%max_level == 5)

  ! calls [[file:mlt_pic_2d_group.F90::sll_mlt_pic_2d_group_delete]]

  call sll_mlt_pic_2d_group_delete(p_group)
  call sll_delete(m)

  print *, 'PASSED'

end program mlt_pic_2d_unit_test
