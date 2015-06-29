!> @ingroup particle_methods
!> @brief Unit test for particles of type \ref sll_simple_pic_4d_particle
program simple_pic_4d_group_tester

#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"

  use sll_simple_pic_4d_group_module
  use sll_cartesian_meshes

  implicit none

type(sll_cartesian_mesh_2d),    pointer     :: mesh_2d ! [[selalib:src/meshes/sll_cartesian_meshes.F90::sll_cartesian_mesh_2d]]
type(sll_simple_pic_4d_group),  pointer     :: simple_pic_particle_group

sll_int32 :: particle_group_id
sll_int32 :: NC_X
sll_int32 :: NC_Y
sll_real64 :: SPECIES_CHARGE
sll_real64 :: SPECIES_MASS
sll_real64 :: XMIN, XMAX
sll_real64 :: YMIN, YMAX

logical :: DOMAIN_IS_X_PERIODIC
logical :: DOMAIN_IS_Y_PERIODIC

particle_group_id = 0
SPECIES_CHARGE = 1.
SPECIES_MASS = 1.
DOMAIN_IS_X_PERIODIC = .true.
DOMAIN_IS_Y_PERIODIC = .true.

NC_X = 10
NC_Y = 10

XMIN = 0.
XMAX = 1.

YMIN = 0.
YMAX = 1.

mesh_2d =>  new_cartesian_mesh_2d( NC_X, NC_Y, XMIN, XMAX, YMIN, YMAX )

simple_pic_particle_group => sll_simple_pic_4d_group_new(       &
    SPECIES_CHARGE,                                             &
    SPECIES_MASS,                                               &
    particle_group_id,                                          &
    DOMAIN_IS_X_PERIODIC,                                       &
    DOMAIN_IS_Y_PERIODIC,                                       &
    mesh_2d )

  print *, 'PASSED'

end program simple_pic_4d_group_tester
