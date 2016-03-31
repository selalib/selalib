  !> @ingroup particle_methods

  !> @brief Unit test for particles of type \ref sll_m_simple_pic_4d_particle::sll_t_simple_pic_4d_particle

  ! [[file:../sll_m_simple_pic_4d_particle.F90::sll_simple_pic_4d_particle]]

program simple_pic_4d_group_tester

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_memory.h"
#include "sll_working_precision.h"
#include "sll_assert.h"

  use sll_m_cartesian_meshes, only: &
    sll_f_new_cartesian_mesh_2d, &
    sll_t_cartesian_mesh_2d

  use sll_m_simple_pic_4d_group, only: &
    sll_t_simple_pic_4d_group, &
    sll_f_simple_pic_4d_group_new

  implicit none
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

! [[selalib:src/mesh/meshes/sll_m_cartesian_meshes.F90::sll_cartesian_mesh_2d]]
type(sll_t_cartesian_mesh_2d),    pointer     :: mesh_2d
type(sll_t_simple_pic_4d_group),  pointer     :: simple_pic_particle_group

sll_int32 :: particle_group_id
sll_int32 :: NC_X
sll_int32 :: NC_Y
sll_int32 :: number_particles
sll_real64 :: SPECIES_CHARGE
sll_real64 :: SPECIES_MASS
sll_real64 :: XMIN, XMAX
sll_real64 :: YMIN, YMAX

logical :: DOMAIN_IS_X_PERIODIC
logical :: DOMAIN_IS_Y_PERIODIC

particle_group_id = 0
SPECIES_CHARGE = 1.0_f64
SPECIES_MASS = 1.0_f64
DOMAIN_IS_X_PERIODIC = .true.
DOMAIN_IS_Y_PERIODIC = .true.

NC_X = 10
NC_Y = 10

number_particles = 10000

XMIN = 0.0_f64
XMAX = 1.0_f64

YMIN = 0.0_f64
YMAX = 1.0_f64

mesh_2d =>  sll_f_new_cartesian_mesh_2d( NC_X, NC_Y, XMIN, XMAX, YMIN, YMAX )

simple_pic_particle_group => sll_f_simple_pic_4d_group_new(       &
    SPECIES_CHARGE,                                             &
    SPECIES_MASS,                                               &
    particle_group_id,                                          &
    DOMAIN_IS_X_PERIODIC,                                       &
    DOMAIN_IS_Y_PERIODIC,                                       &
    number_particles,                                           &
    mesh_2d )

  print *, 'PASSED'

end program simple_pic_4d_group_tester
