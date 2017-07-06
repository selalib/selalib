program test_particle_group_6d

#include "sll_assert.h"
#include "sll_memory.h"
#include "sll_working_precision.h"
#include "sll_accumulators.h"

use sll_m_cartesian_meshes
use sll_m_particle_group_6d

implicit none

type(sll_t_cartesian_mesh_3d) :: mesh_3d
type(sll_t_particle_group_6d) :: p_group_6d
sll_real64                    :: QoverM = 1.0_f64

#define NC_X 10
#define NC_Y 10
#define NC_Z 10
#define XMIN -1.0_f64
#define YMIN -1.0_f64
#define ZMIN -1.0_f64
#define XMAX +1.0_f64
#define YMAX +1.0_f64
#define ZMAX +1.0_f64

call sll_s_cartesian_mesh_3d_init( mesh_3d, NC_X, NC_Y, NC_Z, &
  XMIN, XMAX, YMIN, YMAX, ZMIN, ZMAX )

#define NUM_PARTICLES 10
#define PARTICLE_ARRAY_SIZE 10
#define GUARD_SIZE 1

call sll_s_particle_6d_group_init( p_group_6d, &
         NUM_PARTICLES, &
         PARTICLE_ARRAY_SIZE, &
         GUARD_SIZE, &
         QoverM,     &
         mesh_3d )

call sll_s_particle_6d_group_free( p_group_6d)
print*, "PASSED"

end program test_particle_group_6d
