program test_pic_utilities
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_assert.h"
#include "sll_memory.h"
#include "sll_working_precision.h"
#include "sll_accumulators.h"

use sll_m_cartesian_meshes
use sll_m_particle_group_2d
use sll_m_particle_group_4d
use sll_m_accumulators
use sll_m_pic_utilities

implicit none

type(sll_t_cartesian_mesh_2d)            :: mesh_2d
type(sll_t_particle_group_2d), pointer   :: p_group_2d
type(sll_t_particle_group_4d), pointer   :: p_group_4d
type(sll_t_charge_accumulator_2d_ptr),    pointer :: q_accumulator(:)
type(sll_t_charge_accumulator_2d_cs_ptr), pointer :: q_accumulator_cs(:)
type(sll_t_charge_accumulator_2d_ptr),    pointer :: q_gc_accumulator(:)
type(sll_t_charge_accumulator_2d_cs_ptr), pointer :: q_gc_accumulator_cs(:)
sll_real64 :: QoverM = 1.0_f64
sll_int32  :: i

#define NC_X 10
#define NC_Y 10
#define XMIN -1.0_f64
#define YMIN -1.0_f64
#define XMAX +1.0_f64
#define YMAX +1.0_f64

call sll_s_cartesian_mesh_2d_init( mesh_2d, NC_X, NC_Y, XMIN, XMAX, YMIN, YMAX )

#define NUM_PARTICLES 10
#define PARTICLE_ARRAY_SIZE 10
#define GUARD_SIZE 1


allocate(p_group_4d)
call sll_s_particle_4d_group_init( p_group_4d, &
         NUM_PARTICLES, &
         PARTICLE_ARRAY_SIZE, &
         GUARD_SIZE, &
         QoverM,     &
         mesh_2d )

do i = 1, p_group_4d%number_particles
  p_group_4d%p_list(i)%ic = i
  p_group_4d%p_list(i)%q  = 1.0_f32
end do

allocate(q_accumulator(1))
allocate(q_accumulator(1)%q)
call sll_s_charge_accumulator_2d_init(q_accumulator(1)%q, mesh_2d )
call sll_s_first_charge_accumulation_2d( p_group_4d, q_accumulator )
if (sum(q_accumulator(1)%q%q_acc(:)%q_sw) /= 10.0_f64) stop "FAILED"
if (sum(q_accumulator(1)%q%q_acc(:)%q_se) /= 00.0_f64) stop "FAILED"
if (sum(q_accumulator(1)%q%q_acc(:)%q_nw) /= 00.0_f64) stop "FAILED"
if (sum(q_accumulator(1)%q%q_acc(:)%q_ne) /= 00.0_f64) stop "FAILED"


allocate(q_accumulator_cs(1))
allocate(q_accumulator_cs(1)%q)
call sll_s_charge_accumulator_2d_cs_init(q_accumulator_cs(1)%q, mesh_2d )
call sll_s_first_charge_accumulation_2d_cs( p_group_4d, q_accumulator_cs )

allocate(p_group_2d)
call sll_s_particle_2d_group_init( p_group_2d, &
         NUM_PARTICLES, &
         PARTICLE_ARRAY_SIZE, &
         GUARD_SIZE, &
         QoverM,     &
         mesh_2d )

do i = 1, p_group_2d%number_particles
  p_group_2d%p_list(i)%ic = i
  p_group_2d%p_list(i)%q  = 1.0_f32
end do

allocate(q_gc_accumulator(1))
allocate(q_gc_accumulator(1)%q)
call sll_s_charge_accumulator_2d_init(q_gc_accumulator(1)%q, mesh_2d )
call sll_s_first_gc_charge_accumulation_2d( p_group_2d, q_gc_accumulator )
if (sum(q_gc_accumulator(1)%q%q_acc(:)%q_sw) /= 10.0_f64) stop "FAILED"
if (sum(q_gc_accumulator(1)%q%q_acc(:)%q_se) /= 00.0_f64) stop "FAILED"
if (sum(q_gc_accumulator(1)%q%q_acc(:)%q_nw) /= 00.0_f64) stop "FAILED"
if (sum(q_gc_accumulator(1)%q%q_acc(:)%q_ne) /= 00.0_f64) stop "FAILED"

allocate(q_gc_accumulator_cs(1))
allocate(q_gc_accumulator_cs(1)%q)
call sll_s_charge_accumulator_2d_cs_init(q_gc_accumulator_cs(1)%q, mesh_2d )
call sll_s_first_gc_charge_accumulation_2d_cs( p_group_2d, q_gc_accumulator_cs )


print*, "PASSED"

end program test_pic_utilities
