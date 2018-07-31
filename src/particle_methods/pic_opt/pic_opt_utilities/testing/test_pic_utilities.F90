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
!$ use omp_lib

implicit none

type(sll_t_cartesian_mesh_2d)            :: mesh_2d
type(sll_t_particle_group_2d), pointer   :: p_group_2d
type(sll_t_particle_group_4d), pointer   :: p_group_4d
type(sll_t_charge_accumulator_2d_ptr),    pointer :: q_accumulator(:)
type(sll_t_charge_accumulator_2d_cs_ptr), pointer :: q_accumulator_cs(:)
type(sll_t_charge_accumulator_2d_ptr),    pointer :: q_gc_accumulator(:)
type(sll_t_charge_accumulator_2d_cs_ptr), pointer :: q_gc_accumulator_cs(:)
sll_real64 :: QoverM = 1.0_f64
sll_real32 :: s
sll_int32  :: i, tid

sll_int32  :: nthreads, thread_id

nthreads  = 1
thread_id = 1

!$omp parallel
!$ nthreads =  OMP_GET_NUM_THREADS()
!$omp end parallel


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

allocate(q_accumulator(nthreads))
!$omp parallel PRIVATE(thread_id)
!$    thread_id = OMP_GET_THREAD_NUM()+1
allocate(q_accumulator(thread_id)%q)
call sll_s_charge_accumulator_2d_init(q_accumulator(thread_id)%q, mesh_2d )


!$omp end parallel

call sll_s_first_charge_accumulation_2d( p_group_4d, q_accumulator )

s = 0.0_f32
do tid = 1, nthreads
 s = s + real(sum(q_accumulator(tid)%q%q_acc(:)%q_sw),f32)
end do

print*, s
if ( s /= 10.0_f32 ) stop 'FAILED'

allocate(q_accumulator_cs(nthreads))
!$omp parallel PRIVATE(thread_id)
!$    thread_id = OMP_GET_THREAD_NUM()+1
allocate(q_accumulator_cs(thread_id)%q)
call sll_s_charge_accumulator_2d_cs_init(q_accumulator_cs(thread_id)%q, mesh_2d )
!$omp end parallel
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

allocate(q_gc_accumulator(nthreads))
!$omp parallel PRIVATE(thread_id)
!$    thread_id = OMP_GET_THREAD_NUM()+1
allocate(q_gc_accumulator(thread_id)%q)
call sll_s_charge_accumulator_2d_init(q_gc_accumulator(thread_id)%q, mesh_2d )
!$omp end parallel
call sll_s_first_gc_charge_accumulation_2d( p_group_2d, q_gc_accumulator )

allocate(q_gc_accumulator_cs(nthreads))
!$omp parallel PRIVATE(thread_id)
!$    thread_id = OMP_GET_THREAD_NUM()+1
allocate(q_gc_accumulator_cs(thread_id)%q)
call sll_s_charge_accumulator_2d_cs_init(q_gc_accumulator_cs(thread_id)%q, mesh_2d )
!$omp end parallel
call sll_s_first_gc_charge_accumulation_2d_cs( p_group_2d, q_gc_accumulator_cs )

s = 0.0_f32
do tid = 1, nthreads
 s = s + real(sum(  q_accumulator_cs(tid)%q%q_acc(:)%q_im1j   &
                  + q_accumulator_cs(tid)%q%q_acc(:)%q_ij     &
                  + q_accumulator_cs(tid)%q%q_acc(:)%q_ip1j   &
                  + q_accumulator_cs(tid)%q%q_acc(:)%q_im1jm1 &
                  + q_accumulator_cs(tid)%q%q_acc(:)%q_ijm1   &
                  + q_accumulator_cs(tid)%q%q_acc(:)%q_ip1jm1 &
                  + q_accumulator_cs(tid)%q%q_acc(:)%q_im1jp1 &
                  + q_accumulator_cs(tid)%q%q_acc(:)%q_ijp1   &
                  + q_accumulator_cs(tid)%q%q_acc(:)%q_ip1jp1), f32)
end do

print*, s
if ( abs(s - 10.0_f32) > 1e-6 ) stop 'FAILED'

print*, "PASSED"

end program test_pic_utilities
