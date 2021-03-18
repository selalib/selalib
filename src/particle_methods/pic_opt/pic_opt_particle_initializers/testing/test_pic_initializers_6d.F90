program test_pic_initializers_6d
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_memory.h"
#include "sll_working_precision.h"

   use sll_m_cartesian_meshes
   use sll_m_constants, only: sll_p_pi
   use sll_m_particle_group_6d
   use sll_m_particle_initializers_6d

   implicit none
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#define THERM_SPEED 1._f64
#define NUM_PARTICLES 50000_i32
#define GUARD_SIZE 10000_i32
#define PARTICLE_ARRAY_SIZE 15000000_i32
#define ALPHA 0.5_f64
#define NC_X 128_i32
#define XMIN 0._f64
#define KX   0.5_f64
#define XMAX 2._f64*sll_p_pi/KX
#define NC_Y 32_i32
#define YMIN 0._f64
#define YMAX 1._f64
#define NC_Z 32_i32
#define ZMIN 0._f64
#define ZMAX 1._f64
#define QoverM 1._f64

#define GET_PARTICLE_POSITION(p,m,x,y,z) \
   x = m%eta1_min + m%delta_eta1*(real(p%dx, f64) + real(mod(p%ic - 1, m%num_cells1), f64)); \
   y = m%eta2_min + m%delta_eta2*(real(p%dy, f64) + real(int((p%ic - 1)/(m%num_cells1)), f64)); \
   z = m%eta3_min + m%delta_eta3*(real(p%dz, f64) + real(int((p%ic - 1)/(m%num_cells1*m%num_cells2)), f64))

   type(sll_t_particle_group_6d) :: group
   type(sll_t_cartesian_mesh_3d) :: mesh
   sll_real64 :: mean_ref_v(1:3), variance_ref_v(1:3)
   sll_real64 :: mean_ref_x, variance_ref_x
   sll_real64 :: len_x, len_v
   logical    :: passed
   integer, allocatable, dimension(:) :: iseed
   integer :: nseed

   call random_seed(size=nseed)
   allocate (iseed(nseed))
   iseed = 42
   call random_seed(put=iseed)

   passed = .true.

   call sll_s_cartesian_mesh_3d_init(mesh, NC_X, NC_Y, NC_Z, &
                                     XMIN, XMAX, YMIN, YMAX, ZMIN, ZMAX)

   call sll_s_particle_6d_group_init(group, &
                                     NUM_PARTICLES, &
                                     PARTICLE_ARRAY_SIZE, &
                                     GUARD_SIZE, QoverM, mesh)

   mean_ref_v = [0._f64, 0._f64, 0._f64]! the 2d Gaussian for velocity
   variance_ref_v = [1._f64, 1._f64, 1._f64]!
   mean_ref_x = sll_p_pi/KX! the Landau 1d perturbation for positions
   variance_ref_x = (sll_p_pi**2/3._f64 + 2._f64*ALPHA)/(KX**2)

! the confidence interval for the mean

   len_x = 3._f64*sqrt(variance_ref_x)/sqrt(real(NUM_PARTICLES, f64))
   len_v = 3._f64*sqrt(variance_ref_v(1))/sqrt(real(NUM_PARTICLES, f64))

   print *, 'the Random initialization for the Landau damping'
   call sll_s_initial_random_particles_6d(THERM_SPEED, &
                                          ALPHA, KX, mesh, NUM_PARTICLES, group)
   call test_mean_variance_xv(group, len_x, len_v)

   call sll_s_particle_6d_group_free(group)
   call sll_s_cartesian_mesh_3d_free(mesh)

   if (passed .eqv. .true.) then
      print *, "PASSED"
   else
      print *, "FAILED"
   end if

contains

   subroutine test_mean_variance_xv(part_group, tol_x, tol_v)

      type(sll_t_particle_group_6d), intent(in) :: part_group
      sll_real64, intent(in)                    :: tol_x
      sll_real64, intent(in)                    :: tol_v

      sll_real64 :: mean_v(1:3), variance_v(1:3)! estimators in velocity
      sll_real64 :: mean_x, variance_x!           estimators in position
      sll_int32  :: i
      sll_real64 :: out_v(1:3), out_x
      sll_real64 :: x, y, z

      mean_x = 0._f64
      mean_v = 0._f64
      do i = 1, NUM_PARTICLES
         GET_PARTICLE_POSITION(part_group%p_list(i), mesh, x, y, z)
         mean_x = mean_x + x
         mean_v(1) = mean_v(1) + part_group%p_list(i)%vx
         mean_v(2) = mean_v(2) + part_group%p_list(i)%vy
         mean_v(3) = mean_v(3) + part_group%p_list(i)%vz
      end do
      mean_v = mean_v/real(NUM_PARTICLES, f64)
      mean_x = mean_x/real(NUM_PARTICLES, f64)

      out_v = mean_v - mean_ref_v
      out_x = mean_x - mean_ref_x

      if ((abs(out_v(1)) > tol_v) .or. &
          (abs(out_v(2)) > tol_v) .or. &
          (abs(out_v(3)) > tol_v)) then
         stop 'Error in the expected value in velocity'
      end if
      if (abs(out_x) > tol_x) then
         stop 'Error in the expected value in position'
      end if

      variance_v = 0._f64
      variance_x = 0._f64
      do i = 1, NUM_PARTICLES
         variance_v(1) = variance_v(1) + (part_group%p_list(i)%vx - mean_v(1))**2
         variance_v(2) = variance_v(2) + (part_group%p_list(i)%vy - mean_v(2))**2
         variance_v(3) = variance_v(3) + (part_group%p_list(i)%vz - mean_v(3))**2
         GET_PARTICLE_POSITION(part_group%p_list(i), mesh, x, y, z)
         variance_x = variance_x + (x - mean_x)**2
      end do
      variance_v = variance_v/real(NUM_PARTICLES - 1, f64)
      variance_x = variance_x/real(NUM_PARTICLES - 1, f64)

      print *, 'the error of the variance approximation in X=', &
         variance_x - variance_ref_x
      print *, 'the error of the variance approximation in V=', &
         variance_v - variance_ref_v

   end subroutine test_mean_variance_xv

end program test_pic_initializers_6d
