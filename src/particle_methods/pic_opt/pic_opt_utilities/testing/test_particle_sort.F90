program particle_sorter
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_memory.h"
#include "sll_working_precision.h"
#include "sll_particle_representation.h"

   use sll_m_cartesian_meshes, only: &
      sll_f_new_cartesian_mesh_2d, &
      sll_t_cartesian_mesh_2d

   use sll_m_constants, only: &
      sll_p_pi

   use sll_m_hammersley, only: &
      sll_f_suite_hamm

   use sll_m_particle_group_4d, only: &
      sll_f_new_particle_4d_group, &
      sll_t_particle_group_4d

   use sll_m_particle_sort, only: &
      sll_f_new_particle_sorter_2d, &
      sll_t_particle_sorter_2d, &
      sll_s_sort_particles_2d

   implicit none
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#define THERM_SPEED 1._f64
#define num_particles 100000_i32
#define KX    0.5_f64
#define XMIN  0._f64
#define XMAX (2._f64*sll_p_pi/KX)
#define NC_X  256_i32
#define YMIN  0._f64
#define YMAX  1._f64
#define NC_Y  64_i32
#define ALPHA  0.1_f64
#define QoverM 1._f64

   type(sll_t_particle_group_4d), pointer  :: group
   type(sll_t_cartesian_mesh_2d), pointer    :: m
   type(sll_t_particle_sorter_2d), pointer :: sorter
!  sll_int32  :: i
   sll_real64 :: x, y, vx, vy, nu, xmin, ymin, rdx, rdy
   sll_int32  :: j
   sll_real32 :: weight, off_x, off_y
   sll_int32  :: ncx, ic_x, ic_y
   sll_real64 :: tmp1, tmp2

   m => sll_f_new_cartesian_mesh_2d(NC_X, NC_Y, XMIN, XMAX, YMIN, YMAX)

   group => sll_f_new_particle_4d_group(int(num_particles, i32), &
                                        int(num_particles, i32), int(num_particles/4, i32), QoverM, m)

   sorter => sll_f_new_particle_sorter_2d(m)
! the arguments to new_particle_group should be 32bit ints... change

   weight = real(m%eta1_max - m%eta1_min, f32)* &
            real(m%eta2_max - m%eta2_min, f32)/real(num_particles, f32)
   rdx = 1._f64/m%delta_eta1
   rdy = 1._f64/m%delta_eta2
   xmin = m%eta1_min
   ymin = m%eta2_min
   ncx = m%num_cells1

   open (90, file='initialparticles.dat')
   j = 1
   do while (j <= num_particles)
      call random_number(x)
      x = (m%eta1_max - m%eta1_min)*x + m%eta1_min
      call random_number(y)
      y = 2._f64*y! (2._f64*alpha)*y + 1._f64 - alpha
      if (eval_landau(ALPHA, KX, x) >= y) then
         y = (m%eta2_max - m%eta2_min)*sll_f_suite_hamm(j, 3) + m%eta2_min
!
         nu = THERM_SPEED*sqrt(-2.0_f64*log(1.0_f64 - &
                                            (real(j, f64) - 0.5_f64)/real(num_particles, f64)))
         vx = nu*cos(sll_f_suite_hamm(j, 2)*2.0_f64*sll_p_pi)
         vy = nu*sin(sll_f_suite_hamm(j, 2)*2.0_f64*sll_p_pi)
         write (90, *) x, y, vx, vy
!        call set_group_particle_values( group, j, x, y, vx, vy, weight)
         SET_PARTICLE_VALUES(group%p_list(j), x, y, vx, vy, weight, xmin, ymin, ncx, ic_x, ic_y, off_x, off_y, rdx, rdy, tmp1, tmp2)
         j = j + 1
      end if
   end do; close (90)

   call sll_s_sort_particles_2d(sorter, group)

   do j = 1, num_particles - 1
      if (group%p_list(j)%ic > group%p_list(j + 1)%ic) then
         print *, 'BAD order of ic:', 'j=', j, 'j+1=', j + 1
         stop
      end if
   end do

!  call sll_o_delete (sorter)
   print *, 'PASSED'

contains

   function eval_landau(alp, kx, x)
      sll_real64 :: alp, kx, x
      sll_real64 :: eval_landau
      eval_landau = 1._f64 + alp*cos(kx*x)
   end function eval_landau

end program particle_sorter
