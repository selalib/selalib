program test_pic_poisson_2d
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_memory.h"
#include "sll_working_precision.h"

  use sll_m_constants, only: &
    sll_p_twopi

  use sll_m_particle_group_1d2v, only: &
    sll_s_new_particle_group_1d2v, &
    sll_t_particle_group_1d2v

  use sll_m_particle_group_base, only: &
    sll_c_particle_group_base

  use sll_m_particle_initializer

  use sll_m_pic_poisson_base, only : &
    sll_c_pic_poisson

  use  sll_m_pic_poisson_2d, only: &
    sll_s_new_pic_poisson_2d, &
    sll_t_pic_poisson_2d

  use sll_m_poisson_2d_periodic, only: &
    sll_f_new_poisson_2d_periodic

  use sll_m_poisson_2d_base, only: &
    sll_c_poisson_2d_base

  use sll_m_kernel_smoother_base, only: &
    sll_p_collocation, &
    sll_c_kernel_smoother

  use sll_m_kernel_smoother_spline_2d, only: &
    sll_t_kernel_smoother_spline_2d, &
    sll_s_new_kernel_smoother_spline_2d_ptr

  implicit none
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  sll_int32 :: num_cells(2)
  sll_int32 :: degree_smoother
  sll_real64 :: domain(2,2)
  sll_int32  :: n_particles

  ! Poisson solver
  class(sll_c_poisson_2d_base), pointer :: pic_poisson 
  
  ! PIC Poisson solver
  class(sll_c_pic_poisson), pointer :: poisson

  ! Abstract kernel smoother
  class(sll_c_kernel_smoother), pointer :: kernel_smoother


  num_cells = [32, 32]
  domain(1,:) = [0.0_f64, sll_p_twopi]
  domain(2,:) = [0.0_f64, sll_p_twopi]
  degree_smoother = 3
  n_particles = 1000

  ! Just test allocation and deallocation since functionality is tested in the building blocks.

  ! Initialize the field solver
  pic_poisson => sll_f_new_poisson_2d_periodic( &
       domain(1,1), domain(1,2), num_cells(1), &
       domain(2,1), domain(2,2), num_cells(2) )
  
  ! Initialize the kernel smoother
  call sll_s_new_kernel_smoother_spline_2d_ptr(kernel_smoother, &
       domain, num_cells, n_particles, &
       degree_smoother, sll_p_collocation)
  
  ! Initialize the PIC field solver
  call sll_s_new_pic_poisson_2d( poisson, &
       num_cells, &
       pic_poisson, kernel_smoother)

  call poisson%free()
  call kernel_smoother%free()
  call pic_poisson%free()
  deallocate(poisson)
  deallocate(kernel_smoother)
  deallocate(pic_poisson)

  ! Test passed if init and free successfull.
  print*, 'PASSED'


end program test_pic_poisson_2d
