! Main test program for
! [[file:simulation_4d_vp_generic_pic_cartesian.F90::sll_m_sim_4d_vp_generic_pic_cartesian]]    [link to update]

!> @ingroup particle_methods

!> @author MCP ALH

!> @brief simulation program for the pic_lbfr method

! Program listed by Doxygen as part of the particle_methods in
! [[selalib:doc/build/html/doxygen/html/group__particle__methods.html]]

program sim_pic_vp_2d2v_cart_lbfr

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_working_precision.h"

  use sll_m_collective, only: &
    sll_s_boot_collective, &
    sll_f_get_collective_rank, &
    sll_f_get_collective_size, &
    sll_s_halt_collective, &
    sll_v_world_collective

  use sll_m_sim_pic_vp_2d2v_cart_lbfr, only: &
    sll_t_sim_pic_vp_2d2v_cart_lbfr, &
    sll_p_lbfr_particles

  use sll_m_timer, only: &
    sll_s_set_time_mark, &
    sll_f_time_elapsed_since, &
    sll_t_time_mark

  implicit none
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  type(sll_t_sim_pic_vp_2d2v_cart_lbfr) :: sim
  character(len=256)                          :: filename
  integer                                     :: rank, size
  type(sll_t_time_mark)  ::  t1
  sll_real64           :: time

  call sll_s_boot_collective()
  size = sll_f_get_collective_size(sll_v_world_collective)
  rank = sll_f_get_collective_rank(sll_v_world_collective)

  if (rank==0) call sll_s_set_time_mark(t1)

  call get_command_argument(1, filename)
  call sim%init_from_file(trim(filename))

  if (rank==0) then
    print*, size, 'mpi nodes X', sim%n_particles, 'particles', &
         sim%mesh%num_cells1, 'X',sim%mesh%num_cells2,'cells'
    if( sim%particle_type == sll_p_lbfr_particles )then
      print*, (real(size,f64)/real(sim%mesh%num_cells1 * sim%mesh%num_cells2,f64)) &
           * real(sim%n_particles,f64), 'pushed particles (flow markers and deposition particles) per cell'
    else
      print*, (real(size,f64)/real(sim%mesh%num_cells1 * sim%mesh%num_cells2,f64)) &
           * real(sim%n_particles,f64), 'particles per cell'
    end if
  endif

  call sim%run()

  if (rank==0) then
     time = sll_f_time_elapsed_since(t1)
     print*, 'PASSED', ' Total time=', time
  endif

  call sll_s_halt_collective()

  ! call sim%delete()
end program sim_pic_vp_2d2v_cart_lbfr
