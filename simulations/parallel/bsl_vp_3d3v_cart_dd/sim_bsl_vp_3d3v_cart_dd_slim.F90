program sim_bsl_vp_3d3v_cart_dd_slim
#include "sll_working_precision.h"

  use sll_m_collective, only: &
    sll_s_boot_collective, &
    sll_s_halt_collective

  use sll_m_sim_bsl_vp_3d3v_cart_dd_slim, only: &
       sll_t_sim_bsl_vp_3d3v_cart_dd_slim

  use mpi, only: mpi_thread_single, mpi_thread_funneled, mpi_thread_multiple

  use sll_m_utilities, only: sll_f_query_environment

#ifdef _OPENMP
  use omp_lib
#endif

  implicit none

  type(sll_t_sim_bsl_vp_3d3v_cart_dd_slim) :: sim
  character(len=256) :: filename
  integer :: mpi_mode

#ifdef _OPENMP
  mpi_mode = mpi_thread_multiple
#else
  mpi_mode = mpi_thread_single
#endif
  call sll_s_boot_collective(mpi_mode)

  ! Read in the simulation parameters from file specified in command line
  call get_command_argument(1, filename)
  call sim%init_from_file(trim(filename))

  call sim%run()

  call sim%delete()

  print *, 'Reached end of sim_bsl_vp_3d3v_cart_dd_slim.'

  call sll_s_halt_collective()

end program sim_bsl_vp_3d3v_cart_dd_slim
