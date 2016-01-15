program sim_sl_vp_3d3v_cart_sparsegrid
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_working_precision.h"

  use sll_m_collective, only: &
    sll_s_boot_collective, &
    sll_f_get_collective_rank, &
    sll_f_get_collective_size, &
    sll_s_halt_collective, &
    sll_v_world_collective

  use sll_m_sim_sl_vp_3d3v_cart_sparsegrid, only: &
    sll_t_sim_sl_vp_3d3v_cart_sparsegrid

  implicit none
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  type(sll_t_sim_sl_vp_3d3v_cart_sparsegrid) :: sim
  character(len=256)                         :: filename
  sll_int32                                  :: rank, size

  call sll_s_boot_collective()
  size = sll_f_get_collective_size(sll_v_world_collective)
  rank = sll_f_get_collective_rank(sll_v_world_collective)
  
  ! Read in the simulation parameters from file specified in command line
  call get_command_argument(1, filename)
  call sim%init_from_file(trim(filename))
  
  call sim%run()

  call sll_s_halt_collective()


end program sim_sl_vp_3d3v_cart_sparsegrid
