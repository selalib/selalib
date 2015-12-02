program sim_pic_vp_2d2v_cart

#include "sll_working_precision.h"
  use sll_m_sim_pic_vp_2d2v_cart
  
  type(sll_t_sim_pic_vp_2d2v_cart)  :: sim
  character(len=256)                               :: filename
  integer                                          :: rank, size
 

  call sll_boot_collective()
  size = sll_get_collective_size(sll_world_collective)
  rank = sll_get_collective_rank(sll_world_collective)
  
  ! Read in the simulation parameters from file specified in command line
  call get_command_argument(1, filename)
  call sim%init_from_file(trim(filename))
  
  call sim%run()

  call sll_halt_collective()


end program sim_pic_vp_2d2v_cart
