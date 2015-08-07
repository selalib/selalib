program unit_test_pic_1d2v_vm_cart

#include "sll_working_precision.h"
  use sll_m_sim_pic_1d2v_vm_cart
  
  type(sll_sim_pic_1d2v_vm_cart)  :: sim
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


end program unit_test_pic_1d2v_vm_cart
