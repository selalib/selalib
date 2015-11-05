program sim_bsl_vp_2d2d_cart_fullsparsegrid
#include "sll_working_precision.h"
  use sll_m_sim_bsl_vp_2d2d_cart_fullsparsegrid

  type(sll_t_sim_sl_vp_2d2v_cart_fullsparsegrid) :: sim
  character(len=256)                         :: filename

  ! Read in the simulation parameters from file specified in command line
  call get_command_argument(1, filename)
  call sim%init_from_file(trim(filename))
  
  call sim%run()

end program sim_bsl_vp_2d2d_cart_fullsparsegrid
