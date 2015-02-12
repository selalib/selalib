program vlasov_poisson_2d_no_split
#include "sll_working_precision.h"
#include "sll_constants.h"
  use sll_simulation_2d_vlasov_poisson_no_splitting
  use sll_timer
  implicit none

  class(sll_simulation_base_class), pointer :: sim
  character(len=256) :: filename
  character(len=256) :: filename_local
  type(sll_time_mark)  :: t0
  sll_real64 :: time 
  
  print *, '#Start time mark t0'
  call sll_set_time_mark(t0)

  call get_command_argument(1, filename)
  if (len_trim(filename) == 0)then
    sim => new_vp2d_no_split( )
  else
    filename_local = trim(filename)
    sim => new_vp2d_no_split( filename_local )
  endif
  call sim%run( )
  time = sll_time_elapsed_since(t0)
  print *, '#time elapsed since t0 : ',time
  print *,'#PASSED'

end program vlasov_poisson_2d_no_split
