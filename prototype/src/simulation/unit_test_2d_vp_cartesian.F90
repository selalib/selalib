! Sample computation with the following characteristics:
! - vlasov-poisson
! - 1Dx1D cartesian: x1=x, x2=vx
! - parallel

program vlasov_poisson_2d
#include "sll_working_precision.h"
  use sll_simulation_2d_vlasov_poisson_cartesian
  use sll_collective
  use sll_timer
  use sll_constants
  implicit none
  
  class(sll_simulation_base_class), pointer :: sim
  character(len=256) :: filename
  character(len=256) :: filename_local
  type(sll_time_mark)  :: t0
  sll_real64 :: time 

  call sll_boot_collective()
  if(sll_get_collective_rank(sll_world_collective)==0)then
    print *, '#Start time mark t0'
    call sll_set_time_mark(t0)
    print *, '#Booting parallel environment...'
  endif


  call get_command_argument(1, filename)
  
  print *,'#filename=',filename
  
  if (len_trim(filename) == 0)then
    sim => new_vp2d_par_cart( )
  else
    filename_local = trim(filename)
    sim => new_vp2d_par_cart( filename_local )
  endif
  call sim%run( )

  if(sll_get_collective_rank(sll_world_collective)==0)then
    print *, '#reached end of vp2d test'
    time = sll_time_elapsed_since(t0)
    print *, '#time elapsed since t0 : ',time
    print *, '#PASSED'
  endif

  call sll_halt_collective()


end program vlasov_poisson_2d


