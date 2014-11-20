program guiding_center_2d_curvilinear
#include "sll_working_precision.h"
#include "sll_memory.h"
 use sll_simulation_2d_guiding_center_curvilinear_module
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

  if (len_trim(filename) == 0)then
    sim => new_guiding_center_2d_curvilinear( )
  else
    filename_local = trim(filename)
    sim => new_guiding_center_2d_curvilinear( filename_local )
  endif
print *,'#NEW GUIDING CENTER 2D CURVILINEAR'


call sim%run( )

  if(sll_get_collective_rank(sll_world_collective)==0)then
    print *, '#reached end of guiding_center_2d_curvilinear test'
    time = sll_time_elapsed_since(t0)
    print *, '#time elapsed since t0 : ',time
    print *, '#PASSED'
  endif

end program guiding_center_2d_curvilinear