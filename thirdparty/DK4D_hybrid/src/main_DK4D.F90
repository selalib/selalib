program main_DK4D
#include "sll_working_precision.h"


  use sll_m_collective
  use simulation_DK4D_module
  use sll_m_timer

  implicit none

  character(len=256) :: input_filename
  type(sll_simulation_DK4D) :: simu_DK4D

  !-> For CPU time computation
  type(sll_t_time_mark)  :: t0
  sll_real64 :: time 

  !-> For MPI parallelization
  sll_int32 :: world_size
  sll_int32 :: my_rank

  !*** Parallelization initialization ***
  call sll_s_boot_collective()
  world_size = sll_f_get_collective_size(sll_v_world_collective)
  my_rank    = sll_f_get_collective_rank(sll_v_world_collective)
  
  if ( my_rank==0 )then
    print *, '# *** Booting parallel environment...'
    print *, '# *** Start time mark t0'
    call sll_s_set_time_mark(t0)
  end if

  !*** Read the input data file for the 4D drift-kinetic simulation ***
  call get_command_argument(1,input_filename)
  if ( my_rank==0 )then
    print *, '# *** Read the input data file '
  end if
  call simu_DK4D%init_from_file(trim(input_filename))

  !*** Running ot the 4D drift-kinetic simulation ***
  if ( my_rank==0 )then
    print *, '# *** Run DK4D: START '
  end if
  call simu_DK4D%run( )

  if ( my_rank==0 )then
    print *, '# *** Run DK4D: END'
    time = sll_f_time_elapsed_since(t0)
    print *, '# *** Elapsed time since t0 : ',time
  end if

  call sll_s_halt_collective( )

end program main_DK4D
