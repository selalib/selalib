program sim_pic_vp_1d1v_cart

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_working_precision.h"

  use sll_m_collective, only: &
    sll_s_boot_collective, &
    sll_s_collective_barrier, &
    sll_f_get_collective_rank, &
    sll_s_halt_collective, &
    sll_v_world_collective

  use sll_m_sim_pic_vp_1d1v_cart, only: &
    sll_o_delete, &
    sll_t_simulation_pic1d1v_vp_periodic

  implicit none
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!==============================================================================

  character(len=256) :: filename
  character(len=256) :: filename_local
  sll_int32          :: coll_rank

  type( sll_t_simulation_pic1d1v_vp_periodic ) :: simulation

!==============================================================================

  call sll_s_boot_collective()
  coll_rank = sll_f_get_collective_rank( sll_v_world_collective )
  if( coll_rank == 0 ) print *, '#Booting parallel environment...'

  ! In this test, the name of the file to open is provided as a command line
  ! argument.
  call get_command_argument(1, filename)
  filename_local = trim(filename)

  call simulation%init_from_file( filename_local )
  if( coll_rank == 0 ) print *, 'simulation initialized from file'

  call simulation%new_pic()
  if( coll_rank == 0 ) print *, 'parallel PIC ready to run'

  call simulation%run( )
  if( coll_rank == 0 ) print *, 'run completed'

!  call sll_o_delete( simulation )

  if( coll_rank == 0 ) print *, 'reached end of pic1d1v_vp_periodic test'
  if( coll_rank == 0 ) print *, 'PASSED'

  call sll_s_halt_collective()

end program sim_pic_vp_1d1v_cart
