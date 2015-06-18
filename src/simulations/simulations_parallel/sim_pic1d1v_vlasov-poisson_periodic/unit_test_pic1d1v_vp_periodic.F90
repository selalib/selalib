program pic1d1v_vp_periodic

#include "sll_working_precision.h"
    
    use sll_module_simulation_pic1d1v_vp_periodic, only: &
      sll_simulation_pic1d1v_vp_periodic, &
      sll_delete
    
    use sll_collective, only: &
      sll_world_collective,   &
      sll_collective_barrier, &
      sll_boot_collective,    &
      sll_get_collective_rank,&
      sll_halt_collective
    
    implicit none
!==============================================================================

  character(len=256) :: filename
  character(len=256) :: filename_local
  sll_int32          :: coll_rank

  type( sll_simulation_pic1d1v_vp_periodic ) :: simulation

!==============================================================================

  call sll_boot_collective()
  coll_rank = sll_get_collective_rank( sll_world_collective )
  if( coll_rank == 0 ) then
    print *, '#Booting parallel environment...'
  endif

  ! In this test, the name of the file to open is provided as a command line
  ! argument.
  call getarg(1, filename)
  filename_local = trim(filename)

  call simulation%init_from_file( filename_local )
  call simulation%run( )
  call sll_delete( simulation )

  print *, 'reached end of pic1d1v_vp_periodic test'
  print *, 'PASSED'

  call sll_halt_collective()

end program pic1d1v_vp_periodic
