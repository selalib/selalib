program sim_pic_gc_2d0v_cart_optim
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_working_precision.h"

  use sll_m_collective, only: &
    sll_boot_collective, &
    sll_get_collective_rank, &
    sll_get_collective_size, &
    sll_halt_collective, &
    sll_world_collective

  use sll_m_sim_pic_gc_2d0v_cart_optim, only: &
    sll_pic_simulation_2d_gc_cartesian

  use sll_m_timer, only: &
    sll_set_time_mark, &
    sll_time_elapsed_since, &
    sll_time_mark

  implicit none
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  type(sll_pic_simulation_2d_gc_cartesian) :: sim
  character(len=256)                       :: filename
  integer                                  :: rank, size
  type(sll_time_mark)  ::  t1
  sll_real64           :: time
  
  call sll_boot_collective()
  size = sll_get_collective_size(sll_world_collective)
  rank = sll_get_collective_rank(sll_world_collective)

  if (rank==0) call sll_set_time_mark(t1)

  call get_command_argument(1, filename)
  call sim%init_from_file(trim(filename))

!!$   if (rank==0) then
!!$      print*, size, 'mpi nodes X',sim%parts_number, 'particles', &
!!$           sim%m2d%num_cells1,'x',sim%m2d%num_cells2, 'cells' 
!!$      print*, sim%parts_number/real(sim%m2d%num_cells1* &
!!$           sim%m2d%num_cells2,f64), 'particles per cell'
!!$   endif

  call sim%run()

  if (rank==0)  then
     time = sll_time_elapsed_since(t1)
     print*, 'PASSED', ' Total time of GC sim=', time
  endif
  
  call sll_halt_collective()

  ! call sim%delete()
end program sim_pic_gc_2d0v_cart_optim
