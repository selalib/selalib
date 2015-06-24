! Main program for
! [[file:simulation_4d_vp_generic_pic_cartesian.F90::sll_simulation_4d_vp_generic_pic_cartesian_module]]

program unit_test_4d_vp_generic_pic_cartesian

#include "sll_working_precision.h"
  use sll_simulation_4d_vp_generic_pic_cartesian_module
  use sll_collective 
  use sll_timer

  type(sll_simulation_4d_vp_generic_pic_cartesian) :: sim
  character(len=256)                          :: filename
  integer                                     :: rank, size
  type(sll_time_mark)  ::  t1
  sll_real64           :: time

  call sll_boot_collective()
  size = sll_get_collective_size(sll_world_collective)
  rank = sll_get_collective_rank(sll_world_collective)

  if (rank==0) call sll_set_time_mark(t1)

  call get_command_argument(1, filename)
  call sim%init_from_file(trim(filename))

  if (rank==0) then
     print*, size, 'mpi nodes X', sim%number_particles, 'particles', &
          sim%mesh_2d%num_cells1, 'X',sim%mesh_2d%num_cells2,'cells'
     if( sim%use_lt_pic_scheme )then
         print*, (real(size,f64)/real(sim%mesh_2d%num_cells1 * sim%mesh_2d%num_cells2, f64)) &
              * sim%number_particles, 'transport markers (pushed particles) per cell'

         print*, (real(size,f64)/real(sim%mesh_2d%num_cells1 * sim%mesh_2d%num_cells2, f64)) &
              * sim%virtual_particle_number, 'virtual particles (deposited particles) per cell'
     else
         print*, (real(size,f64)/real(sim%mesh_2d%num_cells1 * sim%mesh_2d%num_cells2, f64)) &
              * sim%number_particles, 'particles per cell'
     end if
  endif

  call sim%run()

  if (rank==0) then
     time = sll_time_elapsed_since(t1)
     print*, 'PASSED', ' Total time=', time
  endif

  call sll_halt_collective()

  ! call sim%delete()
end program unit_test_4d_vp_generic_pic_cartesian
