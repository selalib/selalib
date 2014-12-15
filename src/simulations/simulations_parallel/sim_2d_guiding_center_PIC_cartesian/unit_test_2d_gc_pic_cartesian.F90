program pic_2d_cartesian
  use sll_pic_simulation_2d_cartesian_module
  use sll_collective 

  type(sll_pic_simulation_2d_gc_cartesian) :: sim
  character(len=256)                       :: filename
  integer                                  :: rank, size

  call sll_boot_collective()
  size = sll_get_collective_size(sll_world_collective)
  rank = sll_get_collective_rank(sll_world_collective)

  call get_command_argument(1, filename)
  call sim%init_from_file(trim(filename))
! #        OMEGA  = 1.323
! #        GAMMA  = -0.151

!!$   if (rank==0) then
!!$      print*, sim%ions_number, 'particles',sim%m2d%num_cells1, &
!!$           'x',sim%m2d%num_cells2,'cells' 
!!$      print*, sim%ions_number/real(sim%m2d%num_cells1* &
!!$           sim%m2d%num_cells2,f64), 'particles per cell'
!!$   endif

  call sim%run()

  if (rank==0) print*, 'PASSED'
  call sll_halt_collective()

  ! call sim%delete()
end program pic_2d_cartesian
