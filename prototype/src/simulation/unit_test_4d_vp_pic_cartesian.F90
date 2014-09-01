program pic_4d_cartesian
  use sll_pic_simulation_4d_cartesian_module
  use sll_collective 

  type(sll_pic_simulation_4d_cartesian) :: sim

  integer :: rank, size

  call sll_boot_collective()
  size = sll_get_collective_size(sll_world_collective)
  rank = sll_get_collective_rank(sll_world_collective)

  call sim%init_from_file('params_pic_4d.nml')
! #        OMEGA  = 1.323
! #        GAMMA  = -0.151

  if (rank==0) then
     print*, sim%ions_number, 'particles',sim%m2d%num_cells1, &
       'x',sim%m2d%num_cells2,'cells'
     print*, sim%ions_number/real(sim%m2d%num_cells1* &
          sim%m2d%num_cells2,f64), 'particles per cell'
  endif

  call sim%run()

  if (rank==0) print*, 'PASSED'
  call sll_halt_collective()

  ! call sim%init_from_file("some_name")
  ! call sim%run()
  ! call sim%delete()
end program pic_4d_cartesian
