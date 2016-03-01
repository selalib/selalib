program sim_pic_vp_2d2v_cart

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  use sll_m_collective, only: &
    sll_s_boot_collective, &
    sll_f_get_collective_rank, &
    sll_f_get_collective_size, &
    sll_s_halt_collective, &
    sll_v_world_collective

  use sll_m_sim_pic_vp_2d2v_cart, only: &
    sll_t_sim_pic_vp_2d2v_cart

  implicit none
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
  type(sll_t_sim_pic_vp_2d2v_cart)  :: sim
  character(len=256)                :: filename
  integer                           :: rank, size
  logical                           :: passed
 

  call sll_s_boot_collective()
  size = sll_f_get_collective_size(sll_v_world_collective)
  rank = sll_f_get_collective_rank(sll_v_world_collective)
  
  ! Read in the simulation parameters from file specified in command line
  call get_command_argument(1, filename)
  call sim%init_from_file(trim(filename))
  
  call sim%run()

  if (rank == 0) then
     passed = sim%ctest_passed
  end if

  call sim%delete()

  call sll_s_halt_collective()

  if (rank == 0) then
     if ( passed .eqv. .true. ) then
        print*, 'PASSED.'
     else
        print*, 'FAILED.'
     end if
  end if

end program sim_pic_vp_2d2v_cart
