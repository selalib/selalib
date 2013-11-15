program guiding_center_2d_polar
  use sll_simulation_2d_guiding_center_polar_module
  implicit none
!
!  character(len=256) :: filename
!  character(len=256) :: filename_local
  class(sll_simulation_base_class), pointer :: sim
  
!  call sll_boot_collective()
!  if(sll_get_collective_rank(sll_world_collective)==0)then
!    print *, '#Booting parallel environment...'
!  endif
!
!  ! In this test, the name of the file to open is provided as a command line
!  ! argument.
!  call getarg(1, filename)
!  filename_local = trim(filename)
!  call simulation%init_from_file(filename_local)
  
  sim => new_guiding_center_2d_polar()
  
  call sim%run( )
!  call delete_vp2d_par_cart(simulation)
!  if(sll_get_collective_rank(sll_world_collective)==0)then
!    print *, '#reached end of vp2d test'
!    print *, '#PASSED'
!  endif
!  call sll_halt_collective()
!

  print *,'#PASSED'

end program guiding_center_2d_polar
