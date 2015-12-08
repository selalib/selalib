! Sample computation with the following characteristics:
! - vlasov-poisson
! - 6D cartesian: x, y, z, vx, vy, vz (or x1, x2, x3, x4, x5, x6)
! - parallel

program sim_pif_vp_ndnv_cart
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  use iso_fortran_env, only: &
    output_unit

  use sll_m_collective, only: &
    sll_boot_collective, &
    sll_get_collective_rank, &
    sll_halt_collective, &
    sll_world_collective

  use sll_m_sim_pif_vp_ndnv_cart, only: &
    sll_simulation_general_vlasov_poisson_pif

  implicit none
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! 
   character(len=256) :: filename
   character(len=256) :: filename_local
   integer :: coll_rank
   integer :: idx
   
   !print *, 'Booting parallel environment for simulation...'
   call sll_boot_collective() ! Wrap this up in something else
  
   coll_rank=sll_get_collective_rank( sll_world_collective )
 
   if (coll_rank == 0) then
       print *, 'Proceed to run simulation.'
       print *, 'There are ', COMMAND_ARGUMENT_COUNT(), ' simulation files to be run.'
   endif
   flush( output_unit )
   
  ! Provide files to open as a command line argument, seprated by spaces
  ! argument.
  
  do idx=1, COMMAND_ARGUMENT_COUNT()
  call get_command_argument(idx, filename)
  filename_local = trim(filename)
  call run_from_file(filename)
  end do
  
!   call delete_vp6d_par_cart(simulation)
   print *, 'reached end of PIF test'
   print *, 'PASSED'
   call sll_halt_collective()

   contains
   
   subroutine run_from_file(filename)
      type(sll_simulation_general_vlasov_poisson_pif) :: simulation
      character(len=256), intent(in) :: filename
      call simulation%init_from_file(filename)
      call simulation%run()    
   end subroutine run_from_file
   
   
end program sim_pif_vp_ndnv_cart


