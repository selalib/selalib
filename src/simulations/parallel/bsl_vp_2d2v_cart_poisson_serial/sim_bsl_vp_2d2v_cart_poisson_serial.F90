! Sample computation with the following characteristics:
! - vlasov-poisson
! - 4D cartesian: x,y, vx, vy (or x1, x2)
! - parallel for vlasov
! - sequential for poisson 

program sim_bsl_vp_2d2v_cart_poisson_serial
#include "sll_working_precision.h"
  use sll_m_sim_bsl_vp_2d2v_cart_poisson_serial
  use sll_m_collective
  use sll_m_timer
  implicit none 

!  character(len=256) :: filename
!  character(len=256) :: filename_local
!  class(sll_simulation_base_class), pointer :: sim
!  call sll_boot_collective()
!  if(sll_get_collective_rank(sll_world_collective)==0)then
!    print *, '#Booting parallel environment...'
!  endif

  !class(sll_simulation_base_class), pointer :: sim
  class(sll_simulation_4d_vlasov_par_poisson_seq_cart), pointer :: sim  
  
  character(len=256) :: filename
  character(len=256) :: filename_local
  type(sll_time_mark)  :: t0
  sll_int32 :: count
  sll_int32 :: i
  sll_int32 :: num_min
  sll_int32 :: num_max
  character(len=256) :: str
  call sll_boot_collective()
  if(sll_get_collective_rank(sll_world_collective)==0)then
    print *, '#Start time mark t0'
    call sll_set_time_mark(t0)
    print *, '#Booting parallel environment...'
  endif

  count = command_argument_count()
  if(sll_get_collective_rank(sll_world_collective)==0)then
    print *, '#count=',count
  endif
              
  call get_command_argument(1, filename)
  if (len_trim(filename) == 0)then
    sim => new_vlasov_par_poisson_seq_cart( )
    call sim%run( )
  else
    filename_local = trim(filename)
    call get_command_argument(2, str)
    if(len_trim(str) == 0)then
      sim => new_vlasov_par_poisson_seq_cart( filename_local )
      call sim%run( )
    else
      read(str , *) num_max
      num_min = 0
      call get_command_argument(3, str)
      if(len_trim(str) .ne. 0)then
        num_min = num_max
        read(str , *) num_max
      endif
      !print *,'#num=',num_min,num_max
      do i=num_min,num_max
        sim => new_vlasov_par_poisson_seq_cart( filename_local, i)
        call sim%run( )
        call delete_vp4d_par_cart( sim )
      enddo  
    endif    
  endif
  

  if(sll_get_collective_rank(sll_world_collective)==0)then
    print *, '#reached end of sim4d_vp_cart test'
    print *, '#time elapsed since t0 : ', sll_time_elapsed_since(t0)
    print *, '#PASSED'
  endif

  call sll_halt_collective()




  
  ! In this test, the name of the file to open is provided as a command line
  ! argument.

  
!  call get_command_argument(1, filename)
!  filename_local = trim(filename)
!  
!  sim => new_vlasov_par_poisson_seq_cart( filename_local )
!  
!  call sim%run()
!  
!  !call simulation%init_from_file(filename_local)
!  !call simulation%run( )
!  !call delete_vp2d_par_cart(simulation)
!  
!  if(sll_get_collective_rank(sll_world_collective)==0)then
!    print *, '#reached end of sim4d_vp_cart test'
!    print *, '#PASSED'
!  endif
!  call sll_halt_collective()


end program sim_bsl_vp_2d2v_cart_poisson_serial
