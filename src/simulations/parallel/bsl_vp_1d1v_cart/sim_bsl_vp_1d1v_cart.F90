! Sample computation with the following characteristics:
! - vlasov-poisson
! - 1Dx1D cartesian: x1=x, x2=vx
! - parallel

program sim_bsl_vp_1d1v_cart

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_memory.h"
#include "sll_working_precision.h"

  use sll_m_collective, only: &
    sll_boot_collective, &
    sll_get_collective_rank, &
    sll_halt_collective, &
    sll_world_collective

  use sll_m_common_array_initializers, only: &
    sll_landau_initializer_2d, &
    sll_scalar_initializer_2d

  use sll_m_sim_bsl_vp_1d1v_cart, only: &
    change_initial_function_vp2d_par_cart, &
    delete_vp2d_par_cart, &
    new_vp2d_par_cart, &
    sll_simulation_2d_vlasov_poisson_cart

  use sll_m_timer, only: &
    sll_set_time_mark, &
    sll_time_elapsed_since, &
    sll_time_mark

  implicit none
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class(sll_simulation_2d_vlasov_poisson_cart), pointer :: sim

character(len=256)  :: filename
character(len=256)  :: filename_local
type(sll_time_mark) :: t0
sll_real64          :: time
sll_int32           :: ierr
sll_int32           :: i
sll_int32           :: num_min
sll_int32           :: num_max
character(len=256)  :: str

procedure(sll_scalar_initializer_2d), pointer :: init_func

sll_real64, dimension(:), pointer :: params
sll_int32                         :: num_params
logical                           :: init_from_unit_test  



init_from_unit_test = .false.

call sll_boot_collective()

if(sll_get_collective_rank(sll_world_collective)==0)then

  print *, '#Start time mark t0'
  call sll_set_time_mark(t0)
  print *, '#Booting parallel environment...'

endif


  call get_command_argument(1, filename)
  if (len_trim(filename) == 0)then
    sim => new_vp2d_par_cart( )
    call sim%run( )
  else
    filename_local = trim(filename)
    call get_command_argument(2, str)
    if(len_trim(str) == 0)then
      sim => new_vp2d_par_cart( filename_local )
      
      if (init_from_unit_test) then

        print *,'#Warning: init_function is redefined form unit_test'
        init_func => sll_landau_initializer_2d
        num_params = 2
        SLL_ALLOCATE(params(num_params),ierr)  
        params(1) = 0.26_f64
        params(2) = 100._f64  

        call change_initial_function_vp2d_par_cart( &
          sim,                                      &
          init_func,                                &
          params,                                   &
          num_params)
        
      endif

      
      
      
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
        sim => new_vp2d_par_cart( filename_local, i)
        call sim%run( )
        call delete_vp2d_par_cart( sim )
        nullify( sim )
      enddo  
    endif    
  endif





if(sll_get_collective_rank(sll_world_collective)==0)then

  print *, '#reached end of vp2d test'
  time = sll_time_elapsed_since(t0)
  print *, '#time elapsed since t0 : ',time
  print *, '#PASSED'

endif

call sll_halt_collective()

end program sim_bsl_vp_1d1v_cart
