! Sample computation with the following characteristics:
!   - Vlasov-Poisson
!   - 1Dx1D cartesian: x1=x, x2=vx
!   - Parallel
!   - Two sll_t_species

program sim_bsl_vp_1d1v_cart_multi_species

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_memory.h"
#include "sll_working_precision.h"

  use sll_m_collective, only: &
    sll_s_boot_collective, &
    sll_f_get_collective_rank, &
    sll_s_halt_collective, &
    sll_v_world_collective

  use sll_m_common_array_initializers, only: &
    sll_f_landau_initializer_2d, &
    sll_i_scalar_initializer_2d

  use sll_m_sim_bsl_vp_1d1v_cart_multi_species, only: &
    sll_s_delete_vp2d_par_cart_multi_species, &
    sll_f_new_vp2d_par_cart_multi_species, &
    sll_t_simulation_2d_vlasov_poisson_cart_multi_species

  use sll_m_species, only: &
    sll_s_change_initial_function_species

  use sll_m_timer, only: &
    sll_s_set_time_mark, &
    sll_f_time_elapsed_since, &
    sll_t_time_mark

  implicit none
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class(sll_t_simulation_2d_vlasov_poisson_cart_multi_species), pointer :: sim

character(len=256)  :: filename
character(len=256)  :: filename_local
type(sll_t_time_mark) :: t0
sll_real64          :: time
sll_int32           :: ierr
sll_int32           :: i
sll_int32           :: num_min
sll_int32           :: num_max
sll_int32           :: istep
character(len=256)  :: str

procedure(sll_i_scalar_initializer_2d), pointer :: init_func

sll_real64, dimension(:), pointer :: params
sll_int32                         :: num_params
logical                           :: init_from_unit_test  


init_from_unit_test = .false.

call sll_s_boot_collective()

if(sll_f_get_collective_rank(sll_v_world_collective)==0)then

  print *, '#Start time mark t0'
  call sll_s_set_time_mark(t0)
  print *, '#Booting parallel environment...'

endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

call get_command_argument(1, filename)

if (len_trim(filename) == 0)then

  sim => sll_f_new_vp2d_par_cart_multi_species( )
  do istep = 1, sim%num_iterations
    sim%istep = istep
    call sim%run( )
  end do

else

  filename_local = trim(filename)
  call get_command_argument(2, str)

  if(len_trim(str) == 0)then

    sim => sll_f_new_vp2d_par_cart_multi_species( filename_local )
    
    if (init_from_unit_test) then

      print *,'#Warning: init_function is redefined form unit_test'
      init_func => sll_f_landau_initializer_2d
      num_params = 2
      SLL_ALLOCATE(params(num_params),ierr)  
      params(1) = 0.26_f64
      params(2) = 100._f64  

      call sll_s_change_initial_function_species(sim%sp(1),init_func,params,num_params)
      call sll_s_change_initial_function_species(sim%sp(2),init_func,params,num_params)
      
    endif

    do istep = 1, sim%num_iterations
      sim%istep = istep
      call sim%run( )
    end do

  else

    read(str , *) num_max
    num_min = 0
    call get_command_argument(3, str)

    if(len_trim(str) .ne. 0)then

      num_min = num_max
      read(str , *) num_max

    endif

    do i=num_min,num_max
      sim => sll_f_new_vp2d_par_cart_multi_species( filename_local, i)
      do istep = 1, sim%num_iterations
        sim%istep = istep
        call sim%run( )
      end do
      call sll_s_delete_vp2d_par_cart_multi_species( sim )
      nullify( sim )
    enddo  

  endif    

endif

if(sll_f_get_collective_rank(sll_v_world_collective)==0)then

  print *, '#reached end of vp1d1v multi sll_t_species test'
  time = sll_f_time_elapsed_since(t0)
  print *, '#time elapsed since t0 : ',time
  print *, '#PASSED'

endif

call sll_s_halt_collective()

end program sim_bsl_vp_1d1v_cart_multi_species
