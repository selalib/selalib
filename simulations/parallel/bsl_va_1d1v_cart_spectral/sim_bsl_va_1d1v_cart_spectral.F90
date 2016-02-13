! Sample computation with the following characteristics:
! - vlasov-poisson
! - 1Dx1D cartesian: x1=x, x2=vx
! - parallel

program sim_bsl_va_1d1v_cart_spectral
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_working_precision.h"

  use sll_m_collective, only: &
    sll_s_boot_collective, &
    sll_f_get_collective_rank, &
    sll_s_halt_collective, &
    sll_v_world_collective

  use sll_m_common_array_initializers, only: &
    sll_i_scalar_initializer_2d

  use sll_m_sim_bsl_va_1d1v_cart_spectral, only: &
    sll_f_new_va2d_par_cart, &
    sll_t_simulation_2d_vlasov_ampere_cart

  use sll_m_timer, only: &
    sll_s_set_time_mark, &
    sll_f_time_elapsed_since, &
    sll_t_time_mark

  implicit none
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class(sll_t_simulation_2d_vlasov_ampere_cart), pointer :: sim
character(len=256) :: filename
character(len=256) :: filename_local
type(sll_t_time_mark)  :: t0
sll_real64 :: time
!sll_int32 :: ierr 
procedure(sll_i_scalar_initializer_2d), pointer :: init_func, equil_func
!sll_real64, dimension(:), pointer :: params, equil_params
!sll_int32 :: num_params

call sll_s_boot_collective()
if(sll_f_get_collective_rank(sll_v_world_collective)==0)then
  print *, '#Start time mark t0'
  call sll_s_set_time_mark(t0)
  print *, '#Booting parallel environment...'
endif


call get_command_argument(1, filename)
filename_local = trim(filename)
sim => sll_f_new_va2d_par_cart( filename_local )

call sim%run( )

if(sll_f_get_collective_rank(sll_v_world_collective)==0)then
  print *, '#reached end of vp2d test'
  time = sll_f_time_elapsed_since(t0)
  print *, '#time elapsed since t0 : ',time
  print *, '#PASSED'
endif

call sll_s_halt_collective()

end program sim_bsl_va_1d1v_cart_spectral


