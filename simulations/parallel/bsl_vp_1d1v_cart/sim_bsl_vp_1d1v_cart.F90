! Sample computation with the following characteristics:
! - vlasov-poisson
! - 1Dx1D cartesian: x1=x, x2=vx
! - parallel

program sim_bsl_vp_1d1v_cart

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

   use sll_m_sim_bsl_vp_1d1v_cart, only: &
      sll_s_change_initial_function_vp2d_par_cart, &
      sll_s_delete_vp2d_par_cart, &
      sll_f_new_vp2d_par_cart, &
      sll_t_simulation_2d_vlasov_poisson_cart

   use sll_m_timer, only: &
      sll_s_set_time_mark, &
      sll_f_time_elapsed_since, &
      sll_t_time_mark

   implicit none
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   type(sll_t_simulation_2d_vlasov_poisson_cart) :: sim

   character(len=256)    :: filename
   character(len=256)    :: filename_local
   type(sll_t_time_mark) :: t0
   sll_int32             :: ierr
   sll_int32             :: i
   sll_int32             :: num_min
   sll_int32             :: num_max
   character(len=256)    :: str

   sll_real64, dimension(:), pointer :: params
   sll_int32                         :: num_params
   logical                           :: init_from_unit_test
   procedure(sll_i_scalar_initializer_2d), pointer :: init_func

   init_from_unit_test = .false.

   call sll_s_boot_collective()

   if (sll_f_get_collective_rank(sll_v_world_collective) == 0) then

      print *, '#Start time mark t0'
      call sll_s_set_time_mark(t0)
      print *, '#Booting parallel environment...'

   end if

   call get_command_argument(1, filename)
   if (len_trim(filename) == 0) then
      call sim%init()
      call sim%run()
   else
      filename_local = trim(filename)
      call get_command_argument(2, str)
      if (len_trim(str) == 0) then
         call sim%init(filename_local)

         if (init_from_unit_test) then

            print *, '#Warning: init_function is redefined form unit_test'
            num_params = 2
            SLL_ALLOCATE(params(num_params), ierr)
            params(1) = 0.26_f64
            params(2) = 100._f64
            init_func => sll_f_landau_initializer_2d

            call sll_s_change_initial_function_vp2d_par_cart( &
               sim, init_func, params, num_params)

         end if

         call sim%run()
      else
         read (str, *) num_max
         num_min = 0
         call get_command_argument(3, str)
         if (len_trim(str) .ne. 0) then
            num_min = num_max
            read (str, *) num_max
         end if
         !print *,'#num=',num_min,num_max
         do i = num_min, num_max
            call sim%init(filename_local, i)
            call sim%run()
            call sim%free()
         end do
      end if
   end if

   if (sll_f_get_collective_rank(sll_v_world_collective) == 0) then

      print *, '#reached end of vp2d test'
      print *, '#time elapsed since t0 : ', sll_f_time_elapsed_since(t0)
      print *, '#PASSED'

   end if

   call sll_s_halt_collective()

end program sim_bsl_vp_1d1v_cart
