program test_poisson_1d_periodic
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_memory.h"
#include "sll_working_precision.h"

   use sll_m_poisson_1d_base, only: &
      sll_c_poisson_1d_base

   use sll_m_constants, only: &
      sll_p_pi

   use sll_m_poisson_1d_periodic, only: &
      sll_o_initialize, &
      sll_t_poisson_1d_periodic, &
      sll_f_new_poisson_1d_periodic, &
      sll_o_solve

   implicit none
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   sll_real64, dimension(:), allocatable :: ex
   sll_real64, dimension(:), allocatable :: ex_exact
   sll_real64, dimension(:), allocatable :: rho

   type(sll_t_poisson_1d_periodic)       :: poisson
   class(sll_c_poisson_1d_base), pointer :: poisson_class

   sll_int32   :: nc_eta1
   sll_real64  :: eta1_min
   sll_real64  :: eta1_max
   sll_real64  :: delta_eta1
   sll_real64  :: x
   sll_int32   :: error
   sll_int32   :: mode
   sll_int32   :: i
   sll_real64  :: err

   nc_eta1 = 128

   SLL_ALLOCATE(rho(nc_eta1 + 1), error)
   SLL_ALLOCATE(ex(nc_eta1 + 1), error)
   SLL_ALLOCATE(ex_exact(nc_eta1 + 1), error)

   eta1_min = 0.0_f64
   eta1_max = 2*sll_p_pi
   delta_eta1 = (eta1_max - eta1_min)/nc_eta1
   mode = 4
   do i = 1, nc_eta1 + 1
      x = (i - 1)*delta_eta1
      rho(i) = real(mode*mode, f64)*sin(mode*x)
      ex_exact(i) = -real(mode, f64)*cos(mode*x)
   end do

   call sll_o_initialize(poisson, eta1_min, eta1_max, nc_eta1, error)

   call sll_o_solve(poisson, ex, rho)

   err = maxval(abs(ex - ex_exact))
   print *, 'mode=', mode, '   error=', err

   if (err > 1.e-14) stop 'FAILED'

   poisson_class => sll_f_new_poisson_1d_periodic(eta1_min, &
                                                  eta1_max, &
                                                  nc_eta1)
   do i = 1, nc_eta1 + 1
      x = (i - 1)*delta_eta1
      rho(i) = real(mode*mode, f64)*sin(mode*x)
   end do

   call poisson_class%compute_e_from_rho(ex, rho)

   err = maxval(abs(ex - ex_exact))
   print *, 'mode=', mode, '   error=', err

   if (err > 1.e-14) stop 'FAILED'

   print *, '#PASSED'

end program test_poisson_1d_periodic
