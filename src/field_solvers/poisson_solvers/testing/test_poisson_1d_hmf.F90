program test_poisson_1d_hmf
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_memory.h"
#include "sll_working_precision.h"

  use sll_m_poisson_1d_base, only: &
    sll_c_poisson_1d_base

  use sll_m_constants, only: &
    sll_p_pi

  use sll_m_poisson_1d_hmf, only: &
    sll_o_initialize, &
    sll_t_poisson_1d_hmf, &
    sll_f_new_poisson_1d_hmf, &
    sll_o_solve

  implicit none
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  sll_real64, dimension(:),   allocatable :: ex
  sll_real64, dimension(:),   allocatable :: ex_exact
  sll_real64, dimension(:),   allocatable :: rho
  sll_real64, dimension(:,:), allocatable :: f

  type(sll_t_poisson_1d_hmf)       :: poisson
  class(sll_c_poisson_1d_base), pointer :: poisson_class

  sll_int32   :: nc_eta1, nc_eta2
  sll_real64  :: eta1_min, eta2_min
  sll_real64  :: eta1_max, eta2_max
  sll_real64  :: delta_eta1, delta_eta2
  sll_real64  :: x, v
  sll_int32   :: error
  sll_int32   :: i, j
  sll_real64, parameter :: alpha = 0.000119_f64
  sll_real64, parameter :: beta  = 10.0_f64
  sll_real64, parameter :: mass  = 1.0_f64
  sll_real64, parameter :: m     = 0.946_f64

  nc_eta1 = 128
  nc_eta2 = 128

  SLL_ALLOCATE(rho(nc_eta1+1),error)
  SLL_ALLOCATE(ex(nc_eta1+1),error)
  SLL_ALLOCATE(ex_exact(nc_eta1+1),error)
  SLL_ALLOCATE(f(nc_eta1+1, nc_eta2+1),error)

  eta1_min = -sll_p_pi
  eta1_max = +sll_p_pi
  delta_eta1 = (eta1_max-eta1_min) / nc_eta1
  eta2_min = -3.0_f64
  eta2_max = +3.0_f64
  delta_eta2 = (eta2_max-eta2_min) / nc_eta2

  do j = 1, nc_eta2+1
    v = eta2_min+(j-1)*delta_eta2
    do i = 1, nc_eta1+1
       x = eta1_min+(i-1)*delta_eta1
       f(i,j) = alpha * exp(-beta*(((v*v)*0.5_f64) - m * cos(x)))
    end do
  end do

  do i=1,nc_eta1+1
    x = eta1_min+(i-1)*delta_eta1
    rho(i)      =  sum(f(i,:)) * delta_eta2
    ex_exact(i) = -m * sin(x)
  end do

  call sll_o_initialize(poisson, eta1_min, eta1_max, nc_eta1, error) 

  call sll_o_solve(poisson, ex, rho)

  do i = 1, nc_eta1
    x = eta1_min+(i-1)*delta_eta1
    write(11,*) x, ex(i), ex_exact(i)
  end do

  print*,'   error=',maxval(abs(ex-ex_exact))

  if (error > 1.e-14) stop 'FAILED'
 
  poisson_class => sll_f_new_poisson_1d_hmf( eta1_min, &
                                                  eta1_max, &
                                                  nc_eta1   )
  do i=1,nc_eta1+1
    x = eta1_min+(i-1)*delta_eta1
    rho(i)      =  sum(f(i,:)) * delta_eta2
    ex_exact(i) = -m * sin(x)
  end do

  call poisson_class%compute_e_from_rho( ex, rho )

  print*,'   error=',maxval(abs(ex-ex_exact))

  if (error > 1.e-14) stop 'FAILED'

  print*, '#PASSED'

end program test_poisson_1d_hmf
