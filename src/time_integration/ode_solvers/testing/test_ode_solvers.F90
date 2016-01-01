program test_ode_solvers
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_memory.h"
#include "sll_working_precision.h"

  use sll_m_constants, only: &
    sll_p_pi

  use sll_m_ode_solvers, only: &
    sll_p_compact_ode, &
    sll_s_implicit_ode_nonuniform, &
    sll_p_periodic_ode

  implicit none
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
  sll_int32 :: i, ncells, ierr, order
  sll_real64 :: xmin, xmax, deltax, deltat, x, error, a
  sll_real64, dimension(:), pointer :: a_n, a_np1, xin, xout

  print*, 'checking sll_s_implicit_ode_nonuniform'
  print*, ' '
  ncells = 100
  xmin = 0.0_f64
  xmax = sll_p_pi
  deltax = (xmax-xmin) / ncells
  
  deltat = 0.01_f64

  SLL_ALLOCATE(a_n(ncells+1), ierr)
  SLL_ALLOCATE(a_np1(ncells+1), ierr)
  SLL_ALLOCATE(xout(ncells+1), ierr)
  SLL_ALLOCATE(xin(ncells+1), ierr)

  call random_number(xin)
  do i=1, ncells
     xin(i+1) = xin(i) + xin(i+1)
  end do
  xin = (xin-xin(1))/(xin(ncells+1)-xin(1))*sll_p_pi

  print*, 'checking compact boundary conditions'
  print*, '------------------------------------'

    print*,'testing time dependent field. ' 
  ! consider the ode dx/dt = a (t,x) = t*x
  do i = 1, ncells+1
     x = xin(i)
     a_n(i) = 0.0_f64
     a_np1(i) = deltat*x
  end do

  print*, 'testing order 1:'
  order = 1
  call sll_s_implicit_ode_nonuniform( order, deltat, xin, ncells, sll_p_compact_ode, xout,  a_n ) 
print *, 'do I see this?'
  error = 0.0_f64
  do i = 1, ncells
     x = xin(i)
     error = max(error, abs(modulo(x*exp(-deltat**2/2), xmax-xmin) - xout(i)))
     !print*, i, x*exp(-deltat**2/2), xout(i)
  end do
  print*,'     error=', error

    print*, 'testing order 2:'
  order = 2
  call sll_s_implicit_ode_nonuniform( order, deltat, xin, ncells, sll_p_compact_ode, xout,  a_n, a_np1 ) 
  error = 0.0_f64
  do i = 1, ncells
     x = xin(i)
     error = max(error, abs(modulo(x*exp(-deltat**2/2), xmax-xmin) - xout(i)))
     !print*, i, modulo(xmin+(i-1)*deltax - deltat, xmax-xmin), xout(i)
  end do
  print*,'     error=', error
  ! test the ode dx / dt = sin(x) -> solution is x(t) = x(0) * exp(t)
  do i = 1, ncells+1
     x = xin(i)
     a_n(i) = sin(x)
     a_np1(i) = sin(x)
  end do

  print*, 'testing order 1 on a(x) = sin(x)'
  order = 1
  call sll_s_implicit_ode_nonuniform( order, deltat, xin, ncells, sll_p_compact_ode, xout,  a_n ) 
  error = 0.0_f64
  do i = 1, ncells
     x = xin(i)
     error = max(error,abs(xmin+modulo(2.0_f64*atan(tan(x/2)*exp(-deltat))-xmin, xmax-xmin)-xout(i)))
     !print*, i, x*exp(-deltat**2/2), xout(i)
  end do
  print*,'     error=', error

    print*, 'testing order 2 on a(x) = sin(x)'
  order = 2
  call sll_s_implicit_ode_nonuniform( order, deltat, xin, ncells,  sll_p_compact_ode, xout,  a_n, a_np1 ) 
  error = 0.0_f64
  do i = 1, ncells
     x = xin(i)
     error = max(error,abs(xmin+modulo(2.0_f64*atan(tan(x/2)*exp(-deltat))-xmin, xmax-xmin)-xout(i)))
     !print*, i, 2*atan(tan(x/2)*exp(-deltat)), xout(i)
  end do
  print*,'     error=', error


  print*, 'checking periodic boundary conditions'
  print*, '-------------------------------------'

  print*,'testing constant field. ' 
  ! consider the ode dx/dt = a
  print*, '    advection to the right:'
  a = 10.0_f64
  a_n(:) = a
  a_np1(:) = a 
  order = 2
  call sll_s_implicit_ode_nonuniform( order, deltat, xin, ncells, sll_p_periodic_ode, xout,  a_n, a_np1 ) 
  error = 0.0_f64
  do i = 1, ncells
     x = xin(i)
     error = max(error, abs(modulo(x-a*deltat, xmax-xmin) - xout(i)))
!     print*, i, abs(modulo(x-a*deltat, xmax-xmin)), xout(i)
  end do
  print*,'     error=', error

  print*, '    advection to the left:'
  a = -512.0_f64
  a_n(:) = a
  a_np1(:) = a
  order = 2
  call sll_s_implicit_ode_nonuniform( order, deltat, xin, ncells, sll_p_periodic_ode, xout,  a_n, a_np1 ) 
  error = 0.0_f64
  do i = 1, ncells
     x = xin(i)
     error = max(error, abs(modulo(x-a*deltat, xmax-xmin) - xout(i)))
     !print*, i, abs(modulo(x-a*deltat, xmax-xmin) - xout(i))
  end do
  print*,'     error=', error

    print*,'testing time dependent field. ' 
  ! consider the ode dx/dt = a (t,x) = t*x
  do i = 1, ncells+1
     x = xin(i)
     a_n(i) = 0.0_f64
     a_np1(i) = deltat*x
  end do

  print*, 'testing order 1:'
  order = 1
  call sll_s_implicit_ode_nonuniform( order, deltat, xin, ncells, sll_p_periodic_ode, xout,  a_n ) 
  error = 0.0_f64
  do i = 1, ncells
     x = xin(i)
     error = max(error, abs(modulo(x*exp(-deltat**2/2), xmax-xmin) - xout(i)))
     !print*, i, x*exp(-deltat**2/2), xout(i)
  end do
  print*,'     error=', error

    print*, 'testing order 2:'
  order = 2
  call sll_s_implicit_ode_nonuniform( order, deltat, xin, ncells, sll_p_periodic_ode, xout,  a_n, a_np1 ) 
  error = 0.0_f64
  do i = 1, ncells
     x = xin(i)
     error = max(error, abs(modulo(x*exp(-deltat**2/2.0_f64), xmax-xmin) - xout(i)))
     !print*, i, modulo(xmin+(i-1)*deltax - deltat, xmax-xmin), xout(i)
  end do
  print*,'     error=', error
  ! test the ode dx / dt = sin(x) -> solution is x(t) = x(0) * exp(t)
  do i = 1, ncells+1
     x = xin(i)
     a_n(i) = sin(x)
     a_np1(i) = sin(x)
  end do

  print*, 'testing order 1 on a(x) = sin(x)'
  order = 1
  call sll_s_implicit_ode_nonuniform( order, deltat, xin, ncells, sll_p_periodic_ode, xout,  a_n ) 
  error = 0.0_f64
  do i = 1, ncells
     x = xin(i)
     error = max(error,abs(2.0_f64*atan(tan(x/2)*exp(-deltat))-xout(i)))
  end do
  print*,'     error=', error

    print*, 'testing order 2 on a(x) = sin(x)'
  order = 2
  call sll_s_implicit_ode_nonuniform( order, deltat, xin, ncells,  sll_p_periodic_ode, xout,  a_n, a_np1 ) 
  error = 0.0_f64
  do i = 1, ncells
     x = xin(i)
     error = max(error,abs(2.0_f64*atan(tan(x/2.0_f64)*exp(-deltat))-xout(i)))
     !print*, i, 2*atan(tan(x/2)*exp(-deltat)), xout(i)
  end do
  print*,'     error=', error

  print *, 'Successful, exiting program.'
  
end program test_ode_solvers
