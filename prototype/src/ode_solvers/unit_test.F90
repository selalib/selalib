program unit_test
#include "sll_working_precision.h"
#include "sll_memory.h"
  use ode_solvers
  implicit none
  
  sll_int32 :: i, ncells, ierr, order, nsubsteps
  sll_real64 :: xmin, xmax, deltax, deltat, x, error
  sll_real64, dimension(:), pointer :: a_n, a_np1, xout

  print*, 'checking implicit_ode'
  ncells = 100
  xmin = 0.0_f64
  xmax = sll_pi
  deltax = (xmax-xmin) / ncells
  
  deltat = 0.01_f64

  SLL_ALLOCATE(a_n(ncells+1), ierr)
  SLL_ALLOCATE(a_np1(ncells+1), ierr)
  SLL_ALLOCATE(xout(ncells+1), ierr)
  print*,'testing constant field. ' 
  print*,'     deltax = ',deltax, 'deltat = ',deltat

  ! consider the ode dx/dt = a (t,x) = t*x
  x = xmin
  do i = 1, ncells+1
     a_n(i) = 0.0_f64
     a_np1(i) = deltat*x
     x = x + deltax
  end do

  print*, 'testing order 1:'
  order = 1
  call implicit_ode( order, deltat, xmin, ncells, deltax, PERIODIC_ODE, xout,  a_n ) 
  error = 0.0_f64
  do i = 1, ncells
     x = xmin+(i-1)*deltax
     error = max(error, abs(modulo(x*exp(-deltat**2/2), xmax-xmin) - xout(i)))
     !print*, i, x*exp(-deltat**2/2), xout(i)
  end do
  print*,'     error=', error

    print*, 'testing order 2:'
  order = 2
  call implicit_ode( order, deltat, xmin, ncells, deltax, PERIODIC_ODE, xout,  a_n, a_np1 ) 
  error = 0.0_f64
  do i = 1, ncells
     x = xmin+(i-1)*deltax
     error = max(error, abs(modulo(x*exp(-deltat**2/2), xmax-xmin) - xout(i)))
     !print*, i, modulo(xmin+(i-1)*deltax - deltat, xmax-xmin), xout(i)
  end do
  print*,'     error=', error
  ! test the ode dx / dt = sin(x) -> solution is x(t) = x(0) * exp(t)
  x = xmin
  do i = 1, ncells+1
     a_n(i) = sin(x)
     a_np1(i) = sin(x)
     x = x + deltax
  end do

  print*, 'testing order 1 on a(x) = sin(x)'
  order = 1
  call implicit_ode( order, deltat, xmin, ncells, deltax, PERIODIC_ODE, xout,  a_n ) 
  error = 0.0_f64
  do i = 1, ncells
     x = xmin+(i-1)*deltax
     error = max(error,abs(xmin+modulo(2*atan(tan(x/2)*exp(-deltat))-xmin, xmax-xmin)-xout(i)))
     !print*, i, x*exp(-deltat**2/2), xout(i)
  end do
  print*,'     error=', error

    print*, 'testing order 2 on a(x) = sin(x)'
  order = 2
  call implicit_ode( order, deltat, xmin, ncells, deltax, PERIODIC_ODE, xout,  a_n, a_np1 ) 
  error = 0.0_f64
  do i = 1, ncells
     x = xmin+(i-1)*deltax
     error = max(error,abs(xmin+modulo(2*atan(tan(x/2)*exp(-deltat))-xmin, xmax-xmin)-xout(i)))
     !print*, i, 2*atan(tan(x/2)*exp(-deltat)), xout(i)
  end do
  print*,'     error=', error

  print*, 'testing RK2:'
  nsubsteps = 1
  call rk2( nsubsteps, deltat, xmin, ncells, deltax, PERIODIC_ODE, xout,  a_n ) 
  ! compute error
  error = 0.0_f64
  x = xmin
  do i = 1, ncells
     error = max(error,abs(xmin+modulo(2*atan(tan(x/2)*exp(deltat))-xmin, xmax-xmin)-xout(i)))
     !print*, i, 2*atan(tan(x/2)*exp(deltat)), xout(i), abs(2*atan(tan(x/2)*exp(deltat))-xout(i))
     x = x + deltax
  end do
  print*, '     error=', error

  ! test the ode dx / dt = 1 -> solution is x(t) = x(0) + t
  do i = 1, ncells+1
     a_n(i) = 1.0_f64/(1+0.2*cos(2*x))
  end do

  print*, 'testing order 1 on a(x) = 1'
  order = 1
  call implicit_ode( order, deltat, xmin, ncells, deltax, PERIODIC_ODE, xout,  a_n ) 
  error = 0.0_f64
  do i = 1, ncells+1
     x = xmin+(i-1)*deltax
     error = max(error,abs(xmin+modulo(x-deltat-xmin, xmax-xmin)-xout(i)))
     !print*, i, x*exp(-deltat**2/2), xout(i)
  end do
  print*,'     error=', error

  print *, 'Successful, exiting program.'
  
end program unit_test
