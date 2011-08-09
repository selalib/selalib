program unit_test
#include "sll_working_precision.h"
#include "sll_memory.h"
  use ode_solvers
  implicit none
  
  sll_int32 :: i, ncells, ierr, order
  sll_real64 :: xmin, xmax, deltax, deltat, x, error
  sll_real64, dimension(:), pointer :: a_n, a_np1, xout

  print*, 'checking implicit_ode'
  ncells = 100
  xmin = 0.0_f64
  xmax = 1.0_f64
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


  print *, 'Successful, exiting program.'
  
end program unit_test
