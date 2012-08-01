!***************************************************************************
!
! Selalib 2012     
! Module: test_sll_odd_order_splines_hermite.F90
!
!> @brief 
!> Selalib odd order Hermite splines interpolator tester
!
!> Start date: July 26, 2012
!> Last modification: August 1, 2012
!   
!> @authors                    
!> Aliou DIOUF (aliou.l.diouf@inria.fr)
!                                  
!***************************************************************************

program test_sll_odd_order_splines_hermite

#include "sll_memory.h"
#include "sll_working_precision.h"
use numeric_constants
use sll_odd_order_splines_hermite

  type(odd_order_splines_hermite_plan), pointer :: plan
  sll_real64, dimension(:), allocatable         :: f
  sll_real64                                    :: xmin, xmax, h
  sll_real64                                    :: x, x_rand, y
  sll_real64                                    :: err1, err2
  sll_real64                                    :: maxi, bound
  sll_int32                                     :: n, nmax, i
  sll_int32                                     :: order, ierr, ok

  nmax = 1000
  order = 5

  print*,'Testing quintic_spline_hermite...'

  do n=1,nmax

     call random_number(xmin)
     call random_number(xmax)
     xmax = xmin + xmax + 1.d0
     h = (xmax-xmin)/n
     bound = 8.*sll_pi**4 / ( 2.*real(n,f64)**3*sqrt(h**3) )

     SLL_ALLOCATE( f(n+order+1), ierr)

     do i=-2,n+3
        x = xmin + i*h
        f(i+3) = 2.*sll_pi*x/(n*h) - sin(2.*sll_pi*(x-xmin)/(n*h)) 
     enddo

     plan => new_odd_order_splines_hermite(n, order, xmin, xmax, f)

     err1 = 0.d0
     err2 = 0.d0
     ok = 1

     do i=0,n

        x = xmin + i*h
        err1 = err1 + abs( f(i+3) - spline_hermite(x, plan) )

        call random_number(x_rand)
        x_rand = xmin + n*h*x_rand
        y = ( 2.*sll_pi*x_rand/(n*h) - sin( 2.*sll_pi*(x_rand-xmin)/(n*h) ) )

        maxi = abs( y - spline_hermite(x_rand, plan) )
        if (err2 < maxi) then
           err2 = maxi
        endif

        if ( err2 > bound ) then
           print*, 'Program stopped by iteration number', n
           print*, y, spline_hermite(x_rand, plan), err2, bound
           ok = 0
           print*, 'Exciting...'
           stop
        endif

     enddo

     print*, 'Average absolute errors: ', err1/(n+1), err2, ' bound:', bound

    if ( (err1/(n+1) > 1./10**9) ) then    
       print*, 'Program stopped by iteration number', n
       ok = 0
       print*, 'Exciting...'
       stop
    endif

    SLL_DEALLOCATE_ARRAY(f, ierr)
    call delete_odd_order_splines_hermite(plan)

  enddo

  if (ok==1) then
     print*, 'sll_odd_order_splines_hermite: PASS'
  endif

end program test_sll_odd_order_splines_hermite
