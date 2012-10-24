!***************************************************************************
!
! Selalib 2012     
! Module: test_sll_odd_degree_splines.F90
!
!> @brief 
!> Selalib odd degree splines interpolator tester
!
!> Start date: July 26, 2012
!> Last modification: September 20, 2012
!   
!> @authors                    
!> Aliou DIOUF (aliou.l.diouf@inria.fr)
!                                  
!***************************************************************************

program test_sll_odd_degree_splines

#include "sll_memory.h"
#include "sll_working_precision.h"
use numeric_constants
use sll_odd_degree_splines
use arbitrary_degree_splines
  implicit none

  type(odd_degree_splines_plan), pointer :: plan
  sll_real64, dimension(:), allocatable  :: f
  sll_real64                             :: xmin, xmax, h
  sll_real64                             :: x, y, mu
  sll_real64                             :: err1, err2, norm
  ! err1 is the relative error in the nodes
  ! err2 is the relative error in random points
  sll_int32                              :: n, i, m, pow, pow_max
  sll_int32                              :: degree, degree_max, ierr, ok

  print*,' '
  print*,'Testing odd degree splines module...'
  print*,' '

  xmin = -10.d0
  xmax = 10.d0
  degree_max = 11
  pow_max = 2
  ok = 1

  do degree=1,degree_max,2

     print*, 'DEGREE = ', degree
     m = degree/2+1

     do pow=0,pow_max

        n = degree*10**pow
        h = (xmax-xmin)/n
        mu = (xmin+xmax)/2

        SLL_ALLOCATE( f(n+degree+1), ierr)
        f = 0.d0

        do i=0,n
           x = xmin + i*h
           f(i+m) = exp( - ( x - mu )**2  )
        enddo

        plan => new_odd_degree_splines(n, degree, xmin, xmax, f)

        err1 = 0.d0
        err2 = 0.d0
        norm = 0.d0

        do i=0,n

           x = xmin + i*h
           err1 = err1 + ( f(i+m) - odd_degree_spline_interpolator(x, plan) )**2

           call random_number(x)
           x = xmin + x*(xmax-xmin) ! generate randomly x in [xmin, xmax]            
           y = exp( - ( x - mu )**2  )
           err2 = err2 + ( y - odd_degree_spline_interpolator(x, plan) )**2
           norm = norm + y*y; 

        enddo

        err1 = sqrt( err1/sum( f(m:n+m)**2 ) )
        err2 = sqrt( err2/norm )
        print*, 'Nb_points =', n+1, ', err1 = ', err1, ', err2 =', err2

        if ( (err1 >= 1.e-13) ) then
           print*, 'sll_odd_degree_splines: FAILED'
           ok = 0
           print*, 'Exiting...'
           stop
        endif
            
        SLL_DEALLOCATE_ARRAY(f, ierr)
        call delete_odd_degree_splines(plan)

     enddo

     print*, ' '

  enddo

  if (ok==1) then
     print*, 'sll_odd_degree_splines: PASSED'
    print*, ' '
  endif


end program test_sll_odd_degree_splines
