!*******************************************************************
!
! Selalib 2012     
! Module: test_sll_quintic_splines.F90
!
!> @brief 
!> Selalib quintic splines interpolator
!
!> Last modification: September 20, 2012
!   
!> @authors                    
!> Aliou DIOUF (aliou.l.diouf@inria.fr)
!                                  
!*******************************************************************

program test_sll_quintic_splines

#include "sll_memory.h"
#include "sll_working_precision.h"
use numeric_constants
use sll_quintic_splines
use arbitrary_degree_splines
  implicit none
  
  type(quintic_splines_plan_uniform), pointer   :: plan
  sll_real64, dimension(:), allocatable :: f
  sll_real64                            :: xmin, xmax, h
  sll_real64                            :: x, y, mu
  sll_real64                            :: err1, err2, norm
  ! err1 is the relative error in the nodes
  ! err2 is the relative error in random points
  sll_int32                             :: n, i, pow, pow_max
  sll_int32                             :: ierr, ok
  
  print*,' '
  print*,'Testing quintic splines module...'
  print*,' '

  xmin = -10.d0
  xmax = 10.d0
  pow_max = 4
  ok = 1

  do pow=0,pow_max

     n = 5*10**pow
     h = (xmax-xmin)/n
     mu = (xmin+xmax)/2

     SLL_ALLOCATE( f(n+6), ierr)
     f = 0.d0

     do i=0,n
        x = xmin + i*h
        f(i+3) = exp( - ( x - mu )**2  )
     enddo

     plan => new_quintic_splines(n, xmin, xmax, f)

     err1 = 0.d0
     err2 = 0.d0
     norm = 0.d0

     do i=0,n

        x = xmin + i*h
        err1 = err1 + ( f(i+3) - quintic_splines(x, plan) )**2

        call random_number(x)
        x = xmin + x*(xmax-xmin) ! generate randomly x in [xmin, xmax]            
        y = exp( - ( x - mu )**2  )
        err2 = err2 + ( y - quintic_splines(x, plan) )**2
        norm = norm + y*y; 

     enddo

     err1 = sqrt( err1/sum( f(3:n+3)**2 ) )
     err2 = sqrt( err2/norm )
     print*, 'Nb_points =', n+1, ', err1 = ', err1, ', err2 =', err2

     if ( (err1 >= 1.e-13) ) then
        print*, 'sll_quintic_splines: FAIL'
        ok = 0
        print*, 'Exiting...'
        stop
     endif
            
     SLL_DEALLOCATE_ARRAY(f, ierr)
     call delete_quintic_splines(plan)

  enddo

  print*, ' '

  if (ok==1) then
     print*, 'sll_quintic_splines: PASSED'
    print*, ' '
  endif

end program test_sll_quintic_splines
