!***************************************************************************
!
! Selalib 2012     
! Module: test_sll_odd_degree_splines.F90
!
!> @brief 
!> Selalib odd degree splines interpolator tester
!
!> Start date: July 26, 2012
!> Last modification: September 11, 2012
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

  type(odd_degree_splines_plan), pointer :: plan1, plan2
  sll_real64, dimension(:), allocatable  :: f1, f2
  sll_real64                             :: xmin, xmax, h
  sll_real64                             :: x, y, mu, t0
  sll_real64                             :: err1, err2, norm
  sll_int32                              :: n, nmax
  sll_int32                              :: i, left, m
  sll_int32                              :: degree, degree_max, ierr, ok
    sll_real64, dimension(8)   :: b=0.d0
  print*,'Testing odd degree splines module...'

  nmax = 4
  degree_max = 3
  ok = 1
  xmin = -10.d0
  xmax = 10.d0

  do degree=3,degree_max,2

     m = degree/2+1

     do n=nmax,nmax

        h = (xmax-xmin)/n
        mu = (xmin+xmax)/2

        SLL_ALLOCATE( f1(n+degree+1), ierr)
        SLL_ALLOCATE( f2(n+degree+1), ierr)
        f1 = 0.d0
        f2 = 0.d0

        do i=0,n

           x = xmin + i*h
           f1(i+m) = exp( - ( x - mu )**2  )     
           f2(i+m) = x*(x-xmin)*(xmax-x)

        enddo

        plan1 => new_odd_degree_splines(n, degree, xmin, xmax, f1)
        plan2 => new_odd_degree_splines(n, degree, xmin, xmax, f2)

        err1 = 0.d0
        err2 = 0.d0
        norm = 0.d0

        do i=0,n

           x = xmin + i*h
           err1 = err1 + ( f1(i+m) - spline(x, plan1) )**2

           !call random_number(x)
           !x = x + 0.5d0*h!x*(xmax-xmin) ! generate randomly x in [xmin, xmax]            
   !t0 = (x-xmin)/h
    !left = int(t0) ! Determine the leftmost support index 'i' of x
   ! t0 = t0 - left ! compute normalized_offset

    b(i+1:i+degree+1) = uniform_b_splines_at_x( degree, 0.d0 )
           y = x*(x-xmin)*(xmax-x)
           err2 = err2 + ( y - sum(plan2%coeffs*b) )**2
           norm = norm + y*y; 
print*, y, sum(plan2%coeffs*b),spline(x, plan2)
        enddo

        err1 = sqrt( err1/sum( f1(m:n+m)**2 ) )
        err2 = sqrt( err2/norm )
        print*, 'Relative errors:', err1, err2

        if ( (err1 >= 1.e-12) .or. (err2 >= 1.e-12) ) then
           print*, 'sll_odd_degree_splines: FAIL'
           print*, 'Degree =', degree
           print*, 'nb_points =', n+1
           ok = 0
           print*, 'Exiting...'
           stop
        endif
            
        SLL_DEALLOCATE_ARRAY(f2, ierr)
        SLL_DEALLOCATE_ARRAY(f1, ierr)
        call delete_odd_degree_splines(plan1)
        call delete_odd_degree_splines(plan2)

     enddo

  enddo

  if (ok==1) then
     print*, 'sll_odd_degree_splines: PASS'
  endif


end program test_sll_odd_degree_splines
