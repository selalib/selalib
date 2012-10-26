!***************************************************************************
!
! Selalib 2012     
! Module: test_sll_odd_degree_splines.F90
!
!> @brief 
!> Selalib odd degree splines interpolator tester
!
!> Start date: July 26, 2012
!> Last modification: October 25, 2012
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

  type(odd_degree_splines_uniform_plan), pointer :: plan
  sll_real64, dimension(:), allocatable          :: f
  sll_real64                                     :: xmin, xmax, h
  sll_real64                                     :: x, y, mu
  sll_real64                                     :: err1, err2, norm
  ! err1 is the relative error in the nodes
  ! err2 is the relative error in random points
  sll_int32                                      :: num_pts, i, pow, pow_max
  sll_int32                                      :: degree, degree_max, ierr, ok

  print*,' '
  print*,'Testing odd degree splines module...'
  print*,' '

  xmin = -10.d0
  xmax = 10.d0
  mu = (xmin+xmax)/2
  degree_max = 11
  pow_max = 2
  ok = 1

  do degree=1,degree_max,2

     print*, 'DEGREE = ', degree

     do pow=0,pow_max

        num_pts = degree*10**pow + 1
        h = (xmax-xmin)/(num_pts-1)

        SLL_ALLOCATE( f(num_pts), ierr)

        do i=1,num_pts
           x = xmin + (i-1)*h
           f(i) = exp( - ( x - mu )**2  )
        enddo

        plan => new_odd_degree_splines_uniform(num_pts, degree, xmin, xmax)
        call compute_coeffs_uniform(f, plan)

        err1 = 0.d0
        err2 = 0.d0
        norm = 0.d0

        do i=1,num_pts

           x = xmin + (i-1)*h
           err1 = err1 + ( f(i) - odd_degree_splines_interpolator_uniform_value( &
                                                                     x, plan) )**2

           call random_number(x)
           x = xmin + x*(xmax-xmin) ! generate randomly x in [xmin, xmax]            
           y = exp( - ( x - mu )**2  )
           err2 = err2 + ( y - odd_degree_splines_interpolator_uniform_value( &
                                                                 x, plan) )**2
           norm = norm + y*y; 

        enddo

        err1 = sqrt( err1/sum( f(1:num_pts)**2 ) )
        err2 = sqrt( err2/norm )
        print*, 'Nb_points =', num_pts, ', err1 = ', err1, ', err2 =', err2

        if ( (err1 >= 1.e-13) ) then
           print*, 'sll_odd_degree_splines: FAILED'
           ok = 0
           print*, 'Exiting...'
           stop
        endif
            
        SLL_DEALLOCATE_ARRAY(f, ierr)
        call delete_odd_degree_splines_uniform(plan)

     enddo

     print*, ' '

  enddo

  if (ok==1) then
     print*, 'sll_odd_degree_splines: PASSED'
    print*, ' '
  endif


end program test_sll_odd_degree_splines
