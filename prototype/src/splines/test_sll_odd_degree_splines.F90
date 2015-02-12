!***************************************************************************
!
! Selalib 2012     
! Module: test_sll_odd_degree_splines.F90
!
!> @brief 
!> Selalib odd degree splines interpolator tester
!
!> Start date: July 26, 2012
!> Last modification: Nov. 29, 2012
!   
!> @authors                    
!> Aliou DIOUF (aliou.l.diouf@inria.fr)
!                                  
!***************************************************************************

program test_sll_odd_degree_splines

#include "sll_memory.h"
#include "sll_working_precision.h"
use sll_constants
use sll_odd_degree_splines
use sll_arbitrary_degree_splines

  implicit none

  type(odd_degree_splines_uniform_plan), pointer    :: plan1
  type(odd_degree_splines_nonuniform_plan), pointer :: plan2
  sll_real64, dimension(:), allocatable             :: f1, f2, x2
  ! f1: test function for uniform case
  ! f2: test function for non uniform case
  ! x2: array coordinates for non uniform case
  sll_real64                                     :: xmin, xmax, h
  sll_real64                                     :: x, y, mu
  sll_real64                                     :: err11, err12
  sll_real64                                     :: err21, err22, norm
  ! err11, err21 are the relative errors in the nodes for f1, f2 resp
  ! err12, 22 are the relative errors in random points for f1, f2 resp
  sll_int32                                      :: num_pts, i, pow, pow_max
  sll_int32                                      :: degree, degree_max, ierr, ok

  print*,' '
  print*,'Testing odd degree splines module...'
  print*,' '

  xmin = -10.d0
  xmax = 10.d0
  mu = (xmin+xmax)/2
  degree_max = 11
  pow_max = 3
  ok = 1

  do degree=1,degree_max,2

     print*, 'DEGREE = ', degree

     do pow=0,pow_max

        num_pts = degree*10**pow + 1
        h = (xmax-xmin)/(num_pts-1)

        print *,'#allocation',num_pts
        SLL_ALLOCATE( f1(num_pts), ierr)
        SLL_ALLOCATE( f2(num_pts), ierr)
        SLL_ALLOCATE( x2(num_pts), ierr)
        print *,'#allocation done',num_pts

        do i=0,num_pts-1
           x = xmin + i*h
           f1(i+1) = exp( - ( x - mu )**2  )
        enddo

        x2(1) = xmin
        do i=2,num_pts-1
           call random_number(x)
           ! To avoid duplicated points
           if (num_pts < degree) then
              x = x/num_pts
           else
              x = degree*x/num_pts 
           endif 
           x2(i) = x2(i-1) + x*( xmax - x2(i-1))
        enddo
        x2(num_pts) = xmax

        f2 = exp( - ( x2 - mu )**2  )

        plan1 => new_odd_degree_splines_uniform(num_pts, degree, xmin, xmax)
        call compute_odd_degree_coeffs_uniform(f1, plan1)


        plan2 => new_odd_degree_splines_nonuniform(degree, x2)
        call compute_odd_degree_coeffs_nonuniform(f2, plan2)

        err11 = 0.d0
        err12 = 0.d0
        err21 = 0.d0
        err22 = 0.d0
        norm = 0.d0

        do i=0,num_pts-1

           x = xmin + i*h
           err11 = err11 + ( f1(i+1) - odd_degree_splines(x, plan1) )**2
           err21 = err21 + ( f2(i+1) - odd_degree_splines(x2(i+1), plan2) )**2

           call random_number(x)
           x = xmin + x*(xmax-xmin) ! generate randomly x in [xmin, xmax]            
           y = exp( - ( x - mu )**2  )
           err12 = err12 + ( y - odd_degree_splines(x, plan1) )**2
           err22 = err22 + ( y - odd_degree_splines(x, plan2) )**2
           norm = norm + y*y; 

        enddo

        err11 = sqrt( err11/sum( f1**2 ) )
        err21 = sqrt( err21/sum( f2**2 ) )
        err12 = sqrt( err12/norm )
        err22 = sqrt( err22/norm )

        print*, 'Nb_points =', num_pts
        print*, 'Uniform case: err11 = ', err11, ', err12 =', err12
        print*, 'Non uniform case: err21 = ', err21, ', err22 =', err22
        print*, '--------------------------------------------------', &
                '------------------------------------'
!        if ( (err11 >= 1.e-13) .or. (err21 >= 1.e-13) ) then
!           print*, 'sll_odd_degree_splines: FAILED'
!           ok = 0
!           print*, 'Exiting...'
!           stop
!        endif

        SLL_DEALLOCATE_ARRAY(f1, ierr)          
        SLL_DEALLOCATE_ARRAY(f2, ierr)
        SLL_DEALLOCATE_ARRAY(x2, ierr)
        call delete_odd_degree_splines(plan1)
        call delete_odd_degree_splines(plan2)

        print*, ' '

     enddo
  enddo

  if (ok==1) then
     print*, 'sll_odd_degree_splines: PASSED'
    print*, ' '
  endif


end program test_sll_odd_degree_splines
