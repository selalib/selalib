!*******************************************************************
!
! Selalib 2012     
! Program: test_sll_quintic_splines.F90
!
!> @brief 
!> Selalib quintic splines interpolator
!
!> Last modification: October 03, 2012
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
  
  type(quintic_splines_uniform_plan), pointer :: plan1
  !type(quintic_splines_non_uni_plan), pointer :: plan2
  sll_real64, dimension(:), allocatable       :: f1!, f2, x2
  ! f1: test function for uniform case
  ! f2: test function for non uniform case
  ! x2: array coordinates for non uniform case
  sll_real64                                  :: xmin, xmax, h
  sll_real64                                  :: x, y, mu, norm
  sll_real64                                  :: err11, err12!, err21, err22
  ! err11, err21 are the relative errors in the nodes for f1, f2 resp
  ! err12, 22 are the relative errors in random points for f1, f2 resp
  sll_int32                                   :: num_pts, i, pow, pow_max
  sll_int32                                   :: ierr, ok
  
  print*,' '
  print*,'Testing quintic splines module...'
  print*,' '

  xmin = -10.d0
  xmax = 10.d0
  mu = (xmin+xmax)/2
  pow_max = 4
  ok = 1

  do pow=0,pow_max

     num_pts = 5*10**pow + 1
     h = (xmax-xmin) / (num_pts-1)

     SLL_ALLOCATE( f1(num_pts), ierr)
     !SLL_ALLOCATE( f2(num_pts), ierr)
     !SLL_ALLOCATE( x2(num_pts), ierr)

     do i=0,num_pts-1
        x = xmin + i*h
        f1(i+1) = exp( - ( x - mu )**2  )
     enddo

     !x2(1) = xmin
     !f2(1) = exp( - ( x2(1) - mu )**2  )
     !do i=2,num_pts-1
        !call random_number(x)
        !x2(i) = x2(i-1) + x*( xmax - x2(i-1))
        !f2(i) = exp( - ( x2(i) - mu )**2  )
     !enddo
     !x2(num_pts) = xmax

     plan1 => new_quintic_splines_uniform(num_pts, xmin, xmax)
     call compute_quintic_coeffs_uniform(f1, plan1)

     !plan2 => new_quintic_splines_non_uni(x2)
     !compute_quintic_coeffs_non_uni(f2, plan2)

     err11 = 0.d0
     err12 = 0.d0
     !err21 = 0.d0
     !err22 = 0.d0
     norm = 0.d0

     do i=0,num_pts-1

        x = xmin + i*h
        err11 = err11 + ( f1(i+1) - quintic_splines(x, plan1) )**2
        !err21 = err21 + ( f2(i+1) - quintic_splines(x2(i+1), plan2) )**2
!print*, quintic_splines(x, plan2)
        call random_number(x)
        x = xmin + x*(xmax-xmin) ! generate randomly x in [xmin, xmax]            
        y = exp( - ( x - mu )**2  )
        err12 = err12 + ( y - quintic_splines(x, plan1) )**2
        !err22 = err22 + ( y - quintic_splines(x, plan2) )**2
        norm = norm + y*y; 

     enddo

     err11 = sqrt( err11/sum( f1**2 ) )
     !err21 = sqrt( err21/sum( f2**2 ) )
     err12 = sqrt( err12/norm )
     !err22 = sqrt( err22/norm )

     print*, 'Nb_points =', num_pts
     print*, 'Uniform case: err11 = ', err11, ', err12 =', err12
     !print*, 'Non uniform case: err21 = ', err21, ', err22 =', err22
     print*, '--------------------------------------------------', &
             '------------------------------------'

     !if ( (err11 >= 1.e-13) .or. (err21 >= 1.e-13) ) then
     if ( (err11 >= 1.e-13) ) then
        print*, 'sll_quintic_splines: FAIL'
        ok = 0
        print*, 'Exiting...'
        stop
     endif

     SLL_DEALLOCATE_ARRAY(f1, ierr)          
     !SLL_DEALLOCATE_ARRAY(f2, ierr)
     !SLL_DEALLOCATE_ARRAY(x2, ierr)
     call delete_quintic_splines(plan1)
     !call delete_quintic_splines(plan2)

  enddo

  if (ok==1) then
    print*, ' '
     print*, 'sll_quintic_splines: PASSED'
    print*, ' '
  endif

end program test_sll_quintic_splines
