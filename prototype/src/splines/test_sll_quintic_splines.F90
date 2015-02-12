!*******************************************************************
!
! Selalib 2012     
! Program: test_sll_quintic_splines.F90
!
!> @brief 
!> Selalib quintic splines interpolator
!
!> Last modification: Nov. 29, 2012
!   
!> @authors                    
!> Aliou DIOUF (aliou.l.diouf@inria.fr)
!                                  
!*******************************************************************

program test_sll_quintic_splines

#include "sll_memory.h"
#include "sll_working_precision.h"
use sll_constants
use sll_quintic_splines
use sll_arbitrary_degree_splines

implicit none
  
type(quintic_splines_uniform_plan), pointer    :: plan1
type(quintic_splines_nonuniform_plan), pointer :: plan2
sll_real64, dimension(:), allocatable          :: x1, f1 ! uniform
sll_real64, dimension(:), allocatable          :: x2, f2 ! non uniform

! f1: test function for uniform case
! f2: test function for non uniform case
! x2: array coordinates for non uniform case
sll_real64 :: xmin, xmax, h
sll_real64 :: x, y, norm
sll_real64 :: err11, err12, err21, err22
! err11, err21 are the relative errors in the nodes for f1, f2 resp
! err12, 22 are the relative errors in random points for f1, f2 resp
sll_int32  :: num_pts, i
sll_int32  :: ierr, ok
  
print*,' '
print*,'Testing quintic splines module...'
print*,' '

xmin    = -10.d0
xmax    = 10.d0
ok      = 1

num_pts = 501

h = (xmax-xmin) / (num_pts-1)

SLL_CLEAR_ALLOCATE( f1(1:num_pts), ierr)
SLL_CLEAR_ALLOCATE( f2(1:num_pts), ierr)
SLL_CLEAR_ALLOCATE( x1(1:num_pts), ierr)
SLL_CLEAR_ALLOCATE( x2(1:num_pts), ierr)

do i=1,num_pts
   x1(i) = xmin + (i-1)*h
enddo

f1 = exp( - x1*x1  )

x2(1) = xmin
do i=2,num_pts-1
   call random_number(x)
   x = 5*x/num_pts 
   x2(i) = x2(i-1) + x*( xmax - x2(i-1))
enddo
x2(num_pts) = xmax

f2 = exp( - x2*x2  )

plan1 => new_quintic_splines_uniform(num_pts, xmin, xmax)
call compute_quintic_coeffs_uniform(f1, plan1)




plan2 => new_quintic_splines_nonuniform(x2)
call compute_quintic_coeffs_nonuniform(f2, plan2)

err11 = 0.d0
err12 = 0.d0
err21 = 0.d0
err22 = 0.d0
norm  = 0.d0

do i=1,num_pts

   err11 = err11 + ( f1(i) - quintic_splines(x1(i), plan1) )**2
   err21 = err21 + ( f2(i) - quintic_splines(x2(i), plan2) )**2

   call random_number(x)
   x = xmin + x*(xmax-xmin) ! generate randomly x in [xmin, xmax]            
   y = exp( - x*x  )
   err12 = err12 + ( y - quintic_splines(x, plan1) )**2
   err22 = err22 + ( y - quintic_splines(x, plan2) )**2
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
if ( (err11 >= 1.e-13) .or. (err21 >= 1.e-13) ) then
   print*, 'sll_quintic_splines: FAILED'
   ok = 0
   print*, 'Exiting...'
   stop
endif

SLL_DEALLOCATE_ARRAY(f1, ierr)          
SLL_DEALLOCATE_ARRAY(f2, ierr)
SLL_DEALLOCATE_ARRAY(x1, ierr)
SLL_DEALLOCATE_ARRAY(x2, ierr)
call delete_quintic_splines(plan1)
call delete_quintic_splines(plan2)

if (ok==1) then
   print*, ' '
   print*, 'sll_quintic_splines: PASSED'
   print*, ' '
endif

end program test_sll_quintic_splines
