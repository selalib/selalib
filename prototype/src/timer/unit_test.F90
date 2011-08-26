program main
#include "sll_timer.h"
!  use sll_timer
!  use iso_c_binding
  implicit none
  integer(8) :: i, j,k
  ! etime stuff
  real, dimension(2) :: tarray1, tarray2
  real :: result
  ! sll stuff
  real(c_double), dimension(:), allocatable :: dt
  real(c_double) :: acc
  real(c_double) :: min, max
  time_mark      :: tmark
!  type(c_ptr)    :: tmark

  min =  50.0
  max = -50.0

#define ITERATIONS 10

  allocate(dt(ITERATIONS))
  tmark = set_time_mark()
  print *, 'initialized time marker...'
  print *, 'is tmark associated?:', c_associated(tmark)
!  print *, 'printing tmark: ', tmark
  ! main testing loop 
  do i=1,ITERATIONS
     tmark = reset_time_mark(tmark)
     do k=1,100000    ! Just a delay
        j = i * i - i
     end do
!     print *, 'tmark = ', tmark
!     call flush()
     dt(i) = time_elapsed_since(tmark) 
  end do

  do i=1,ITERATIONS
     write (*, '(a,es20.12)') 'delta = ', dt(i)
   end do

  do i=1,ITERATIONS
     if( dt(i) < min ) then
        min = dt(i)
     end if
     if( dt(i) > max ) then
        max = dt(i)
     end if
  end do
  print *, 'min = ', min
  print *, 'max = ', max
  acc = 0.0
  do i=1, ITERATIONS
     acc = acc + dt(i)
  end do
  print *, 'average measurement = ', acc/real(ITERATIONS)

  call ETIME(tarray1, result)
!  do i=1,100000000    ! Just a delay
!     j = i * i - i
!  end do
  call ETIME(tarray2, result)

  print *, 'after first call to ETIME: ', result
  print *, tarray1(1)
  print *, tarray1(2)
  print *, 'after second call... ', result
  print *, tarray2(1)
  print *, tarray2(2)
end program main
     
