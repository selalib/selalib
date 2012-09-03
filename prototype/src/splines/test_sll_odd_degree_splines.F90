!***************************************************************************
!
! Selalib 2012     
! Module: test_sll_odd_degree_splines.F90
!
!> @brief 
!> Selalib odd degree splines interpolator tester
!
!> Start date: July 26, 2012
!> Last modification: September 04, 2012
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

  type(odd_degree_splines_plan), pointer :: plan1, plan2
  sll_real64, dimension(:), allocatable :: f1, f2
  sll_real64                            :: xmin, xmax, h
  sll_real64                            :: x, y, mu
  sll_real64                            :: err1, err2, norm
  sll_int32                             :: n, nmax
  sll_int32                             :: i, left, m
  sll_int32                             :: degree, degree_max, ierr, ok

  nmax = 10
  degree_max = 19
  ok = 1

  print*,'Testing quintic_spline...'

  do degree=1,degree_max,2

     m = degree/2+1

     do n=1,nmax

        xmin = 0.d0
        xmax = 10.d0
        h = (xmax-xmin)/n
        mu = ( (xmin-degree*h) + (xmin+(n+degree+1)*h) ) / 2

        SLL_ALLOCATE( f1(n+degree+1), ierr)
        SLL_ALLOCATE( f2(n+degree+1), ierr)

        do i=degree/2+1-degree,n+degree/2+1

           x = xmin + i*h
           f1(i+m) = exp( - .5*( x - mu )**2  )
      
           f2(i+m) = 0.d0
           do j=i-degree,i
              if ( (j>=-degree) .and. (j<=n) ) then
                 f2(i+m) = f2(i+m) + B_test(degree, j, x, xmin, h)
              endif
           enddo

        enddo

        plan1 => new_odd_degree_splines(n, degree, xmin, xmax, f1)
        plan2 => new_odd_degree_splines(n, degree, xmin, xmax, f2)

        err1 = 0.d0
        err2 = 0.d0
        norm = 0.d0

        do i=0,n

           x = xmin + i*h
           err1 = err1 + ( f1(i+m) - spline(x, plan1) )**2

           call random_number(x)
           x = xmin + x*(xmax-xmin) ! x in [xmin, xmax]
           left = int((x-xmin)/h)!Determine the leftmost support index 'i' of x

           y = 0.d0
           do j=left-degree,left
              if ( (j>=-degree) .and. (j<=n) ) then
                 y = y + B_test(degree, j, x, xmin, h)
              endif
           enddo

           norm = norm + y*y
           err2 = err2 + ( y - spline(x, plan2) )**2

        enddo

        err1 = sqrt( err1/sum( f1(m:n+m)**2 ) )
        err2 = sqrt( err2/norm )
        print*, 'Relative errors:', err1, err2

        if ( (err1 >= 1.e-12) .or. (err2 >= 1.e-12) ) then
           print*, 'Program stopped by iteration number', n
           ok = 0
           print*, 'Exiting...'
           stop
        endif
    
        SLL_DEALLOCATE_ARRAY(f1, ierr)
        SLL_DEALLOCATE_ARRAY(f2, ierr)
        call delete_odd_degree_splines(plan1)
        call delete_odd_degree_splines(plan2)

     enddo

  enddo

  if (ok==1) then
     print*, 'sll_odd_degree_splines: PASS'
  endif


contains

  sll_real64 recursive function B_test(j, i, x, xmin, h) result(res) ! j:=degree

    sll_real64                                    :: x, xmin, h
    sll_int32                                     :: j, i

    !             x-t(i)                       t(i+j+1)-x
    ! B[j,i](x) = ----------- * B[j-1,i](x) + ----------------* B[j-1,i+1](x)
    !           t(i+j)-t(i)                  t(i+j+1)-t(i+1)
    
    ! And
  
    ! B[0,i] = 1 if t(i) <= x < t(i+1), and 0 otherwise.
    
    ! t(i) = xmin + i*h

    if (j/=0) then
                                            
      res = ((x-xmin-i*h)/(j*h)) * B_test(j-1, i, x, xmin, h) + &
        ((xmin+(i+j+1)*h-x)/(j*h)) * B_test(j-1, i+1, x, xmin, h)
   
    else

      if ( ( xmin + i*h <= x ) .and. ( x < xmin + (i+1)*h ) ) then
        res = 1.d0
      else
        res = 0.d0
      endif

    endif

  end function B_test


end program test_sll_odd_degree_splines
