!***************************************************************************
!
! Selalib 2012     
! Module: sll_odd_order_splines_hermite.F90
!
!> @brief 
!> Selalib odd order Hermite splines interpolator
!
!> Start date: July 26, 2012
!> Last modification: July 27, 2012
!   
!> @authors                    
!> Aliou DIOUF (aliou.l.diouf@inria.fr)
!                                  
!***************************************************************************

module sll_odd_order_splines_hermite

#include "sll_memory.h"
#include "sll_working_precision.h"
  implicit none

  type odd_order_splines_hermite_plan
    sll_int32                             :: n
    sll_int32                             :: order
    sll_real64                            :: xmin
    sll_real64                            :: xmax
    sll_real64, dimension(:), allocatable :: coeffs
  end type odd_order_splines_hermite_plan


contains 


  function new_odd_order_splines_hermite(n,order,xmin,xmax,f)result(plan)

    sll_real64, dimension(:)                      :: f
    type(odd_order_splines_hermite_plan), pointer :: plan
    sll_int32                                     :: ierr, n, order
    sll_real64                                    :: xmin
    sll_real64                                    :: xmax

    if (mod(order,2) == 0) then
       print*, 'This needs to run with odd spline order'
       print*, 'Exiting...'
       stop
    endif

    SLL_ALLOCATE(plan, ierr)
    SLL_ALLOCATE(plan%coeffs(size(f)), ierr)

    plan%n = n
    plan%order = order
    plan%xmin = xmin
    plan%xmax = xmax
    call compute_coeffs_hermite(f, plan)

  end function new_odd_order_splines_hermite


  subroutine compute_coeffs_hermite(f, plan)

    ! f is the vector of the values of the function 
    !  in the nodes of the mesh

    sll_real64, dimension(:)                      :: f
    type(odd_order_splines_hermite_plan), pointer :: plan
    sll_real64                                    :: xmin, xmax, h
    sll_int32                                     :: n, order, i, j
    sll_real64, dimension(plan%order/2+1)         :: v
    sll_real64, dimension(size(f),size(f))        :: A, AB
    sll_int32                                     :: KD, LDAB, ierr
    
    order = plan%order
    xmin = plan%xmin
    xmax = plan%xmax
    h = (xmax-xmin)/plan%n

    do j=1,order/2+1
        ! For a given line i, we have:
        ! 0 ... 0 A(i,i-order) ... A(i,i) ... A(i+order) 0 ... 0
        ! As A is symmetric, we just need to storage A(i,i)...A(i+order)
        v(j) = B(order, plan%n-( order/2 + 2 - j ), xmax, plan)
        ! v(i) = B(order, n-(order/2+1-j+1), xmax, plan)
    enddo
print*,v(j)
    ! Solve the linear system

    n = size(f)
    KD = n - 1
    LDAB = KD + 1

    A = 0.
    do i=1,n
       do j= 0, order/2
          if ( i+j<=n ) then
             A(i,i+j) = v(j+1) 
          endif
          if ( i-j>0 ) then
            A(i,i-j) = A(i-j,i)
          endif
       enddo
    enddo

    do j=1,n
       do i=j,min(n,j+KD)
          AB(1+i-j,j) = A(i,j)
       enddo
    enddo
  
    ! Cholesky factorization
    call DPBTRF( 'L', n, KD, AB, LDAB, ierr )
    ! Solve the linear system with Cholesky factorization
    plan%coeffs = f
    call DPBTRS( 'L', n, KD, 1, AB, LDAB, plan%coeffs, n, ierr )

  end subroutine compute_coeffs_hermite


  sll_real64 recursive function B(j, i, x, plan) result(res)

    sll_real64                                    :: x, xmin, xmax, h
    sll_int32                                     :: n, j, i
    type(odd_order_splines_hermite_plan), pointer :: plan

    xmin = plan%xmin
    xmax = plan%xmax
    n    = plan%n
    h    = (xmax-xmin)/n

    !             x-t(i)                       t(i+j+1)-x
    ! B[j,i](x) = ----------- * B[j-1,i](x) + ----------------* B[j-1,i+1](x)
    !           t(i+j)-t(i)                  t(i+j+1)-t(i+1)
    
    ! And
  
    ! B[0,i] = 1 if t(i) <= x < t(i+1), and 0 otherwise.
    
    ! t(i) = xmin + i*h

    if (j/=0) then
                                            
      res = ((x-xmin-i*h)/(j*h)) * B(j-1, i, x, plan) + &
          ((xmin+(i+j+1)*h-x)/(j*h)) * B(j-1, i+1, x, plan)
   
    else

      if ( ( xmin + i*h <= x ) .and. ( x < xmin + (i+1)*h ) ) then
        res = 1.d0
      else
        res = 0.d0
      endif

    endif

  end function B


  function spline_hermite( x, plan) result(s) ! The interpolator spline function

    sll_real64                                    :: x, xmin, xmax
    sll_real64                                    :: h, s
    type(odd_order_splines_hermite_plan), pointer :: plan
    sll_int32                                     :: n, j, left, order

    xmin = plan%xmin
    xmax = plan%xmax
    n = plan%n
    order = plan%order
    h = (xmax-xmin)/n
    left = int((x-xmin)/h)!Determine the leftmost support index 'i' of x

    s = 0.d0
    do j=left-order,left
       if ( (j>=-order) .and. (j<=n) ) then
          s = s + plan%coeffs(j+order) * B(order, j, x, plan)
       endif
    enddo

  end function spline_hermite


  subroutine delete_odd_order_splines_hermite(plan)

    type(odd_order_splines_hermite_plan), pointer :: plan
    sll_int32                                     :: ierr

    SLL_DEALLOCATE_ARRAY(plan%coeffs, ierr)
    SLL_DEALLOCATE_ARRAY(plan, ierr)

  end subroutine delete_odd_order_splines_hermite

end module sll_odd_order_splines_hermite
